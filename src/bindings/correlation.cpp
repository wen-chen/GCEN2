#include <algorithm>
#include <cstdio>
#include <limits>
#include <string>
#include <thread>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "build_file.hpp"
#include "correlation.hpp"
#include "edge_format.hpp"
#include "network_stream.hpp"
#include "tom.hpp"

namespace py = pybind11;

namespace {

// Edge-list writer with GCEN v1 number formatting (see edge_format.hpp).
class EdgeFormatter {
public:
    explicit EdgeFormatter(std::vector<std::string> names) : names_(std::move(names)) {}

    py::bytes format(py::array_t<int, py::array::c_style | py::array::forcecast> idx,
                     py::array_t<double, py::array::c_style | py::array::forcecast> cor,
                     py::array_t<double, py::array::c_style | py::array::forcecast> pval,
                     py::object fdr, int n_threads, py::object tom) const {
        auto ib = idx.request();
        py::ssize_t n = cor.size();
        if (ib.ndim != 2 || ib.shape[1] != 2 || ib.shape[0] != n || pval.size() != n)
            throw std::invalid_argument("idx must be (n, 2); cor and pval must have length n.");
        const int *ip = static_cast<const int *>(ib.ptr);
        const double *rp = cor.data();
        const double *pp = pval.data();
        py::array_t<double, py::array::c_style | py::array::forcecast> q;
        const double *qp = nullptr;
        if (!fdr.is_none()) {
            q = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(fdr);
            if (!q || q.size() != n)
                throw std::invalid_argument("fdr must have length n.");
            qp = q.data();
        }
        py::array_t<double, py::array::c_style | py::array::forcecast> tv;
        const double *tp = nullptr;
        if (!tom.is_none()) {
            tv = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(tom);
            if (!tv || tv.size() != n)
                throw std::invalid_argument("tom must have length n.");
            tp = tv.data();
        }
        const int n_names = static_cast<int>(names_.size());
        for (py::ssize_t k = 0; k < 2 * n; ++k)
            if (ip[k] < 0 || ip[k] >= n_names)
                throw std::out_of_range("gene index out of range.");

        std::string out;
        {
            py::gil_scoped_release release;
            out = format_edge_lines(names_, static_cast<size_t>(n),
                                    (qp ? kEdgeFdr : 0) | (tp ? kEdgeTom : 0), n_threads,
                                    [&](size_t k) {
                                        return EdgeRow{ip[2 * k], ip[2 * k + 1], rp[k], pp[k],
                                                       qp ? qp[k] : 0.0, tp ? tp[k] : 0.0};
                                    });
        }
        return py::bytes(out);
    }

private:
    std::vector<std::string> names_;
};

py::buffer_info check_matrix(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &matrix) {
    py::buffer_info buf = matrix.request();
    if (buf.ndim != 2)
        throw std::invalid_argument("Input matrix must be 2-dimensional (genes x samples).");
    if (buf.shape[0] > std::numeric_limits<int>::max() ||
        buf.shape[1] > std::numeric_limits<int>::max())
        throw std::invalid_argument("matrix is too large.");
    return buf;
}

NetworkOptions options(double cor, double pval, bool positive_only, bool fdr, bool all_edges,
                       int64_t max_pairs, int64_t keep_limit) {
    NetworkOptions o;
    o.cor = cor;
    o.pval = pval;
    o.positive_only = positive_only;
    o.fdr = fdr;
    o.all_edges = all_edges;
    o.max_pairs_per_block = max_pairs;
    o.fdr_keep_limit = static_cast<size_t>(std::max<int64_t>(0, keep_limit));
    return o;
}

TomSpec tom_spec(const py::object &power, double tom_min, const std::vector<double> &powers) {
    TomSpec t;
    if (power.is_none()) return t;
    t.enabled = true;
    t.power = power.cast<double>();
    if (!(t.power >= 0)) throw std::invalid_argument("tom_power must be >= 0 (0: pick it).");
    t.tom_min = tom_min;
    t.powers = powers;
    for (const double p : powers)
        if (!(p > 0)) throw std::invalid_argument("tom_powers must be positive.");
    return t;
}

py::list fits_list(const std::vector<SoftThresholdFit> &fits) {
    py::list out;
    for (const auto &f : fits)
        out.append(py::make_tuple(f.power, f.r2, f.slope, f.mean_k, f.median_k, f.max_k));
    return out;
}

py::dict tom_dict(const TomInfo &info) {
    py::dict d;
    d["power"] = info.power;
    d["estimated"] = info.estimated;
    d["fits"] = fits_list(info.fits);
    return d;
}

}  // namespace

void bind_correlation(py::module_ &m) {
    py::class_<EdgeFormatter>(m, "EdgeFormatter", R"doc(
Format edge rows as GCEN v1 does: ``gene_a<TAB>gene_b<TAB>r<TAB>p[<TAB>fdr][<TAB>tom]``
with ``r`` printed as ``%f`` and ``p`` / ``fdr`` / ``tom`` as ``%e`` (6 decimals).
)doc")
        .def(py::init<std::vector<std::string>>(), py::arg("names"))
        .def("format", &EdgeFormatter::format, py::arg("idx"), py::arg("cor"),
             py::arg("pval"), py::arg("fdr") = py::none(), py::arg("n_threads") = 1,
             py::arg("tom") = py::none(),
             "Return the formatted lines for gene-index pairs ``idx`` (n x 2) as bytes.");

    m.def(
        "pairwise_correlation",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> matrix,
           const std::string &method,
           int n_threads,
           double min_abs_cor,
           double max_p,
           bool positive_only,
           int row_begin,
           int row_end) -> std::tuple<py::array_t<int>, py::array_t<double>, py::array_t<double>>
        {
            auto buf = matrix.request();
            if (buf.ndim != 2)
                throw std::runtime_error("Input matrix must be 2-dimensional (genes x samples).");

            int n_genes = static_cast<int>(buf.shape[0]);
            int n_samples = static_cast<int>(buf.shape[1]);
            double *data = static_cast<double *>(buf.ptr);

            if (row_begin < 0 || (row_end >= 0 && row_end < row_begin) || row_end > n_genes)
                throw std::invalid_argument("invalid row range.");
            const bool nonzero = method == "pearson_nonzero";
            const bool complete = method == "pearson_pairwise_complete";
            if (method != "pearson" && method != "spearman" && !nonzero && !complete)
                throw std::invalid_argument(
                    "method must be 'pearson', 'spearman', 'pearson_nonzero' or "
                    "'pearson_pairwise_complete'.");
            bool use_spearman = (method == "spearman");
            std::vector<std::tuple<int, int, double, double>> result;
            {
                py::gil_scoped_release release;
                if (nonzero || complete) {
                    const NonzeroRows rows = prepare_nonzero_rows(
                        data, n_genes, n_samples,
                        nonzero ? NonzeroMode::BothNonzero : NonzeroMode::AnyNonzero);
                    for (const CorEdge &e : correlate_rows_nonzero(
                             rows, n_threads, min_abs_cor, max_p, positive_only, row_begin,
                             row_end))
                        result.emplace_back(e.i, e.j, e.corr, e.pval);
                } else {
                    result = compute_pairwise_correlation(
                        data, n_genes, n_samples, use_spearman, n_threads, min_abs_cor, max_p,
                        positive_only, row_begin, row_end);
                }
            }

            // result is a vector of (i, j, cor, pval) tuples
            size_t n_edges = result.size();
            py::array_t<int> indices({static_cast<py::ssize_t>(n_edges), static_cast<py::ssize_t>(2)});
            py::array_t<double> correlations(static_cast<py::ssize_t>(n_edges));
            py::array_t<double> p_values(static_cast<py::ssize_t>(n_edges));

            auto idx_buf = indices.mutable_unchecked<2>();
            auto cor_buf = correlations.mutable_unchecked<1>();
            auto pval_buf = p_values.mutable_unchecked<1>();

            for (size_t k = 0; k < n_edges; ++k) {
                idx_buf(k, 0) = std::get<0>(result[k]);
                idx_buf(k, 1) = std::get<1>(result[k]);
                cor_buf(k) = std::get<2>(result[k]);
                pval_buf(k) = std::get<3>(result[k]);
            }

            return {indices, correlations, p_values};
        },
        py::arg("matrix"),
        py::arg("method") = "spearman",
        py::arg("n_threads") = 2,
        py::arg("min_abs_cor") = 0.0,
        py::arg("max_p") = std::numeric_limits<double>::infinity(),
        py::arg("positive_only") = false,
        py::arg("row_begin") = 0,
        py::arg("row_end") = -1,
        R"doc(
Compute all pairwise Pearson or Spearman correlations.

Parameters
----------
matrix : np.ndarray, shape (n_genes, n_samples), dtype float64
    Gene expression matrix (genes as rows, samples as columns).
method : str
    'pearson', 'spearman', 'pearson_nonzero' (Pearson over the samples
    where both genes are non-zero) or 'pearson_pairwise_complete' (over the
    samples where at least one is non-zero); p-values on that number of
    samples (see src/correlation.hpp).
n_threads : int
    Number of parallel threads.
min_abs_cor : float, default 0.0
    Only return pairs with |r| >= min_abs_cor (0 disables this filter).
max_p : float, default +inf
    Only return pairs with p <= max_p (+inf disables this filter).
positive_only : bool, default False
    Only return pairs with r >= 0.
row_begin, row_end : int, default 0, -1
    Only compute pairs (i, j) with row_begin <= i < row_end and j > i
    (row_end = -1: up to the last gene). Used to stream the pairs in blocks.

Filtering is done inside the worker threads, so memory use scales with the
number of returned edges rather than with n_genes^2.

Returns
-------
indices : np.ndarray, shape (n_edges, 2), dtype int32
    Gene pair indices (i, j), i < j, sorted by (i, j).
correlations : np.ndarray, shape (n_edges,), dtype float64
    Correlation coefficients.
p_values : np.ndarray, shape (n_edges,), dtype float64
    Two-tailed p-values.
)doc"
    );

    m.def(
        "build_network_file",
        [](const std::string &in_path, const std::string &out_path, bool has_header,
           const std::string &method, double cor, double pval, bool positive_only,
           bool all_edges, int n_threads, int64_t max_pairs_per_block, bool fdr,
           py::object tom_power, double tom_min, std::vector<double> tom_powers,
           const std::string &kernel) -> py::object {
            use_gemm_kernel(kernel, 0);  // validate
            if (method != "pearson" && method != "spearman")
                throw std::invalid_argument("method must be 'pearson' or 'spearman'.");
            if (max_pairs_per_block < 1)
                throw std::invalid_argument("max_pairs_per_block must be positive.");
            const TomSpec tom = tom_spec(tom_power, tom_min, tom_powers);
            BuildFileResult res;
            {
                py::gil_scoped_release release;
                res = build_network_file(in_path, out_path, has_header, method == "spearman",
                                         cor, pval, positive_only, all_edges, n_threads,
                                         max_pairs_per_block, fdr, tom, kernel);
            }
            if (!res.supported) return py::none();
            return py::make_tuple(res.node_count, res.edge_count,
                                  tom.enabled ? tom_dict(res.tom) : py::object(py::none()));
        },
        py::arg("in_path"), py::arg("out_path"), py::arg("has_header"), py::arg("method"),
        py::arg("cor"), py::arg("pval"), py::arg("positive_only"), py::arg("all_edges"),
        py::arg("n_threads"), py::arg("max_pairs_per_block"), py::arg("fdr") = false,
        py::arg("tom_power") = py::none(), py::arg("tom_min") = 0.0,
        py::arg("tom_powers") = std::vector<double>(), py::arg("kernel") = "auto",
        R"doc(
Build a network from a plain expression matrix file and write the edge list,
without numpy: same file as gcen2._matrix_io.read_matrix followed by
gcen2.network.build_to_file (no log transform). Returns
(node_count, edge_count, tom_info or None), or None, writing nothing, when
the input is outside the supported plain subset or the output cannot be
opened (the caller then uses the Python path). TOM options as in
build_network_array.
)doc");

    m.def(
        "build_network_array",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> matrix,
           std::vector<std::string> genes, const std::string &out_path, const std::string &method,
           double cor, double pval, bool positive_only, bool fdr, bool all_edges, int n_threads,
           int64_t max_pairs_per_block, int64_t fdr_keep_limit, py::object tom_power,
           double tom_min, std::vector<double> tom_powers, const std::string &kernel) -> py::tuple {
            const auto buf = check_matrix(matrix);
            use_gemm_kernel(kernel, 0);  // validate before the output file is opened
            if (genes.size() != static_cast<size_t>(buf.shape[0]))
                throw std::invalid_argument("genes must give one name per row of matrix.");
            if (max_pairs_per_block < 1)
                throw std::invalid_argument("max_pairs_per_block must be positive.");
            const TomSpec tom = tom_spec(tom_power, tom_min, tom_powers);
            std::FILE *fh = std::fopen(out_path.c_str(), "wb");
            if (!fh) {
                PyErr_SetFromErrnoWithFilename(PyExc_OSError, out_path.c_str());
                throw py::error_already_set();
            }
            int64_t nodes, edges;
            TomInfo info;
            {
                py::gil_scoped_release release;
                WorkerPool pool(n_threads);
                const double *data = static_cast<const double *>(buf.ptr);
                const int ng = static_cast<int>(buf.shape[0]), ns = static_cast<int>(buf.shape[1]);
                Correlate correlate;
                std::unique_ptr<TomScorer> scorer;
                try {
                    correlate = make_correlate(method, data, ng, ns, pool, kernel);
                    scorer = make_tom_scorer(tom, method, data, ng, ns, positive_only, pool);
                } catch (...) {
                    std::fclose(fh);
                    throw;
                }
                EdgeFileWriter writer(fh, out_path, genes, fdr, pool, scorer != nullptr);
                NetworkOptions opt = options(cor, pval, positive_only, fdr, all_edges,
                                             max_pairs_per_block, fdr_keep_limit);
                opt.tom = scorer.get();
                opt.tom_min = tom.tom_min;
                stream_network(ng, opt, pool.size(), correlate,
                               [&](const std::vector<CorEdge> &e, const std::vector<double> *q,
                                   const std::vector<double> *t) { writer.write(e, q, t); });
                writer.close();
                nodes = writer.node_count();
                edges = writer.edge_count();
                if (scorer) info = tom_info(*scorer);
            }
            return py::make_tuple(nodes, edges,
                                  tom.enabled ? tom_dict(info) : py::object(py::none()));
        },
        py::arg("matrix"), py::arg("genes"), py::arg("out_path"), py::arg("method"),
        py::arg("cor"), py::arg("pval"), py::arg("positive_only"), py::arg("fdr"),
        py::arg("all_edges"), py::arg("n_threads"), py::arg("max_pairs_per_block"),
        py::arg("fdr_keep_limit") = int64_t(1) << 24, py::arg("tom_power") = py::none(),
        py::arg("tom_min") = 0.0, py::arg("tom_powers") = std::vector<double>(),
        py::arg("kernel") = "auto",
        R"doc(
gcen2.network.build_to_file in C++: correlate the genes x samples matrix in
blocks of at most max_pairs_per_block pairs, filter (and adjust with
Benjamini-Hochberg when fdr) and write the edge list to out_path. With fdr,
up to fdr_keep_limit edges are kept between the p-value pass and the output
(24 bytes each); beyond that the blocks are correlated a second time.

tom_power: None, no TOM; a power > 0, or 0 to pick it from tom_powers
(default WGCNA's 1..10, 12..20) as pickSoftThreshold does (src/tom.hpp). A
"tom" column is added and, unless all_edges, edges with TOM < tom_min are
dropped. Signed (positive_only) networks use ((1 + r) / 2)^power, others
|r|^power. Returns (node_count, edge_count, tom_info or None), tom_info a
dict: power, estimated, fits (power, r2, slope, mean_k, median_k, max_k).
)doc");

    m.attr("GEMM_MIN_SAMPLES") = kGemmMinSamples;

    m.def(
        "network_edges",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> matrix,
           const std::string &method, int n_threads, double cor, double pval,
           bool positive_only, bool fdr, bool all_edges, py::object tom_power, double tom_min,
           std::vector<double> tom_powers, const std::string &kernel) -> py::tuple {
            const auto buf = check_matrix(matrix);
            const TomSpec tom = tom_spec(tom_power, tom_min, tom_powers);
            std::vector<CorEdge> all;
            std::vector<double> all_q, all_t;
            TomInfo info;
            {
                py::gil_scoped_release release;
                WorkerPool pool(n_threads);
                const double *data = static_cast<const double *>(buf.ptr);
                const int ng = static_cast<int>(buf.shape[0]), ns = static_cast<int>(buf.shape[1]);
                const Correlate correlate = make_correlate(method, data, ng, ns, pool, kernel);
                std::unique_ptr<TomScorer> scorer =
                    make_tom_scorer(tom, method, data, ng, ns, positive_only, pool);
                // one block: the edges are all returned anyway
                NetworkOptions opt = options(cor, pval, positive_only, fdr, all_edges,
                                             std::numeric_limits<int64_t>::max(),
                                             std::numeric_limits<int64_t>::max());
                opt.tom = scorer.get();
                opt.tom_min = tom.tom_min;
                stream_network(ng, opt, pool.size(), correlate,
                               [&](const std::vector<CorEdge> &e, const std::vector<double> *q,
                                   const std::vector<double> *t) {
                                   all.insert(all.end(), e.begin(), e.end());
                                   if (q) all_q.insert(all_q.end(), q->begin(), q->end());
                                   if (t) all_t.insert(all_t.end(), t->begin(), t->end());
                               });
                if (scorer) info = tom_info(*scorer);
            }
            const auto n = static_cast<py::ssize_t>(all.size());
            py::array_t<int> idx({n, static_cast<py::ssize_t>(2)});
            py::array_t<double> r(n), p(n);
            int *ip = idx.mutable_data();
            double *rp = r.mutable_data(), *pp = p.mutable_data();
            for (py::ssize_t k = 0; k < n; ++k) {
                ip[2 * k] = all[k].i;
                ip[2 * k + 1] = all[k].j;
                rp[k] = all[k].corr;
                pp[k] = all[k].pval;
            }
            py::object q = py::none(), t = py::none(), ti = py::none();
            if (fdr) q = py::array_t<double>(n, all_q.data());
            if (tom.enabled) {
                t = py::array_t<double>(n, all_t.data());
                ti = tom_dict(info);
            }
            return py::make_tuple(idx, r, p, q, t, ti);
        },
        py::arg("matrix"), py::arg("method"), py::arg("n_threads"), py::arg("cor"),
        py::arg("pval"), py::arg("positive_only"), py::arg("fdr"), py::arg("all_edges"),
        py::arg("tom_power") = py::none(), py::arg("tom_min") = 0.0,
        py::arg("tom_powers") = std::vector<double>(), py::arg("kernel") = "auto",
        R"doc(
The edges of gcen2.network.build: (indices (n, 2) int32, correlations,
p-values, Benjamini-Hochberg values or None without fdr, TOM or None,
tom_info or None), filtered as build filters them (nothing filtered with
all_edges), sorted by (i, j). TOM options as in build_network_array.
)doc");

    m.def(
        "soft_threshold",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> matrix,
           const std::string &method, bool signed_net, std::vector<double> powers,
           int n_threads) -> py::list {
            const auto buf = check_matrix(matrix);
            if (method != "pearson" && method != "spearman" && method != "bicor")
                throw std::invalid_argument("method must be 'pearson', 'spearman' or 'bicor'.");
            if (powers.empty()) powers = default_soft_powers();
            std::vector<SoftThresholdFit> fits;
            {
                py::gil_scoped_release release;
                WorkerPool pool(n_threads);
                const int ng = static_cast<int>(buf.shape[0]), ns = static_cast<int>(buf.shape[1]);
                fits = soft_threshold_fits(tom_rows(method, static_cast<const double *>(buf.ptr),
                                                    ng, ns),
                                           ns, signed_net, powers, pool, method != "bicor");
            }
            return fits_list(fits);
        },
        py::arg("matrix"), py::arg("method"), py::arg("signed"),
        py::arg("powers") = std::vector<double>(), py::arg("n_threads") = 2,
        R"doc(
WGCNA 1.74 pickSoftThreshold fit table: one (power, r2, slope, mean_k,
median_k, max_k) tuple per power (default 1..10, 12..20), r2 the plain R^2 of
scaleFreeFitIndex. Signed: ((1 + r) / 2)^power, else |r|^power. method
"bicor": WGCNA's biweight midcorrelation (corFnc = bicor).
)doc");
}
