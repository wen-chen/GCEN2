#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>
#if __has_include(<charconv>)
#include <charconv>
#endif
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "enrichment.hpp"
#include "py_repr.hpp"

namespace py = pybind11;

void bind_enrichment(py::module_ &m) {
    m.def("format_py_float", [](double v) { std::string s; append_py_repr(s, v); return s; },
          py::arg("value"), "Python repr() of a float, as the annotation writer formats it.");

    using i64arr = py::array_t<long long, py::array::c_style | py::array::forcecast>;
    using f64arr = py::array_t<double, py::array::c_style | py::array::forcecast>;
    m.def(
        "format_annotation_rows",
        [](const std::vector<std::string> &term_prefix, i64arr terms, f64arr enrichment,
           i64arr study_count, i64arr study_n, i64arr pop_count, long long pop_n, f64arr p,
           f64arr fdr, i64arr bounds, int n_threads) {
            py::ssize_t n = terms.size();
            if (enrichment.size() != n || study_count.size() != n || study_n.size() != n ||
                pop_count.size() != n || p.size() != n || fdr.size() != n)
                throw std::invalid_argument("row arrays must have the same length.");
            py::ssize_t ng = bounds.size() - 1;
            const long long *bp = bounds.data();
            if (ng < 0 || bp[0] != 0 || bp[ng] != n)
                throw std::invalid_argument("bounds must start at 0 and end at the row count.");
            for (py::ssize_t g = 0; g < ng; ++g)
                if (bp[g + 1] < bp[g]) throw std::invalid_argument("bounds must be sorted.");
            const long long *tp = terms.data();
            for (py::ssize_t i = 0; i < n; ++i)
                if (tp[i] < 0 || tp[i] >= static_cast<long long>(term_prefix.size()))
                    throw std::out_of_range("term index out of range.");
            const double *ep = enrichment.data(), *pp = p.data(), *fp = fdr.data();
            const long long *sc = study_count.data(), *sn = study_n.data(), *pc = pop_count.data();
            std::vector<std::string> out(static_cast<size_t>(ng));
            {
                py::gil_scoped_release release;
                int nt = std::max(1, std::min<int>(n_threads, static_cast<int>(n / 20000) + 1));
                auto work = [&](int t) {
                    char nb[96];
                    for (py::ssize_t g = ng * t / nt; g < ng * (t + 1) / nt; ++g) {
                        std::string &s = out[static_cast<size_t>(g)];
                        for (long long i = bp[g]; i < bp[g + 1]; ++i) {
                            s += term_prefix[static_cast<size_t>(tp[i])];
                            s += '\t';
                            append_py_repr(s, ep[i]);
                            int len = std::snprintf(nb, sizeof nb, "\t%lld\t%lld\t%lld\t%lld\t",
                                                    sc[i], sn[i], pc[i], pop_n);
                            s.append(nb, static_cast<size_t>(len));
                            append_py_repr(s, pp[i]);
                            s += '\t';
                            append_py_repr(s, fp[i]);
                            s += '\n';
                        }
                    }
                };
                std::vector<std::thread> th;
                for (int t = 1; t < nt; ++t) th.emplace_back(work, t);
                work(0);
                for (auto &x : th) x.join();
            }
            py::list res(ng);
            for (py::ssize_t g = 0; g < ng; ++g)
                res[g] = py::bytes(out[static_cast<size_t>(g)]);
            return res;
        },
        py::arg("term_prefix"), py::arg("terms"), py::arg("enrichment"), py::arg("study_count"),
        py::arg("study_n"), py::arg("pop_count"), py::arg("pop_n"), py::arg("p_value"),
        py::arg("fdr"), py::arg("bounds"), py::arg("n_threads") = 1,
        "Format annotation result rows (one bytes object per query group "
        "rows[bounds[g]:bounds[g+1]]) exactly like pandas to_csv with sep='\\t': "
        "term_prefix[term] + floats as Python repr().");

    m.def(
        "hypergeometric_pvalue",
        [](int pop_n, int pop_count, int study_n, int study_count) -> double {
            return hypergeometric_pvalue(pop_n, pop_count, study_n, study_count);
        },
        py::arg("pop_n"),
        py::arg("pop_count"),
        py::arg("study_n"),
        py::arg("study_count"),
        "Compute a single hypergeometric p-value (Fisher exact test, upper tail)."
    );

    m.def(
        "batch_hypergeometric",
        [](int pop_n,
           py::array_t<int, py::array::c_style | py::array::forcecast> pop_counts,
           int study_n,
           py::array_t<int, py::array::c_style | py::array::forcecast> study_counts) -> py::array_t<double>
        {
            auto pc_buf = pop_counts.request();
            auto sc_buf = study_counts.request();
            if (pc_buf.ndim != 1 || sc_buf.ndim != 1)
                throw std::runtime_error("pop_counts and study_counts must be 1-D arrays.");
            if (pc_buf.size != sc_buf.size)
                throw std::runtime_error("pop_counts and study_counts must have the same length.");

            size_t n = static_cast<size_t>(pc_buf.size);
            const int *pc = static_cast<const int *>(pc_buf.ptr);
            const int *sc = static_cast<const int *>(sc_buf.ptr);

            py::array_t<double> result(static_cast<py::ssize_t>(n));
            double *out = result.mutable_data();

            {
                py::gil_scoped_release release;
                batch_hypergeometric(pop_n, pc, study_n, sc, out, n);
            }
            return result;
        },
        py::arg("pop_n"),
        py::arg("pop_counts"),
        py::arg("study_n"),
        py::arg("study_counts"),
        "Vectorized hypergeometric p-value computation over arrays."
    );

    using darr = py::array_t<double, py::array::c_style | py::array::forcecast>;
    m.def(
        "hypergeom_sf",
        [](darr k, darr M, darr n, darr N, int n_threads) {
            py::ssize_t len = k.size();
            if (M.size() != len || n.size() != len || N.size() != len)
                throw std::invalid_argument("k, M, n, N must have the same length.");
            for (py::ssize_t i = 0; i < len; ++i) {
                double kk = k.data()[i], mm = M.data()[i], nn = n.data()[i], NN = N.data()[i];
                if (!(mm >= 0 && nn >= 0 && NN >= 0 && nn <= mm && NN <= mm && mm < 2e9) ||
                    kk != std::floor(kk) || mm != std::floor(mm) || nn != std::floor(nn) ||
                    NN != std::floor(NN))
                    throw std::invalid_argument("hypergeom_sf needs integers 0 <= n, N <= M.");
            }
            py::array_t<double> out(len);
            const double *kp = k.data(), *mp = M.data(), *np_ = n.data(), *Np = N.data();
            double *op = out.mutable_data();
            {
                py::gil_scoped_release release;
                hypergeom_sf(kp, mp, np_, Np, op, static_cast<std::size_t>(len), n_threads);
            }
            return out;
        },
        py::arg("k"), py::arg("M"), py::arg("n"), py::arg("N"), py::arg("n_threads") = 1,
        "One-sided hypergeometric upper tail P(X >= k), X ~ Hypergeom(M, n, N), elementwise "
        "(1-D float arrays of equal length holding integers).");

    using s64arr = py::array_t<int64_t, py::array::c_style | py::array::forcecast>;
    using i32arr = py::array_t<int32_t, py::array::c_style | py::array::forcecast>;
    auto check_csr = [](s64arr &q_ptr, s64arr &q_idx, s64arr &g_ptr, i32arr &g_idx,
                        int32_t n_terms) {
        const py::ssize_t nq = q_ptr.size() - 1, ng = g_ptr.size() - 1;
        if (nq < 0 || ng < 0) throw std::invalid_argument("q_ptr and g_ptr need n + 1 entries.");
        const int64_t *qp = q_ptr.data(), *qi = q_idx.data(), *gp = g_ptr.data();
        const int32_t *gi = g_idx.data();
        if (qp[0] != 0 || qp[nq] != q_idx.size() || gp[0] != 0 || gp[ng] != g_idx.size())
            throw std::invalid_argument("pointer arrays do not match index arrays.");
        for (py::ssize_t r = 0; r < nq; ++r)
            if (qp[r + 1] < qp[r]) throw std::invalid_argument("q_ptr must be non-decreasing.");
        for (py::ssize_t g = 0; g < ng; ++g)
            if (gp[g + 1] < gp[g]) throw std::invalid_argument("g_ptr must be non-decreasing.");
        for (py::ssize_t p = 0; p < q_idx.size(); ++p)
            if (qi[p] < 0 || qi[p] >= ng) throw std::invalid_argument("gene index out of range.");
        for (py::ssize_t p = 0; p < g_idx.size(); ++p)
            if (gi[p] < 0 || gi[p] >= n_terms) throw std::invalid_argument("term index out of range.");
    };
    m.def(
        "study_counts",
        [check_csr](s64arr q_ptr, s64arr q_idx, s64arr g_ptr, i32arr g_idx, int32_t n_terms,
           int n_threads) {
            check_csr(q_ptr, q_idx, g_ptr, g_idx, n_terms);
            const py::ssize_t nq = q_ptr.size() - 1;
            const int64_t *qp = q_ptr.data(), *qi = q_idx.data(), *gp = g_ptr.data();
            const int32_t *gi = g_idx.data();
            std::vector<int64_t> rows, terms, counts;
            {
                py::gil_scoped_release release;
                study_counts(qp, qi, nq, gp, gi, n_terms, n_threads, rows, terms, counts);
            }
            auto to_np = [](std::vector<int64_t> &v) {
                py::array_t<int64_t> a(static_cast<py::ssize_t>(v.size()));
                std::copy(v.begin(), v.end(), a.mutable_data());
                v = {};
                return a;
            };
            auto r = to_np(rows);
            auto t = to_np(terms);
            auto c = to_np(counts);
            return py::make_tuple(r, t, c);
        },
        py::arg("q_ptr"), py::arg("q_idx"), py::arg("g_ptr"), py::arg("g_idx"),
        py::arg("n_terms"), py::arg("n_threads") = 1,
        "Study counts of many queries (CSR queries x genes, CSR genes x terms with unique "
        "entries per row): returns (rows, terms, counts) int64 arrays of the non-zero "
        "entries of Q @ G, grouped by query.");

    m.def(
        "significant_counts",
        [check_csr](s64arr q_ptr, s64arr q_idx, s64arr g_ptr, i32arr g_idx, int32_t n_terms,
                    s64arr study_n, s64arr pop_counts, int64_t pop_n, double pval,
                    int n_threads) {
            check_csr(q_ptr, q_idx, g_ptr, g_idx, n_terms);
            const py::ssize_t nq = q_ptr.size() - 1;
            if (study_n.size() != nq) throw std::invalid_argument("study_n needs one entry per query.");
            if (pop_counts.size() != n_terms)
                throw std::invalid_argument("pop_counts needs one entry per term.");
            if (!(pop_n >= 0 && pop_n < 2000000000))
                throw std::invalid_argument("pop_n out of range.");
            const int64_t *sn = study_n.data(), *pc = pop_counts.data();
            for (py::ssize_t r = 0; r < nq; ++r)
                if (sn[r] < 0 || sn[r] > pop_n)
                    throw std::invalid_argument("study_n must be in [0, pop_n].");
            for (int32_t t = 0; t < n_terms; ++t)
                if (pc[t] < 0 || pc[t] > pop_n)
                    throw std::invalid_argument("pop_counts must be in [0, pop_n].");
            std::vector<int64_t> rows, terms, counts;
            std::vector<double> p, fdr;
            std::vector<uint8_t> hits;
            {
                py::gil_scoped_release release;
                significant_counts(q_ptr.data(), q_idx.data(), nq, g_ptr.data(), g_idx.data(),
                                   n_terms, sn, pc, pop_n, pval, n_threads, rows, terms,
                                   counts, p, fdr, hits);
            }
            auto to_np = [](auto &v) {
                using T = typename std::decay_t<decltype(v)>::value_type;
                py::array_t<T> a(static_cast<py::ssize_t>(v.size()));
                std::copy(v.begin(), v.end(), a.mutable_data());
                std::decay_t<decltype(v)>().swap(v);
                return a;
            };
            auto r = to_np(rows);
            auto t = to_np(terms);
            auto c = to_np(counts);
            auto pv = to_np(p);
            auto q = to_np(fdr);
            py::array_t<bool> h(nq);
            std::copy(hits.begin(), hits.end(), h.mutable_data());
            return py::make_tuple(r, t, c, pv, q, h);
        },
        py::arg("q_ptr"), py::arg("q_idx"), py::arg("g_ptr"), py::arg("g_idx"),
        py::arg("n_terms"), py::arg("study_n"), py::arg("pop_counts"), py::arg("pop_n"),
        py::arg("pval"), py::arg("n_threads") = 1,
        "study_counts fused with the one-sided hypergeometric test (same p as hypergeom_sf): "
        "returns (rows, terms, counts, p_value, fdr, has_hits) for entries with "
        "pop_counts > 0 and p <= pval; fdr is the Benjamini-Hochberg adjusted p over the "
        "terms with pop_counts > 0 (all p-values of the query, untested terms p = 1); "
        "has_hits[r] tells whether query r shared any term with the background.");
}
