#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "binding_util.hpp"
#include "edge_files.hpp"
#include "edge_shuffle.hpp"
#include "network_stats.hpp"
#include "pcg64.hpp"
#include "edge_table.hpp"
#include "table_writer.hpp"
#include "layer_file.hpp"
#include "layer_fuse.hpp"

namespace py = pybind11;

using binding_util::vector_to_numpy;

namespace {

// numpy bit_generator.state dict <-> Pcg64
uint64_t py_u64(py::handle v, int shift) {
    return py::reinterpret_borrow<py::int_>(v).attr("__rshift__")(shift)
        .attr("__and__")(py::int_(UINT64_MAX)).cast<uint64_t>();
}

Pcg64 pcg64_from_state(py::dict st) {
    if (st["bit_generator"].cast<std::string>() != "PCG64")
        throw std::invalid_argument("not a PCG64 state.");
    py::dict inner = st["state"];
    return Pcg64(py_u64(inner["state"], 64), py_u64(inner["state"], 0),
                 py_u64(inner["inc"], 64), py_u64(inner["inc"], 0),
                 st["has_uint32"].cast<int>() != 0, st["uinteger"].cast<uint32_t>());
}

py::dict pcg64_state(const Pcg64 &g) {
    py::dict st, inner;
    auto u128 = [](uint64_t hi, uint64_t lo) {
        return py::int_(hi).attr("__lshift__")(64).attr("__or__")(py::int_(lo));
    };
    inner["state"] = u128(g.state_hi(), g.state_lo());
    inner["inc"] = u128(g.inc_hi(), g.inc_lo());
    st["bit_generator"] = "PCG64";
    st["state"] = inner;
    st["has_uint32"] = g.has_uint32() ? 1 : 0;
    st["uinteger"] = g.uinteger();
    return st;
}

// Read-only 1-D view of `v`, keeping `owner` (the Python object holding v) alive.
template <typename T>
py::array_t<T> readonly_view(const std::vector<T> &v, py::handle owner) {
    py::array_t<T> a(static_cast<py::ssize_t>(v.size()), v.data(), owner);
    a.attr("setflags")(py::arg("write") = false);
    return a;
}

}  // namespace

void bind_network(py::module_ &m) {
    py::class_<EdgeSwapper>(m, "EdgeSwapper", R"doc(
Degree-preserving edge swaps (gcen2.network.shuffle): the swap candidates
are passed in, so the caller's random stream decides the result.
)doc")
        .def(py::init([](py::array_t<int64_t, py::array::c_style | py::array::forcecast> a,
                         py::array_t<int64_t, py::array::c_style | py::array::forcecast> b,
                         int64_t base) {
                 if (a.ndim() != 1 || b.ndim() != 1)
                     throw std::invalid_argument("a and b must be 1-D.");
                 return new EdgeSwapper(std::vector<int64_t>(a.data(), a.data() + a.size()),
                                        std::vector<int64_t>(b.data(), b.data() + b.size()),
                                        base);
             }),
             py::arg("a"), py::arg("b"), py::arg("base"))
        .def("run",
             [](EdgeSwapper &self,
                py::array_t<int64_t, py::array::c_style | py::array::forcecast> i,
                py::array_t<int64_t, py::array::c_style | py::array::forcecast> j) {
                 if (i.ndim() != 1 || j.ndim() != 1 || i.size() != j.size())
                     throw std::invalid_argument("i and j must be 1-D of equal length.");
                 py::gil_scoped_release release;
                 return self.run(i.data(), j.data(), i.size());
             },
             py::arg("i"), py::arg("j"),
             "Try the swaps (i[t], j[t]) in order; returns how many were applied.")
        .def("shuffle",
             [](EdgeSwapper &self, py::dict state, int64_t n_swaps) {
                 Pcg64 rng = pcg64_from_state(state);
                 int64_t applied;
                 {
                     py::gil_scoped_release release;
                     applied = self.shuffle(rng, n_swaps);
                 }
                 return py::make_tuple(applied, pcg64_state(rng));
             },
             py::arg("state"), py::arg("n_swaps"),
             R"doc(
``n_swaps`` swaps drawn from numpy's PCG64 stream reproduced in C++,
started from ``state`` (``bit_generator.state``): the same swaps as calling
``run(rng.integers(0, n, m), rng.integers(0, n, m))`` per batch of up to
2**20. Returns (swaps applied, final state).
)doc")
        .def("b", [](const EdgeSwapper &self) { return vector_to_numpy(self.b()); },
             "Current second endpoint of every edge (the first never changes).");

    py::class_<Pcg64>(m, "Pcg64", R"doc(
numpy's PCG64 bit generator, started from ``bit_generator.state`` (see
gcen2.network._pcg64); ``integers(high, n)`` draws as
``Generator.integers(0, high, n)``.
)doc")
        .def(py::init([](py::dict st) {
                 return new Pcg64(pcg64_from_state(st));
             }),
             py::arg("state"))
        .def("integers", [](Pcg64 &g, int64_t high, int64_t n) {
            if (n < 0) throw std::invalid_argument("n must be >= 0.");
            std::vector<int64_t> v(static_cast<size_t>(n));
            g.integers(high, v.data(), n);
            return vector_to_numpy(std::move(v));
        }, py::arg("high"), py::arg("n"))
        .def_property_readonly("state", [](const Pcg64 &g) { return pcg64_state(g); },
                               "The state, as numpy's bit_generator.state.");

    m.def(
        "graph_stats",
        [](py::array_t<int64_t, py::array::c_style | py::array::forcecast> indptr,
           py::array_t<int32_t, py::array::c_style | py::array::forcecast> indices,
           int n_threads) -> py::dict {
            if (indptr.ndim() != 1 || indices.ndim() != 1 || indptr.shape(0) < 1)
                throw std::invalid_argument("indptr (n + 1) and indices must be 1-D.");
            const py::ssize_t n = indptr.shape(0) - 1;
            if (n > INT32_MAX) throw std::invalid_argument("too many nodes.");
            const int64_t *ip = indptr.data();
            const int32_t *ix = indices.data();
            if (ip[0] != 0 || ip[n] != indices.shape(0))
                throw std::invalid_argument("indptr does not match indices.");
            for (py::ssize_t i = 0; i < n; ++i) {
                if (ip[i + 1] < ip[i]) throw std::invalid_argument("indptr must be non-decreasing.");
                for (int64_t k = ip[i]; k < ip[i + 1]; ++k) {
                    if (ix[k] < 0 || ix[k] >= n) throw std::out_of_range("index out of range.");
                    if (k > ip[i] && ix[k] <= ix[k - 1])
                        throw std::invalid_argument("rows must be sorted and unique.");
                }
            }
            GraphStatsArrays s;
            {
                py::gil_scoped_release release;
                s = graph_stats(ip, ix, static_cast<int32_t>(n), n_threads);
            }
            py::dict d;
            d["degree_no_loop"] = vector_to_numpy(std::move(s.degree_no_loop));
            d["self_loop"] = vector_to_numpy(std::move(s.self_loop));
            d["triangles"] = vector_to_numpy(std::move(s.triangles));
            d["edge_count"] = s.edge_count;
            d["components"] = s.components;
            d["giant_component"] = s.giant_component;
            return d;
        },
        py::arg("indptr"), py::arg("indices"), py::arg("n_threads") = 2,
        R"doc(
Per-node quantities for NetworkStats from a symmetric CSR adjacency (rows
sorted and unique, self-loop = diagonal entry; the matrix must be
symmetric): degree_no_loop, self_loop, triangles (int64 arrays) and
edge_count, components, giant_component.
)doc");
    py::class_<DiffTable>(m, "DiffTable", R"doc(
Joined edges of two edge-list files for gcen2.differential (see
diff_edge_files); cor_a / cor_b / flags are read-only views.
)doc")
        .def_property_readonly("size", [](const DiffTable &t) { return t.flags.size(); })
        .def("cor_a", [](py::object self) {
            return readonly_view(self.cast<const DiffTable &>().cor_a, self);
        }, "Correlation in A per edge (NaN when not in A or NA).")
        .def("cor_b", [](py::object self) {
            return readonly_view(self.cast<const DiffTable &>().cor_b, self);
        }, "Correlation in B per edge (NaN when not in B or NA).")
        .def("flags", [](py::object self) {
            return readonly_view(self.cast<const DiffTable &>().flags, self);
        }, "Per edge: bit 0 set when in A, bit 1 when in B.")
        .def("lo", [](py::object self) {
            return readonly_view(self.cast<const DiffTable &>().lo, self);
        }, "Per edge: index into names() of the first gene (lo <= hi).")
        .def("hi", [](py::object self) {
            return readonly_view(self.cast<const DiffTable &>().hi, self);
        }, "Per edge: index into names() of the second gene.")
        .def("names", [](const DiffTable &t) { return t.names; },
             "Gene names, sorted (list of str).")
        .def("write",
             [](const DiffTable &t, const std::string &path,
                py::array_t<double, py::array::c_style | py::array::forcecast> z,
                py::array_t<double, py::array::c_style | py::array::forcecast> p,
                py::array_t<double, py::array::c_style | py::array::forcecast> q, double pval,
                bool use_fdr, int n_threads) {
                 if (z.ndim() != 1 || p.ndim() != 1 || q.ndim() != 1 || z.size() != p.size() ||
                     z.size() != q.size())
                     throw std::invalid_argument("z, p and q must be 1-D of equal length.");
                 py::gil_scoped_release release;
                 return t.write(path, z.data(), p.data(), q.data(), static_cast<size_t>(z.size()),
                                pval, use_fdr, n_threads);
             },
             py::arg("path"), py::arg("z"), py::arg("p"), py::arg("q"), py::arg("pval"),
             py::arg("use_fdr") = false, py::arg("n_threads") = 1,
             R"doc(
Write the result table (gene_a, gene_b, cor_a, cor_b, z_score, p_value,
fdr, category) to ``path``; ``z`` / ``p`` / ``q`` are the Fisher z test
and its FDR for the shared edges (flags == 3) in row order. Returns the
counts of A_only, B_only, shared_changed (p, or q with use_fdr, <= pval)
and shared_stable. Rows are formatted by n_threads threads.
)doc");

    m.def(
        "bh_adjust",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> p, int n_threads) {
            if (p.ndim() != 1) throw std::invalid_argument("p must be 1-D.");
            std::vector<double> q;
            {
                py::gil_scoped_release release;
                q = bh_adjust(p.data(), static_cast<size_t>(p.size()), n_threads);
            }
            return vector_to_numpy(std::move(q));
        },
        py::arg("p"), py::arg("n_threads") = 1,
        R"doc(
Benjamini-Hochberg adjusted p-values of the finite entries of p (NaN
elsewhere), bit for bit as gcen2.differential._bh; sorted by n_threads
threads.
)doc");

    m.def(
        "diff_edge_files",
        [](const std::tuple<std::string, bool, int> &a, const std::tuple<std::string, bool, int> &b,
           int n_threads, size_t block_bytes) {
            const EdgeFileSpec fa{std::get<0>(a), std::get<1>(a), std::get<2>(a)};
            const EdgeFileSpec fb{std::get<0>(b), std::get<1>(b), std::get<2>(b)};
            py::gil_scoped_release release;
            return diff_join_edge_files(fa, fb, n_threads, block_bytes);
        },
        py::arg("a"), py::arg("b"), py::arg("n_threads") = 1,
        py::arg("block_bytes") = size_t(1) << 26,
        R"doc(
gcen2.differential.compare on two edge-list files, each a (path,
skip_header, correlation column) tuple: their undirected edges (last row
wins), joined and sorted by (gene_a, gene_b), as a DiffTable. The files
are read in blocks of about block_bytes (smaller only in tests), parsed and
the edges sorted by n_threads threads.
)doc");

    auto table_spec = [](const std::string &path, bool skip_header, int n_columns,
                         bool v1_layout, const std::string &header, const std::string &output,
                         int n_threads) {
        EdgeTableSpec spec;
        spec.path = path;
        spec.skip_header = skip_header;
        spec.n_columns = n_columns;
        spec.v1_layout = v1_layout;
        spec.header = header;
        spec.output = output;
        spec.n_threads = n_threads;
        return spec;
    };
    auto result_dict = [](const EdgeTableResult &r) {
        py::dict d;
        d["ok"] = r.ok;
        d["reason"] = r.reason;
        d["node_count"] = r.node_count;
        d["edge_count"] = r.edge_count;
        d["swaps_applied"] = r.swaps_applied;
        return d;
    };

    m.def(
        "shuffle_edge_file",
        [=](const std::string &path, bool skip_header, int n_columns, bool v1_layout,
            const std::string &header, const std::string &output, int64_t swap_factor,
            py::dict state, int n_threads) {
            const EdgeTableSpec spec =
                table_spec(path, skip_header, n_columns, v1_layout, header, output, n_threads);
            Pcg64 rng = pcg64_from_state(state);
            EdgeTableResult r;
            {
                py::gil_scoped_release release;
                r = shuffle_edge_file(spec, swap_factor, rng);
            }
            return result_dict(r);
        },
        py::arg("path"), py::arg("skip_header"), py::arg("n_columns"), py::arg("v1_layout"),
        py::arg("header"), py::arg("output"), py::arg("swap_factor"), py::arg("state"),
        py::arg("n_threads") = 4,
        R"doc(
gcen2.network.shuffle on a file: read the edge list (``n_columns`` columns,
gene_a / gene_b first), swap with numpy's PCG64 stream started from
``state``, and write ``header`` and the rows as Network.to_edgelist would
(``v1_layout``: gene_a gene_b correlation p_value [fdr] as %f / %e; else as
DataFrame.to_csv of the dtypes pandas would infer). Returns a dict with
``ok`` (False, with ``reason``, when pandas' reading of the file is not
certain; nothing is written then), ``node_count``, ``edge_count`` and
``swaps_applied``.
)doc");

    m.def(
        "extract_edge_file",
        [=](const std::string &path, bool skip_header, int n_columns, bool v1_layout,
            const std::string &header, const std::string &output,
            const std::vector<std::string> &genes, int n_threads) {
            const EdgeTableSpec spec =
                table_spec(path, skip_header, n_columns, v1_layout, header, output, n_threads);
            EdgeTableResult r;
            {
                py::gil_scoped_release release;
                r = extract_edge_file(spec, genes);
            }
            return result_dict(r);
        },
        py::arg("path"), py::arg("skip_header"), py::arg("n_columns"), py::arg("v1_layout"),
        py::arg("header"), py::arg("output"), py::arg("genes"), py::arg("n_threads") = 4,
        R"doc(
gcen2.network.extract on a file: the rows whose first or second gene is in
``genes``, written as Network.to_edgelist would write them (see
shuffle_edge_file; only the kept rows decide the column dtypes and whether
the file is handled). Returns a dict with ``ok``, ``reason``,
``node_count`` and ``edge_count``.
)doc");

    m.def(
        "merge_edge_files",
        [](const std::vector<std::tuple<std::string, bool, int>> &files, double cor_cutoff) {
            std::vector<EdgeFileSpec> specs;
            for (const auto &f : files)
                specs.push_back({std::get<0>(f), std::get<1>(f), std::get<2>(f)});
            MergedEdges r;
            {
                py::gil_scoped_release release;
                r = merge_edge_files(specs, cor_cutoff);
            }
            return py::make_tuple(py::bytes(r.text), r.node_count, r.edge_count);
        },
        py::arg("files"), py::arg("cor_cutoff"),
        R"doc(
gcen2.network.merge on edge-list files. ``files`` holds (path, skip_header,
correlation column) tuples; returns (rows, node_count, edge_count) with
``rows`` the tab-separated gene_a, gene_b, correlation, direction lines
(no header) as bytes.
)doc");
    m.def(
        "write_coded_table",
        [](const std::string &path, const std::string &header,
           const std::vector<std::string> &names,
           py::array_t<int64_t, py::array::c_style | py::array::forcecast> a,
           py::array_t<int64_t, py::array::c_style | py::array::forcecast> b,
           const std::vector<py::tuple> &columns, int n_threads) {
            const size_t n = static_cast<size_t>(a.size());
            if (a.ndim() != 1 || b.ndim() != 1 || static_cast<size_t>(b.size()) != n)
                throw std::invalid_argument("a and b must be 1-D of equal length.");
            std::vector<TableColumn> cols(columns.size());
            std::vector<py::array_t<double, py::array::c_style | py::array::forcecast>> fs;
            std::vector<py::array_t<int64_t, py::array::c_style | py::array::forcecast>> is;
            std::vector<std::vector<std::string>> labels(columns.size());
            fs.reserve(columns.size());
            is.reserve(columns.size());
            for (size_t c = 0; c < columns.size(); ++c) {
                const py::tuple &t = columns[c];
                if (t.size() < 2)
                    throw std::invalid_argument("a column is (kind, values[, labels]).");
                const std::string kind = t[0].cast<std::string>();
                if (kind == "float") {
                    fs.push_back(py::array_t<double, py::array::c_style |
                                                     py::array::forcecast>::ensure(t[1]));
                    if (!fs.back() || fs.back().ndim() != 1 ||
                        static_cast<size_t>(fs.back().size()) != n)
                        throw std::invalid_argument("float column of the wrong shape.");
                    cols[c].kind = TableColumn::Float;
                    cols[c].f = fs.back().data();
                } else if (kind == "int" || kind == "label") {
                    is.push_back(py::array_t<int64_t, py::array::c_style |
                                                      py::array::forcecast>::ensure(t[1]));
                    if (!is.back() || is.back().ndim() != 1 ||
                        static_cast<size_t>(is.back().size()) != n)
                        throw std::invalid_argument(kind + " column of the wrong shape.");
                    cols[c].i = is.back().data();
                    if (kind == "int") {
                        cols[c].kind = TableColumn::Int;
                    } else {
                        if (t.size() < 3)
                            throw std::invalid_argument("a label column needs its labels.");
                        labels[c] = t[2].cast<std::vector<std::string>>();
                        cols[c].kind = TableColumn::Label;
                        cols[c].labels = &labels[c];
                    }
                } else {
                    throw std::invalid_argument("column kind must be float, int or label.");
                }
            }
            py::gil_scoped_release release;
            write_coded_table(path, header, names, a.data(), b.data(), n, cols, n_threads);
        },
        py::arg("path"), py::arg("header"), py::arg("names"), py::arg("a"), py::arg("b"),
        py::arg("columns"), py::arg("n_threads") = 4,
        R"doc(
Write a tab-separated table whose first two columns are names[a] and
names[b], followed by ``columns``: ("float", values), ("int", values) or
("label", codes, labels). Same bytes as pandas ``to_csv(sep="\t",
index=False)`` of the decoded table; ``header`` is written as given.
)doc");
    m.def(
        "read_layer_file",
        [](const std::string &path, const std::string &sep, bool skip_header, int n_cols,
           int weight_col) {
            py::dict out;
            if (sep.size() != 1) {
                out["ok"] = false;
                out["reason"] = "separator is not one character";
                return out;
            }
            LayerFileSpec spec{path, sep[0], skip_header, n_cols, weight_col};
            LayerFileResult r;
            {
                py::gil_scoped_release release;
                r = read_layer_file(spec);
            }
            out["ok"] = r.ok;
            out["reason"] = r.reason;
            if (!r.ok) return out;
            py::list names(r.names.size());
            try {
                for (size_t i = 0; i < r.names.size(); ++i) names[i] = py::str(r.names[i]);
            } catch (py::error_already_set &) {
                out["ok"] = false;
                out["reason"] = "gene names are not UTF-8";
                return out;
            }
            out["names"] = names;
            out["a"] = vector_to_numpy(std::move(r.a));
            out["b"] = vector_to_numpy(std::move(r.b));
            out["w"] = vector_to_numpy(std::move(r.w));
            return out;
        },
        py::arg("path"), py::arg("sep"), py::arg("skip_header"), py::arg("n_cols"),
        py::arg("weight_col"),
        R"doc(
Read an edge-list file for gcen2.multilayer.read_layer: returns a dict with
``ok`` and ``reason`` and, when ok, ``names`` (genes in first-seen order),
``a`` / ``b`` (int32 codes of the first two columns) and ``w`` (|weight| of
column ``weight_col``, or 1 when it is -1). ``ok`` is false for inputs whose
pandas reading is not certain (see src/layer_file.hpp); the caller then
reads the file with pandas.
)doc");
    m.def(
        "fuse_layer_codes",
        [](const std::vector<py::tuple> &layers, int64_t n_genes, const std::string &normalize,
           const std::vector<double> &alpha, int min_layers, double min_score, int n_threads) {
            using I64 = py::array_t<int64_t, py::array::c_style | py::array::forcecast>;
            using F64 = py::array_t<double, py::array::c_style | py::array::forcecast>;
            std::vector<I64> keep_i;
            std::vector<F64> keep_f;
            keep_i.reserve(2 * layers.size());
            keep_f.reserve(layers.size());
            std::vector<FuseLayerCodes> codes(layers.size());
            for (size_t k = 0; k < layers.size(); ++k) {
                if (layers[k].size() != 3)
                    throw std::invalid_argument("a layer is (a, b, w).");
                keep_i.push_back(I64::ensure(layers[k][0]));
                keep_i.push_back(I64::ensure(layers[k][1]));
                keep_f.push_back(F64::ensure(layers[k][2]));
                const I64 &a = keep_i[2 * k], &b = keep_i[2 * k + 1];
                const F64 &w = keep_f[k];
                if (!a || !b || !w || a.ndim() != 1 || b.ndim() != 1 || w.ndim() != 1 ||
                    a.size() != b.size() || a.size() != w.size())
                    throw std::invalid_argument("a, b and w must be 1-D of equal length.");
                codes[k] = FuseLayerCodes{a.data(), b.data(), w.data(),
                                          static_cast<size_t>(a.size())};
            }
            FuseSpec spec;
            spec.n_genes = n_genes;
            if (normalize == "max") spec.normalize = FuseSpec::Max;
            else if (normalize == "rank") spec.normalize = FuseSpec::Rank;
            else if (normalize == "none") spec.normalize = FuseSpec::None;
            else throw std::invalid_argument("normalize must be max, rank or none.");
            spec.alpha = alpha;
            spec.min_layers = min_layers;
            spec.min_score = min_score;
            spec.n_threads = n_threads;
            FuseResult r;
            {
                py::gil_scoped_release release;
                r = fuse_layer_codes(codes, spec);
            }
            py::dict out;
            out["ok"] = r.ok;
            out["reason"] = r.reason;
            if (!r.ok) return out;
            out["a"] = vector_to_numpy(std::move(r.a));
            out["b"] = vector_to_numpy(std::move(r.b));
            out["score"] = vector_to_numpy(std::move(r.score));
            out["n_layers"] = vector_to_numpy(std::move(r.n_layers));
            out["label_idx"] = vector_to_numpy(std::move(r.label_idx));
            out["patterns"] = r.patterns;
            py::list wide;
            for (auto &col : r.wide) wide.append(vector_to_numpy(std::move(col)));
            out["wide"] = wide;
            out["node_count"] = r.node_count;
            return out;
        },
        py::arg("layers"), py::arg("n_genes"), py::arg("normalize"), py::arg("alpha"),
        py::arg("min_layers"), py::arg("min_score"), py::arg("n_threads") = 4,
        R"doc(
gcen2.multilayer.fuse on gene codes. ``layers`` holds per layer (a, b, w):
int64 codes into the n_genes genes (sorted by name) and the |weights|.
Returns a dict with ``ok`` and ``reason`` and, when ok, the kept rows in
output order: ``a`` < ``b`` (codes), ``score``, ``n_layers``, ``label_idx``
into ``patterns`` (layer-presence bit patterns, ascending), ``wide`` (one
weight array per layer) and ``node_count``. ``ok`` is false for NaN weights
or more than 62 layers (see src/layer_fuse.hpp).
)doc");
}
