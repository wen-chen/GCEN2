#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>

#include "binding_util.hpp"
#include "module_identify.hpp"

namespace py = pybind11;

using binding_util::check_csr;
using binding_util::vector_to_numpy;

void bind_module_identify(py::module_ &m) {
    m.def(
        "jaccard_modules",
        [](py::array_t<int64_t, py::array::c_style | py::array::forcecast> indptr,
           py::array_t<int32_t, py::array::c_style | py::array::forcecast> indices,
           double similarity_cutoff,
           int n_threads) -> py::array_t<int32_t>
        {
            if (indptr.ndim() != 1 || indices.ndim() != 1)
                throw std::invalid_argument("indptr and indices must be 1-D.");
            if (indptr.shape(0) < 1)
                throw std::invalid_argument("indptr must have n + 1 entries.");
            const py::ssize_t n = indptr.shape(0) - 1;
            if (n > INT32_MAX)
                throw std::invalid_argument("too many nodes.");
            const int64_t *ip = indptr.data();
            const int32_t *ix = indices.data();
            check_csr(ip, ix, n, indices.shape(0));

            std::vector<int32_t> labels;
            {
                py::gil_scoped_release release;
                labels = jaccard_modules_csr(ip, ix, static_cast<int32_t>(n),
                                             similarity_cutoff, n_threads);
            }
            py::array_t<int32_t> out(static_cast<py::ssize_t>(labels.size()));
            std::copy(labels.begin(), labels.end(), out.mutable_data());
            return out;
        },
        py::arg("indptr"),
        py::arg("indices"),
        py::arg("similarity_cutoff") = 0.5,
        py::arg("n_threads") = 2,
        R"doc(
Jaccard-similarity module identification on a CSR adjacency matrix
(GCEN v1 algorithm, visiting only 1- and 2-hop node pairs).

Parameters
----------
indptr : np.ndarray[int64], shape (n + 1,)
indices : np.ndarray[int32]
    CSR structure of a symmetric adjacency matrix with unique entries per
    row (a diagonal entry marks a self-loop).
similarity_cutoff : float
    Pairs with Jaccard similarity > cutoff are linked.
n_threads : int
    Worker threads for the similarity phase.

Returns
-------
np.ndarray[int32], shape (n,)
    Module label per node, -1 for nodes without any similarity link.
)doc"
    );

    m.def(
        "read_edge_adjacency",
        [](const std::string &path, bool skip_header, bool with_counts,
           bool descending) -> py::tuple {
            EdgeAdjacency adj;
            {
                py::gil_scoped_release release;
                adj = read_edge_adjacency(path, skip_header, with_counts, descending);
            }
            // hand the vectors to numpy without copying them
            auto indptr = vector_to_numpy(std::move(adj.indptr));
            auto indices = vector_to_numpy(std::move(adj.indices));
            py::list names(adj.names.size());
            for (size_t i = 0; i < adj.names.size(); ++i)
                names[i] = py::str(adj.names[i]);
            adj.names = {};
            if (!with_counts) return py::make_tuple(indptr, indices, names);
            auto counts = vector_to_numpy(std::move(adj.counts));
            return py::make_tuple(indptr, indices, names, counts);
        },
        py::arg("path"),
        py::arg("skip_header"),
        py::arg("with_counts") = false,
        py::arg("descending") = false,
        R"doc(
Read an edge-list file (first two tab-separated columns) into a symmetric
0/1 CSR adjacency without pandas: duplicate and reversed edges collapse,
a self-loop is a diagonal entry, nodes are sorted by name.

Returns
-------
(indptr : np.ndarray[int64], indices : np.ndarray[int32], names : list[str])
plus ``counts : np.ndarray[int32]`` (entry multiplicities, a self-loop line
counts 2) when ``with_counts`` is true. ``descending`` stores each row's
columns in descending order (the order of a scipy ``W @ D`` product).
)doc"
    );
}
