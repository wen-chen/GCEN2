#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdexcept>

#include "binding_util.hpp"
#include "rwr.hpp"

namespace py = pybind11;
using binding_util::check_csr;

void bind_rwr(py::module_ &m) {
    m.def(
        "random_walk_restart",
        [](py::array_t<int64_t, py::array::c_style | py::array::forcecast> indptr,
           py::array_t<int32_t, py::array::c_style | py::array::forcecast> indices,
           py::array_t<double, py::array::c_style | py::array::forcecast> data,
           py::array_t<double, py::array::c_style | py::array::forcecast> P0,
           double gamma,
           int max_iter,
           double tol) -> std::pair<py::array_t<double>, int>
        {
            if (indptr.ndim() != 1 || indices.ndim() != 1 || data.ndim() != 1 || P0.ndim() != 1)
                throw std::invalid_argument("all inputs must be 1-D arrays.");
            const py::ssize_t n = P0.shape(0);
            if (indptr.shape(0) != n + 1)
                throw std::invalid_argument("indptr must have len(P0) + 1 entries.");
            if (indices.shape(0) != data.shape(0))
                throw std::invalid_argument("indices and data must have equal length.");
            const int64_t *ip = indptr.data();
            const int32_t *ix = indices.data();
            check_csr(ip, ix, n, indices.shape(0));

            std::pair<std::vector<double>, int> result;
            {
                py::gil_scoped_release release;
                result = random_walk_restart_csr(ip, ix, data.data(), P0.data(),
                                                 static_cast<int32_t>(n), gamma,
                                                 max_iter, tol);
            }
            py::array_t<double> out(n);
            std::copy(result.first.begin(), result.first.end(), out.mutable_data());
            return {out, result.second};
        },
        py::arg("indptr"),
        py::arg("indices"),
        py::arg("data"),
        py::arg("P0"),
        py::arg("gamma")    = 0.5,
        py::arg("max_iter") = 1000000,
        py::arg("tol")      = 1e-10,
        R"doc(
Random Walk with Restart on a sparse (CSR) transition matrix.

PT+1 = (1 - gamma) * W * PT + gamma * P0, until ||PT+1 - PT||_1 < tol.

Parameters
----------
indptr, indices, data : CSR arrays (int64, int32, float64) of the
    column-normalised adjacency matrix W, shape (n, n).
P0 : np.ndarray, shape (n,), float64 — seed probability vector (sums to 1).
gamma : float — restart probability.
max_iter : int — maximum iterations.
tol : float — L1 convergence tolerance.

Returns
-------
(PT, steps) : steady-state probability vector and iteration count.
)doc"
    );

    m.def(
        "random_walk_restart_counts",
        [](py::array_t<int64_t, py::array::c_style | py::array::forcecast> indptr,
           py::array_t<int32_t, py::array::c_style | py::array::forcecast> indices,
           py::array_t<int32_t, py::array::c_style | py::array::forcecast> counts,
           py::array_t<double, py::array::c_style | py::array::forcecast> P0,
           double gamma,
           int max_iter,
           double tol) -> std::pair<py::array_t<double>, int>
        {
            if (indptr.ndim() != 1 || indices.ndim() != 1 || counts.ndim() != 1 || P0.ndim() != 1)
                throw std::invalid_argument("all inputs must be 1-D arrays.");
            const py::ssize_t n = P0.shape(0);
            if (indptr.shape(0) != n + 1)
                throw std::invalid_argument("indptr must have len(P0) + 1 entries.");
            if (indices.shape(0) != counts.shape(0))
                throw std::invalid_argument("indices and counts must have equal length.");
            const int64_t *ip = indptr.data();
            const int32_t *ix = indices.data();
            const int32_t *cp = counts.data();
            check_csr(ip, ix, n, indices.shape(0));
            for (py::ssize_t p = 0; p < counts.shape(0); ++p)
                if (cp[p] < 0) throw std::invalid_argument("counts must be non-negative.");
            std::pair<std::vector<double>, int> result;
            {
                py::gil_scoped_release release;
                result = random_walk_restart_counts(ip, ix, cp, P0.data(),
                                                    static_cast<int32_t>(n), gamma,
                                                    max_iter, tol);
            }
            py::array_t<double> out(n);
            std::copy(result.first.begin(), result.first.end(), out.mutable_data());
            return {out, result.second};
        },
        py::arg("indptr"), py::arg("indices"), py::arg("counts"), py::arg("P0"),
        py::arg("gamma") = 0.5, py::arg("max_iter") = 1000000, py::arg("tol") = 1e-10,
        "random_walk_restart on the unnormalised adjacency (entry multiplicities "
        "``counts``, int32): W = A / column sums is formed on the fly, giving the same "
        "result as passing ``(1.0 / col_sums)[indices] * counts`` as ``data``.");
}
