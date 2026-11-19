#pragma once
// Helpers shared by the pybind11 binding files.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace binding_util {

namespace py = pybind11;

// A numpy array that takes ownership of a std::vector's buffer (no copy).
template <typename T>
py::array_t<T> vector_to_numpy(std::vector<T> &&v) {
    auto *owner = new std::vector<T>(std::move(v));
    py::capsule free_when_done(owner, [](void *p) { delete static_cast<std::vector<T> *>(p); });
    return py::array_t<T>(static_cast<py::ssize_t>(owner->size()), owner->data(),
                          free_when_done);
}

// Checks a CSR structure of n rows and nnz entries: indptr starts at 0, ends
// at nnz and never decreases; every column index is in [0, n).
inline void check_csr(const int64_t *ip, const int32_t *ix, py::ssize_t n, py::ssize_t nnz) {
    if (ip[0] != 0 || ip[n] != nnz)
        throw std::invalid_argument("indptr does not match indices.");
    for (py::ssize_t i = 0; i < n; ++i)
        if (ip[i + 1] < ip[i])
            throw std::invalid_argument("indptr must be non-decreasing.");
    for (py::ssize_t p = 0; p < nnz; ++p)
        if (ix[p] < 0 || ix[p] >= n)
            throw std::invalid_argument("column index out of range.");
}

}  // namespace binding_util
