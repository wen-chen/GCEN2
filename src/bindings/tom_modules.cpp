#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "hclust.hpp"
#include "network_stream.hpp"
#include "tom.hpp"
#include "tom_gemm.hpp"

namespace py = pybind11;

namespace {

py::tuple dendro_tuple(const Dendrogram &d, int64_t n) {
    const auto m = static_cast<py::ssize_t>(n > 0 ? n - 1 : 0);
    py::array_t<int> merge({m, static_cast<py::ssize_t>(2)});
    py::array_t<double> height(m);
    std::copy(d.merge.begin(), d.merge.end(), merge.mutable_data());
    std::copy(d.height.begin(), d.height.end(), height.mutable_data());
    return py::make_tuple(merge, height);
}

template <class T>
py::tuple hclust_impl(py::array_t<T, py::array::c_style | py::array::forcecast> diss) {
    auto buf = diss.request();
    if (buf.ndim != 2 || buf.shape[0] != buf.shape[1])
        throw std::invalid_argument("diss must be a square matrix.");
    const int64_t n = buf.shape[0];
    Dendrogram d;
    {
        py::gil_scoped_release release;
        d = hclust_average(static_cast<const T *>(buf.ptr), n);
    }
    return dendro_tuple(d, n);
}

}  // namespace

void bind_tom_modules(py::module_ &m) {
    m.def("tom_kernel", &gemm_nt_kernel,
          "Kernel of the TOM products: 'avx' or 'baseline' (same results).");
    m.def("hclust_average", &hclust_impl<float>, py::arg("diss"));
    m.def("hclust_average", &hclust_impl<double>, py::arg("diss"),
          R"doc(
Average-linkage clustering of a square dissimilarity matrix, as R's
fastcluster::hclust(as.dist(diss), method = "average"): (merge (n-1, 2)
int32 in R's convention, height (n-1)). Only the upper triangle is read.
)doc");
    m.def(
        "hclust_average_inplace",
        [](py::array_t<double, py::array::c_style> diss) -> py::tuple {
            auto buf = diss.request(true);
            if (buf.ndim != 2 || buf.shape[0] != buf.shape[1])
                throw std::invalid_argument("diss must be a square matrix.");
            const int64_t n = buf.shape[0];
            Dendrogram d;
            {
                py::gil_scoped_release release;
                d = hclust_average_inplace(static_cast<double *>(buf.ptr), n);
            }
            return dendro_tuple(d, n);
        },
        py::arg("diss"),
        R"doc(
hclust_average on a writable C-ordered float64 square matrix without copying
it: the lower triangle is set to the upper one, the upper one is the working
storage and is restored at the end (the matrix ends up symmetric).
)doc");

    m.def(
        "tom_dissimilarity",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> matrix,
           const std::string &method, bool signed_net, double power,
           std::vector<double> powers, int n_threads, bool similarity) -> py::tuple {
            auto buf = matrix.request();
            if (buf.ndim != 2) throw std::invalid_argument("matrix must be 2-dimensional.");
            const int ng = static_cast<int>(buf.shape[0]), ns = static_cast<int>(buf.shape[1]);
            if (method != "pearson" && method != "spearman" && method != "bicor")
                throw std::invalid_argument("TOM needs method 'pearson', 'spearman' or 'bicor'.");
            if (powers.empty()) powers = default_soft_powers();
            py::array_t<double> out({static_cast<py::ssize_t>(ng), static_cast<py::ssize_t>(ng)});
            double *op = out.mutable_data();
            TomInfo info;
            {
                py::gil_scoped_release release;
                WorkerPool pool(n_threads);
                info = tom_dissimilarity(tom_rows(method, static_cast<const double *>(buf.ptr),
                                                  ng, ns),
                                         ns, signed_net, power, powers, pool, op,
                                         method != "bicor", similarity);
            }
            py::list fits;
            for (const auto &f : info.fits)
                fits.append(py::make_tuple(f.power, f.r2, f.slope, f.mean_k, f.median_k, f.max_k));
            py::dict d;
            d["power"] = info.power;
            d["estimated"] = info.estimated;
            d["fits"] = fits;
            return py::make_tuple(out, d);
        },
        py::arg("matrix"), py::arg("method"), py::arg("signed"), py::arg("power"),
        py::arg("powers") = std::vector<double>(), py::arg("n_threads") = 2,
        py::arg("similarity") = false,
        R"doc(
1 - TOM of all gene pairs (float64 n x n, diagonal 0; computed in double,
unlike the float "tom" edge column) and the TOM info dict (power,
estimated, fits); power 0: picked as pickSoftThreshold. similarity: the TOM
itself (diagonal 1) instead of 1 - TOM.
)doc");

    m.def(
        "row_quantile",
        [](std::vector<py::array_t<double, py::array::c_style | py::array::forcecast>> sets,
           double q, int n_threads) -> py::array_t<double> {
            if (sets.empty()) throw std::invalid_argument("no sets.");
            if (!(q >= 0 && q <= 1)) throw std::invalid_argument("q must be in [0, 1].");
            const int64_t n = sets[0].size();
            std::vector<const double *> ptrs;
            for (auto &a : sets) {
                if (a.size() != n) throw std::invalid_argument("sets differ in length.");
                ptrs.push_back(a.data());
            }
            py::array_t<double> out(n);
            double *op = out.mutable_data();
            {
                py::gil_scoped_release release;
                WorkerPool pool(n_threads);
                wgcna_row_quantile(ptrs, n, q, op, pool);
            }
            return out;
        },
        py::arg("sets"), py::arg("q"), py::arg("n_threads") = 2,
        R"doc(
WGCNA rowQuantileC over sets: element i of the result is the q quantile of
element i of each array (NaN left out), by WGCNA's pivot(), so with its
rounding (the consensus TOM of blockwiseConsensusModules).
)doc");
}
