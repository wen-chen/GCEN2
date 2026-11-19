#include <pybind11/pybind11.h>

namespace py = pybind11;

// Forward declarations for each submodule
void bind_correlation(py::module_ &m);
void bind_enrichment(py::module_ &m);
void bind_module_identify(py::module_ &m);
void bind_rwr(py::module_ &m);
void bind_network(py::module_ &m);
void bind_tom_modules(py::module_ &m);

PYBIND11_MODULE(_gcen_core, m) {
    m.doc() = "GCEN2 C++ computation engine (pybind11 bindings)";

    auto m_cor = m.def_submodule("correlation", "Pairwise correlation kernels");
    bind_correlation(m_cor);

    auto m_enr = m.def_submodule("enrichment", "Hypergeometric enrichment kernels");
    bind_enrichment(m_enr);

    auto m_mod = m.def_submodule("module", "Module identification (Jaccard + DFS)");
    bind_module_identify(m_mod);
    bind_tom_modules(m_mod);

    auto m_rwr = m.def_submodule("rwr", "Random Walk with Restart");
    bind_rwr(m_rwr);

    auto m_net = m.def_submodule("network", "Network statistics");
    bind_network(m_net);
}
