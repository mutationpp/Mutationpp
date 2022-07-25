#include <pybind11/pybind11.h>

/**
 * This is the module definition for the Mutation++ Python wrapper.
 * The objects interfaces are kept as close as possible to the C++ interface
 * form minimal maintenance effort.
 */

namespace py = pybind11;

void py_export_MixtureOptions(py::module &);
void py_export_Mixture(py::module &);

PYBIND11_MODULE(_mutationpp, m) {
    m.doc() = "Mutation++ Python bindings";
    py_export_MixtureOptions(m);
    py_export_Mixture(m);
}
