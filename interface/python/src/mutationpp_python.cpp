#include <pybind11/pybind11.h>

/**
 * This is the module definition for the Mutation++ Python wrapper.
 * The objects interfaces are kept as close as possible to the C++ interface
 * form minimal maintenance effort.
 */

namespace py = pybind11;

void python_export_MixtureOptions(py::module &);
void python_export_Mixture(py::module &);

PYBIND11_MODULE(mutationpp, m) {
  m.doc() = "Mutation++ Python bindings";
  python_export_MixtureOptions(m);
  python_export_Mixture(m);
}
