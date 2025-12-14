#include <nanobind/nanobind.h>

#define PY_MODULE NB_MODULE
namespace py = nanobind;

/**
 * This is the module definition for the Mutation++ Python wrapper.
 * The objects interfaces are kept as close as possible to the C++ interface
 * form minimal maintenance effort.
 */

void py_export_GlobalOptions(py::module_ &);
void py_export_MixtureOptions(py::module_ &);
void py_export_Mixture(py::module_ &);

PY_MODULE(_mutationpp, m) {
    m.doc() = "Mutation++ Python bindings";
    py_export_GlobalOptions(m);
    py_export_MixtureOptions(m);
    py_export_Mixture(m);
}
