#include <pybind11/pybind11.h>

/**
 * This is the module definition for the Mutation++ Python wrapper.
 * The objects interfaces are kept as close as possible to the C++ interface
 * form minimal maintenance effort.
 */

namespace py = pybind11;

void init_MixtureOptions(py::module &);
//void export_Mixture();

PYBIND11_MODULE(mutationpp_python, m)
{
	m.doc() = "Mutation++ Python bindings";
	init_MixtureOptions(m);
}
