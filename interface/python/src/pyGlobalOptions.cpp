#include <GlobalOptions.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

/**
 * Python wrapper definition for the GlobalOptions class.
 */

void py_export_GlobalOptions(py::module &m) {

  /**
   * Overloaded member functions wrappers
   */

  /**
   * Python class definition
   */
  py::class_<Mutation::GlobalOptions>(m, "GlobalOptions")
      .def_static("dataDirectory", [](const std::string& datadir) {
        Mutation::GlobalOptions::dataDirectory(datadir);
      })
      .def_static("dataDirectory", []() {
        return Mutation::GlobalOptions::dataDirectory();
      })
      .def_static("workingDirectory", [](const std::string& workdir) {
        Mutation::GlobalOptions::workingDirectory(workdir);
      })
      .def_static("workingDirectory", []() {
        return Mutation::GlobalOptions::workingDirectory();
      })
      .def_static("reset", []() {
        Mutation::GlobalOptions::reset();
      })
    ;
}
