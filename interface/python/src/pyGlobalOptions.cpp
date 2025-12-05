#include <GlobalOptions.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

/**
 * Python wrapper definition for the GlobalOptions class.
 */

namespace {

  inline void setDefaultDataDirectory(const std::string& datadir) {
    if (std::getenv("MPP_DATA_DIRECTORY") == nullptr) {
      Mutation::GlobalOptions::dataDirectory(datadir);
    }
  };

}

void py_export_GlobalOptions(py::module &m) {

  const py::object resources = py::module_::import("importlib.resources");
  const py::object pkgname = m.attr("__spec__").attr("parent");
  const py::object rootdir = resources.attr("files")(pkgname);
  const py::object datadir = rootdir / py::str("data");
  const std::string data_directory = py::str(datadir).cast<std::string>();

  setDefaultDataDirectory(data_directory);

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
      .def_static("reset", [data_directory]() {
        Mutation::GlobalOptions::reset();
        setDefaultDataDirectory(data_directory);
      })
    ;
}
