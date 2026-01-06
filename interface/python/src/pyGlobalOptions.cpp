#include <GlobalOptions.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

namespace py = nanobind;

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

void py_export_GlobalOptions(py::module_ &m) {

  const py::object resources = py::module_::import_("importlib.resources");
  const py::object pkgname = m.attr("__spec__").attr("parent");
  const py::object rootdir = resources.attr("files")(pkgname);
  const py::object datadir = rootdir / py::str("data");
  const std::string data_directory = py::cast<std::string>(py::str(datadir));

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
