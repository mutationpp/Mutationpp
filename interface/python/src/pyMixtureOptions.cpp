#include <MixtureOptions.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

/**
 * Python wrapper definition for the MixtureOptions class. All member
 * functions using base data types, strings or classes already exposed
 * to Python as arguments or return values are exposed. The
 * ones which do not are kept unexposed.
 */

void py_export_MixtureOptions(py::module &m) {
  /**
   * Overloaded member functions wrappers
   */

  /**
   * Python class definition
   */
  py::class_<Mutation::MixtureOptions>(m, "MixtureOptions")
      .def(py::init<std::string &>())
      .def(py::init<>())
      .def("loadFromFile", &Mutation::MixtureOptions::loadFromFile)
      .def("loadFromXmlElement", &Mutation::MixtureOptions::loadFromXmlElement)
      .def("getSource", &Mutation::MixtureOptions::getSource)
      .def("setDefaultOptions", &Mutation::MixtureOptions::setDefaultOptions)
      .def("getSpeciesDescriptor",
           &Mutation::MixtureOptions::getSpeciesDescriptor)
      .def("setSpeciesDescriptor",
           &Mutation::MixtureOptions::setSpeciesDescriptor)
      .def("getStateModel", &Mutation::MixtureOptions::getStateModel)
      .def("setStateModel", &Mutation::MixtureOptions::setStateModel)
      .def("getThermodynamicDatabase",
           &Mutation::MixtureOptions::getThermodynamicDatabase)
      .def("setThermodynamicDatabase",
           &Mutation::MixtureOptions::setThermodynamicDatabase)
      .def("getMechanism", &Mutation::MixtureOptions::getMechanism)
      .def("setMechanism", &Mutation::MixtureOptions::setMechanism)
      .def("getViscosityAlgorithm",
           &Mutation::MixtureOptions::getViscosityAlgorithm)
      .def("setViscosityAlgorithm",
           &Mutation::MixtureOptions::setViscosityAlgorithm)
      .def("getThermalConductivityAlgorithm",
           &Mutation::MixtureOptions::getThermalConductivityAlgorithm)
      .def("setThermalConductivityAlgorithm",
           &Mutation::MixtureOptions::setThermalConductivityAlgorithm)
      .def("setDefaultComposition",
           &Mutation::MixtureOptions::setDefaultComposition)
      .def("hasDefaultComposition",
           &Mutation::MixtureOptions::hasDefaultComposition);
}