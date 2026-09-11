#include <MillikanWhite.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/**
 * Python wrapper for the MillikanWhiteModelData class, which carries
 * the vibration-translation relaxation parameters of one vibrator.
 */
void py_export_MillikanWhite(py::module &m) {
  py::class_<Mutation::Transfer::MillikanWhiteModelData>(
      m, "MillikanWhiteModelData")
    .def_property_readonly("species_index",
			   &Mutation::Transfer::MillikanWhiteModelData::
			     speciesIndex,
			   "Returns the index of the vibrator in the "
			   "mixture.")

    .def_property_readonly("molecular_weight",
			   &Mutation::Transfer::MillikanWhiteModelData::
			     molecularWeight,
			   "Returns the molecular weight of the vibrator "
			   "in kg/mol.")

    .def_property_readonly("reference_cross_section",
			   &Mutation::Transfer::MillikanWhiteModelData::
			     referenceCrossSection,
			   "Returns Park's reference cross section in m^2.")

    .def_property_readonly("a",
			   [](const Mutation::Transfer::
			      MillikanWhiteModelData &self) {
			     return std::vector<double>(
			       self.a().data(),
			       self.a().data() + self.a().size());
			   },
			   "Returns the Millikan-White a constant for each "
			   "heavy collision partner.")

    .def_property_readonly("b",
			   [](const Mutation::Transfer::
			      MillikanWhiteModelData &self) {
			     return std::vector<double>(
			       self.b().data(),
			       self.b().data() + self.b().size());
			   },
			   "Returns the Millikan-White b constant for each "
			   "heavy collision partner.");
}
