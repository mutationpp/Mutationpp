#include <MillikanWhite.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <vector>

namespace py = nanobind;

/**
 * Python wrapper for the MillikanWhiteModelData class, which carries the
 * vibration-translation relaxation parameters of one vibrator.
 */
void py_export_MillikanWhite(py::module_ &m) {
  using Mutation::Transfer::MillikanWhiteModelData;

  py::class_<MillikanWhiteModelData>(m, "MillikanWhiteModelData")
    .def_prop_ro("species_index", &MillikanWhiteModelData::speciesIndex,
                 "Returns the index of the vibrator in the mixture.")

    .def_prop_ro("molecular_weight", &MillikanWhiteModelData::molecularWeight,
                 "Returns the molecular weight of the vibrator in kg/mol.")

    .def_prop_ro("reference_cross_section",
                 &MillikanWhiteModelData::referenceCrossSection,
                 "Returns Park's reference cross section in m^2.")

    .def_prop_ro("a",
                 [](const MillikanWhiteModelData &self) {
                   return std::vector<double>(
                       self.a().data(), self.a().data() + self.a().size());
                 },
                 "Returns the Millikan-White a constant for each heavy "
                 "collision partner.")

    .def_prop_ro("b",
                 [](const MillikanWhiteModelData &self) {
                   return std::vector<double>(
                       self.b().data(), self.b().data() + self.b().size());
                 },
                 "Returns the Millikan-White b constant for each heavy "
                 "collision partner.");
}
