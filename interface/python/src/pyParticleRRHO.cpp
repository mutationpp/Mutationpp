#include <ParticleRRHO.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>

#include <utility>
#include <vector>

namespace py = nanobind;

/**
 * Python wrapper for the ParticleRRHO class, which carries the
 * rigid-rotor harmonic-oscillator parameters of a single species.
 */
void py_export_ParticleRRHO(py::module_ &m) {
  using Mutation::Thermodynamics::ParticleRRHO;

  py::class_<ParticleRRHO>(m, "ParticleRRHO")
    .def_prop_ro("formation_enthalpy", &ParticleRRHO::formationEnthalpy,
                 "Returns the formation enthalpy in J/mol.")

    .def_prop_ro("rotational_temperature",
                 &ParticleRRHO::rotationalTemperature,
                 "Returns the characteristic rotational temperature in K.")

    .def_prop_ro("steric_factor", &ParticleRRHO::stericFactor,
                 "Returns the steric factor.")

    .def_prop_ro("linearity", &ParticleRRHO::linearity,
                 "Returns the number of rotational degrees of freedom: 2 if "
                 "linear, 3 if nonlinear, 0 for an atom.")

    .def_prop_ro("electronic_levels",
                 [](const ParticleRRHO &self) {
                   std::vector<std::pair<int, double> > levels;
                   for (int i = 0; i < self.nElectronicLevels(); ++i)
                     levels.push_back(self.electronicEnergy(i));
                   return levels;
                 },
                 "Returns the (degeneracy, characteristic temperature in K) "
                 "pairs of the electronic levels.")

    .def_prop_ro("vibrational_temperatures",
                 [](const ParticleRRHO &self) {
                   std::vector<double> temperatures;
                   for (int i = 0; i < self.nVibrationalLevels(); ++i)
                     temperatures.push_back(self.vibrationalEnergy(i));
                   return temperatures;
                 },
                 "Returns the characteristic vibrational temperatures in K.");
}
