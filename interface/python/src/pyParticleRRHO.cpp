#include <Mixture.h>
#include <ParticleRRHO.h>
#include <Utilities.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/**
 * Python wrapper for the ParticleRRHO class, which carries the
 * rigid-rotor harmonic-oscillator parameters of a single species.
 */
void py_export_ParticleRRHO(py::module &m) {
  py::class_<Mutation::Thermodynamics::ParticleRRHO>(m, "ParticleRRHO")
    .def_property_readonly("formation_enthalpy",
			   &Mutation::Thermodynamics::ParticleRRHO::
			     formationEnthalpy,
			   "Returns the formation enthalpy in J/mol.")

    .def_property_readonly("rotational_temperature",
			   &Mutation::Thermodynamics::ParticleRRHO::
			     rotationalTemperature,
			   "Returns the characteristic rotational "
			   "temperature in K.")

    .def_property_readonly("steric_factor",
			   &Mutation::Thermodynamics::ParticleRRHO::
			     stericFactor,
			   "Returns the steric factor.")

    .def_property_readonly("linearity",
			   &Mutation::Thermodynamics::ParticleRRHO::linearity,
			   "Returns the number of rotational degrees of "
			   "freedom: 2 if linear, 3 if nonlinear, 0 for an "
			   "atom.")

    .def_property_readonly("electronic_levels",
			   [](const Mutation::Thermodynamics::ParticleRRHO
			      &self) {
			     std::vector<std::pair<int, double> > levels;
			     for (int i = 0; i < self.nElectronicLevels(); ++i)
			       levels.push_back(self.electronicEnergy(i));
			     return levels;
			   },
			   "Returns the (degeneracy, characteristic "
			   "temperature in K) pairs of the electronic "
			   "levels.")

    .def_property_readonly("vibrational_temperatures",
			   [](const Mutation::Thermodynamics::ParticleRRHO
			      &self) {
			     std::vector<double> temperatures;
			     for (int i = 0; i < self.nVibrationalLevels(); ++i)
			       temperatures.push_back(
				 self.vibrationalEnergy(i));
			     return temperatures;
			   },
			   "Returns the characteristic vibrational "
			   "temperatures in K.");
}
