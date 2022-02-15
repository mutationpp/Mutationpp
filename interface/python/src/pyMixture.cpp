#include <Mixture.h>
#include <MixtureOptions.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// Todo return the values with py::array_t<double> so it returns as numpy type

/**
 * Python wrapper definition for the Mixture class. All member
 * functions using base data types, strings or classes already exposed
 * to Python as arguments or return values are exposed. The
 * ones which do not are kept unexposed.
 */

void py_export_Mixture(py::module &m) {
  py::class_<Mutation::Mixture>(m, "Mixture")
      .def(py::init<Mutation::MixtureOptions>())
      .def(py::init<std::string &>())
      .def("nElements", &Mutation::Mixture::nElements,
           "Returns the number of elements considered in the mixture.")

      .def("nAtoms", &Mutation::Mixture::nAtoms,
           "Returns the number of atomic species in the mixture.")

      .def("nMolecules", &Mutation::Mixture::nMolecules,
           "Returns the number of molecules in the mixture.")

      .def("nHeavy", &Mutation::Mixture::nHeavy,
           "Returns the number of heavy particles (non electrons) in the "
           "mixture.")

      .def("nSpecies", &Mutation::Mixture::nSpecies,
           "Returns the number of species considered in the mixture.")

      .def("nPhases", &Mutation::Mixture::nPhases,
           "Returns the number of phases belonging to this mixture.")

      .def("nGas", &Mutation::Mixture::nGas,
           "Returns number of gas species in the mixture.")

      .def("nCondensed", &Mutation::Mixture::nCondensed,
           "Returns the number of condensed phase species in the mixture.")

      .def("nEnergyEqns", &Mutation::Mixture::nEnergyEqns,
           "Returns the number of energy equations associated with the mixture "
           "StateModel.")

      .def("nMassEqns", &Mutation::Mixture::nMassEqns,
           "Returns the number of mass equations associated with the mixture "
           "StateModel.")

      .def("hasElectrons", &Mutation::Mixture::hasElectrons,
           "Returns true if this mixture includes electrons, false otherwise.")

      .def("setState", &Mutation::Mixture::setState,
           "Sets the state of the mixture using the StateModel belonging to "
           "the mixture."
           "The input variables depend on the type of StateModel being used.")

      .def("setState",
           [](Mutation::Mixture &self,
                   std::vector<double> rho_i,
                   const double T,
                   const int vars)
          {
               self.setState(rho_i.data(),&T,vars);
           },
           "Sets the state of the mixture using the StateModel belonging to "
           "the mixture."
           "The input variables depend on the type of StateModel being used.")

     .def("setState",
           [](Mutation::Mixture &self,
                   std::vector<double> rho_i,
                   std::vector<double> T,
                   const int vars)
          {
               self.setState(rho_i.data(),T.data(),vars);
           },
           "Sets the state of the mixture using the StateModel belonging to "
           "the mixture."
           "The input variables depend on the type of StateModel being used.")

      .def("T", &Mutation::Mixture::T,
           "Returns the mixture translational temperature.")

      .def("Tr", &Mutation::Mixture::Tr,
           "Returns the mixture rotational temperature.")

      .def("Tv", &Mutation::Mixture::Tv,
           "Returns the mixture vibrational temperature.")

      .def("Te", &Mutation::Mixture::Te,
           "Returns the mixture electron temperature.")

      .def("Tel", &Mutation::Mixture::Tel,
           "Returns the mixture electronic temperature.")

      .def("P", &Mutation::Mixture::P,
           "Returns the current mixture static pressure in Pa.")

      .def("speciesName", &Mutation::Mixture::speciesName,
           "Returns the name of the species with the given index in the "
           "species array.")

      .def("elementName", &Mutation::Mixture::elementName,
           "Returns the name of the element with the given index in the "
           "element array.")

      .def("standardStateT", &Mutation::Mixture::standardStateT,
           "Returns the standard state temperature for the thermodynamic "
           "database begin used in K.")

      .def("standardStateP", &Mutation::Mixture::standardStateP,
           "Returns the standard state pressure for the thermodynamic database "
           "being used in Pa.")

      .def(
          "mixtureMwMass",
          [](const Mutation::Mixture &self, std::vector<double> y) {
            return self.mixtureMwMass(y.data());
          },
          "Returns the average molecular weight of the mixture given the "
          "species mass fractions.")

      .def(
          "mixtureMwMole",
          [](const Mutation::Mixture &self, std::vector<double> x) {
            return self.mixtureMwMole(x.data());
          },
          "Returns the average molecular weight of the mixture given the "
          "species mole fractions.")

      .def("mixtureMw", &Mutation::Mixture::mixtureMw,
           "Returns the average molecular weight of the mixture.")

      .def(
          "speciesMw",
          [](const Mutation::Mixture &self, int i) {
            return self.speciesMw(i);
          },
          "Returns the molecular weight in kg/mol of the species with the "
          "given index in the species array.")

      .def(
          "speciesMw",
          [](const Mutation::Mixture &self) { return self.speciesMw(); },
          "Returns an array of the species molecular weight in kg/mol.")

      .def("atomicMass", &Mutation::Mixture::atomicMass,
           "Returns the atomic weight in kg/mol of the element with the given "
           "index in the element array.")

      .def(
          "elementMatrix",
          [](const Mutation::Mixture &self) { return self.elementMatrix(); },
          "* Returns the element matrix \f$ E_{i,j} \f$ in which the (i,j) "
          "component is equal to the"
          "number of atoms of the jth element in one molecule of the ith "
          "species.")

      .def(
          "addComposition",
          [](Mutation::Mixture &self, std::string mixture,
             bool add_composition) {
            self.addComposition(mixture.c_str(), add_composition);
          },
          "Add a named element composition to the mixture which may be "
          "retrieved with getComposition().")

      .def(
          "getComposition_mass",
          [](Mutation::Mixture &self, std::string mixture_composition) {
            std::vector<double> x_e(self.nElements());
            self.getComposition(mixture_composition, x_e.data(),
                                Mutation::Thermodynamics::Composition::MASS);
            return py::array(py::cast(x_e));
          },
          "Gets the element mass fractions associated with a named composition "
          "in the mixture.")

      .def(
          "getComposition_mole",
          [](Mutation::Mixture &self, std::string mixture_composition) {
            std::vector<double> x_e(self.nElements());
            self.getComposition(mixture_composition, x_e.data(),
                                Mutation::Thermodynamics::Composition::MOLE);
            return py::array(py::cast(x_e));
          },
          "Gets the element mole fractions associated with a named composition "
          "in the mixture.")

      .def("numberDensity",
           static_cast<double (Mutation::Mixture::*)(const double, const double)
                           const>(&Mutation::Mixture::numberDensity),
           "Returns the number density of the mixture given the mixture "
           "temperature and pressure.")

      .def("numberDensity",
           static_cast<double (Mutation::Mixture::*)(void) const>(
               &Mutation::Mixture::numberDensity),
           "Returns the current number density of the mixture.")

      .def(
          "pressure",
          [](const Mutation::Mixture &self, const double T, const double P,
             const std::vector<double> Y) {
            return self.pressure(T, P, Y.data());
          },
          "Returns the pressure of the mixture in Pa given the mixture "
          "temperature and density and species mass fractions.")

      .def(
          "density",
          [](const Mutation::Mixture &self, const double T, const double P,
             const std::vector<double> X) {
            return self.density(T, P, X.data());
          },
          "Returns the density of the mixture in kg/m^3 given the mixture "
          "temperature and pressure and species mole fractions.")

      .def(
          "density",
          [](const Mutation::Mixture &self) { return self.density(); },
          "Returns the current mixture density in kg/m^3.")

      .def(
          "densities",
          [](const Mutation::Mixture &self) {
            std::vector<double> rho_i(self.nSpecies());
            self.densities(rho_i.data());
            return py::array(py::cast(rho_i));
          },
          "Returns the current species densities.")

      .def(
          "speciesCpOverR",
          [](const Mutation::Mixture &self) {
            std::vector<double> cp_i(self.nSpecies());
            self.speciesCpOverR(cp_i.data());
            return py::array(py::cast(cp_i));
          },
          "Returns the unitless vector of species specific heats at constant "
          "pressure \f$ C_{p,i} / R_u \f$.")

      .def(
          "speciesCpOverR",
          [](const Mutation::Mixture &self, const double T) {
            std::vector<double> cp_i(self.nSpecies());
            self.speciesCpOverR(T, cp_i.data());
            return py::array(py::cast(cp_i));
          },
          "Returns the unitless vector of species specific heats at constant "
          "pressure "
          "\f$ C_{p,i} / R_u \f$ given temperature in K.")

      .def("mixtureFrozenCpMole", &Mutation::Mixture::mixtureFrozenCpMole,
           "Returns the frozen mixture averaged specific heat at constant "
           "pressure in J/mol-K.")

      .def("mixtureFrozenCpMass", &Mutation::Mixture::mixtureFrozenCpMass,
           "Returns the frozen mixture averaged specific heat at constant "
           "pressure in J/kg-K.")

      .def("mixtureEquilibriumCpMole",
           &Mutation::Mixture::mixtureEquilibriumCpMole,
           "Returns the current equilibrium mixture averaged specific heat at"
           " constant pressure in J/mol-K. This method assumes that the "
           "current state"
           " of the mixture was already set using the equilibrate() method.")

      .def("mixtureEquilibriumCpMass",
           &Mutation::Mixture::mixtureEquilibriumCpMass,
           "Returns the current equilibrium mixture averaged specific heat at"
           " constant pressure in J/kg-K.  This method assumes that the "
           "current state"
           " of the mixture was already set using the equilibrate() method.")

      .def("mixtureFrozenCvMole", &Mutation::Mixture::mixtureFrozenCvMole,
           "Returns the mixture averaged specific heat at constant volume in "
           "J/mol-K.")

      .def("mixtureFrozenCvMass", &Mutation::Mixture::mixtureFrozenCvMass,
           "Returns the mixture averaged specific heat at constant volume in "
           "J/kg-K.")

      .def("mixtureEquilibriumCvMass",
           &Mutation::Mixture::mixtureEquilibriumCvMass,
           "Returns the current equilibrium mixture averaged specific heat at"
           " constant volume in J/kg-K.  This method assumes that the current "
           "state of"
           " the mixture was already set using the equilibrate() method.")

      .def("mixtureFrozenGamma", &Mutation::Mixture::mixtureFrozenGamma,
           "Returns the ratio of mixture averaged specific heats which is a "
           "unitless quantity.")

      .def("mixtureEquilibriumGamma", &Mutation::Mixture::mixtureFrozenGamma,
           "Returns the current equilibrium ratio of mixture averaged specific "
           "heats which is a unitless quantity.")

      .def(
          "X",
          [](const Mutation::Mixture &self) {
            std::vector<double> x(self.X(), self.X() + self.nSpecies());
            return py::array(py::cast(x));
          },
          "Returns the current species mole fractions.")

      .def(
          "Y",
          [](const Mutation::Mixture &self) {
            std::vector<double> y(self.Y(), self.Y() + self.nSpecies());
            return py::array(py::cast(y));
          },
          "Returns the current species mass fractions.")

      .def("setBField", &Mutation::Mixture::setBField,
           "Sets the current magnitude of the magnetic field in teslas.")

      .def("getBField", &Mutation::Mixture::getBField,
           "Gets the magnitude of the magnetic field in teslas.")

      .def("speciesIndex", &Mutation::Mixture::speciesIndex,
           "Returns the index in the species array assigned to the "
           "species with the given name. If no species exists with "
           "that name, then the index is < 0.")

      .def("elementIndex", &Mutation::Mixture::elementIndex,
           "Returns the index in the element array assigned to the "
           "element with the given name. If no element exists with "
           "that name, then the index is < 0.")

      .def("speciesCharge", &Mutation::Mixture::speciesCharge,
           "Returns the charge of the species with the given index.")

      .def("speciesThermoValidAtT", &Mutation::Mixture::speciesThermoValidAtT,
           "Returns true if the thermodynamic data in this database is valid "
           "at the"
           "given temperature for the given species index.")

      .def(
          "equilibriumComposition",
          [](const Mutation::Mixture &self, double T, double P) {
            std::vector<double> x(self.nSpecies());
            self.equilibriumComposition(T, P, x.data());
            return py::array(py::cast(x));
          },
          "Computes the equilibrium composition of the mixture at the given "
          "fixed"
          " temperature and pressure using the default elemental composition.")

      .def(
          "equilibriumComposition",
          [](const Mutation::Mixture &self, double T, double P,
             std::vector<double> xe) {
            std::vector<double> x(self.nSpecies());
            self.equilibriumComposition(T, P, xe.data(), x.data());
            return py::array(py::cast(x));
          },
          "Computes the equilibrium composition of the mixture at the given "
          "fixed"
          " temperature, pressure, and elemental moles/mole fractions.")

      .def(
          "equilibrate",
          [](const Mutation::Mixture &self, double T, double P) {
            self.equilibrate(T, P);
          },
          "Equilibrates the mixture to the given temperature, pressure, and "
          "elemental mole fractions. "
          "If the mole fractions are not given, then the default fractions are "
          "used.")

      .def(
          "speciesHOverRT",
          [](const Mutation::Mixture &self) {
            std::vector<double> h_i(self.nSpecies());
            self.speciesCpOverR(h_i.data());
            return py::array(py::cast(h_i));
          },
          "Returns the unitless vector of species enthalpies \f$ H_i / R_u T "
          "\f$.")

      .def("mixtureHMole",
           static_cast<double (Mutation::Mixture::*)(void) const>(
               &Mutation::Mixture::mixtureHMole),
           "Returns the mixture averaged enthalpy in J/mol.")

      .def("mixtureHMole",
           static_cast<double (Mutation::Mixture::*)(double) const>(
               &Mutation::Mixture::mixtureHMole),
           "Returns the mixture averaged enthalpy in J/mol at the given "
           "temperature.")

      .def("mixtureHMass",
           static_cast<double (Mutation::Mixture::*)(void) const>(
               &Mutation::Mixture::mixtureHMass),
           "Returns the mixture averaged enthalpy in J/kg.")

      .def("mixtureHMass",
           static_cast<double (Mutation::Mixture::*)(double) const>(
               &Mutation::Mixture::mixtureHMass),
           "Returns the mixture averaged enthalpy in J/kg at the given "
           "temperature.")

      .def("mixtureHMinusH0Mass", &Mutation::Mixture::mixtureHMinusH0Mass,
           "Returns the mixture averaged enthalpy minus the enthalpy of the "
           "mixture at 0K with.")

      .def("mixtureEnergyMole", &Mutation::Mixture::mixtureEnergyMole,
           "Returns the mixture energy in J/mol.")

      .def("mixtureEnergyMass", &Mutation::Mixture::mixtureEnergyMass,
           "Returns the mixture energy in J/kg.")

      .def("frozenSoundSpeed", &Mutation::Mixture::frozenSoundSpeed,
           "Returns the frozen sound speed of the mixture in m/s.")

      .def("equilibriumSoundSpeed", &Mutation::Mixture::equilibriumSoundSpeed,
           "Returns the current equilibrium sound speed of the mixture in m/s.")

      .def(
          "convert_rho_to_y",
          [](const Mutation::Mixture &self, std::vector<double> x) {
            std::vector<double> y(self.nSpecies());
            self.Mutation::Mixture::Thermodynamics::convert<
                Mutation::Thermodynamics::RHO_TO_Y>(x.data(), y.data());
            return py::array(py::cast(y));
          },
          "Converts species densities to mass fraction.")

      .def(
          "convert_rho_to_x",
          [](const Mutation::Mixture &self, std::vector<double> x) {
            std::vector<double> y(self.nSpecies());
            self.Mutation::Mixture::Thermodynamics::convert<
                Mutation::Thermodynamics::RHO_TO_X>(x.data(), y.data());
            return py::array(py::cast(y));
          },
          "Converts species densities to mole fraction.")

      .def(
          "convert_y_to_x",
          [](const Mutation::Mixture &self, std::vector<double> x) {
            std::vector<double> y(self.nSpecies());
            self.Mutation::Mixture::Thermodynamics::convert<
                Mutation::Thermodynamics::Y_TO_X>(x.data(), y.data());
            return py::array(py::cast(y));
          },
          "Converts species mass to mole fraction.")

      .def(
          "convert_x_to_y",
          [](const Mutation::Mixture &self, std::vector<double> x) {
            std::vector<double> y(self.nSpecies());
            self.Mutation::Mixture::Thermodynamics::convert<
                Mutation::Thermodynamics::X_TO_Y>(x.data(), y.data());
            return py::array(py::cast(y));
          },
          "Converts species mole to mass fraction.")

      .def(
          "convert_x_to_xe",
          [](const Mutation::Mixture &self, std::vector<double> x) {
            std::vector<double> y(self.nElements());
            self.Mutation::Mixture::Thermodynamics::convert<
                Mutation::Thermodynamics::X_TO_XE>(x.data(), y.data());
            return py::array(py::cast(y));
          },
          "Converts species mole to elemental mole fraction.")

      .def(
          "convert_y_to_ye",
          [](const Mutation::Mixture &self, std::vector<double> x) {
            std::vector<double> y(self.nElements());
            self.Mutation::Mixture::Thermodynamics::convert<
                Mutation::Thermodynamics::Y_TO_YE>(x.data(), y.data());
            return py::array(py::cast(y));
          },
          "Converts species mass to elemental mass fraction.")

      .def("dRhodP", &Mutation::Mixture::dRhodP,
           "Returns the density derivative with respect to pressure for the "
           "current equilibrium state.")

      .def(
          "dXidT",
          [](const Mutation::Mixture &self) {
            std::vector<double> value(self.nSpecies());
            self.Mutation::Mixture::dXidT(value.data());
            return py::array(py::cast(value));
          },
          "Returns the derivative of mole fraction w.t.r temperature for the "
          "given equilibrium mixture."
          " dxdt - on return, the dX_i/dT derivatives")

      .def(
          "dXidP",
          [](const Mutation::Mixture &self) {
            std::vector<double> value(self.nSpecies());
            self.Mutation::Mixture::dXidP(value.data());
            return py::array(py::cast(value));
          },
          "Returns the species derivatives of mole fraction w.r.t. pressure "
          "for the"
          " given equilibrium mixture.  Note that it is assumed the state "
          "model is an"
          " equilibrium one.")

      .def("nReactions", &Mutation::Mixture::nReactions,
           "Returns the number of reactions in the mechanism.")

      .def(
          "forwardRateCoefficients",
          [](Mutation::Mixture &self) {
            std::vector<double> value(self.nReactions());
            self.Mutation::Mixture::forwardRateCoefficients(value.data());
            return py::array(py::cast(value));
          },
          "Fills the vector kf with the forward rate coefficients \f$ k_{f,j} "
          "\f$"
          " for each reaction.  The coefficients are computed based on the "
          "rate laws"
          " assigned to each individual reaction.")

      .def(
          "forwardRatesOfProgress",
          [](Mutation::Mixture &self) {
            std::vector<double> value(self.nReactions());
            self.Mutation::Mixture::forwardRatesOfProgress(value.data());
            return py::array(py::cast(value));
          },
          "Fills the vector ropf with the forward rate of progress variables "
          "for each reaction.")

      .def(
          "backwardRateCoefficients",
          [](Mutation::Mixture &self) {
            std::vector<double> value(self.nReactions());
            self.Mutation::Mixture::backwardRateCoefficients(value.data());
            return py::array(py::cast(value));
          },
          "Fills the vector kb with the backward rate coefficients for each "
          "reaction.")

      .def(
          "backwardRatesOfProgress",
          [](Mutation::Mixture &self) {
            std::vector<double> value(self.nReactions());
            self.Mutation::Mixture::backwardRatesOfProgress(value.data());
            return py::array(py::cast(value));
          },
          "Fills the vector ropb with the backward rates of progress variables "
          "for each reaction.")

      .def(
          "netRatesOfProgress",
          [](Mutation::Mixture &self) {
            std::vector<double> value(self.nReactions());
            self.Mutation::Mixture::netRatesOfProgress(value.data());
            return py::array(py::cast(value));
          },
          "Fills the vector rop with the net rates of progress for each "
          "reaction")

      .def(
          "netProductionRates",
          [](Mutation::Mixture &self) {
            std::vector<double> value(self.nSpecies());
            self.Mutation::Mixture::netProductionRates(value.data());
            return py::array(py::cast(value));
          },
          "Fills the vector wdot with the net species production rates due to "
          "the chemical reactions.")

      .def("setViscosityAlgo", &Mutation::Mixture::setViscosityAlgo,
           "Sets the viscosity algorithm.")

      .def("setThermalConductivityAlgo",
           &Mutation::Mixture::setThermalConductivityAlgo,
           "Sets the heavy particle, thermal conductivity algorithm.")

      .def("setDiffusionMatrixAlgo", &Mutation::Mixture::setDiffusionMatrixAlgo,
           "Sets the diffusion matrix algorithm.")

      .def("nCollisionPairs", &Mutation::Mixture::nCollisionPairs,
           "Returns the number of collision pairs accounted for in this "
           "mixture.")

      .def("viscosity", &Mutation::Mixture::viscosity,
           "Returns the mixture viscosity.")

      .def("frozenThermalConductivity",
           &Mutation::Mixture::frozenThermalConductivity,
           "Returns the mixture thermal conductivity for a frozen mixture."
           " To be used only at thermal equilibrium.")

      .def(
          "frozenThermalConductivityVector",
          [](Mutation::Mixture &self) {
            std::vector<double> lambda_i(self.nEnergyEqns());
            self.frozenThermalConductivityVector(lambda_i.data());
            return py::array(py::cast(lambda_i));
          },
          "Returns the mixture thermal conductivity vector for a frozen "
          "mixture according to the state model.")

      .def("heavyThermalConductivity",
           &Mutation::Mixture::heavyThermalConductivity,
           "Returns the heavy particle translational thermal conductivity "
           "using the"
           " set algorithm.")

      .def("equilibriumThermalConductivity",
           &Mutation::Mixture::equilibriumThermalConductivity,
           "Returns the mixture thermal conductivity for a mixture in "
           "thermochemical equilibrium")

      .def("internalThermalConductivity",
           &Mutation::Mixture::internalThermalConductivity,
           "Returns the internal energy thermal conductivity using Euken's "
           "formulas.")

      .def("rotationalThermalConductivity",
           &Mutation::Mixture::rotationalThermalConductivity,
           "Returns the rotational energy thermal conductivity using Euken's "
           "formulas.")

      .def("vibrationalThermalConductivity",
           &Mutation::Mixture::vibrationalThermalConductivity,
           "Returns the vibrattest_mixtureEquilibriumCvMassional energy "
           "thermal conductivity using Euken's formulas.")

      .def("electronicThermalConductivity",
           &Mutation::Mixture::electronicThermalConductivity,
           "Returns the electronic energy thermal conductivity using Euken's "
           "formulas.")

      .def("reactiveThermalConductivity",
           &Mutation::Mixture::reactiveThermalConductivity,
           "Returns the reactive thermal conductivity which accounts for "
           "reactions"
           " for mixtures in thermochemical equilibrium.")

      .def("butlerBrokawThermalConductivity",
           &Mutation::Mixture::reactiveThermalConductivity,
           "Returns the reactive thermal conductivity which accounts for "
           "reactions"
           " for mixtures in thermochemical equilibrium using the "
           "Butler-Brokaw formula.")

      .def(
          "soretThermalConductivity",
          &Mutation::Mixture::soretThermalConductivity,
          "Returns the Soret thermal conductivity the mixture in thermochemical"
          " equilibrium.")

      .def(
          "heavyThermalDiffusionRatios",
          [](Mutation::Mixture &self) {
            std::vector<double> diff_i(self.nSpecies());
            self.heavyThermalDiffusionRatios(diff_i.data());
            return py::array(py::cast(diff_i));
          },
          "Returns the heavy thermal diffusion ratios for each species.")

      .def(
          "diffusionMatrix",
          [](Mutation::Mixture &self) { return self.diffusionMatrix(); },
          "Returns the multicomponent diffusion coefficient matrix.")

      .def("meanFreePath", &Mutation::Mixture::meanFreePath,
           "Mean free path of the mixture in m.")

      .def("electronMeanFreePath", &Mutation::Mixture::electronMeanFreePath,
           "Mean free path of electrons in m.")

      .def("speciesThermalSpeed", &Mutation::Mixture::speciesThermalSpeed,
           "Thermal speed of species i in m/s.")

      .def("averageHeavyThermalSpeed",
           &Mutation::Mixture::averageHeavyThermalSpeed,
           "Average heavy particle thermal speed of mixture in m/s.")

      .def("electronThermalSpeed", &Mutation::Mixture::electronThermalSpeed,
           "Electron thermal speed of mixture in m/s.")

      .def("electronHeavyCollisionFreq",
           &Mutation::Mixture::electronHeavyCollisionFreq,
           "Electron-heavy collision frequency in 1/s.")

      .def("averageHeavyCollisionFreq",
           &Mutation::Mixture::electronHeavyCollisionFreq,
           "Average collision frequency of heavy particles in mixture in 1/s.")

      .def("electricConductivity", &Mutation::Mixture::electricConductivity,
           "Isotropic electric conductivity in S/m (no magnetic field).")

      .def("electronThermalConductivity",
           &Mutation::Mixture::electronThermalConductivity,
           py::arg("order") = 2,
           "Returns the electron thermal conductivity in W/m-K.")

      .def("electronDiffusionCoefficient",
           &Mutation::Mixture::electronDiffusionCoefficient,
           "Isotropic electron diffusion coefficient.")

      .def("electronThermalDiffusionRatio",
           &Mutation::Mixture::electronThermalDiffusionRatio,
           "Isotropic electron thermal diffusion ratio.")

      .def(
          "stefanMaxwell",
          [](Mutation::Mixture &self, std::vector<double> grad_x) {
            std::vector<double> diff_velocities(self.nSpecies());
            double electricField;

            self.stefanMaxwell(grad_x.data(), diff_velocities.data(),
                               electricField);
            return py::make_tuple(py::array(py::cast(diff_velocities)), electricField);
          },
          "Computes the species diffusion velocities and ambipolar electric "
          "field"
          " using the Ramshaw approximation of the generalized Stefan-Maxwell"
          " equations and the supplied modified driving forces.")

      .def(
          "averageDiffusionCoeffs",
          [](Mutation::Mixture &self) {
            std::vector<double> array(self.nSpecies());
            self.averageDiffusionCoeffs(array.data());
            return py::array(py::cast(array));
          },
          "Returns the average diffusion coefficients.");
}