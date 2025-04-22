#include <Reaction.h>
#include <RateManager.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/**
 * Python wrapper for the Reaction class.
 */
void py_export_RateLaw(py::module &m) {
  py::class_<Mutation::Kinetics::Arrhenius>(m, "Arrhenius")
    .def_property_readonly("log_pre_exponential",
			   &Mutation::Kinetics::Arrhenius::logA,
			   "Returns the logarithm of the pre-exponential for this rate law.")
    .def_property_readonly("exponent",
			   &Mutation::Kinetics::Arrhenius::n,
			   "Returns the temperature exponent for this rate law.")
    .def_property_readonly("activation_temperature",
			   &Mutation::Kinetics::Arrhenius::T,
			   "Returns the activation temperature for this rate law.");
}

/**
 * Python wrapper for the Reaction class.
 */
void py_export_Reaction(py::module &m) {
  py::class_<Mutation::Kinetics::Reaction>(m, "Reaction")
    .def_property_readonly("reactants",
			   &Mutation::Kinetics::Reaction::reactants,
			   "Returns the reactants for this reaction.")

    .def_property_readonly("products",
			   &Mutation::Kinetics::Reaction::products,
			   "Returns the products for this reaction.")

    .def("rate_law",
	 &Mutation::Kinetics::Reaction::rateLaw,
	 py::return_value_policy::reference_internal,
	 "Returns the rate law for this reaction.")

    .def_property_readonly("fwd_rate_coeff_temperature",
			   &Mutation::Kinetics::Reaction::fwdRateTemperature,
			   "Returns the temperature at which the fwd "
			   "rate coeff. for this reaction is evaluated.")
    .def_property_readonly("rev_rate_coeff_temperature",
			   &Mutation::Kinetics::Reaction::revRateTemperature,
			   "Returns the temperature at which the rev "
			   "rate coeff. for this reaction is evaluated.");
}
