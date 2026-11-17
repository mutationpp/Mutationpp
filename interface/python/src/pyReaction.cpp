#include <Reaction.h>
#include <RateLaws.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <stdexcept>

namespace py = nanobind;

/**
 * Python wrapper for the Arrhenius rate law.
 */
void py_export_RateLaw(py::module_ &m) {
  py::class_<Mutation::Kinetics::Arrhenius>(m, "Arrhenius")
    .def_prop_ro("log_pre_exponential",
                 &Mutation::Kinetics::Arrhenius::logA,
                 "Returns the logarithm of the pre-exponential for this "
                 "rate law.")
    .def_prop_ro("exponent",
                 &Mutation::Kinetics::Arrhenius::n,
                 "Returns the temperature exponent for this rate law.")
    .def_prop_ro("activation_temperature",
                 &Mutation::Kinetics::Arrhenius::T,
                 "Returns the activation temperature for this rate law.");
}

/**
 * Python wrapper for the Reaction class.
 */
void py_export_Reaction(py::module_ &m) {
  py::class_<Mutation::Kinetics::Reaction>(m, "Reaction")
    .def_prop_ro("reactants",
                 &Mutation::Kinetics::Reaction::reactants,
                 "Returns the reactants for this reaction.")

    .def_prop_ro("products",
                 &Mutation::Kinetics::Reaction::products,
                 "Returns the products for this reaction.")

    .def("rate_law",
         [](const Mutation::Kinetics::Reaction &self)
             -> const Mutation::Kinetics::Arrhenius & {
           const auto* p_arrhenius =
               dynamic_cast<const Mutation::Kinetics::Arrhenius*>(
                   self.rateLaw());
           if (p_arrhenius == nullptr)
             throw std::runtime_error(
                 "only Arrhenius rate laws are exposed to Python");
           return *p_arrhenius;
         },
         py::rv_policy::reference_internal,
         "Returns the rate law for this reaction.")

    .def_prop_ro("fwd_rate_coeff_temperature",
                 &Mutation::Kinetics::Reaction::fwdRateTemperature,
                 "Returns the temperature at which the forward rate "
                 "coefficient of this reaction is evaluated.")

    .def_prop_ro("rev_rate_coeff_temperature",
                 &Mutation::Kinetics::Reaction::revRateTemperature,
                 "Returns the temperature at which the reverse rate "
                 "coefficient of this reaction is evaluated.")

    .def_prop_ro("formula",
                 &Mutation::Kinetics::Reaction::formula,
                 "Returns the formula for this reaction.")

    .def_prop_ro("order",
                 &Mutation::Kinetics::Reaction::order,
                 "Returns the order of this reaction.")

    .def_prop_ro("is_reversible",
                 &Mutation::Kinetics::Reaction::isReversible,
                 "Returns true if this reaction is reversible.")

    .def_prop_ro("is_thirdbody",
                 &Mutation::Kinetics::Reaction::isThirdbody,
                 "Returns true if this reaction has a third body.")

    .def_prop_ro("efficiencies",
                 &Mutation::Kinetics::Reaction::efficiencies,
                 "Returns the (species index, efficiency) pairs of the "
                 "third bodies in this reaction.")

    .def_prop_ro("reaction_type",
                 [](const Mutation::Kinetics::Reaction &self) {
                   return std::string(
                       Mutation::Kinetics::reactionTypeString(self.type()));
                 },
                 "Returns the type of this reaction.");
}
