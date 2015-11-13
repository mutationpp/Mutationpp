#include "../general/MixtureOptions.h"
#include <boost/python.hpp>

using namespace boost::python;

/**
 * Python wrapper definition for the MixtureOptions class. All member
 * functions using base data types, strings or classes already exposed
 * to Python as arguments or return values are exposed. The
 * ones which do not are kept unexposed.
 */

void export_MixtureOptions()
{
	/**
	 * Overloaded member functions wrappers
	 */
	bool (Mutation::MixtureOptions::*getLoadTransport)() const =
		&Mutation::MixtureOptions::loadTransport;
	void (Mutation::MixtureOptions::*setLoadTransport)(bool load) =
		&Mutation::MixtureOptions::loadTransport;

	/**
	 * Python class definition
	 */
	class_<Mutation::MixtureOptions>("MixtureOptions", init<>())
		.def(init<std::string>())
		.def(
			"getSource",
			&Mutation::MixtureOptions::getSource,
			return_value_policy< copy_const_reference >())
		.def(
			"setDefaultOptions",
			&Mutation::MixtureOptions::setDefaultOptions)
		.def(
			"getSpeciesDescriptor",
			&Mutation::MixtureOptions::getSpeciesDescriptor,
			return_value_policy< copy_const_reference >())
		.def(
			"setSpeciesDescriptor",
			&Mutation::MixtureOptions::setSpeciesDescriptor)
		.def(
			"getStateModel",
			&Mutation::MixtureOptions::getStateModel,
			return_value_policy< copy_const_reference >())
		.def(
			"setStateModel",
			&Mutation::MixtureOptions::setStateModel)
		.def(
			"getThermodynamicDatabase",
			&Mutation::MixtureOptions::getThermodynamicDatabase,
			return_value_policy< copy_const_reference >())
		.def(
			"setThermodynamicDatabase",
			&Mutation::MixtureOptions::setThermodynamicDatabase)
		.def(
			"getMechanism",
			&Mutation::MixtureOptions::getMechanism,
			return_value_policy< copy_const_reference >())
		.def(
			"setMechanism",
			&Mutation::MixtureOptions::setMechanism)
		.def(
			"getLoadTransport",
			getLoadTransport)
		.def(
			"setLoadTransport",
			setLoadTransport)
		.def(
			"getViscosityAlgorithm",
			&Mutation::MixtureOptions::getViscosityAlgorithm,
			return_value_policy< copy_const_reference >())
		.def(
			"setViscosityAlgorithm",
			&Mutation::MixtureOptions::setViscosityAlgorithm)
		.def(
			"getThermalConductivityAlgorithm",
			&Mutation::MixtureOptions::getThermalConductivityAlgorithm,
			return_value_policy< copy_const_reference >())
		.def(
			"setThermalConductivityAlgorithm",
			&Mutation::MixtureOptions::setThermalConductivityAlgorithm)
		.def(
			"setDefaultComposition",
			&Mutation::MixtureOptions::setDefaultComposition)
		;
}
