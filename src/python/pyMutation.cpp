#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "../thermo/Thermodynamics.h"
#include "../general/Mixture.h"
#include <boost/python.hpp>
#include <boost/numpy.hpp>

using namespace boost::python;
namespace np = boost::numpy;

/**
 * This is the module definition for the Mutation++ Python wrapper.
 * The objects interfaces are kept as close as possible to the C++ interface
 * form minimal maintenance effort.
 */

void export_Mixture();
void export_MixtureOptions();

BOOST_PYTHON_MODULE(libpyMutation)
{
	np::initialize();

	export_Mixture();
	export_MixtureOptions();
}
