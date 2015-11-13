#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "../general/Mixture.h"
#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <iostream>

namespace bp = boost::python;
namespace np = boost::numpy;

/**
 * Converter functions
 */

//double* convert_nparray_cdouble(PyObject* p_pyobj) {
//	/* Todo:
//	 * - check for array dimension
//	 * - check for array contiguity
//	 * - check for array data type
//	 */
//
//	if (p_pyobj == NULL)
//		return NULL;
//	else
//		return (double*)PyArray_DATA((PyArrayObject*)p_pyobj);
//}

template <typename T>
const np::ndarray wrap_ndarray(const T* const a, int n) {
	return np::from_data
		(
			a,
			np::dtype::get_builtin<T>(),
			bp::make_tuple(n),
			bp::make_tuple(sizeof(T)),
			bp::object()
		);
}

double* convert_nparray_cdouble(np::ndarray &array) {
	/* Todo:
	 * - check for array dimension
	 * - check for array contiguity
	 * - check for array data type
	 */

	return reinterpret_cast<double*>(array.get_data());
}

//boost::python::numeric::array convert_cdouble_nparray(double* data, int size)
//{
//	// Inspired by:
//	// http://stackoverflow.com/questions/10701514/how-to-return-numpy-array-from-boostpython
//	// This is very relevant as well:
//	// http://stackoverflow.com/questions/12957492/writing-python-bindings-for-c-code-that-use-opencv
//	npy_intp pysize = size;
//
//	// create a PyObject * from pointer and data
//	PyObject* pyObj = PyArray_SimpleNewFromData(1, &pysize, NPY_DOUBLE, data);
//	boost::python::handle<> handle(pyObj);
//	boost::python::numeric::array arr(handle);
//
//	/* The problem of returning arr is twofold:
//	 * - firstly the user can modify the data which will betray
//	 *   the const-correctness
//	 * - secondly, the lifetime of the data is managed by the C++ API
//	 *   and not the lifetime of the numpy array whatsoever
//	 * But we have a simple solution:
//	*/
//
//	return arr; // copy the object. numpy owns the copy now.
//}

/**
 * Python wrapper definition for the Mixture class. The interface
 * is kept as close as possible to the C++ interface, but might be
 * missing some public members.
 */

namespace Mutation {

class MixtureWrapper : public Mixture
{
public:
	MixtureWrapper(const Mutation::MixtureOptions& options) : Mixture::Mixture(options) {};

	void averageDiffusionCoeffsWrapper(np::ndarray &Di)
	{
		double* p_Di = convert_nparray_cdouble(Di);
		this->averageDiffusionCoeffs(p_Di);
	}

	void backwardRateCoefficientsWrapper(
		np::ndarray &kb)
	{
		double* p_kb = convert_nparray_cdouble(kb);
		this->backwardRateCoefficients(p_kb);
	}

	np::ndarray diffusionMatrixWrapper()
	{
		// Todo: check bullshit happening wtf??.
		      int         n = this->nSpecies();
		const double* p_arr = this->diffusionMatrix().data();

		np::ndarray py_array =
			np::from_data
			(
				p_arr,
				np::dtype::get_builtin<double>() ,
				bp::make_tuple(n,n),
				bp::make_tuple(sizeof(double),n*sizeof(double)),
				bp::object()
			);

		return py_array;
	}

	void dXidTWrapper(np::ndarray& dXdT)
	{
		double* p_dXdT = convert_nparray_cdouble(dXdT);
		this->dXidT(p_dXdT);
	}

	void elementPotentialsWrapper(np::ndarray &lambda)
	{
		double* p_lambda = convert_nparray_cdouble(lambda);
		this->elementPotentials(p_lambda);
	}

	void equilDiffFluxFacsPWrapper(np::ndarray &F)
	{
		double* p_F = convert_nparray_cdouble(F);
		this->equilDiffFluxFacsP(p_F);
	}

	void equilibrateWrapper1(
		const double T,
		const double P,
		np::ndarray &Xe)
	{
		double* p_Xe = convert_nparray_cdouble(Xe);
		this->equilibrate(T, P, p_Xe);
	}

	void equilibrateWrapper2(
		const double T,
		const double P)
	{
		this->equilibrate(T, P);
	}

	void forwardRateCoefficientsWrapper(
		np::ndarray &kf)
	{
		double* p_kf = convert_nparray_cdouble(kf);
		this->forwardRateCoefficients(p_kf);
	}

	void netProductionRatesWrapper(
		np::ndarray &wdot)
	{
		double* p_wdot = convert_nparray_cdouble(wdot);
		this->netProductionRates(p_wdot);
	}

	void omega11iiWrapper(
		np::ndarray &omega)
	{
		double* p_omega = convert_nparray_cdouble(omega);
		this->omega11ii(p_omega);
	}

	void omega22iiWrapper(
		np::ndarray &omega)
	{
		double* p_omega = convert_nparray_cdouble(omega);
		this->omega22ii(p_omega);
	}

	void phaseMolesWrapper(np::ndarray &moles)
	{
		double* p_moles = convert_nparray_cdouble(moles);
		this->phaseMoles(p_moles);
	}

	void setStateWrapper(
		np::ndarray& mass,
		np::ndarray& energy,
		const int vars = 0) 
	{
		const double* p_mass = convert_nparray_cdouble(mass);
		const double* p_energy = convert_nparray_cdouble(energy);
		this->setState(p_mass, p_energy, vars);
	}
	
	void speciesCpOverRWrapper(
		np::ndarray &cp)
	{
		double* p_cp = convert_nparray_cdouble(cp);
		this->speciesCpOverR(p_cp);
	}

	void speciesHOverRTWrapper(
		double T,
		np::ndarray &h)
	{
		double* p_h = convert_nparray_cdouble(h);
		this->speciesHOverRT(T, p_h);
	}

	void speciesSOverRWrapper(
		np::ndarray &s)
	{
		double* p_s = convert_nparray_cdouble(s);
		this->speciesSOverR(p_s);
	}

	void speciesGOverRTWrapper(
		np::ndarray &g)
	{
		double* p_g = convert_nparray_cdouble(g);
		this->speciesGOverRT(p_g);
	}

	void stefanMaxwellWrapper(
		np::ndarray dp,
		np::ndarray& p_V,
		double &E) 
	{
		const double* p_dp = convert_nparray_cdouble(dp);
		double* p_p_V = convert_nparray_cdouble(p_V);
		this->stefanMaxwell(p_dp, p_p_V, E);
	}

	void thermalDiffusionRatiosWrapper(
		np::ndarray &k)
	{
		double* p_k = convert_nparray_cdouble(k);
		this->thermalDiffusionRatios(p_k);
	}

	const np::ndarray XWrapper() const {
		return wrap_ndarray(X(),nSpecies());
	}

	const np::ndarray YWrapper() const {
		return wrap_ndarray(Y(),nSpecies());;
	}
};

} // namespace Mutation

void export_Mixture()
{
	/**
	 * Overloaded member functions wrappers (i.e. the C++ function is
	 * overloaded, and the version of the function in which we are interested
	 * takes basic types as arguments). More complex overloading requires
	 * the definition of a wrapper function.
	 */
	double (Mutation::MixtureWrapper::*density)() const =
		&Mutation::MixtureWrapper::density;
	double (Mutation::MixtureWrapper::*numberDensity)() const =
		&Mutation::MixtureWrapper::numberDensity;
	double (Mutation::MixtureWrapper::*mixtureHMole)() const =
		&Mutation::MixtureWrapper::mixtureHMole;
	double (Mutation::MixtureWrapper::*mixtureHMass)() const =
		&Mutation::MixtureWrapper::mixtureHMass;
	double (Mutation::MixtureWrapper::*mixtureEquilibriumGamma)() =
		&Mutation::MixtureWrapper::mixtureEquilibriumGamma;
	
    bp::class_<Mutation::MixtureWrapper>("Mixture", bp::init<Mutation::MixtureOptions>())
		.def(bp::init<std::string>())
		.def(
			"averageDiffusionCoeffs",
			&Mutation::MixtureWrapper::averageDiffusionCoeffsWrapper)
		.def(
			"averageHeavyCollisionFreq",
			&Mutation::MixtureWrapper::averageHeavyCollisionFreq)
		.def(
			"averageHeavyThermalSpeed",
			&Mutation::MixtureWrapper::averageHeavyThermalSpeed)
		.def(
			"backwardRateCoefficients",
			&Mutation::MixtureWrapper::backwardRateCoefficientsWrapper)
		.def(
			"density",
			density)
		.def(
			"diffusionMatrix",
			&Mutation::MixtureWrapper::diffusionMatrixWrapper)
		.def(
			"dRhodP",
			&Mutation::MixtureWrapper::dRhodP)
		.def(
			// Todo: assess validity of call policy management
			"dXidT",
			&Mutation::MixtureWrapper::dXidTWrapper)
		.def(
			"electronHeavyCollisionFreq",
			&Mutation::MixtureWrapper::electronHeavyCollisionFreq)
		.def(
			"electronMeanFreePath",
			&Mutation::MixtureWrapper::electronMeanFreePath)
		.def(
			"electronThermalConductivity",
			&Mutation::MixtureWrapper::electronThermalConductivity)
		.def(
			"electronThermalSpeed",
			&Mutation::MixtureWrapper::electronThermalSpeed)
		.def(
			"elementPotentials",
			&Mutation::MixtureWrapper::elementPotentialsWrapper)
		.def("equilibrate",
			&Mutation::MixtureWrapper::equilibrateWrapper1)
		.def("equilibrate",
			&Mutation::MixtureWrapper::equilibrateWrapper2)
		.def(
			"equilibriumSoundSpeed",
			&Mutation::MixtureWrapper::equilibriumSoundSpeed)
		.def(
			"equilibriumThermalConductivity",
			&Mutation::MixtureWrapper::equilibriumThermalConductivity)
		.def(
			"forwardRateCoefficients",
			&Mutation::MixtureWrapper::frozenSoundSpeed)
		.def(
			"frozenSoundSpeed",
			&Mutation::MixtureWrapper::frozenSoundSpeed)
		.def(
			"getBField",
			&Mutation::MixtureWrapper::getBField)
		.def(
			"hasElectrons",
			&Mutation::MixtureWrapper::hasElectrons)
		.def(
			"heavyThermalConductivity",
			&Mutation::MixtureWrapper::heavyThermalConductivity)
		.def(
			"internalThermalConductivity",
			&Mutation::MixtureWrapper::internalThermalConductivity)
		.def(
			"meanFreePath",
			&Mutation::MixtureWrapper::meanFreePath)
		.def(
			"mixtureEnergyMass",
			&Mutation::MixtureWrapper::mixtureEnergyMass)
		.def(
			"mixtureEnergyMole",
			&Mutation::MixtureWrapper::mixtureEnergyMole)
		.def(
			"mixtureEquilibriumCpMass",
			&Mutation::MixtureWrapper::mixtureEquilibriumCpMass)
		.def(
			"mixtureEquilibriumCpMole",
			&Mutation::MixtureWrapper::mixtureEquilibriumCpMole)
		.def(
			"mixtureEquilibriumGamma",
			mixtureEquilibriumGamma)
		.def(
			"mixtureEquilibriumCvMass",
			&Mutation::MixtureWrapper::mixtureEquilibriumCvMass)
		.def(
			"mixtureFrozenCpMass",
			&Mutation::MixtureWrapper::mixtureFrozenCpMass)
		.def(
			"mixtureFrozenCpMole",
			&Mutation::MixtureWrapper::mixtureFrozenCpMole)
		.def(
			"mixtureFrozenCvMass",
			&Mutation::MixtureWrapper::mixtureFrozenCvMass)
		.def(
			"mixtureFrozenCvMole",
			&Mutation::MixtureWrapper::mixtureFrozenCvMole)
		.def(
			"mixtureFrozenGamma",
			&Mutation::MixtureWrapper::mixtureFrozenGamma)
		.def(
			"mixtureHMass",
			mixtureHMass)
		.def(
			"mixtureHMole",
			mixtureHMole)
		.def(
			"mixtureMw",
			&Mutation::MixtureWrapper::mixtureMw)
		.def(
			"mixtureSMass",
			&Mutation::MixtureWrapper::mixtureSMass)
		.def(
			"mixtureSMole",
			&Mutation::MixtureWrapper::mixtureSMole)
		.def(
    		"netProductionRates",
    		&Mutation::MixtureWrapper::netProductionRatesWrapper)
		.def(
    		"nAtoms",
    		&Mutation::MixtureWrapper::nAtoms)
		.def(
			"nCondensed",
			&Mutation::MixtureWrapper::nCondensed)
		.def(
			"nElements",
			&Mutation::MixtureWrapper::nElements)
		.def(
			"nEnergyEqns",
			&Mutation::MixtureWrapper::nEnergyEqns)
		.def(
			"nEquilibriumSteps",
			&Mutation::MixtureWrapper::nEquilibriumSteps)
		.def(
			"nEquilibriumNewtons",
			&Mutation::MixtureWrapper::nEquilibriumNewtons)
		.def(
			"nGas",
			&Mutation::MixtureWrapper::nGas)
		.def(
			"nHeavy",
    		&Mutation::MixtureWrapper::nHeavy)
    	.def(
			"nMassEqns",
			&Mutation::MixtureWrapper::nMassEqns)
		.def(
			"nMolecules",
			&Mutation::MixtureWrapper::nMolecules)
		.def(
			"nPhases",
			&Mutation::MixtureWrapper::nPhases)
		.def(
			"nSpecies",
			&Mutation::MixtureWrapper::nSpecies)
		.def(
			"numberDensity",
			numberDensity)
		.def(
			"omega11ii",
			&Mutation::MixtureWrapper::omega11iiWrapper)
		.def(
			"omega22ii",
			&Mutation::MixtureWrapper::omega22iiWrapper)
		.def(
			"phaseMoles",
			&Mutation::MixtureWrapper::phaseMolesWrapper)
		.def(
			"P",
			&Mutation::MixtureWrapper::P)
		.def(
			"reactiveThermalConductivity",
			&Mutation::MixtureWrapper::reactiveThermalConductivity)
		.def(
			"setBField",
			&Mutation::MixtureWrapper::setBField)
		.def(
			// Todo: assess validity of call policy management
			"setState",
			&Mutation::MixtureWrapper::setStateWrapper,
			bp::with_custodian_and_ward<1,2,
				bp::with_custodian_and_ward<1,3> >())
		.def(
			"sigma",
			&Mutation::MixtureWrapper::sigma)
		.def(
			"sigmaParallel",
			&Mutation::MixtureWrapper::sigmaParallel)
		.def(
			"sigmaPerpendicular",
			&Mutation::MixtureWrapper::sigmaPerpendicular)
		.def(
			"sigmaTransverse",
			&Mutation::MixtureWrapper::sigmaTransverse)
		.def(
			"soretThermalConductivity",
			&Mutation::MixtureWrapper::soretThermalConductivity)
		.def(
			"speciesCpOverR",
			&Mutation::MixtureWrapper::speciesCpOverRWrapper)
		.def(
			"speciesHOverRT",
			&Mutation::MixtureWrapper::speciesHOverRTWrapper)
		.def(
			"speciesMw",
			&Mutation::MixtureWrapper::speciesMw)
		.def(
			"speciesSOverR",
			&Mutation::MixtureWrapper::speciesSOverRWrapper)
		.def(
			"speciesGOverRT",
			&Mutation::MixtureWrapper::speciesGOverRTWrapper)
		.def(
			// Todo: assess validity of call policy management
			"stefanMaxwell",
			&Mutation::MixtureWrapper::stefanMaxwellWrapper)
		.def(
			"thermalDiffusionRatios",
			&Mutation::MixtureWrapper::thermalDiffusionRatiosWrapper)
		.def(
			"T",
			&Mutation::MixtureWrapper::T)
		.def(
			"viscosity",
			&Mutation::MixtureWrapper::viscosity)
		.def(
			"X",
			&Mutation::MixtureWrapper::XWrapper)
		.def(
			"Y",
			&Mutation::MixtureWrapper::YWrapper)
		;
}
