#include "Thermodynamics.h"
#include "Constants.h"
#include "GfcEquilSolver.h"
#include "MultiPhaseEquilSolver.h"
#include "ThermoDB.h"
#include "StateModel.h"
#include "Utilities.h"

#include <set>

using namespace std;
using namespace Mutation::Numerics;
using namespace Mutation::Utilities;

//==============================================================================

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

Thermodynamics::Thermodynamics(
    const vector<string> &species_names, const string& thermo_db,
    const string& state_model )
    : mp_work1(NULL), mp_work2(NULL), mp_default_composition(NULL),
      m_has_electrons(false)
{
    // Load the species and element objects for the specified species
    loadSpeciesFromList(species_names);
    
    // Now we can load the relevant thermodynamic database
    mp_thermodb = Config::Factory<ThermoDB>::create(thermo_db, m_species);
    
    // Build the composition matrix
    m_element_matrix = RealMatrix(nSpecies(), nElements());
    
    for (int i = 0; i < nSpecies(); ++i)
        for (int j = 0; j < nElements(); ++j)
            m_element_matrix(i,j) = 
                m_species[i].nAtoms(m_elements[j].name());
    
    // Store the species molecular weights for faster access
    m_species_mw = RealVector(nSpecies());
    for (int i = 0; i < nSpecies(); ++i)
        m_species_mw(i) = m_species[i].molecularWeight();
    
    // Allocate storage for the work array
    mp_work1 = new double [nSpecies()];
    mp_work2 = new double [nSpecies()];
    mp_y     = new double [nSpecies()];
    
    // Default composition (every element has equal parts)
    mp_default_composition = new double [nElements()];
    std::fill(
        mp_default_composition, mp_default_composition + nElements(), 
        1.0 / nElements());
    
    // Allocate a new equilibrium solver
    mp_equil = new MultiPhaseEquilSolver(*this);
    
    // Allocate a new state model
    mp_state = Config::Factory<StateModel>::create(state_model, nSpecies());
    //mp_state->notifyOnUpdate(this);
}

//==============================================================================

Thermodynamics::~Thermodynamics()
{
    delete [] mp_work1;
    delete [] mp_work2;
    delete [] mp_default_composition;
    delete [] mp_y;
    
    delete mp_thermodb;
    delete mp_equil;
    delete mp_state;
}

//==============================================================================

void Thermodynamics::loadSpeciesFromList(
    const std::vector<std::string> &species_names)
{
    // Determine file paths
    string thermo_directory = 
        getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/thermo";
    string elements_path    = thermo_directory + "/elements.xml";
    string species_path     = thermo_directory + "/species.xml";
    
    // First we need to load the entire element database for use in constructing
    // our species list
    IO::XmlDocument element_doc(elements_path);
    
    IO::XmlElement::Iterator element_iter = element_doc.root().begin();
    IO::XmlElement::Iterator element_end  = element_doc.root().end();
    
    vector<Element> elements;
    
    for ( ; element_iter != element_end; ++element_iter)
        elements.push_back(Element(*element_iter));
    
    // Load the species XML database
    IO::XmlDocument species_doc(species_path);    
    string species_name;
    
    // Use a set for the species names to ensure that species are listed only
    // once and so that the find() function can be used
    set<string> species_set(species_names.begin(), species_names.end());
    set<int>    used_elements;
    
    // Iterate over all species in the database and pull out the ones that are
    // needed from the list
    IO::XmlElement::Iterator species_iter = species_doc.root().begin();
    IO::XmlElement::Iterator species_end  = species_doc.root().end();
    
    for ( ; species_iter != species_end; ++species_iter) {
        // Get the name of current species
        species_iter->getAttribute("name", species_name);
        
        // Do we need this one?
        if (species_set.find(species_name) == species_set.end())
            continue;
        
        // We do, so store it
        m_species.push_back(Species(*species_iter, elements, used_elements));
        
        // Make sure that the electron (if it exists) is stored at the beginning
        // to make life easier.
        if (species_name == "e-") {
            if (m_species.size() > 1) 
                std::swap(m_species[0], m_species.back());
            m_has_electrons = true;
        }
        
        // Clear this species from the list of species that we still need to
        // find
        species_set.erase(species_name);
        
        // Break out early when they are all loaded
        if (species_set.size() == 0)
            break;
    }
    
    // Make sure all species were loaded (todo: better error message with file
    // name and list of missing species...)
    if (species_set.size() > 0) {
        cout << "Could not find all the species in listed in the mixture!";
        cout << endl << "Missing species:" << endl;
        
        set<string>::const_iterator missing_iter = species_set.begin();
        set<string>::const_iterator missing_end  = species_set.end();
        
        for ( ; missing_iter != missing_end; ++missing_iter)
            cout << setw(10) << *missing_iter << endl;
        
        exit(1);
    }
    
    
    // Now the species are loaded and the corresponding elements are determined
    // so store only the necessary elements in the class (note the ordering of
    // elements in the database is preserved here because of set)
    set<int>::const_iterator iter = used_elements.begin();
    set<int>::const_iterator end  = used_elements.end();
    
    for ( ; iter != end; ++iter)
        m_elements.push_back(elements[*iter]);
    
    // Finally store the species and element order information for easy access
    for (int i = 0; i < m_elements.size(); ++i)
        m_element_indices[m_elements[i].name()] = i;
    
    for (int i = 0; i < m_species.size(); ++i)
        m_species_indices[m_species[i].name()] = i;
}

//==============================================================================


void Thermodynamics::setDefaultComposition(
        const std::vector<std::pair<std::string, double> >& composition)
{
    // Make sure all elements are included exactly once and 
    bool set_element [nElements()];
    std::fill(set_element, set_element+nElements(), false);
    
    vector< pair<string, double> >::const_iterator iter = 
        composition.begin();
    
    for ( ; iter != composition.end(); ++iter) {
        int index = elementIndex(iter->first);
        if (index >= 0) {
            if (set_element[index]) {
                cerr << "Error: trying to set the default elemental"
                     << " composition for element " << iter->first
                     << " more than once!" << endl;
                exit(1);
            } else {
                mp_default_composition[index] = iter->second;
                set_element[index] = true;
            }
        } else {
            cerr << "Error: trying to set the default elemental"
                 << " composition for element " << iter->first
                 << " which is not in this mixture!" << endl;
            exit(1);
        }
    }
    
    for (int i = 0; i < nElements(); ++i) {
        if (!set_element[i]) {
            cerr << "Error: did not include element " << elementName(i)
                 << " while setting the default elemental compsotion"
                 << " of the mixture!" << endl;
            exit(1);
        }
    }
    
    //cout << "Default Composition Set:" << endl;
    //for (int i = 0; i < nElements(); ++i)
    //    cout << elementName(i) << " " << mp_default_composition[i] << endl;
        
    
    // Scale the fractions to sum to one
    RealVecWrapper wrapper(mp_default_composition, nElements());
    wrapper = wrapper / wrapper.sum();
}

//==============================================================================

void Thermodynamics::setStateTPX(
    const double* const T, const double* const P, const double* const X)
{
    mp_state->setStateTPX(T, P, X);
    convert<X_TO_Y>(X, mp_y);
}

//==============================================================================

void Thermodynamics::setStateTPY(
        const double* const T, const double* const P, const double* const Y)
{
    std::copy(Y, Y+nSpecies(), mp_y);
    convert<Y_TO_X>(Y, mp_work1);
    mp_state->setStateTPX(T, P, mp_work1);
}

//==============================================================================

double Thermodynamics::T() const {
    return mp_state->T();
}

//==============================================================================

double Thermodynamics::Tr() const {
    return mp_state->Tr();
}

//==============================================================================

double Thermodynamics::Tv() const {
    return mp_state->Tv();
}

//==============================================================================

double Thermodynamics::Te() const {
    return mp_state->Te();
}

//==============================================================================

double Thermodynamics::Tel() const {
    return mp_state->Tel();
}

//==============================================================================

double Thermodynamics::P() const {
    return mp_state->P();
}

//==============================================================================

const double* const Thermodynamics::X() const {
    return mp_state->X();
}

//==============================================================================

const double* const Thermodynamics::Y() const {
    return mp_y;
}

//==============================================================================

double Thermodynamics::standardStateT() const {
    return mp_thermodb->standardTemperature();
}

//==============================================================================

double Thermodynamics::standardStateP() const {
    return mp_thermodb->standardPressure();
}

//==============================================================================

double Thermodynamics::mixtureMw() const 
{
    return dot(Numerics::asVector(mp_state->X(), nSpecies()), m_species_mw);
}

//==============================================================================

void Thermodynamics::equilibrate(
    double T, double P, const double* const p_c, double* const p_X, 
    bool set_state)
{
    
    /// @todo Add an internal mp_X array in EquilibriumSolver class instead of
    /// needing the work array on the outside (dangerous to need work array's
    /// because other functions could be using them and they get clobbered when
    /// equilibrate is called.
    
    mp_equil->equilibrate(T, P, p_c, p_X);
    //convert<CONC_TO_X>(p_X, p_X);
    
    //cout << "After Return:" << endl;
    //for (int i = 0; i < nSpecies(); ++i)
    //    cout << setw(14) << p_X[i];
    //cout << endl;
    
    if (set_state) {
        convert<X_TO_Y>(p_X, mp_y);
        mp_state->setStateTPX(T, P, p_X);
    }

    
}

//==============================================================================

void Thermodynamics::equilibrate(double T, double P)
{
    equilibrate(T, P, mp_default_composition, mp_work1);
}

//==============================================================================

double Thermodynamics::numberDensity(const double T, const double P) const
{
    return P / (KB * T);
}

//==============================================================================

double Thermodynamics::numberDensity() const 
{
    double Xe = (hasElectrons() ? mp_state->X()[0] : 0.0);
    return mp_state->P() / KB * 
        ((1.0 - Xe) / mp_state->T() + Xe / mp_state->Te());
}

//==============================================================================

double Thermodynamics::pressure(
    const double T, const double rho, const double *const Y) const
{
    double pressure = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        pressure += Y[i] / speciesMw(i);
    pressure *= rho * T * RU;
    return pressure;
}

//==============================================================================
    
double Thermodynamics::density(
    const double T, const double P, const double *const X) const
{
    double density = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        density += X[i] * speciesMw(i);
    density *= P / (RU * T);
    return density;
}

//==============================================================================

double Thermodynamics::density() const
{
    //return numberDensity() * mixtureMw() / NA;
    return density(mp_state->T(), mp_state->P(), mp_state->X());
}

//==============================================================================

void Thermodynamics::speciesCpOverR(double *const p_cp) const
{
    mp_thermodb->cp(
        mp_state->T(), mp_state->Te(), mp_state->Tr(), mp_state->Tv(),
        mp_state->Tel(), p_cp, NULL, NULL, NULL, NULL);
}

//==============================================================================

void Thermodynamics::speciesCpOverR(double T, double* const p_cp) const
{
    mp_thermodb->cp(
        T, T, T, T, T, p_cp, NULL, NULL, NULL, NULL);
}

//==============================================================================

void Thermodynamics::speciesCpOverR(
    double Th, double Te, double Tr, double Tv, double Tel, double *const p_cp, 
    double *const p_cpt, double *const p_cpr, double *const p_cpv, 
    double *const p_cpel) const
{
    mp_thermodb->cp(
        Th, Te, Tr, Tv, Tel, p_cp, p_cpt, p_cpr, p_cpv, p_cpel);
}

//==============================================================================

void Thermodynamics::speciesCvOverR(
    double Th, double Te, double Tr, double Tv, double Tel, double *const p_cv, 
    double *const p_cvt, double *const p_cvr, double *const p_cvv, 
    double *const p_cvel) const
{
    mp_thermodb->cp(Th, Te, Tr, Tv, Tel, p_cv, p_cvt, p_cvr, p_cvv, p_cvel);
    
    if (p_cv != NULL) {
        for (int i = 0; i < nSpecies(); ++i)
            p_cv[i] -= 1.0;
    }
    
    if (p_cvt != NULL) {
        for (int i = 0; i < nSpecies(); ++i)
            p_cvt[i] -= 1.0;
    }
    
    if (p_cvr != NULL) {
        for (int i = 0; i < nSpecies(); ++i)
            p_cvr[i] -= 1.0;
    }
    
    if (p_cvv != NULL) {
        for (int i = 0; i < nSpecies(); ++i)
            p_cvv[i] -= 1.0;
    }
    
    if (p_cvel != NULL) {
        for (int i = 0; i < nSpecies(); ++i)
            p_cvel[i] -= 1.0;
    }
}

//==============================================================================

double Thermodynamics::mixtureFrozenCpMole() const
{
    double cp = 0.0;
    speciesCpOverR(mp_work1);
    for (int i = 0; i < nSpecies(); ++i)
        cp += mp_work1[i] * X()[i];
    return (cp * RU);
}

//==============================================================================

double Thermodynamics::mixtureFrozenCpMass() const 
{
    return mixtureFrozenCpMole() / mixtureMw();
}

//==============================================================================

double Thermodynamics::mixtureEquilibriumCpMole()
{
    mp_equil->dXdT(mp_work1);
    speciesHOverRT(mp_work2);
    double cp = 0.0;
    for (int i = 0; i < nSpecies(); i++)
        cp += mp_work1[i] * mp_work2[i];
    
    return (cp * RU * T() + mixtureFrozenCpMole());
}

//==============================================================================

double Thermodynamics::mixtureEquilibriumCpMass()
{
    mp_equil->dXdT(mp_work1);
    double dMwdT = 0.0;
    double Mwmix = 0.0;
    for (int i = 0; i < nSpecies(); ++i) {
        dMwdT += mp_work1[i] * speciesMw(i);
        Mwmix += X()[i] * speciesMw(i);
    }
    
    for (int i = 0; i < nSpecies(); ++i)
        mp_work2[i] = speciesMw(i)*(mp_work1[i]*Mwmix - X()[i]*dMwdT)/
            (Mwmix*Mwmix);
        
    speciesHOverRT(mp_work1);
    double cp = 0.0;
    for (int i = 0; i < nSpecies(); i++)
        cp += mp_work1[i] * mp_work2[i] / speciesMw(i);
    
    return (cp * RU * T() + mixtureFrozenCpMass());
}

//==============================================================================

void Thermodynamics::dXidT(double* const dxdt)
{
    /*const double eps = 1.0e-6;
    
    // Compute the current elemental fraction
    elementFractions(Xeq, mp_work1);
    
    // Compute equilibrium mole fractions at a perturbed temperature
    equilibrate(T*(1.0+eps), P, mp_work1, mp_work2, false);
    
    for (int i = 0; i < nSpecies(); ++i)
        dxdt[i] = (mp_work2[i] - Xeq[i]) / (T * eps);*/
    
    mp_equil->dXdT(dxdt);
}

//==============================================================================

double Thermodynamics::dRhodP(
    double T, double P, const double* const Xeq)
{
    const double eps = 1.0e-6;
    
    // Compute the current elemental fraction
    elementFractions(Xeq, mp_work1);
    
    // Compute equilibrium mole fractions at a perturbed pressure
    equilibrate(T, P*(1.0+eps), mp_work1, mp_work2, false);
    
    // Compute drho/dP
    double drhodp = 0.0;
    for (int i = 0; i < nSpecies(); ++i)
        drhodp += speciesMw(i) * ((1.0 + eps) * mp_work2[i] - Xeq[i]);
    return (drhodp / (RU * T * eps));
}

//==============================================================================

double Thermodynamics::mixtureFrozenCvMole() const
{
    return mixtureFrozenCpMole() - RU;
}

//==============================================================================

double Thermodynamics::mixtureFrozenCvMass() const
{
    return (mixtureFrozenCpMole() - RU) / mixtureMw();
}

//==============================================================================

double Thermodynamics::mixtureEquilibriumCvMass() 
{
    // We assume that Y is already in equilibrium at T and P
    // Need to perturb current equilibrium conditions and compute dh/dT, dh/dP,
    // drho/T and drho/P
    
    double rho, rhoP, rhoT, T, P;
    double e, eP, eT, dT, dP;
    double drhodt, drhodp, dedt, dedp;
    
    // Get current mixture properties
    rho = density();
    e = mixtureEnergyMass();
    T = this->T();
    P = this->P();
    elementFractions(X(), mp_work2);
    
    // Perturb pressure
    dP = std::sqrt(NumConst<double>::eps)*P;
    equilibrate(T, P+dP, mp_work2, mp_work1);
    rhoP = density();
    eP   = mixtureEnergyMass();
    
    // Perturb temperature
    dT = std::sqrt(NumConst<double>::eps)*T;
    equilibrate(T+dT, P, mp_work2, mp_work1);
    rhoT  = density();
    eT    = mixtureEnergyMass();
    
    // Return to current state
    equilibrate(T, P, mp_work2, mp_work1);
    drhodt = (rhoT-rho)/dT;
    drhodp = (rhoP-rho)/dP;
    dedt = (eT-e)/dT;
    dedp = (eP-e)/dP;
    
    //cout << endl;
    //cout << T << " " << dT << " " << P << " " << dP << endl;
    //cout << drhodt << " " << drhodp << " " << dedt << " " << dedp << endl;
    
    // Compute Cv
    return dedt - dedp*drhodt/drhodp;
}

//==============================================================================
    
double Thermodynamics::mixtureFrozenGamma() const 
{
    double cp = mixtureFrozenCpMole();
    return (cp / (cp - RU));
}

//==============================================================================

double Thermodynamics::mixtureEquilibriumGamma()
{
    return mixtureEquilibriumCpMass() / mixtureEquilibriumCvMass();
}

//==============================================================================

void Thermodynamics::speciesHOverRT(
    double* const h, double* const ht, double* const hr, double* const hv,
    double* const hel, double* const hf) const
{
    mp_thermodb->enthalpy(
        mp_state->T(), mp_state->Te(), mp_state->Tr(), mp_state->Tv(),
        mp_state->Tel(), h, ht, hr, hv, hel, hf);
}

//==============================================================================

void Thermodynamics::speciesHOverRT(double T, double* const h) const 
{
    mp_thermodb->enthalpy(
        T, T, T, T ,T, h, NULL, NULL, NULL, NULL, NULL);
}

//==============================================================================

double Thermodynamics::mixtureHMole() const
{
    double h = 0.0;
    speciesHOverRT(mp_work1);
    for (int i = 0; i < nSpecies(); ++i)
        h += mp_work1[i] * X()[i];
    return (h * RU * mp_state->T());
}

//==============================================================================

double Thermodynamics::mixtureHMass() const 
{
    return mixtureHMole() / mixtureMw();
}

//==============================================================================

void Thermodynamics::speciesSOverR(double *const p_s) const
{
    mp_thermodb->entropy(
        mp_state->T(), mp_state->Te(), mp_state->Tr(), mp_state->Tv(),
        mp_state->Tel(), standardStateP(), p_s, NULL, NULL, NULL, NULL);
    
    double lnp = std::log(mp_state->P() / standardStateP());
    for (int i = 0; i < nSpecies(); ++i)
        if (m_species[i].phase() == GAS) 
            p_s[i] -= lnp;
}

//==============================================================================

double Thermodynamics::mixtureSMole() const 
{
    double s = 0;
    speciesSOverR(mp_work1);
    for (int i = 0; i < nSpecies(); ++i)
        s += (mp_work1[i] - std::log(X()[i])) * X()[i];
    return (s * RU);
}

//==============================================================================

double Thermodynamics::mixtureSMass() const {
    return mixtureSMole() / mixtureMw();
}

//==============================================================================

void Thermodynamics::speciesGOverRT(double* const p_g) const
{
    mp_thermodb->gibbs(
        mp_state->T(), mp_state->Te(), mp_state->Tr(), mp_state->Tv(),
        mp_state->Tel(), standardStateP(), p_g, NULL, NULL, NULL, NULL);
    
    double lnp = std::log(mp_state->P() / standardStateP());
    for (int i = 0; i < nSpecies(); ++i)
        if (m_species[i].phase() == GAS) 
            p_g[i] += lnp;
}

//==============================================================================

void Thermodynamics::speciesGOverRT(double T, double P, double* const p_g) const
{
    mp_thermodb->gibbs(T, T, T, T, T, standardStateP(), p_g, NULL, NULL, NULL, NULL);
    
    double lnp = std::log(P / standardStateP());
    for (int i = 0; i < nSpecies(); ++i)
        if (m_species[i].phase() == GAS)
            p_g[i] += lnp;
}

//==============================================================================

void Thermodynamics::elementMoles(
    const double *const species_N, double *const element_N) const
{
    asVector(element_N, nElements()) = 
        asVector(species_N, nSpecies()) * m_element_matrix;
}

//==============================================================================

void Thermodynamics::elementFractions(
    const double* const Xs, double* const Xe) const
{
    RealVecWrapper wrapper(Xe, nElements());
    wrapper = asVector(Xs, nSpecies()) * m_element_matrix;
    double sum = wrapper.sum();
    for (int i = 0; i < nElements(); ++i)
        Xe[i] /= sum;
    //wrapper = wrapper / wrapper.sum();
}

//==============================================================================

void Thermodynamics::surfaceMassBalance(
    const double *const p_Yke, const double *const p_Ykg, const double T, 
    const double P, const double Bg, double &Bc, double &hw, double *const p_Xs)
{
    const int ne = nElements();
    const int ns = nSpecies();
    
    double* p_Yw = new double [ne];
    double* p_Xw = new double [ne];
    double* p_X  = (p_Xs != NULL ? p_Xs : new double [ns]);
    double* p_h  = new double [ns]; 
    
    // Initialize the wall element fractions to be the pyrolysis gas fractions
    double sum = 0.0;
    for (int i = 0; i < ne; ++i) {
        p_Yw[i] = p_Yke[i] + Bg*p_Ykg[i];
        sum += p_Yw[i];
    }
    
    // Use "large" amount of carbon to simulate infinite char
    int ic = elementIndex("C");
    //double carbon = std::min(1000.0, std::max(100.0,1000.0*Bg));
    double carbon = 10.0*Bg + 10.0;
    p_Yw[ic] += carbon;
    sum += carbon;
    
    for (int i = 0; i < ne; ++i)
        p_Yw[i] /= sum;
    
    // Compute equilibrium
    convert<YE_TO_XE>(p_Yw, p_Xw);    
    equilibrate(T, P, p_Xw, p_X, false);
    
    // Compute the gas mass fractions at the wall
    double mwg = 0.0;
    
    for (int i = 0; i < ne; ++i)
        p_Yw[i] = 0.0;
    
    for (int j = 0; j < ns; ++j) {
        if (species(j).phase() == GAS) {
            mwg += speciesMw(j) * p_X[j];
            for (int i = 0; i < ne; ++i)
                p_Yw[i] += elementMatrix()(j,i) * p_X[j];
        }
    }
    
    for (int i = 0; i < ne; ++i)
        p_Yw[i] *= atomicMass(i) / mwg;
    
    // Compute char mass blowing rate
    Bc = (p_Yke[ic] + Bg*p_Ykg[ic] - p_Yw[ic]*(1.0 + Bg)) / (p_Yw[ic] - 1.0);
    Bc = std::max(Bc, 0.0);
    
    // Compute the gas enthalpy
    speciesHOverRT(T, p_h);
    
    hw = 0.0;
    for (int i = 0; i < ns; ++i)
        if (species(i).phase() == GAS)
            hw += p_X[i] * p_h[i];
    
    hw *= RU * T / mwg;
    
    delete [] p_Yw;
    delete [] p_Xw;
    if (p_X != p_Xs) delete [] p_X;
    delete [] p_h;
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation


