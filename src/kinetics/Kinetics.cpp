#include "Kinetics.h"
#include "Constants.h"
#include "Utilities.h"
#include "MillikanWhite.h"

using namespace std;
using namespace Mutation::Numerics;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;


namespace Mutation {
    namespace Kinetics {

using Mutation::Thermodynamics::Thermodynamics;

//==============================================================================

Kinetics::Kinetics(
    const Thermodynamics& thermo, string mechanism)
    : m_num_rxns(0), m_thermo(thermo), m_jacobian(thermo), mp_g(NULL),
      m_T_last(-1.0)
{
    if (mechanism == "none")
        return;
    
    mechanism = 
        getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/mechanisms/" +
        mechanism + ".xml";

    // Open the mechanism file as an XML document
    IO::XmlDocument doc(mechanism);        
    IO::XmlElement root = doc.root();
    
    if (root.tag() != "mechanism") {
        cout << "Root element in mechanism file " << mechanism
             << " is not of 'mechanism' type!";
        exit(1); 
    }
    
    // Now loop over all of the reaction nodes and add each reaction to the
    // corresponding data structure pieces
    IO::XmlElement::Iterator iter = root.begin();
    for ( ; iter != root.end(); ++iter) {        
        if (iter->tag() == "reaction")
            addReaction(Reaction(*iter, thermo));
        else if (iter->tag() == "arrhenius_units")
            Arrhenius::setUnits(*iter);
    }
    
    // Finally close the reaction mechanism
    closeReactions(true);
}

//==============================================================================

void Kinetics::addReaction(const Reaction& reaction)
{
    // Add reaction to reaction list
    m_reactions.push_back(reaction);
    
    // Insert the reactants
    m_reactants.addReaction(m_num_rxns, reaction.reactants());
    
    // Insert products
    if (reaction.isReversible())
        m_rev_prods.addReaction(m_num_rxns, reaction.products());
    else
        m_irr_prods.addReaction(m_num_rxns, reaction.products());
    
    // Add thirdbodies if necessary
    if (reaction.isThirdbody())
        m_thirdbodies.addReaction(m_num_rxns, reaction.efficiencies());
    
    // Insert new ratelaw to keep track of
    m_rates.addRateCoefficient(
        m_num_rxns, reaction.rateLaw());
    
    // Add the reaction to the jacobian managaer
    m_jacobian.addReaction(reaction);
    
    m_num_rxns++;
}

//==============================================================================

void Kinetics::closeReactions(const bool validate_mechanism) 
{
    const size_t ns = m_thermo.nSpecies();
    RealVector tempsp(ns);
    
    // Validate the mechanism
    if (validate_mechanism) {
        bool is_valid = true;
        //cout << "Validating reaction mechanism..." << endl;
        
        // Check for duplicate reactions
        RealVector stoichi(ns);
        RealVector stoichj(ns);
        for (size_t i = 0; i < m_num_rxns-1; ++i) {
            for (size_t k = 0; k < ns; ++k)
                stoichi(k) =
                    m_reactions[i].product(k) - m_reactions[i].reactant(k);
            stoichi.normalize();
            for (size_t j = i+1; j < m_num_rxns; ++j) {
                for (size_t k = 0; k < ns; ++k)
                    stoichj(k) =
                        m_reactions[j].product(k) - m_reactions[j].reactant(k);
                stoichj.normalize();
                if (stoichi == stoichj) {
                    cerr << "Reactions " << i+1 << " \"" 
                         << m_reactions[i].formula()
                         << "\" and " << j+1 << " \""
                         << m_reactions[j].formula()
                         << "\" are identical." << endl;
                    is_valid = false;
                }
            }
        }
        
        // Check for elemental mass and charge conservation
        RealVector mass(m_num_rxns);
        for (int i = 0; i < m_thermo.nElements(); ++i) {
            for (size_t k = 0; k < ns; ++k)
                stoichi(k) = m_thermo.elementMatrix()(k,i);
            getReactionDelta(stoichi, mass);
            for (size_t j = 0; j < m_num_rxns; ++j) {
                if (mass(j) != 0.0) {
                    cerr << "Reaction " << j+1 << " \""
                         << m_reactions[j].formula()
                         << "\" does not conserve element "
                         << m_thermo.elementName(i) << "." << endl;
                    is_valid = false;
                }
            }
        }
        
        if (!is_valid) {
            cout << "Validation checks failed!" << endl;
            exit(1);
        }
    }
    
    // Compute dnu        
    tempsp = 1.0;
    m_dnu = RealVector(m_num_rxns);
    getReactionDelta(tempsp, m_dnu);
    
    // Allocate work arrays
    mp_g    = new double [m_thermo.nSpecies()];
    m_lnkf  = RealVector(m_num_rxns);
    m_lnkeq = RealVector(m_num_rxns);
    m_ropf  = RealVector(m_num_rxns);
    m_ropb  = RealVector(m_num_rxns);
    m_rop   = RealVector(m_num_rxns);
    
}

//==============================================================================

void Kinetics::getReactionDelta(const RealVector& s, RealVector& r)
{
    r = 0.0;
    m_reactants.decrReactions(s, r);
    m_rev_prods.incrReactions(s, r);
    m_irr_prods.incrReactions(s, r);
}

//==============================================================================

/*vector<size_t> Kinetics::speciesIndices(
    const multiset<string>& set)
{    
    multiset<string>::const_iterator iter = set.begin();
    vector<size_t> indices;
    
    for ( ; iter != set.end(); ++iter)
        indices.push_back(m_thermo.speciesIndex(*iter));
    
    return indices;
}*/

//==============================================================================

/*vector<pair<size_t, double> > Kinetics::thirdbodyEffs(
    const vector<pair<string, double> >& string_effs)
{
    vector<pair<size_t, double> > effs;
    vector<pair<string, double> >::const_iterator iter;
    
    for (iter = string_effs.begin(); iter != string_effs.end(); ++iter)
        effs.push_back(
            make_pair(m_thermo.speciesIndex(iter->first), iter->second));
    
    return effs;
}*/

//==============================================================================

void Kinetics::updateT(const double T)
{
    // We shouldn't care about temperatures that are only different by small
    // amounts
    if (abs(T - m_T_last) < 1.0E-6) return;

    // Update forward rates
    m_rates.lnForwardRateCoefficients(T, m_lnkf);
    
    // Update the equilibrium constants
    m_lnkeq = m_dnu * std::log(101325.0 / (RU * T));
    m_thermo.speciesGOverRT(mp_g);
    
    const RealVecWrapper g(mp_g, m_thermo.nSpecies());
    m_reactants.incrReactions(g, m_lnkeq);
    m_rev_prods.decrReactions(g, m_lnkeq);
    m_irr_prods.decrReactions(g, m_lnkeq);
    
    m_T_last = T;
}

//==============================================================================
    
void Kinetics::equilibriumConstants(const double T, RealVector& keq)
{      
    updateT(T);
    keq = exp(m_lnkeq);       
}

//==============================================================================

void Kinetics::forwardRateCoefficients(const double T, RealVector& kf)
{
    updateT(T);
    kf = exp(m_lnkf);
}

//==============================================================================

void Kinetics::backwardRateCoefficients(const double T, RealVector& kb)
{
    updateT(T);
    kb = exp(m_lnkf - m_lnkeq);
}

//==============================================================================

void Kinetics::forwardRatesOfProgress(
    const double T, const RealVector& conc, RealVector& ropf)
{
    forwardRateCoefficients(T, ropf);
    m_reactants.multReactions(conc, ropf);
    m_thirdbodies.multiplyThirdbodies(conc, ropf);
}

//==============================================================================

void Kinetics::backwardRatesOfProgress(
    const double T, const RealVector& conc, RealVector& ropb)
{
    backwardRateCoefficients(T, ropb);
    m_rev_prods.multReactions(conc, ropb);
    m_thirdbodies.multiplyThirdbodies(conc, ropb);
}

//==============================================================================

void Kinetics::netRatesOfProgress(
    const double T, const RealVector& conc, RealVector& rop)
{
    forwardRateCoefficients(T, m_ropf);
    m_reactants.multReactions(conc, m_ropf);        
    backwardRateCoefficients(T, m_ropb);
    m_rev_prods.multReactions(conc, m_ropb);        
    m_thirdbodies.multiplyThirdbodies(conc, rop = (m_ropf - m_ropb));
}

//==============================================================================

void Kinetics::netProductionRates(
    const double T, const double* const p_conc, double* const p_wdot)
{
    const RealVector conc(p_conc, m_thermo.nSpecies());
    RealVector wdot(p_wdot, m_thermo.nSpecies());
    
    netRatesOfProgress(T, conc, m_rop);
    m_reactants.decrSpecies(m_rop, wdot = 0.0);
    m_rev_prods.incrSpecies(m_rop, wdot);
    m_irr_prods.incrSpecies(m_rop, wdot);
    
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        p_wdot[i] = wdot(i) * m_thermo.speciesMw(i);
}

//==============================================================================

void Kinetics::netProductionRates(double* const p_wdot)
{
    RealVector conc(m_thermo.nSpecies());
    const double rho = m_thermo.density();
    const double T   = m_thermo.T();
    
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        conc(i) = m_thermo.Y()[i] * rho / m_thermo.speciesMw(i);
    
    RealVector wdot(p_wdot, m_thermo.nSpecies());
    
    netRatesOfProgress(T, conc, m_rop);
    m_reactants.decrSpecies(m_rop, wdot = 0.0);
    m_rev_prods.incrSpecies(m_rop, wdot);
    m_irr_prods.incrSpecies(m_rop, wdot);
    
    for (int i = 0; i < m_thermo.nSpecies(); ++i)
        p_wdot[i] = wdot(i) * m_thermo.speciesMw(i);
}

//==============================================================================

void Kinetics::jacobianRho(
    const double T, const double* const p_conc, double* const p_jac)
{
    updateT(T);
    double* p_kf = new double [m_num_rxns];
    double* p_kb = new double [m_num_rxns];
    
    for (int i = 0; i < m_num_rxns; ++i) {
        p_kf[i] = std::exp(m_lnkf(i));
        p_kb[i] = std::exp(m_lnkf(i) - m_lnkeq(i));
    }
    
    m_jacobian.computeJacobian(p_kf, p_kb, p_conc, p_jac);
    
    delete [] p_kf;
    delete [] p_kb;
}

//==============================================================================

double Kinetics::omegaTV()
{
    // Load Millikan-White model data on first call to this method
    static MillikanWhite data(m_thermo);
    
    const int nheavy = m_thermo.nHeavy();
    const int offset = m_thermo.hasElectrons() ? 1 : 0;
    
    for (int i = 0; i < data.nVibrators(); ++i) {
        cout << "Vibrator: " << m_thermo.speciesName(data[i].index())
             << ", omega = " << data[i].omega() << endl;
        for (int j = 0; j < nheavy; ++j) {
            cout << setw(14) << m_thermo.speciesName(data[i].index()) +
                "-" + m_thermo.speciesName(j+offset);
            cout << setw(14) << data[i][j].a();
            cout << setw(14) << data[i][j].b();
            cout << setw(14) << data[i][j].mu();
            cout << endl;
        }
        cout << endl;
    }
    
    return 0.0;
}

    } // namespace Kinetics
} // namespace Mutation

