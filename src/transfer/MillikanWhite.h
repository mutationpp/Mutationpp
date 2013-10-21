
#ifndef MILLIKAN_WHITE_H
#define MILLIKAN_WHITE_H

#include <string>
#include <vector>
#include <cmath>
#include "Thermodynamics.h"

namespace Mutation {
    namespace Transfer {

/** 
 * Data storage for single Millikan-White collision pair.
 */
class MillikanWhitePartner
{
public:

    /**
     * Construct using given values of a, b, and mu.
     */
    MillikanWhitePartner(double a, double b, double mu)
        : m_a(a), m_b(b), m_mu(mu)
    {}
    
    /**
     * Construct using default values of a and b based on mu and theta_v.
     */
    MillikanWhitePartner(double mu, double theta)
        : m_a(1.16E-3*std::sqrt(mu)*std::pow(theta, 4.0/3.0)),
          m_b(0.015*std::pow(mu, 0.25)),
          m_mu(mu)
    { }
    
    double a()  const {return m_a; }
    double b()  const {return m_b; }
    double mu() const {return m_mu;}

private:

    double m_a;
    double m_b;
    double m_mu;
};

/**
 * Data storage for all collision partners of a vibrating species.
 */
class MillikanWhiteVibrator
{
public:

    /**
     * Load available data from VT.xml and use defaults when data does not 
     * exist.
     */
    MillikanWhiteVibrator(
        Mutation::Utilities::IO::XmlElement& node,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Use default data for all collision pairs with this vibrator.
     */
    MillikanWhiteVibrator(
        const std::string& name,
        const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns the (a,b,mu) data for the i'th heavy species colliding with this
     * vibrator.
     */
    const MillikanWhitePartner& operator[](int i) const {
        assert(i >= 0);
        assert(i < m_partners.size());
        return m_partners[i];
    }
    
    /**
     * Returns the limiting cross-section for vibrational excitation at 50,000 K
     * for this vibrator.  Defaults to 3.0E-17 unless given differently in the 
     * data file.
     */
    double omega() const {
        return m_omegav;
    }
    
    /**
     * Returns the index of this vibrator in the mixture species list.
     */
    int index() const {
        return m_index;
    }

private:

    std::vector<MillikanWhitePartner> m_partners;
    double m_omegav;
    int m_index;
};


/**
 * Data storage for all of the information needed to evaluate the Millikan-White
 * vibration-translation energy transfer.
 */
class MillikanWhite
{
public:
    /**
     * Loads all of the data associated with the Millikan-White model.
     */
    MillikanWhite(const Mutation::Thermodynamics::Thermodynamics& thermo);
    
    /**
     * Returns the number of vibrators being kept in this model.
     */
    int nVibrators() const {
        return m_vibrators.size();
    }
    
    /**
     * Returns the i'th vibrator data.
     */
    const MillikanWhiteVibrator& operator[](int i) const {
        assert(i >= 0);
        assert(i < m_vibrators.size());
        return m_vibrators[i];
    }

private:

    std::vector<MillikanWhiteVibrator> m_vibrators;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // MILLIKAN_WHITE_H
