/**
 * @file CollisionDB.h
 *
 * @brief Declaration of CollisionDB type.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef TRANSPORT_COLLISION_DB_H
#define TRANSPORT_COLLISION_DB_H

#include "CollisionPair.h"
#include "CollisionGroup.h"
#include "Constants.h"
#include "Thermodynamics.h"
#include "XMLite.h"

#include <cassert>
#include <complex>
#include <map>
#include <string>
#include <vector>

namespace Mutation {
    namespace Transport {

/**
 * Provides collision integral data loaded from a database file.
 */
class CollisionDB
{
public:

    /**
     * Constructs a new CollisionDB type.
     */
    CollisionDB(
        const std::string& db_name,
        const Thermodynamics::Thermodynamics& thermo);

    /// Returns number of species in the database.
    int nSpecies() const;

    /// Returns number of heavy particles in the database.
    int nHeavy() const;

    /// Returns the number of collision pairs in this database.
    int size() const { return m_pairs.size(); }

    /// Returns an array of species masses in kg.
    const Eigen::ArrayXd& mass() const { return m_mass; }

    /// Just a wrapper around thermo().X().
    Eigen::Map<const Eigen::ArrayXd> X() const;

    /// Just a wrapper around thermo().Y().
    Eigen::Map<const Eigen::ArrayXd> Y() const;

    /// Returns the reference of thermo.
    const Thermodynamics::Thermodynamics& thermo() const { return m_thermo; }

    /// Returns i'th CollisionPair in the database.
    const CollisionPair& operator [] (int i) const {
        assert(i >= 0 && i < size());
        return m_pairs[i];
    }

    /**
     * Returns the collision group corresponding to name, updated to the current
     * state.  If the group is not loaded yet, it is first loaded from the
     * database.
     */
    const CollisionGroup& group(const std::string& name);

    /// Provides Q11 collision integral for the electron-electron interaction.
    double Q11ee() { return group("Q11ee")[0]; }

    /// Provides Q11 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q11ei() { return group("Q11ei").array(); }

    /// Provides Q11 collision integrals for diagonal heavy-heavy interactions.
    const Eigen::ArrayXd& Q11ii() { return group("Q11ii").array(); }

    /// Provides Q11 collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Q11ij() { return group("Q11ij").array(); }

    /// Provides Q12 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q12ei() { return group("Q12ei").array(); }

    /// Provides Q13 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q13ei() { return group("Q13ei").array(); }

    /// Provides Q14 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q14ei() { return group("Q14ei").array(); }

    /// Provides Q15 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q15ei() { return group("Q15ei").array(); }

    /// Provides Q22 collision integral for the electron-electron interaction.
    double Q22ee() { return group("Q22ee")[0]; }

    /// Provides Q22 collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Q22ei() { return group("Q22ei").array(); }

    /// Provides Q22 collision integrals for diagonal heavy-heavy interactions.
    const Eigen::ArrayXd& Q22ii() { return group("Q22ii").array(); }

    /// Provides Q22 collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Q22ij() { return group("Q22ij").array(); }

    /// Provides Q23 collision integrals for electron-electron interaction.
    double Q23ee() { return group("Q23ee")[0]; }

    /// Provides Q24 collision integrals for electron-electron interaction.
    double Q24ee() { return group("Q24ee")[0]; }

    /// Provides A* collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Astei() { return group("Astei").array(); }

    /// Provides A* collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Astij() { return group("Astij").array(); }

    /// Provides B* collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Bstei() { return group("Bstei").array(); }

    /// Provides B* collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Bstij() { return group("Bstij").array(); }

    /// Provides C* collision integrals for electron-heavy interactions.
    const Eigen::ArrayXd& Cstei() { return group("Cstei").array(); }

    /// Provides C* collision integrals for heavy-heavy interactions.
    const Eigen::ArrayXd& Cstij() { return group("Cstij").array(); }

    /// Returns the pure, heavy species shear viscosities array.
    const Eigen::ArrayXd& etai();

    /// Returns binary diffusion coefficients for electron-heavy interactions.
    const Eigen::ArrayXd& nDei();

    /// Returns binary diffusion coefficients for heavy-heavy interactions.
    const Eigen::ArrayXd& nDij();

    /**
     * Returns the average diffusion coefficients.
     * \f[ D_{im} = \frac{(1-x_i)}{\sum_{j\ne i}x_j/\mathscr{D}_{ij}} \f]
     * If include_one_minus_x is false, the \f$1-x_i\f$ term taken to be 1.
     * This is sometimes useful when only the partial evaluation of Dim is 
     * required.
     */
    const Eigen::ArrayXd& Dim(bool include_one_minus_x = true);

    /// Returns the \f$\Lambda^{01}_{ei}\f$ array.
    const Eigen::ArrayXd& L01ei();

    /// Returns the \f$\Lambda^{02}_{ei}\f$ array.
    const Eigen::ArrayXd& L02ei();

    /**
     * Returns the factor which must be multiplied by the result of the
     * Lee() function to retrieve \f$\L_{ee}^{pq}\f$.
     */
    double Leefac() {
        const double Te = m_thermo.Te();
        const double me = m_mass(0);
        const double xe = m_thermo.X()[0];
        return 16.*m_thermo.P()/(3.*KB*Te)*std::sqrt(me/(TWOPI*KB*Te));
    }

    /**
     * Returns the symmetric \f$\Lambda_{ee}\f$ matrix of the given size divided
     * by the \f$\frac{64x_e}{75k_B}\sqrt{\frac{m_e}{2\pi k_B T_e}}\f$.
     */
    template <int SIZE>
    Eigen::Matrix<double, SIZE, SIZE> Lee()
    {
        Eigen::Matrix<double, SIZE, SIZE> L;
        const int nh = nHeavy();

        // L00ee
        const Eigen::ArrayXd& Q11 = Q11ei();
        L(0,0) = dotxh(Q11);

        // Return if done
        if (SIZE == 1) return L;

        const Eigen::ArrayXd& Q12 = Q12ei();
        const Eigen::ArrayXd& Q13 = Q13ei();
        const double Q22 = Q22ee();

        // L01ee, L10ee
        L(0,1) = dotxh(2.5*Q11-3.*Q12);
        L(1,0) = L(0,1);

        // L11ee
        L(1,1) = dotxh(6.25*Q11-15.*Q12+12.*Q13);
        L(1,1) += X()(0)*SQRT2*Q22;

        // Return if done
        if (SIZE == 2) return L;

        const Eigen::ArrayXd& Q14 = Q14ei();
        const Eigen::ArrayXd& Q15 = Q15ei();
        const double Q23 = Q23ee();
        const double Q24 = Q24ee();

        // L02ee, L20ee
        L(0,2) = dotxh(35./8.*Q11-10.5*Q12+6.*Q13);
        L(2,0) = L(0,2);

        // L12ee, L21ee
        L(1,2) = dotxh(175./16.*Q11-315./8.*Q12+57.*Q13-30.*Q14);
        L(1,2) += X()(0)*SQRT2*(1.75*Q22 - 2.*Q23);
        L(2,1) = L(1,2);

        // L22ee
        L(2,2) = dotxh(1225./64.*Q11-735./8.*Q12+199.5*Q13-210.*Q14+90.*Q15);
        L(2,2) += X()(0)*SQRT2*(77./16.*Q22-7.*Q23+5.*Q24);

        return L;
    }

    /**
     * Returns the \f$\L_{ee}^{Bpq}\f$ matrix with the given size.
     */
    template <int SIZE>
    Eigen::Matrix<std::complex<double>, SIZE, SIZE> LBee()
    {
        Eigen::Matrix<std::complex<double>, SIZE, SIZE> LB;
        const double fac = m_thermo.getBField()*QE/(KB*m_thermo.Te());

        LB(0,0).imag(fac);
        if (SIZE == 1) return LB;

        LB(1,1).imag(2.5*fac);
        if (SIZE == 2) return LB;

        LB(2,2).imag(35./8.*fac);
        return LB;
    }

    template <int ORDER>
    Eigen::Matrix<double,ORDER,ORDER> L1inv() {
        return Lee<ORDER>().inverse()/Leefac();
    }

    template <int ORDER>
    Eigen::Matrix<std::complex<double>,ORDER,ORDER> L2inv()
    {
        Eigen::Matrix<std::complex<double>,ORDER,ORDER> L2 = LBee<ORDER>();
        L2.real() += Leefac()*Lee<ORDER>();
        return L2.inverse();
    }

    template <int ORDER>
    Eigen::Matrix<double,ORDER,1> L1sol(const Eigen::Matrix<double,ORDER,1>& b)
    {
        return L1inv<ORDER>()*b;
    }

    template <int ORDER>
    Eigen::Matrix<std::complex<double>,ORDER,1> L2sol(
        const Eigen::Matrix<double,ORDER,1>& b)
    {
        return L2inv<ORDER>()*b;
    }

private:

    /// Helper function to evaluate \f$ sum_{i\in\mathcal{H}} x_i a_i \f$.
    template <typename E>
    double dotxh(const Eigen::ArrayBase<E>& a) {
        return (X()*a).tail(nHeavy()).sum();
    }

    /// Types of collision integral groups (ie: electron-heavy, ...).
    enum GroupType {
        EE = 0, EI, II, IJ, BAD_TYPE
    };

    /// Determines the type of group from the group name.
    GroupType groupType(const std::string& name);

private:

    Mutation::Utilities::IO::XmlDocument m_database;
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    // Storing gaseous species (all and heavy)
    const int m_ng;
    const int m_nh;

    // Tabulation parameters
    bool m_tabulate;
    double m_table_min;
    double m_table_max;
    double m_table_del;

    // List of collision pairs
    std::vector<CollisionPair> m_pairs;

    // CollisionGroup container
    std::map<std::string, CollisionGroup> m_groups;

    // Derived data and helper arrays
    Eigen::ArrayXd m_mass;
    Eigen::ArrayXd m_etai;
    Eigen::ArrayXd m_etafac;
    Eigen::ArrayXd m_nDei;
    Eigen::ArrayXd m_Deifac;
    Eigen::ArrayXd m_nDij;
    Eigen::ArrayXd m_Dijfac;
    Eigen::ArrayXd m_Dim;
    Eigen::ArrayXd m_L01ei;
    Eigen::ArrayXd m_L02ei;
};

    } // namespace Transport
} // namespace Mutation


#endif // TRANSPORT_COLLISION_DB_H
