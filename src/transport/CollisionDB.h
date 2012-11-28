#ifndef TRANSPORT_COLLISIONDB_H
#define TRANSPORT_COLLISIONDB_H

#include <cmath>
#include <vector>
#include <string>
#include <iostream>

#include "Thermodynamics.h"
#include "utilities.h"
#include "Numerics.h"

/**
 * Implements the function exp(A + B*ln(x) + C*ln(x)^2 + ... ) for an arbitrary
 * number of coefficients.
 */
template <int N>
class CollisionFunc
{
public:

    CollisionFunc()
    {
        for (int i = 0; i < N; ++i)
            m_c[i] = 0.0;
    }
    
    inline double operator()(const double lnT) const {
        double val = m_c[N-1];
        for (int i = N-2; i >= 0; --i)
            val = lnT * val + m_c[i];
        return (val == 0.0 ? 0.0 : std::exp(val));
    }
    
    template <int S> friend std::istream& operator >> (
        std::istream& in, CollisionFunc<S>& f);
    
protected:

    double m_c[N];
    
};

template <int N>
std::istream& operator >> (std::istream& in, CollisionFunc<N>& f)
{
    for (int i = N-1; i >= 0; --i)
        in >> f.m_c[i];
    return in;
}

class CollisionFunc4 : public CollisionFunc<4>
{
public:
    CollisionFunc4() { }

    CollisionFunc4(double c1, double c2, double c3, double c4)
    {
        m_c[0] = c1;
        m_c[1] = c2;
        m_c[2] = c3;
        m_c[3] = c4;
    }
};

class CollisionFunc5 : public CollisionFunc<5>
{
public:
    CollisionFunc5() { }

    CollisionFunc5(double c1, double c2, double c3, double c4, double c5)
    {
        m_c[0] = c1;
        m_c[1] = c2;
        m_c[2] = c3;
        m_c[3] = c4;
        m_c[4] = c5;
    }
};


/**
 * Simple object representing a collision pair which treats any collision
 * pair and it's reverse as equal collisions.  This is used in the 
 * constructor for building a useful collision pair map.
 */
class CollisionPair
{
public:
    
    /**
     * Constructs a collision pair given the collision pair name.  Assumes the 
     * net charge is neutral.
     */
    CollisionPair(const std::string& pair)
    {
        std::vector<std::string> tokens;
        std::string str = pair;
        
        utils::StringUtils::eraseAll(str, ".");
        utils::StringUtils::tokenize(str, tokens, "-");
        
        initialize(tokens[0], tokens[1]);
        m_charge = 0;
    }
    
    /**
     * Constructs a collision pair given two species objects.
     */
    CollisionPair(const Species& s1, const Species& s2) 
    {
        initialize(s1.name(), s2.name());
        m_charge = s1.charge() * s2.charge();
    }
    
    /**
     * Equality operator.
     */
    bool operator == (const CollisionPair& right) const {
        return m_collision == right.m_collision;
    }
    
    /**
     * Less than operator used to sort collision pairs.
     */
    bool operator < (const CollisionPair& right) const {
        return m_collision < right.m_collision;
    }
    
    /**
     * Returns the name of the collision pair.
     */
    const std::string& name() const { 
        return m_collision; 
    }    
    
    /**
     * Returns the charge of this collision pair. Negative values are attractive
     * positive values are attractive and zero represents neutral collisions.
     */
    int charge() const { 
        return m_charge; 
    }
    
    /**
     * Returns the name of the first species in the pair.
     */
    const std::string& speciesName1() const {
        return m_s1;
    }
    
    /**
     * Returns the name of the second species in the pair.
     */
    const std::string& speciesName2() const {
        return m_s2;
    }
    
private:

    /**
     * Some common initialization stuff.
     */
    void initialize(const std::string &name1, const std::string &name2)
    {
        m_s1 = name1;
        m_s2 = name2;
        
        // Just some simple logic to make sure the names match in the database
        // This will eventually be fased out by using an XML database
        if (m_s1[m_s1.size()-1] == '+')
            m_s1[m_s1.size()-1] = 'p';
        else if (m_s1[m_s1.size()-1] == '-')
            m_s1[m_s1.size()-1] = 'm';
        
        if (m_s2[m_s2.size()-1] == '+')
            m_s2[m_s2.size()-1] = 'p';
        else if (m_s2[m_s2.size()-1] == '-')
            m_s2[m_s2.size()-1] = 'm';
        
        if (m_s1 < m_s2)
            m_collision = "." + m_s1 + "-" + m_s2 + ".";
        else
            m_collision = "." + m_s2 + "-" + m_s1 + ".";
    }

private:

    std::string m_collision;
    std::string m_s1;
    std::string m_s2;
    int         m_charge;
    
}; // CollisionPair


class CollisionDBTester;

/**
 * This class is responsible for managing the computation of collision integrals
 * which are loaded from a database.  The idea is to provide a wrapper class
 * which hides the implementation of the database from user classes which
 * require collision integrals in order to minimize effects of database format
 * changes.
 */
class CollisionDB
{
public:
    /**
     * Construct a CollisionDB object for a given list of species.  Note that
     * all possible species pair interactions will be loaded.
     */
    explicit CollisionDB(const Thermodynamics& thermo);
    
    /**
     * Destructor.
     */
    ~CollisionDB() {}
    
    /**
     * Returns the number of species in the database.
     */
    size_t nSpecies() const {
        return m_ns;
    }
    
    /**
     * Returns the total number of collision pairs loaded from the species list
     * (including ones not found in the database, ie: NS * (NS + 1) / 2).
     */
    size_t nCollisionPairs() const { 
        return m_ncollisions; 
    }
    
    /**
     * Returns the mass of a single atom/molecule of each species as a NS 
     * dimensional Vector.
     */
    const Numerics::RealVector& mass() const {
        return m_mass; 
    }
    
    /**
     * Returns the symmetric matrix \f$ A_{ij} = m_i + m_j \f$ where \f$ m_i \f$
     * is the mass of one molecule of species i.  
     */
    const Numerics::RealSymMat& massSum() const {
        return m_mass_sum; 
    }
    
    /**
     * Returns the symmetric matrix \f$ A_{ij} = m_i * m_j \f$ where \f$ m_i \f$
     * is the mass of one molecule of species i.
     */
    const Numerics::RealSymMat& massProd() const {
        return m_mass_prod; 
    }
    
    /**
     * Returns the reduced molecular mass as a symmetric matrix in kg/molecule.
     */
    const Numerics::RealSymMat& reducedMass() const {
        return m_red_mass;
    }
    
    /**
     * Computes Q_{ij}^{1,1}(T) for i,j \in \set{H}.  The i,j ordering was
     * determined by the constructor based on the order of the species listed in
     * the species vector. 
     */
    const Numerics::RealSymMat& Q11(
        const double Th, const double Te, const double nd, 
        const double *const p_x) {
        updateCollisionData(Th, Te, nd, p_x, Q11IJ);
        return m_Q11;
    }

    /**
     * Computes Q_{ij}^{2,2}(T) for i,j \in \set{H}.  The i,j ordering was
     * determined by the constructor based on the order of the species listed in
     * the species vector. 
     */
    const Numerics::RealSymMat& Q22(
        const double Th, const double Te, const double nd, 
        const double *const p_x) {
        updateCollisionData(Th, Te, nd, p_x, Q22IJ);
        return m_Q22;
    }

    /**
     * Returns the dimensionless ratio \f$ Q_{ij}^{(2,2)} / Q_{ij}^{(1,1)} \f$
     * as a symmetric matrix.
     */
    const Numerics::RealSymMat& Astar(
        const double Th, const double Te, const double nd, 
        const double *const p_x) {
        updateCollisionData(Th, Te, nd, p_x, ASTAR);
        return m_Ast;
    }

    /**
     * Returns the dimensionless ratio \f$ (5Q_{ij}^{(1,2)} - 4Q_{ij}^{(1,3)}) /
     * Q_{ij}^{(1,1)} \f$ as a symmetric matrix.
     */
    const Numerics::RealSymMat& Bstar(
        const double Th, const double Te, const double nd, 
        const double *const p_x) {
        updateCollisionData(Th, Te, nd, p_x, BSTAR);
        return m_Bst;
    }
    
    /**
     * Returns the pure species shear viscosities.
     */
    const Numerics::RealVector& etai(
        const double Th, const double Te, const double nd, 
        const double *const p_x) {
        updateCollisionData(Th, Te, nd, p_x, ETAI);
        return m_eta;
    }
    
    /**
     * Returns binary diffusion coefficients.
     */
    const Numerics::RealSymMat& nDij(
        const double Th, const double Te, const double nd,
        const double *const p_x) {
        updateCollisionData(Th, Te, nd, p_x, NDIJ);
        return m_Dij;
    }
    
    friend class CollisionDBTester;
    
private:

    /*static void computeQ11(const double& T, double* const p_Q11, void* p_data) {
        CollisionDB& collisions = 
            static_cast<CollisionDB&>(*(CollisionDB*)(p_data));
        const int ns = collisions.nSpecies();
        const int nc = collisions.nCollisionPairs();
        
        std::cout << "compute Q11" << ns << nc << std::endl;       
        
    }*/

    enum CollisionData {
        Q11IJ = 0,
        Q22IJ,
        ASTAR,
        BSTAR,
        ETAI,
        NDIJ,
        DATA_SIZE
    };

    /**
     * Loads the heavy particle collision pairs.
     */
    void loadCollisionIntegrals(const std::vector<Species>& species);

    /**
     * Checks to see if a string represents a collision pair (ie: ".Ar-H.")
     */
    static bool isValidCollisionString(const std::string& str);
    
    /**
     * Updates the collision data specified by data.
     */
    void updateCollisionData(
        const double Th, const double Te, const double nd, 
        const double *const p_x, const CollisionData data);

private:

    // Sizing parameters
    int m_ns;
    int m_ncollisions;
    int m_em_index;
    
    // Keeps track of which collisions are neutral and which are charged
    std::vector<int> m_neutral_indices;
    std::vector<int> m_attract_indices;
    std::vector<int> m_repulse_indices;
    
    // Storage for neutral-neutral and neutral-charge collision integral 
    // functions
    std::vector<CollisionFunc4> m_Q11_funcs;
    std::vector<CollisionFunc4> m_Q22_funcs;
    std::vector<CollisionFunc4> m_Bst_funcs;
    
    // Charge-charge collision integrals (defined in .cpp file)
    static const CollisionFunc5 sm_Q11_att;
    static const CollisionFunc5 sm_Q11_rep;
    static const CollisionFunc5 sm_Q22_att;
    static const CollisionFunc5 sm_Q22_rep;
    static const CollisionFunc5 sm_Bst_att;
    static const CollisionFunc5 sm_Bst_rep;
    
    // Mass quantities
    Numerics::RealVector m_mass;
    Numerics::RealSymMat m_mass_sum;
    Numerics::RealSymMat m_mass_prod;
    Numerics::RealSymMat m_red_mass;
    
    // Stores the last computed collision integral data    
    Numerics::RealSymMat m_Q11;
    Numerics::RealSymMat m_Q22;
    Numerics::RealSymMat m_Ast;
    Numerics::RealSymMat m_Bst;
    Numerics::RealVector m_eta;
    Numerics::RealSymMat m_Dij;
    
    // Keeps track of last temperature a particular set of collision data values
    // was updated
    double mp_last_T [DATA_SIZE];
};

#endif // TRANSPORT_COLLISIONDB_H
