
#include "mutation++.h"

#include <string>
#include <fstream>

/**
 * Base class for all Mutation++ test fixtures.  This is used as a framework to
 * setup a testing fixture for use in Boost.Test.
 */
class MppTestFixture
{
public:

    /**
     * Constructor, sets up the test fixture.
     */
    MppTestFixture(std::string name)
        : mp_mix(NULL), mp_sp1(NULL), m_name(name)
    {
        BOOST_TEST_MESSAGE(std::string("Setup fixture: ") + m_name);
    }
    
    /**
     * Destructor, tears down the test fixture.
     */
    virtual ~MppTestFixture()
    {
        BOOST_TEST_MESSAGE(std::string("Tearing down fixture: ") + m_name);
        
        if (mp_mix != NULL) delete mp_mix;
        if (mp_sp1 != NULL) delete [] mp_sp1;
    }
    
    /**
     * Returns a pointer to an instantiated mixture object.
     */
    Mutation::Mixture* const mix() 
    { 
        if (mp_mix == NULL) {
            createMixture();
        }
        
        return mp_mix; 
    }
    
    /**
     * Returns a work array allocated for number of species in the mixture.
     */
    double* const sp1()
    {
        if (mp_sp1 == NULL) {
            if (mp_mix == NULL) createMixture();
            mp_sp1 = new double [mp_mix->nSpecies()];
            BOOST_CHECK_MESSAGE(mp_sp1 != NULL, "sp1 is NULL!");
        }
        
        return mp_sp1;
    }
    
    /**
     * Compares results for equilibrium calculations to a given data file.
     */
    template <typename T>
    void compareEquilibriumValues(std::string filename, T func, int nvals)
    {
        // Open file
        std::ifstream file(filename.c_str());
        BOOST_CHECK_MESSAGE(file.is_open(), 
            "Input file '" + filename + "' could not be opened!");
    
        // Read number of temperatures and pressures from top of file
        int nt = 0;
        int np = 0;
        
        file >> nt;
        file >> np;
        
        double T, P;
        
        double values [nvals];
        double result [nvals];
        
        // Now loop over each pressure and temperature to compare the results
        for (int ip = 0; ip < np; ++ip) {
            // Read the pressure
            file >> P;
            
            for (int it = 0; it < nt; ++it) {
                // Read the temperature and associated values
                file >> T;
                for (int i = 0; i < nvals; ++i)
                    file >> values[i];
                
                // Equilibrate
                mix()->equilibrate(T, P);
                
                // Now compute the given values using the equilibrium mixture
                func(this, result);
                
                // Compare
                for (int i = 0; i < nvals; ++i)
                    BOOST_CHECK_CLOSE(values[i], result[i], 1.0e-3);
            }
        }
        
        file.close();
    }

protected:

    /**
     * Should be overloaded by the child class to set the MixtureOptions which
     * will be used to create a Mixture object suitable for this test fixture.
     */
    virtual Mutation::MixtureOptions mixtureOptions() const = 0;

private:

    /**
     * Creates a mixture object using the MixtureOptions defined by the derived
     * class.
     */
    void createMixture()
    {
        mp_mix = new Mutation::Mixture(mixtureOptions());
        BOOST_CHECK_MESSAGE(mp_mix != NULL, "Mixture is NULL!");
    }

private:
    
    Mutation::Mixture* mp_mix;
    double*            mp_sp1;
    std::string        m_name;

};


/// Simple macro to simplify creating children of MppTest
#define MPP_TEST_FIXTURE(__NAME__,__MIX__,__OPTIONS__)\
class __NAME__ : public MppTestFixture\
{\
public:\
    __NAME__ ()\
        : MppTestFixture( #__NAME__ )\
    { }\
protected:\
    Mutation::MixtureOptions mixtureOptions() const\
    {\
        Mutation::MixtureOptions options( __MIX__ );\
        __OPTIONS__\
        return options;\
    }\
};


