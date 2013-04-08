// Boost testing header and setup
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Thermo
#include <boost/test/unit_test.hpp>

#include "MppTestFixture.h"

// Setup the test fixture options
MPP_TEST_FIXTURE( Air11RRHO, "air11",
    options.setThermodynamicDatabase("RRHO");
)
#define TEST_FIXTURE_NAME Air11RRHO
#define TEST_FIXTURE_DATA_DIRECTORY "data/air11/RRHO"

#define TEST_SET_STATE
//#define TEST_SPECIES_X
//#define TEST_SPECIES_X



//#define TEST_MIXTURE_S
//#define TEST_MIXTURE_E

//#define TEST_MIXTURE_CP
#define TEST_MIXTURE_H

//#define TEST_MIXTURE_HR
//#define TEST_MIXTURE_HF
//#define TEST_MIXTURE_HEL
// #define TEST_MIXTURE_HV


#define TEST_MIXTURE_A_F
#define TEST_MIXTURE_A_EQ
#define TEST_MIXTURE_GAMMA
#define TEST_MIXTURE_GAM_EQ


#define TEST_MIXTURE_S
#define TEST_MIXTURE_E
#include "ThermoTests.h"

