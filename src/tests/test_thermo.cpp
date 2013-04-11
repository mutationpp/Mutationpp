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

//#define TEST_SPECIES_X 	// 180 fails
//#define TEST_SPECIES_Y 	// 180 fails 



#define TEST_MIXTURE_CV  	// 148 fails
//#define TEST_MIXTURE_CV_EQ 	// 148 fails
//#define TEST_MIXTURE_CP_EQ 	// 148 fails 
//#define TEST_MIXTURE_GAM_EQ  	// 54 Fails
//#define TEST_MIXTURE_A_EQ  	// 47 fails 





//#define TEST_MIXTURE_HR  	// fails in compilation 
//#define TEST_MIXTURE_HF  	// fails in compilation
//#define TEST_MIXTURE_HEL 	// fails in compilation
//#define TEST_MIXTURE_HV  	// fails in compilation



#define TEST_MIXTURE_H  	// passes
#define TEST_MIXTURE_CP 	// passes 
#define TEST_MIXTURE_A_F 	// passes
#define TEST_MIXTURE_GAMMA	// passes
#define TEST_MIXTURE_RHO	// passes
#define TEST_MIXTURE_S		// passes
#define TEST_MIXTURE_E		// passes



#include "ThermoTests.h"

