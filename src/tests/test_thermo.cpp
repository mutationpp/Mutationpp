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
#define TEST_FIXTURE_DATA_DIRECTORY "air11/RRHO"

#define TEST_SET_STATE
#define TEST_MIXTURE_CP
#define TEST_MIXTURE_H
#define TEST_SPECIES_X

#include "ThermoTests.h"

