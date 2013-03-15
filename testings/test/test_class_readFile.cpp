/*
 * test_class_De.cpp
 *
 *  Created on: Mar 7, 2013
 *      Author: didi
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test Module for Boost Testings"
#include <iostream>
#include <fstream>
#include <string>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "readFile.hpp"

using namespace std;

BOOST_AUTO_TEST_SUITE( test_class_position)

BOOST_AUTO_TEST_CASE( sample_test) {
	
	int pos = 12 ;
	int somme= 4;

	BOOST_CHECK_EQUAL(pos, 12);
	BOOST_CHECK(pos == 12);
	BOOST_CHECK_EQUAL((pos+2), 14);

	BOOST_CHECK_EQUAL((ReadFile::compareFiles()), 1);

	

}


BOOST_AUTO_TEST_SUITE_END()

