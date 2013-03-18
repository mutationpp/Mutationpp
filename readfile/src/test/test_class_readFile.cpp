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
#include <stdlib.h>
#include <iomanip>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "readFile.hpp"

using namespace std;

BOOST_AUTO_TEST_SUITE( test_class_readfile)

BOOST_AUTO_TEST_CASE( sample_test) {
	
	int pos = 12 ;
	int somme= 4;
	ReadFile r;
	BOOST_CHECK_EQUAL(pos, 12);
	BOOST_CHECK(pos == 12);
	BOOST_CHECK_EQUAL((pos+2), 14);
	string path = boost::filesystem::current_path().generic_string();
//	std::cout << std::endl << std::endl << "PATH ==> "  << path << std::endl << std::endl << std::endl;

//	std::cout << std::endl << std::endl << "SUCCES ? : "  << (r.compareFiles()) << std::endl << std::endl << std::endl;
	BOOST_CHECK_EQUAL((r.compareFiles()),1);

	

}


BOOST_AUTO_TEST_SUITE_END()

