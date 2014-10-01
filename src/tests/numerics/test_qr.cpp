/**
 * @file test_qr.cpp
 * @brief Several test cases which test the QR and QRP decompositions.
 * @author J.B. Scoggins
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE QR
#include <boost/test/unit_test.hpp>

#include "QRP.h"
using namespace Mutation::Numerics;

//==============================================================================
BOOST_AUTO_TEST_CASE( LehmerQRSolve )
{
    for (int i = 1; i <= 10; ++i) {
        Vector<double> b(i,1.0);
        Vector<double> x(i);
        QR<double> qr((Lehmer<double>(i)));
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(Lehmer<double>(i)*x, Vector<double>(i, 1.0));
        b = Vector<double>(i,0.0);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(Lehmer<double>(i)*x, Vector<double>(i, 0.0));
    }
}

//==============================================================================
BOOST_AUTO_TEST_CASE( MinIJQRSolve )
{
    for (int i = 1; i <= 10; ++i) {
        SymmetricMatrix<double> A = MinIJ<double>(i);
        Vector<double> b(i,1.0);
        Vector<double> x(i);
        QR<double> qr(A);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 1.0));
        b = Vector<double>(i,0.0);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 0.0));
    }
}

//==============================================================================
BOOST_AUTO_TEST_CASE( IdentityQRSolve )
{
    for (int i = 1; i <= 10; ++i) {
        SymmetricMatrix<double> A = Identity<double>(i);
        Vector<double> b(i,1.0);
        Vector<double> x(i);
        QR<double> qr(A);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 1.0));
        b = Vector<double>(i,0.0);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 0.0));
    }
}

//==============================================================================
BOOST_AUTO_TEST_CASE( LehmerQRPSolve )
{
    for (int i = 2; i <= 10; ++i) {
        SymmetricMatrix<double> A = Lehmer<double>(i);
        Vector<double> b(i,1.0);
        Vector<double> x(i);
        QRP<double> qr(A);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 1.0));
        b = Vector<double>(i,0.0);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 0.0));
    }
}

//==============================================================================
BOOST_AUTO_TEST_CASE( MinIJQRPSolve )
{
    for (int i = 2; i <= 10; ++i) {
        SymmetricMatrix<double> A = MinIJ<double>(i);
        Vector<double> b(i,1.0);
        Vector<double> x(i);
        QRP<double> qr(A);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 1.0));
        b = Vector<double>(i,0.0);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 0.0));
    }
}

//==============================================================================
BOOST_AUTO_TEST_CASE( IdentityQRPSolve )
{
    for (int i = 2; i <= 10; ++i) {
        SymmetricMatrix<double> A = Identity<double>(i);
        Vector<double> b(i,1.0);
        Vector<double> x(i);
        QRP<double> qr(A);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 1.0));
        b = Vector<double>(i,0.0);
        qr.solve(x, b);
        BOOST_CHECK_EQUAL(A*x, Vector<double>(i, 0.0));
    }
}



