#include "Numerics.h"
#include "LDLT.h"
#include "cg.h"

using namespace std;
using namespace Numerics;

void vec_mat_mult()
{
    RealMatrix m = Lehmer<double>(10,5);
    cout << m << endl;
    
    RealVector v(10,1.0);
    cout << v << endl;
    
    cout << v * m << endl;
}

void qr_test()
{
    RealMatrix m = Lehmer<double>(5,5);
    cout << m << endl;
    
    QRP<double> qr(m);
    cout << qr.R() << endl;
    
    RealVector x(5);
    RealVector b(5,1);
    qr.solve(x, b);
    cout << x << endl;
    cout << m*x << endl;
    
    m = MinIJ<double>(10,5);
    cout << m << endl;
    
    qr = QRP<double>(m);
    cout << qr.R() << endl;    
}

void svd_test()
{
    RealMatrix m = Lehmer<double>(10,5);
    cout << m << endl;
    
    SVD<double> svd(m);
    cout << svd.U() << endl;
    cout << svd.V() << endl;
    cout << svd.U().transpose() * m * svd.V() << endl;
}


void ls_test()
{
    RealMatrix m = Lehmer<double>(10,5);
    m.col(1) = m.col(4);
    cout << m << endl;
    
    RealVector x(m.cols());
    RealVector b1(m.rows(), 1.0);
    RealVector b2(m.cols(), 1.0);
    
    LeastSquares<double> ls(m);
    
    ls.solve(x, b1);
    cout << x << endl;
    
    ls.solveATA(x, b2);
    cout << x << endl;
}

void symmetric_test()
{
    RealSymMat A = Lehmer<double>(5,5);
    cout << A.size() << endl;
    cout << A << endl;
    
    SymmetricMatrix<double> B(5, 2.0);
    cout << (A += B) << endl;
    
    double * p = new double [5];
    for (int i = 0; i < 5; ++i)
        p[i] = i;
    
    cout << A * asVector(p, 5) << endl;
}

void ldlt_test()
{
    RealSymMat A = Lehmer<double>(1000,1000);
    
    //cout << "A:" << endl;
    //cout << A << endl;
    
    RealVector x(A.cols());
    RealVector b(A.rows(), 1.0);
    
    //cout << "b:" << endl;
    //cout << b << endl;
    
    LDLT<double> ldlt(A);
    ldlt.solve(x, b);
    
    cout << "x (LDLT):" << endl;
    cout << x << endl;
    
    x = 0.0;    
    cg(A, x, b, IdentityPreconditioner<double>(), A.cols());
    
    cout << "x (CG):" << endl;
    cout << x << endl;
}

int main()
{
    ldlt_test();
}
