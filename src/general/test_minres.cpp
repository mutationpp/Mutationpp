
#include "Numerics.h"
#include "minres.h"

using namespace Mutation::Numerics;

#include <iostream>

using namespace std;

int main()
{
    RealSymMat A(4, 0.0);
    RealVector b(4, 1.0), x(4);
    
    A(0,0) = 1.0;
    A(1,1) = 2.0;
    A(2,2) = 0.0;
    A(3,3) = 3.0;
    A(1,3) = 1.0;
    b(3) = 0.0;
    
    cout << "Initialization:" << endl;
    cout << "A:" << endl;
    cout << A << endl;
    cout << "b:" << endl;
    cout << b << endl;
    cout << "x:" << endl;
    cout << x << endl;
    
    std::pair<int, double> conv = minres(A, x, b);
    cout << "Solution to Ax = b in " << conv.first << " iterations:" << endl;
    cout << "x:" << endl;
    cout << x << endl;
    cout << "relative residual: " << conv.second << endl;
}


