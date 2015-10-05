
#include "mutation++.h"
#include "CollisionDBNew.h"
using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport;

#include <Eigen/Dense>
using namespace Eigen;

#include <iostream>
using namespace std;

void printIntegrals(CollisionDBNew& db, string name) {
    cout << name << " = \n" << db.group(name).array() << endl;
}

int main(int argc, char* argv[])
{
    // Create the mixture
    Mixture mix(argv[1]);

    CollisionDBNew db("collisions", mix);

    for (int i = 0; i < 100; ++i) {
        mix.equilibrate(i*100.0+1000.0, ONEATM);
        cout << setw(10) << mix.T();
        const ArrayXd& v = db.group(argv[2]).array();
        for (int j = 0; j < v.size(); ++j)
            cout << setw(15) << v[j]*1.0e20;
        cout << endl;
    }

    return 0;
}
