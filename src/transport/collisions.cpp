
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
    const CollisionGroup& Q = db.group(name);
    cout << name << " = \n" << Map<const ArrayXd>(&Q[0], Q.size()) << endl;
}

int main(int argc, char* argv[])
{
    // Create the mixture
    Mixture mix(argv[1]);

    CollisionDBNew db("collisions", mix);

    mix.equilibrate(10000.0, ONEATM);
    printIntegrals(db, argv[2]);

    return 0;
}
