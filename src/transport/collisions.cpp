
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

    mix.equilibrate(12000.0, ONEATM);

    printIntegrals(db, "Q11ee");
    printIntegrals(db, "Q11ei");
    printIntegrals(db, "Q11ii");
    printIntegrals(db, "Q11ij");

    printIntegrals(db, "Bstee");
    printIntegrals(db, "Bstei");
    printIntegrals(db, "Bstii");

    return 0;
}
