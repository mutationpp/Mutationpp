
#include "mutation++.h"
#include "CollisionDBNew.h"
using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport;

#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
    // Create the mixture
    Mixture mix(argv[1]);

    CollisionDBNew db("collisions", mix);

    double T = 2000.0;

    for (int i = 0; i < db.size(); ++i)
        cout << "Q11(" << db[i].species1() << ", "
             << db[i].species2() << ")(" << T << ") = "
             << db[i].Q11()->compute(T) << endl;

    return 0;
}
