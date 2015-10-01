
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

    mix.equilibrate(12000.0, ONEATM);

    const CollisionGroup& Q11 = db.Q11ij();
    for (int i = 0; i < Q11.size(); ++i)
        cout << Q11[i] << endl;

    return 0;
}
