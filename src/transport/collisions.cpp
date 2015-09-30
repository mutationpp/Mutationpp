
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

    for (int i = 0; i < db.size(); ++i) {
    	db[i].Q11()->getOtherParams(mix);
        cout << "Q11(" << db[i].species1() << ", "
             << db[i].species2() << ")(" << mix.T() << ") = "
             << db[i].Q11()->compute(mix.T()) << endl;

        db[i].Q22()->getOtherParams(mix);
        cout << "Q22(" << db[i].species1() << ", "
			 << db[i].species2() << ")(" << mix.T() << ") = "
			 << db[i].Q22()->compute(mix.T()) << endl;

        db[i].Bst()->getOtherParams(mix);
		cout << "Bst(" << db[i].species1() << ", "
					 << db[i].species2() << ")(" << mix.T() << ") = "
					 << db[i].Bst()->compute(mix.T()) << endl;

		db[i].Cst()->getOtherParams(mix);
		cout << "Cst(" << db[i].species1() << ", "
			 << db[i].species2() << ")(" << mix.T() << ") = "
			 << db[i].Cst()->compute(mix.T()) << endl;
    }

    return 0;
}
