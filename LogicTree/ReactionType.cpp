#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "mutation++.h"

using namespace std;


int main(int argc, char** argv)

{ 
  if (argc <= 1) exit(1);

  cout << "-------------------------------------------" << endl;
  cout << "Automatic classification of reaction types." << endl;
  cout << "-------------------------------------------" << endl;
  cout << "Aurelie Bellemans. September 2013." << endl;
  cout << endl;
  
  std:: string mix(argv[argc-1]);
  cout << "The considered mixture is:  " << mix << endl; 

  Mutation::Mixture mixture(mix);

  const int ne = mixture.nElements();
  const int ns = mixture.nSpecies();
  const int nr = mixture.nReactions();

  cout << endl;
  cout << "Number of species:  " << ns << endl;
  cout << "Number of reactions:  " << nr << endl;
  cout << "Number of elements:  " << ne << endl;
  cout << endl;



//Implementation of the Logic Tree

  cout << "Reaction type:" << endl;
  cout << "--------------" << endl;
  cout << endl;

 
  cout << setw(7) << "# " ;
  cout.setf(std::ios::left, std::ios::adjustfield);
  cout << setw(15) << "Formula";
  cout << setw(20) << "Reaction type";
  cout << endl;  





 for (int i = 0; i < nr; i++) {
   const Mutation::Kinetics::Reaction& r = mixture.reactions()[i];
  
   if (r.ionReactants(mixture)==1){
     if (r.ionProducts(mixture)==1) {
     cout << setw(5) << i+1 << setw(10) << r.formula() << " :Charge exchange" << endl;
     }
     else {
       if (r.isThirdbody()) {
       cout << setw(5) << i+1 << setw(10) << r.formula() << " :Ion recombination" << endl;
       }
       else {
         if (r.electronProducts(mixture)==1){
         cout << setw(5) << i+1 << setw(10) << r.formula() << " :Electronic detachment" << endl;
         }
         else {
         cout << setw(5) << i+1 << setw(10) << r.formula() << ": Dissociative recombination" << endl;
         }
       }
     }
   }

   else {
     if (r.isThirdbody()==1){
       if (r.ionProducts(mixture)==1){
       cout << setw(5) << i+1 << setw(10)<< r.formula() << " :Ionization" << endl;
       }
       else {
         if (r.nReactants() > r.nProducts()){
         cout << setw(5) << i+1 << setw(10) << r.formula() << " :Recombination" << endl;
         }
         else {
         cout << setw(5) << i+1 << setw(10) << r.formula() << " :Dissociation" << endl;
         }
       }

     }


     else{
       if (r.electronProducts(mixture)==1){
       cout << setw(5) << i+1 << setw(10) << r.formula() << " :Associative ionization" << endl;
       }
       else if (r.electronReactants(mixture)==1){
       cout << setw(5) << i+1 << setw(10) << r.formula() << " :Electronic attachment" << endl;
       }
       else {
       cout << setw(5) << i+1 << setw(10) << r.formula() << " :Exchange" << endl;
       }
     }   
   }
 }
  cin.get();

}

