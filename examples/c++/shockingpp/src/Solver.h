#ifndef SOLVER_H
#define SOLVER_H

#include <boost/numeric/odeint.hpp>
#include <string>
#include "mutation++.h"

#include "SetupProperties.h"
#include "SetupProblem.h"

#include "Problem.h"
#include "Data.h"
#include "Mesh.h"
#include "BasicOutput.h"

using namespace boost::numeric::odeint;

class Solver {
public:
    Solver(size_t l_n_args, std::string l_file_name);
    ~Solver();

	void initialize();
	void solve();

private:
	size_t m_n_args;
	std::string s_file_name;
    SetupProperties* p_setup_props;
    SetupProblem* p_setup_problem;

    Mutation::Mixture* p_mix;
    Problem* p_problem;
    Mesh* p_mesh;
    ShockRelations* p_shock_relations;
    Data* p_data_pre;
    Data* p_data_post;
    BasicOutput* p_output;

    inline void setMutationpp();
    inline void errorInsufficientNumberOfArguments();

};

#endif // SOLVER_H

//=============================================================================================================

//    int ns = mix->nSpecies();
//    double Yi[ns];
//    double Xi[ns];
//    double rhoi[ns];
//    for (int i=0; i<ns; i++) {
//        Yi[i] = 0.0;
//        Xi[i] = 0.0;
//    }
//    double T = 0.0;
//    double rho  = 0.0;
//    double p    = 0.0;
//    double Us = 0.0;
//    bool pressure = false;
//    bool density  = false;
//    bool temperature = false;
//    bool massFraction = false;
//    bool velocity= false;
//    double timeOptions[5];//timeOptions = {t0,tend,dt,dtsave,dtprint}
//    double meshOptions[5];//meshOptions = {x0,xend,dx,dxsave,dxprint}
//    
//    while(getline(in, line)){
//        if(line[0]=='-')continue;
//        if (line.compare("Pressure [Pa] (Pre-shock Conditions) :")==0) {
//            getline(in,line);
//            stringstream ss(line);
//            ss>>p;
//            pressure = true;
//        }
//        if (line.compare("Density [kg/m3] (Pre-shock Conditions) :")==0) {
//            getline(in,line);
//            stringstream ss(line);
//            ss>>rho;
//            density = true;
//        }
//        if (line.compare("Temperature [K] (Pre-shock Conditions) :")==0) {
//            getline(in,line);
//            stringstream ss(line);
//            ss>>T;
//            temperature = true;
//        }
//        if (line.compare("Shock speed [m/s]:")==0) {
//            getline(in,line);
//            stringstream ss(line);
//            ss>>Us;
//            velocity = true;
//        }
//        if (line.compare("Composition (Mass fraction):")==0) {
//            getline(in,line);
//            int iCol = -1;
//            int indexSpecies[ns];
//            int i = 0;
//            while (line[iCol+1]!='*'){
//                iCol++;
//                char letter  = line[iCol];
//                std::string SpeciesName(1, letter);
//                iCol++;
//                
//                while (line[iCol]!='/') {
//                    SpeciesName+=line[iCol];
//                    iCol++;
//                }
//                indexSpecies[i] = mix->speciesIndex(SpeciesName);
//                i++;
//            }
//            
//            getline(in,line);
//            iCol = -1;
//            i=0;
//            while (line[iCol+1]!='*'){
//                iCol++;
//                char letter  = line[iCol];
//                std::string Massfrac(1, letter);
//                iCol++;
//                
//                while (line[iCol]!='/') {
//                    Massfrac+=line[iCol];
//                    iCol++;
//                }
//                double value = 0.0;
//                stringstream ss(Massfrac);
//                ss >> value;
//                Yi[indexSpecies[i]] = value;
//                i++;
//            }
//            massFraction = true;
//        }
//        if (line.compare("Composition (Mole fraction):")==0) {
//            getline(in,line);
//            int iCol = -1;
//            int indexSpecies[ns];
//            int i = 0;
//            while (line[iCol+1]!='*'){
//                iCol++;
//                char letter  = line[iCol];
//                std::string SpeciesName(1, letter);
//                iCol++;
//                
//                while (line[iCol]!='/') {
//                    SpeciesName+=line[iCol];
//                    iCol++;
//                }
//                indexSpecies[i] = mix->speciesIndex(SpeciesName);
//                i++;
//            }
//            
//            getline(in,line);
//            iCol = -1;
//            i=0;
//            while (line[iCol+1]!='*'){
//                iCol++;
//                char letter  = line[iCol];
//                std::string Massfrac(1, letter);
//                iCol++;
//                
//                while (line[iCol]!='/') {
//                    Massfrac+=line[iCol];
//                    iCol++;
//                }
//                double value = 0.0;
//                stringstream ss(Massfrac);
//                ss >> value;
//                Xi[indexSpecies[i]] = value;
//                i++;
//            }
//        }
//    }
//    
//    double Mm = 0.0;
//    if (~massFraction) {
//        Mm = mix->mixtureMwMole(Xi);
//        for (int i = 0; i<ns; i++) {
//            Yi[i] = (mix->speciesMw(i)/Mm)*Xi[i];
//        }
//    }
//    if (massFraction) {
//        Mm = mix->mixtureMwMass(Yi);
//        for (int i = 0; i<ns; i++) {
//            Xi[i] = (Mm/mix->speciesMw(i))*Yi[i];
//        }
//    }
//    
//    Mm = mix->mixtureMwMass(Yi);
//    double TP [2];
//    if (temperature && density && ~pressure) {
//        TP[0] = T;
//        TP[1] = mix->pressure(TP[0], rho, Yi);
//    }
//    else if(temperature &&pressure && ~density){
//        TP[0] = T;
//        TP[1] = p;
//        rho = p*Mm/(RU*TP[0]);
//    }
//    else if(pressure && density && ~temperature){
//        TP[1] = p;
//        TP[0] = p*Mm/(rho*RU);
//        
//    }
//    else if(pressure &&temperature && density){
//        cout<< "The initial condition are not consistent (P,T, and rho are given)"<<endl;
//        exit(0);
//    }
//    else{
//        cout << "Not enough initial data , the program will terminated"<<endl;
//        exit(0);
//    }
//    if (stateModel.compare("EquilTP")==0) {
//        cout<< "For now only ChemNonEq1T state Model works for RH-LTE"<<endl;
//        exit(0);
//    }
//    if (stateModel.compare("ChemNonEq1T")==0) {
//        for (int i=0; i<ns; i++) {
//            rhoi[i] = rho*Yi[i];
//            
//        }
//        mix->setState(rhoi, &T,1);
//    }
//    else{
//        cout <<"Other state model not yet implemented, this version works only with ChemNonEq1T"<<endl;
//    }
//    
//    cout<<"-----------------------------------------------------------------------------------------------------------"<<endl;
//    cout<<"Initial Condition : "<<endl;
//    cout << setw(15) << "T [K]";
//    cout << setw(15) << "P [Pa]";
//    cout << setw(15) << "rho [kg/m^3]";
//    for (int i = 0; i<ns; i++) {
//        string Yname = "Y_";
//        Yname+=mix->speciesName(i);
//        cout << setw(15) << Yname;
//
//    }
//    cout << endl;
//    cout << setw(15) << TP[0];
//    cout << setw(15) << TP[1];
//    cout << setw(15) << rho;
//    for (int i = 0; i<ns; i++) {
//        cout << setw(15) << Yi[i];
//        
//    }
//    cout << endl;
//    cout<<"-----------------------------------------------------------------------------------------------------------"<<endl;
//    
//
//    //Cold Gas Approximation
//    double T1 = TP[0];
//    double P1 = TP[1];
//    double U1 = Us;
//    
//
//   
//    
//    
////    double Q[3];
////    double* TP_test = new double[2];
////
////    Q[0] = 0.000918567871246755;
////    Q[1] = 5.62285105782305;
////    Q[2] = 400446.041166623;
////    
////    TP_test[0] = 153382.333103596;
////    TP_test[1] = 162415.970331782;
////    
////    Xi[0] = 0.499998097118838;
////    Xi[1] = 3.04034003326572e-06;
////    Xi[2] = 0.394998462936083;
////    Xi[3] = 7.65422288279804e-07;
////    Xi[4] = 0.104999634182754;
////    Xi[5] = 2.9444665323132e-22;
////    Xi[6] = 5.65762970919407e-22;
////    Xi[7] = 2.87251595483064e-15;
////    Xi[8] = 9.74339949122833e-23;
////    Xi[9] = 2.57561698350342e-16;
////    Xi[10] = 1.08030889977981e-15;
////
////    TCEqNewton(mix,Q,TP_test,Xi);
////    exit(0);
////
//    vector_type sol(3.0,0.0);
//    RHLTE(U1,P1,T1,rho,Xi,mix,sol);
//
//    
//
//    
//    return 0;
//}

//=============================================================================================================
