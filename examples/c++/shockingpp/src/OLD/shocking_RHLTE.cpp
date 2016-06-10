//
//  reaction0D.cpp
//
//  Created by Pierre Schrooyen on 7/11/13.
//
//  1D Euler equations with different species with mutation++
/*  The test program integrate the mass conservation law \partial_t \rho_i = w_i, the momemtum equations and the energy equations. A mathematical trick is used to transform the system of equation to an ODE. The Stationnary solution is given. The space integration uses ODEINT which provides different solvers (here a rosenbrock of order 4 is used) that can be changed in the main. The mesh size has to be provided. For an implicit solver the jacobian has to be provided, by default a finite difference jacobian is built.
 
 */


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <utility>
#include <string.h>

#include <boost/numeric/odeint.hpp>


#include "mutation++.h"

using Mutation::RU;
using namespace std;

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
typedef boost::numeric::ublas::permutation_matrix<size_t> pmatrix;

/*
    TPX_equil
    Pre : compute the Residual to find (Temperature and pressure ensuring equilibrium and knowing rho)
          - gasMix = Mutation++ mixture
          - rho    = current density
          - etot   = current internal energy
          - TP     = array of Temperature and Pressure
    Post:
          - R  = residual (energy conservation equations) (overwritten)
          - Xi = Mole fraction at equilibrium (overwritten)
          - Mm = Molar mass (overwritten)
    Rem : both density and energy should be solved but in order to make it converge, the molar mass is lagged to find a guess of the pressure at each iteration. Only the energy equation is solved by the Newton-Raphson algorithm.
*/
void TPX_equil(Mutation::Mixture* gasMix,const double rho,const double etot,const double* TP,double*  Xi, double& Mm,double& R){
    double T = TP[0];
    double p = rho*RU*T/Mm;
    double Yi[gasMix->nSpecies()];
    gasMix->equilibriumComposition(T,p,Xi);
    gasMix->convert<Mutation::Thermodynamics::X_TO_Y>(Xi, Yi);// converts X_i to Y_i and store it under Xi
    double rhoi[gasMix->nSpecies()];
    for (int iSpecies=0; iSpecies<gasMix->nSpecies(); iSpecies++) {
        rhoi[iSpecies] = rho * Yi[iSpecies];
    }
    gasMix->setState(rhoi,&T,1);
    Mm = gasMix->mixtureMw();
    double e  = gasMix->mixtureEnergyMass();
    R = etot - e;
}

/*
    TCEqNewton
    Pre : compute the Temperature, Pressure and composition at equilibrium
        - gasMix = Mutation++ mixture
        - Qtemp  = current conservative variables
        - TP     = array of Temperature and Pressure (Guess)
 
    Post:
        - Xi     = Mole fraction at equilibrium
        - TP     = array of Temperature and Pressure (Output), The guess is overwritten.
 */
void TCEqNewton(Mutation::Mixture* gasMix,const double Q_temp[], double* TP,double*  Xi){
    
    cout<<std::setprecision(15)<<Q_temp[0]<<" "<<Q_temp[1]<<" "<<Q_temp[2]<<" "<<TP[0]<<" "<<TP[1]<<endl;
//    for (int iSpecies=0; iSpecies<gasMix->nSpecies(); iSpecies++) {
//        cout<<Xi[iSpecies]<<endl;
//    }
    
    double rho  = Q_temp[0];
    double U = Q_temp[1]/Q_temp[0];
    double Ekin = 0.5 * U *U;
    double etot    = Q_temp[2]/Q_temp[0] - Ekin;
    const int ns = gasMix->nSpecies();
    double error =1.0;
    int niter = 0;
    
    double Residual = 0.0;
    double J=0.0;
    double TP_temp[2];
    TP_temp[0] =TP[0];
    TP_temp[1] =TP[1];
    double T_temp = TP_temp[0];
    double Yi[ns];
    gasMix->convert<Mutation::Thermodynamics::X_TO_Y>(Xi, Yi);
    double rhoi[ns];
    for (int iSpecies=0; iSpecies<gasMix->nSpecies(); iSpecies++) {
        rhoi[iSpecies] = rho * Yi[iSpecies];
    }
    gasMix->setState(rhoi,&T_temp,1);
    double Mm = gasMix->mixtureMw();
    TPX_equil(gasMix,rho,etot,TP_temp,Xi,Mm,Residual);
    double error0 = std::sqrt(Residual*Residual);
    double error_absolute = std::sqrt(Residual*Residual);
    while ((std::abs(error) > 1.0e-6 && std::abs(error_absolute>1e-5)) && niter<100) {
//        cout<<rho<<" "<<etot<<" "<<TP[0]<<" "<<TP[1]<<endl;
        TPX_equil(gasMix,rho,etot,TP,Xi,Mm,Residual);
        //--- Jacobian
        double Mm_temp = Mm;
        double Xi_temp[ns];
        for (int is=0; is<gasMix->nSpecies(); is++) {
            Xi_temp[is] = Xi[is];
        }
        double epsilon = 1e-2;
        double Rpert = 0.0;
        TP[0] *= (1.0+epsilon);
        TPX_equil(gasMix,rho,etot,TP,Xi_temp,Mm_temp,Rpert);
        TP[0] /=(1.0+epsilon);
        J = -(Rpert-Residual)/(TP[0]*epsilon);
        //--- Jacobian
        double alpha = 1.0;
        double T_new =  TP[0] + alpha*Residual/J;
        int nbunderRelaxation = 0;
        while (T_new<0.0 && nbunderRelaxation<10) {
            alpha *= 0.7;
            T_new = TP[0] + alpha*Residual/J; //under-relaxation if the initial temperature is too high
            cout<< "[TCEqNewton] under-relaxation : alpha = "<<alpha<<endl;
            nbunderRelaxation++;
        }
        TP[0] = T_new;
        if (TP[0]<0) {
            cout<<"ERROR ---[TCEqNewton]--- The data are still negative after 10 iterations  : Temperature = "<<TP[0]<< " and pressure = "<<TP[1]<<endl;
            exit(0);
        }
        TPX_equil(gasMix,rho,etot,TP,Xi,Mm,Residual);
        error_absolute = (std::sqrt(Residual*Residual));
        error = error_absolute/error0;
        cout<< "      Newton-Rahson to find the temperature ---- iteration : "<<niter<<" and the relative error is "<<error<<" and the absolute error is "<< (std::sqrt(Residual*Residual))<<endl;
        niter +=1;
    }
    
    TP[1] = TP[0]*RU*rho/Mm;
    
    if (niter>99){
        cout<< "[TCEqNewton] Newton-Rahson to find the temperature did not converge"<<endl;
        exit(0);
    }
}

/*
 systemToSolve
    Pre : compute the residual to get the jump relation at equilibrium
        - gasMix = Mutation++ mixture
        - Conv1  = conservative flux before the shock
        - Q_temp = current conservative variable
        - Xi     = mole fraction
    Post:
        - TP = array of Temperature and Pressure (Output), The guess is overwritten.
        - F  = residual (overwritten)
 */

void systemToSolve(Mutation::Mixture* gasMix,const double Conv1 [],const double Q_temp[], double* Xi, double* TP,vector_type& F){
    double RHO_temp  = Q_temp[0];
    double RHOU_temp = Q_temp[1];
    double RHOE_temp = Q_temp[2];
    double U_temp    = RHOU_temp/RHO_temp;
    //Compute the Temperature, Pressure and Mole/Mass fraction at equilibrium
    TCEqNewton(gasMix,Q_temp,TP,Xi);
    double P_temp = TP[1];
    
    F[0] = (Conv1[0] -  RHOU_temp);
    F[1] = (Conv1[1] -  RHOU_temp * U_temp - P_temp);
    F[2] = (Conv1[2] - (RHOE_temp + P_temp) * U_temp);
}


/*
 RHLTE
    Pre : solve the jump relation for a mixture assuming local thermal equilibrium
        - Us   = sound speed
        - P1   = pressure pre-shock
        - T1   = temperature pre-shock
        - RHO1 = density pre-shock
        - Xi_initial = mole fraction
        - gasMix = Mutation++ mixture
 Post:
    - sol  = solution
 */

void RHLTE(double Us, double P1, double T1,double RHO1,double Xi_initial[] ,Mutation::Mixture* gasMix, vector_type& sol){
    int ns = gasMix->nSpecies();
    double U1 = Us;
    double TP_CG[2];
    TP_CG[0] = T1;
    TP_CG[1] = P1;
    //Set State
    double Yi_initial[ns];
    double rho_i[ns];
    gasMix->convert<Mutation::Thermodynamics::X_TO_Y>(Xi_initial, Yi_initial);
    for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
        rho_i[iSpecies] = RHO1*Yi_initial[iSpecies];
    }
    gasMix->setState(rho_i,&T1,1);
    cout<<"Compute jump conditions..."<<endl;
    cout<<"---------------------------------------------------------------------------"<<endl;
    cout <<" Pre-shock conditions (given)"<<endl;
    cout <<"P1 : "<<setw(15)<< P1<<" [Pa]"<<endl;
    cout <<"T1 : "<<setw(15)<< T1<<" [K]"<<endl;
    cout <<"US : "<<setw(15)<< U1<<" [m/s]"<<endl;
    cout<<"---------------------------------------------------------------------------"<<endl;
    double Mm1 = gasMix->mixtureMw();
    double Rstar = RU/Mm1;
    double Gamma = gasMix->mixtureFrozenCpMass()/gasMix->mixtureFrozenCvMass();
    
    double GM1 = Gamma-1.0;
    double GP1 = Gamma+1.0;
    
    double C1 = std::sqrt(Gamma*Rstar*T1);
    double Ms = U1/C1;
    
    double P2 = P1 *(1.0+2.0*Gamma/GP1*(Ms*Ms-1.0));
    double T2 = T1 *(1.0+((2.0*GM1/(GP1*GP1))*((Gamma*Ms*Ms+1.0)/(Ms*Ms))*(Ms*Ms-1.0)));
    double U2 = C1 *2.0 /GP1 *(Ms - 1.0 /Ms);
    double RHO2 = P2/(Rstar*T2);
    
    cout <<" Post-shock conditions (cold gas approximation)"<<endl;
    cout <<"P2 :    "<<setw(15)<< P2<<" [Pa]"<<endl;
    cout <<"T2 :    "<<setw(15)<< T2<<" [K]"<<endl;
    cout <<"US-U2 : "<<setw(15)<< Us-U2<<" [m/s]"<<endl;
    cout <<"RHO2  : "<<setw(15)<< RHO2<<" [kg/m^3]"<<endl;
    cout<<"---------------------------------------------------------------------------"<<endl;
    
    //Hot gas approximation
    cout <<" Compute post-shock conditions using LTE ..."<<endl;
    double rho1  = P1/(Rstar*T1);
    double rhoU1 = rho1*U1;
    gasMix->setState(rho_i,&T1,1);
    double h1    = gasMix->mixtureHMass();
    double Ekin1 = 0.5*U1*U1;
    double e1 = gasMix->mixtureEnergyMass();
    double rhoE1 = rho1 * (e1 + Ekin1);
    
    int niter = 0;
    int nbitermax = 100;
    vector_type F (3,0.0);
    matrix_type J(3,3);
    double error = 1.0;

    //Initial guess :
    TP_CG[0] = T2;
    TP_CG[1] = P2;
    for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
        rho_i[iSpecies] = RHO2*Yi_initial[iSpecies];
    }

    gasMix->setState(rho_i,&T2,1);
    double e2    = gasMix->mixtureEnergyMass();
    double Conv1 [3];
    Conv1[0] = rhoU1;
    Conv1[1] = rhoU1*U1 + P1;
    Conv1[2] = (rhoE1 + P1) * U1;
    double TP_pert[2];
    double Q_temp[3];
    Q_temp[0] = RHO2;
    Q_temp[1] = RHO2 * (Us-U2);
    Q_temp[2] = RHO2 * e2 + 0.5*RHO2 * (Us-U2)*(Us-U2);
    
    systemToSolve(gasMix,Conv1, Q_temp,Xi_initial,TP_CG, F);
    double error0 = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);

    while (std::abs(error)>1.0e-8 && niter<nbitermax) {
        cout<<"Begin iteration "<<niter<<" ..."<<endl;
        systemToSolve(gasMix,Conv1, Q_temp,Xi_initial,TP_CG, F);
        cout<<"      ------------------------------------BEGIN Compute Jacobian ------------------------------------ "<<endl;
        TP_pert[0] = TP_CG[0];
        TP_pert[1] = TP_CG[1];
        
        for (int i = 0; i<3; i++) {
            double epsilon = 1e-3;
            vector_type Fpert (3,0.0);
            Q_temp[i] *= (1.+epsilon);
            systemToSolve(gasMix,Conv1, Q_temp,Xi_initial,TP_pert, Fpert);
            Q_temp[i] /=(1.+epsilon);
            for (int j = 0; j<3; j++) {
                J(j,i) = -(Fpert(j)-F(j))/(Q_temp[i]*epsilon);
            }
        }
        cout<<"      ------------------------------------ END Compute Jacobian ------------------------------------ "<<endl;
        
        // Solution of the system
//        cout<<"System to solve"<<endl;
//        cout<<std::setprecision(15)<<J(0,0)<<" "<<J(0,1)<<" "<<J(0,2)<<endl;
//        cout<<std::setprecision(15)<<J(1,0)<<" "<<J(1,1)<<" "<<J(1,2)<<endl;
//        cout<<std::setprecision(15)<<J(2,0)<<" "<<J(2,1)<<" "<<J(2,2)<<endl;
//        cout<<std::setprecision(15)<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
        pmatrix piv(3);
        lu_factorize(J, piv);
        lu_substitute(J, piv, F);
        double alpha = 1.0;
        // Compute Error
//        cout<<"Solution"<<endl;
//        cout<<std::setprecision(15)<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
        Q_temp[0] += alpha*F[0];
        Q_temp[1] += alpha*F[1];
        Q_temp[2] += alpha*F[2];
        systemToSolve(gasMix,Conv1, Q_temp,Xi_initial,TP_CG, F);
        error = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2])/error0;
        
        cout<<"--------------------------------------------------------------------------------------------------"<<endl;
        cout<< "Newton-Rhapson to compute post-shock condition ----- iteration : "<<niter<<" , residual : "<<error<<endl;
        cout<<"--------------------------------------------------------------------------------------------------"<<endl;
        niter++;
        
    }
    if (niter>=(nbitermax-1)) {
        cout<<"Newton-Raphson for the LTE condition did not converge, the program will be ended"<<endl;
        exit(0);
    }
    else{
        cout<<"Newton-Raphson converged in "<< niter<< " iterations"<<endl;
    }
    
    //Final Results
    T2 = TP_CG[0];
    P2 = TP_CG[1];
    double USmU2 = Q_temp[1]/Q_temp[0];
    sol[0] = Q_temp[0];
    sol[1] = Q_temp[1];
    sol[2] = Q_temp[2];
    cout<<"---------------------------------------------------------------------------"<<endl;
    cout <<" Post-shock conditions (hot gas approximation)"<<endl;
    cout <<"P2 :    "<<setw(15)<< P2<<" [Pa]"<<endl;
    cout <<"T2 :    "<<setw(15)<< T2<<" [K]"<<endl;
    cout <<"US-U2 : "<<setw(15)<< USmU2<<" [m/s]"<<endl;
    cout<<"---------------------------------------------------------------------------"<<endl;
    cout<<" Composition"<<endl;
    for (int i = 0; i<ns; i++) {
        string Yname = "Y_";
        Yname+=gasMix->speciesName(i);
        cout << setw(25) << Yname;
        
    }
    cout << endl;
    double Xi[gasMix->nSpecies()];
    double Yi[gasMix->nSpecies()];
    gasMix->equilibriumComposition(T2,P2,Xi);
    gasMix->convert<Mutation::Thermodynamics::X_TO_Y>(Xi, Yi);
    for (int i = 0; i<ns; i++) {
        cout << setw(25) << Yi[i];
    }
    cout<<endl;
    
}


// -------------------------------------------------------------------------------
// Main
// -------------------------------------------------------------------------------
int main( int argc , char **argv ){
    
    //Initialisation (read input file)
    cout<< "Loading Input file ..."<<endl;
    const char* filename;
    if (argv[1]==NULL) {
        cout <<"WARNING : No input file has been provided... The shockair.in input file will be loaded by default!"<<endl;
        filename ="shockair.in";
    }
    else{
        filename =argv[1];
    }
    
    ifstream in;
    in.open(filename);
    if (!in) {
        cout << "ERROR: couldn't open input file." << endl;
    }
    std::string mixtName;
    std::string reacMech;
    std::string stateModel;
    std::string thermoDataBase;
    string line;
    
    while(getline(in, line)){
        if(line[0]=='-')continue;
        if (line.compare("Name of the mixture:")==0) {
            getline(in,line);
            mixtName = line;
        }
        if (line.compare("Reaction mechanism:")==0) {
            getline(in,line);
            reacMech = line;
        }
        if (line.compare("State Model:")==0) {
            getline(in,line);
            stateModel = line;
        }
        if (line.compare("Thermodynamic Database:")==0) {
            getline(in,line);
            thermoDataBase = line;
        }
        if(line.compare("End Mixture Options")==0){
            break;
        }
    }
    cout<<"Creating mixture..."<<endl;
    Mutation::MixtureOptions opts(mixtName.c_str());
    opts.setStateModel(stateModel);//stateModel);
    opts.setThermodynamicDatabase(thermoDataBase);
    Mutation::Mixture* mix = new Mutation::Mixture(opts);

    int ns = mix->nSpecies();
    double Yi[ns];
    double Xi[ns];
    double rhoi[ns];
    for (int i=0; i<ns; i++) {
        Yi[i] = 0.0;
        Xi[i] = 0.0;
    }
    double T = 0.0;
    double rho  = 0.0;
    double p    = 0.0;
    double Us = 0.0;
    bool pressure = false;
    bool density  = false;
    bool temperature = false;
    bool massFraction = false;
    bool velocity= false;
    double timeOptions[5];//timeOptions = {t0,tend,dt,dtsave,dtprint}
    double meshOptions[5];//meshOptions = {x0,xend,dx,dxsave,dxprint}
    
    while(getline(in, line)){
        if(line[0]=='-')continue;
        if (line.compare("Pressure [Pa] (Pre-shock Conditions) :")==0) {
            getline(in,line);
            stringstream ss(line);
            ss>>p;
            pressure = true;
        }
        if (line.compare("Density [kg/m3] (Pre-shock Conditions) :")==0) {
            getline(in,line);
            stringstream ss(line);
            ss>>rho;
            density = true;
        }
        if (line.compare("Temperature [K] (Pre-shock Conditions) :")==0) {
            getline(in,line);
            stringstream ss(line);
            ss>>T;
            temperature = true;
        }
        if (line.compare("Shock speed [m/s]:")==0) {
            getline(in,line);
            stringstream ss(line);
            ss>>Us;
            velocity = true;
        }
        if (line.compare("Composition (Mass fraction):")==0) {
            getline(in,line);
            int iCol = -1;
            int indexSpecies[ns];
            int i = 0;
            while (line[iCol+1]!='*'){
                iCol++;
                char letter  = line[iCol];
                std::string SpeciesName(1, letter);
                iCol++;
                
                while (line[iCol]!='/') {
                    SpeciesName+=line[iCol];
                    iCol++;
                }
                indexSpecies[i] = mix->speciesIndex(SpeciesName);
                i++;
            }
            
            getline(in,line);
            iCol = -1;
            i=0;
            while (line[iCol+1]!='*'){
                iCol++;
                char letter  = line[iCol];
                std::string Massfrac(1, letter);
                iCol++;
                
                while (line[iCol]!='/') {
                    Massfrac+=line[iCol];
                    iCol++;
                }
                double value = 0.0;
                stringstream ss(Massfrac);
                ss >> value;
                Yi[indexSpecies[i]] = value;
                i++;
            }
            massFraction = true;
        }
        if (line.compare("Composition (Mole fraction):")==0) {
            getline(in,line);
            int iCol = -1;
            int indexSpecies[ns];
            int i = 0;
            while (line[iCol+1]!='*'){
                iCol++;
                char letter  = line[iCol];
                std::string SpeciesName(1, letter);
                iCol++;
                
                while (line[iCol]!='/') {
                    SpeciesName+=line[iCol];
                    iCol++;
                }
                indexSpecies[i] = mix->speciesIndex(SpeciesName);
                i++;
            }
            
            getline(in,line);
            iCol = -1;
            i=0;
            while (line[iCol+1]!='*'){
                iCol++;
                char letter  = line[iCol];
                std::string Massfrac(1, letter);
                iCol++;
                
                while (line[iCol]!='/') {
                    Massfrac+=line[iCol];
                    iCol++;
                }
                double value = 0.0;
                stringstream ss(Massfrac);
                ss >> value;
                Xi[indexSpecies[i]] = value;
                i++;
            }
        }
    }
    
    double Mm = 0.0;
    if (~massFraction) {
        Mm = mix->mixtureMwMole(Xi);
        for (int i = 0; i<ns; i++) {
            Yi[i] = (mix->speciesMw(i)/Mm)*Xi[i];
        }
    }
    if (massFraction) {
        Mm = mix->mixtureMwMass(Yi);
        for (int i = 0; i<ns; i++) {
            Xi[i] = (Mm/mix->speciesMw(i))*Yi[i];
        }
    }
    
    Mm = mix->mixtureMwMass(Yi);
    double TP [2];
    if (temperature && density && ~pressure) {
        TP[0] = T;
        TP[1] = mix->pressure(TP[0], rho, Yi);
    }
    else if(temperature &&pressure && ~density){
        TP[0] = T;
        TP[1] = p;
        rho = p*Mm/(RU*TP[0]);
    }
    else if(pressure && density && ~temperature){
        TP[1] = p;
        TP[0] = p*Mm/(rho*RU);
        
    }
    else if(pressure &&temperature && density){
        cout<< "The initial condition are not consistent (P,T, and rho are given)"<<endl;
        exit(0);
    }
    else{
        cout << "Not enough initial data , the program will terminated"<<endl;
        exit(0);
    }
    if (stateModel.compare("EquilTP")==0) {
        cout<< "For now only ChemNonEq1T state Model works for RH-LTE"<<endl;
        exit(0);
    }
    if (stateModel.compare("ChemNonEq1T")==0) {
        for (int i=0; i<ns; i++) {
            rhoi[i] = rho*Yi[i];
            
        }
        mix->setState(rhoi, &T,1);
    }
    else{
        cout <<"Other state model not yet implemented, this version works only with ChemNonEq1T"<<endl;
    }
    
    cout<<"-----------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Initial Condition : "<<endl;
    cout << setw(15) << "T [K]";
    cout << setw(15) << "P [Pa]";
    cout << setw(15) << "rho [kg/m^3]";
    for (int i = 0; i<ns; i++) {
        string Yname = "Y_";
        Yname+=mix->speciesName(i);
        cout << setw(15) << Yname;

    }
    cout << endl;
    cout << setw(15) << TP[0];
    cout << setw(15) << TP[1];
    cout << setw(15) << rho;
    for (int i = 0; i<ns; i++) {
        cout << setw(15) << Yi[i];
        
    }
    cout << endl;
    cout<<"-----------------------------------------------------------------------------------------------------------"<<endl;
    

    //Cold Gas Approximation
    double T1 = TP[0];
    double P1 = TP[1];
    double U1 = Us;
    

   
    
    
//    double Q[3];
//    double* TP_test = new double[2];
//
//    Q[0] = 0.000918567871246755;
//    Q[1] = 5.62285105782305;
//    Q[2] = 400446.041166623;
//    
//    TP_test[0] = 153382.333103596;
//    TP_test[1] = 162415.970331782;
//    
//    Xi[0] = 0.499998097118838;
//    Xi[1] = 3.04034003326572e-06;
//    Xi[2] = 0.394998462936083;
//    Xi[3] = 7.65422288279804e-07;
//    Xi[4] = 0.104999634182754;
//    Xi[5] = 2.9444665323132e-22;
//    Xi[6] = 5.65762970919407e-22;
//    Xi[7] = 2.87251595483064e-15;
//    Xi[8] = 9.74339949122833e-23;
//    Xi[9] = 2.57561698350342e-16;
//    Xi[10] = 1.08030889977981e-15;
//
//    TCEqNewton(mix,Q,TP_test,Xi);
//    exit(0);
//
    vector_type sol(3.0,0.0);
    RHLTE(U1,P1,T1,rho,Xi,mix,sol);

    

    
    return 0;
}


