//
//  shocking.cpp
//
//  Created by Pierre Schrooyen on 7/01/14.
//
//  Shock Relaxation
/*  The test program computes the post shock parameters and integrate the stationary Euler equations behind the shock. The space integration uses ODEINT which provides different solvers (here a rosenbrock of order 4 is used) that can be changed in the main routine. The mesh size has to be provided. For an implicit solver the jacobian has to be provided, by default a finite difference jacobian is built.
 
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
using namespace boost::numeric::odeint;
using namespace Mutation;
using namespace Mutation::Thermodynamics;
// Definition of ublas vector needed for implicit integrator type
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
typedef boost::numeric::ublas::permutation_matrix<size_t> pmatrix;

// types for the lagrangian solver
typedef boost::numeric::ublas::vector< double > vec_type;
typedef boost::numeric::ublas::matrix< double > mat_type;

typedef std::vector< double > state_type;

#include "LagrangianSolver.h"

void problem(const vector_type &x , vector_type &dxdt , double t, Mutation::Mixture* mix, double mdot,double overmdot ){
    size_t ns = mix->nSpecies();
    size_t iUIndex = ns;
    size_t iTIndex = ns+1;

    double* Wdot = new double [ns];
    double* Yi = new double [ns];
    double rhoi [ns];
    double* Cp = new double [ns];
    double* H  = new double [ns];
    
    double T = x[iTIndex];
    double u = x[iUIndex];
    double rho = mdot/u;
    
    for (size_t iSpecies = 0; iSpecies<ns; iSpecies++) {
        Yi[iSpecies] = x[iSpecies];
        rhoi[iSpecies] = rho*Yi[iSpecies];
    }
    double Mm = mix->mixtureMwMass(Yi);
    
    mix->speciesCpOverR(T, Cp);
    mix->speciesHOverRT(T, H);
    
    double p = rho*RU*T/Mm;
    
    mix->setState(rhoi,&T,1);
    mix->netProductionRates(Wdot);
    
    double a = mdot*u/(RU*T) - p/(RU*T);
    double b = 0.0;
    double c = 0.0;
    double d = mdot*u;
    double e = 0.0;
    double f = 0.0;
    
    for (int iSpecies = 0; iSpecies<ns; iSpecies++){
        double Mmi = mix->speciesMw(iSpecies);
        
        b += (mdot/T)*(x[iSpecies]/Mmi);
        c -= Wdot[iSpecies]/Mmi;
        e += mdot*(RU/Mmi)*(x[iSpecies]*(Cp[iSpecies]));
        f -= Wdot[iSpecies]*H[iSpecies]*(RU/Mmi)*T;
        
        dxdt[ iSpecies ] = Wdot[iSpecies]*overmdot;
    }
    
    double overdetA = 1.0/(a*e-b*d);
    dxdt[iUIndex] = (e*c-b*f)*overdetA;
    dxdt[iTIndex] = (a*f-c*d)*overdetA;
    
    delete [] Wdot;
    delete [] Yi;
    delete [] Cp;
    delete [] H;
}

void rankineHugoniot(double Us, double P1, double T1,double RHO1,double Yi[] ,Mutation::Mixture* gasMix, vector_type& sol){
    
    double U1 = Us;
    size_t ns = gasMix->nSpecies();
    double rho_i [ns];
    for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
        rho_i[iSpecies] = RHO1*Yi[iSpecies];
    }
    gasMix->setState(rho_i,&T1,1);
    
    cout<<"Compute jump conditions..."<<endl;
    cout<<"---------------------------------------------------------------------------"<<endl;
    cout <<" Pre-shock conditions (given)"<<endl;
    cout <<"P1 : "<<setw(15)<< P1<<" [Pa]"<<endl;
    cout <<"T1 : "<<setw(15)<< T1<<" [K]"<<endl;
    cout <<"US : "<<setw(15)<< U1<<" [m/s]"<<endl;
    cout<<"---------------------------------------------------------------------------"<<endl;
    double Mm1 = gasMix->mixtureMwMass(Yi);
    double Rstar = RU/Mm1;
    double Gamma = gasMix->mixtureFrozenCpMass()/gasMix->mixtureFrozenCvMass();
    
    double GM1 = Gamma-1.0;
    double GP1 = Gamma+1.0;
    
    double C1 = std::sqrt(Gamma*Rstar*T1);
    double Ms = U1/C1;
    
    double P2 = P1 *(1.0+2.0*Gamma/GP1*(Ms*Ms-1.0));
    double T2 = T1 *(1.0+((2.0*GM1/(GP1*GP1))*((Gamma*Ms*Ms+1.0)/(Ms*Ms))*(Ms*Ms-1.0)));
    double U2 = C1 *2.0 /GP1 *(Ms - 1.0/Ms);
    double RHO2 = P2/(Rstar*T2);
    
    cout <<" Post-shock conditions (cold gas approximation)"<<endl;
    cout <<"P2 :    "<<setw(15)<< P2<<" [Pa]"<<endl;
    cout <<"T2 :    "<<setw(15)<< T2<<" [K]"<<endl;
    cout <<"US-U2 : "<<setw(15)<< Us-U2<<" [m/s]"<<endl;
    cout <<"RHO2  : "<<setw(15)<< RHO2<<" [kg/m^3]"<<endl;

    cout<<"---------------------------------------------------------------------------"<<endl;
    
    //Hot gas approximation
    double rho1  = RHO1;
    double rhou1 = RHO1*U1;
    double h1    = gasMix->mixtureHMass();
    double Ekin1 = 0.5*U1*U1;
    int niter = 0;
    double TP_temp[2];
    
    vector_type F (3,0.0);
    matrix_type J(3,3);
    double error = 1.0;
    vector_type xi(3,0.0);
    
    //Initial guess :
    xi[0] = Us-U2;
    xi[1] = P2;
    xi[2] = T2;
    double Ui = xi[0];
    double Pi = xi[1];
    double Ti = xi[2];
    double hi = 0.0;
    hi = gasMix->mixtureHMass();
    for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
        rho_i[iSpecies] = Yi[iSpecies]*Pi*Mm1/(RU*Ti);
    }
    gasMix->setState(rho_i,&Ti,1);
    F[0] = -(rhou1 - Pi/(Rstar*Ti)*Ui);
    F[1] = -(P1 - Pi - Pi/(Rstar*Ti)*(Ui*Ui)+rhou1*U1);
    F[2] = -(h1 + Ekin1 - hi - 0.5*Ui*Ui);
    double Resini = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
    
    while (std::abs(error)>1.0e-12 && niter<100) {
        Ui = xi[0];
        Pi = xi[1];
        Ti = xi[2];
        TP_temp[0] = Ti;
        TP_temp[1] = Pi;
        double rho = Pi*Mm1/(RU*Ti);
        for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
            rho_i[iSpecies] = Yi[iSpecies]*rho;
        }
        gasMix->setState(rho_i,&Ti,1);
        hi = gasMix->mixtureHMass();
        
        F[0] = -(rhou1 - Pi/(Rstar*Ti)*Ui);
        F[1] = -(P1 - Pi - Pi/(Rstar*Ti)*(Ui*Ui)+rhou1*U1);
        F[2] = -(h1 + Ekin1 - hi - 0.5*Ui*Ui);
        
        J(0,0) = -Pi/(Rstar*Ti);
        J(0,1) = -Ui/(Rstar*Ti);
        J(0,2) = (Ui*Pi)/(Rstar*Ti*Ti);
        J(1,0) = -(2.0*Ui*Pi)/(Rstar*Ti);
        J(1,1) = -1.0-(Ui*Ui)/(Rstar*Ti);
        J(1,2) = (Ui*Ui*Pi)/(Rstar*Ti*Ti);
        J(2,0) = - Ui;
        J(2,1) = 0.0;
        J(2,2) = - gasMix->mixtureFrozenCpMass();
        
        pmatrix piv(3);
        lu_factorize(J, piv);
        lu_substitute(J, piv, F);
        xi[0] += F[0];
        xi[1] += F[1];
        xi[2] += F[2];
        
        //Compute Error
        Ui = xi[0];
        Pi = xi[1];
        Ti = xi[2];
        TP_temp[0] = Ti;
        TP_temp[1] = Pi;
        rho = Pi*Mm1/(RU*Ti);
        for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
            rho_i[iSpecies] = Yi[iSpecies]*rho;
        }
        gasMix->setState(rho_i,&Ti,1);
        hi = gasMix->mixtureHMass();
        
        F[0] = -(rhou1 - Pi/(Rstar*Ti)*Ui);
        F[1] = -(P1 - Pi - Pi/(Rstar*Ti)*(Ui*Ui)+rhou1*U1);
        F[2] = -(h1 + Ekin1 - hi - 0.5*Ui*Ui);
        
        error = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2])/Resini;
        niter++;
    }
    if (niter>=99) {
        cout<<"Newton-Raphson for the post-shock value did not converge, the program will terminate :"<<error<<endl;
        exit(0);
    }
    else{
        cout<<"Newton-Raphson converged in "<< niter<< " iterations and the residual are "<<error<<endl;
    }
    sol[0] = xi[0];
    sol[1] = xi[1];
    sol[2] = xi[2];
    double c2 = gasMix->frozenSoundSpeed();
    double Mach2 =(Us-U2)/c2;
    cout <<" Post-shock conditions (hot gas approximation)"<<endl;
    cout <<"P2 :    "<<setw(15)<< sol[1]<<" [Pa]"<<endl;
    cout <<"T2 :    "<<setw(15)<< sol[2]<<" [K]"<<endl;
    cout <<"US-U2 : "<<setw(15)<< sol[0]<<" [m/s]"<<endl;
        cout<<"Mach number"<<setw(15)<< Mach2<<" [-]"<<endl;
    cout<<"---------------------------------------------------------------------------"<<endl;
}

// -------------------------------------------------------------------------------
// Problem definition (ODE to solve)
// dxdt = S(x,t)
// -------------------------------------------------------------------------------

// System
struct equationsToSolve
{
    Mutation::Mixture* mix;
    size_t ns;
    double mdot;
    double overmdot;
    equationsToSolve(Mutation::Mixture* gasMix,double mdot){
        this->mix = gasMix;
        ns = mix->nSpecies();
        this->mdot = mdot;
        this->overmdot = 1./mdot;
    }
    
    void operator()( const vector_type &x , vector_type &dxdt , double t)
    {
        problem(x , dxdt ,t,mix,mdot,overmdot );
    }
};
// -------------------------------------------------------------------------------

// -------------------------------------------------------------------------------
// Jacobian of the system to solve (implicit solver)
// -------------------------------------------------------------------------------
struct jacobian
{
    Mutation::Mixture* mix;
    size_t ns;
    size_t numberOfVariables;
    double mdot;
    double overmdot;
    jacobian(Mutation::Mixture* gasMix,double mdot){
        this->mix = gasMix;
        ns = mix->nSpecies();
        this->mdot = mdot;
        this->overmdot = 1./mdot;
        this->numberOfVariables = ns+2;
    }
    
    void operator()( const vector_type &x  , matrix_type &J , const double &t , vector_type &dfdt )
    {
        
//        //--------------------------------------------------------------------------------------
//        // Finite Difference Jacobian centered
//        //--------------------------------------------------------------------------------------
//        vector_type dxdt(numberOfVariables,0.0);
//        vector_type x_pert_u(numberOfVariables,0.0);
//        vector_type dxdt_pert_u(numberOfVariables,0.0);
//        vector_type x_pert_d(numberOfVariables,0.0);
//        vector_type dxdt_pert_d(numberOfVariables,0.0);
//        double pert = 1e-8;
//        matrix_type J_fd(numberOfVariables,numberOfVariables);
//        
//        for (int i=0; i<numberOfVariables; i++){
//            x_pert_u[i] = x[i];
//            x_pert_d[i] = x[i];
//        }
//        
//        for (int i=0; i<numberOfVariables; i++){
//            x_pert_u[i]+=pert;
//            problem(x_pert_u,dxdt_pert_u,t,mix,mdot,overmdot);
//            x_pert_d[i]-=pert;
//            problem(x_pert_d,dxdt_pert_d,t,mix,mdot,overmdot);
//            
//            for (int j=0; j<numberOfVariables; j++){
//                J_fd(j,i) = (dxdt_pert_u[j]-dxdt_pert_d[j])/(2*pert);
//                J(j,i) = J_fd(j,i);
//            }
//            x_pert_u[i]-=pert;
//            x_pert_d[i]+=pert;
//            dfdt[i] = 0.0;
//        }
        //--------------------------------------------------------------------------------------
        // Finite Difference Jacobian upwind
        //--------------------------------------------------------------------------------------
//        cout<<"---------------------- FD JACOBIAN ----------------------"<<endl;
        vector_type dxdt(numberOfVariables,0.0);
        vector_type x_pert(numberOfVariables,0.0);
        vector_type dxdt_pert(numberOfVariables,0.0);
        double pert = 1e-8;
        matrix_type J_fd(numberOfVariables,numberOfVariables);
        
        for (int i=0; i<numberOfVariables; i++){
            x_pert[i] = x[i];
        }
        problem(x,dxdt,t,mix,mdot,overmdot);
        
        for (int i=0; i<numberOfVariables; i++){
            x_pert[i]*=(1.0+pert);//*=(1.0+pert);//pert;
            problem(x_pert,dxdt_pert,t,mix,mdot,overmdot);
            x_pert[i]/=(1.0+pert);//pert;//(1.0+pert);
            for (int j=0; j<numberOfVariables; j++){
                J_fd(j,i) = (dxdt_pert[j]-dxdt[j])/((x_pert[i]*pert+1e-16));//pert(x_pert[i]*pert+1e-16);
                J(j,i) = J_fd(j,i);
            }
//            x_pert[i]-=pert;
//            dfdt[i] = 0.0;
        }
//                cout<<"---------------------- FD JACOBIAN ----------------------"<<endl;
//        cout<< "Compare Jacobian ..."<<endl;
//        cout<< "------------------------------------------------"<<endl;
//        cout<< "                   Analytical                   "<<endl;
//        cout<< "------------------------------------------------"<<endl;
//         for (int i=0; i<ns; i++){
//              for (int j=0; j<ns; j++){
//                  cout << J(i,j)<<"     ";
//              }
//             cout<<endl;
//         }
//        cout<< "------------------------------------------------"<<endl;
//        cout<< "                Finite difference               "<<endl;
//        cout<< "------------------------------------------------"<<endl;
//        for (int i=0; i<ns; i++){
//            for (int j=0; j<ns; j++){
//                cout << J_fd(i,j)<<"     ";
//            }
//            cout<<endl;
//        }
//        cout<< "------------------------------------------------"<<endl;
//        exit(0);
    }
};
// -------------------------------------------------------------------------------


// -------------------------------------------------------------------------------
// Intermediate solution (print or save)
// -------------------------------------------------------------------------------

struct intermediateSolution
{
    Mutation::Mixture* mix;
    int nsave;
    int nprint;
    double mdot;
    double dtprint, dtsave;
    intermediateSolution(double mdot,Mutation::Mixture* gasMix,double dtsave,double dtprint){
        this->mix = gasMix;
        nprint = 0;
        nsave = 0;
        this->dtsave  = dtsave;
        this->dtprint = dtprint;
        this->mdot    = mdot;
    }
    
    void operator()( const vector_type &x , const double t)
    {
        int ns = mix->nSpecies();
        size_t iUIndex = ns;
        size_t iTIndex = ns+1;
        
        double u = x[iUIndex];
        double T = x[iTIndex];
        double rho = mdot/u;

        double Yi[ns];
        for (int i=0; i<ns; i++) {
            Yi[i] = x[i];
        }

        double p = mix->pressure(T,rho,Yi);

        // SAVE
        std::ofstream myflux("output.dat", std::ios::app);
        if (t >= dtsave*nsave){
            myflux << t<<" "<<T<<" "<<p<<" "<<rho<<" "<<u<<" ";
            for (int i=0; i<ns; i++) {
                myflux<<Yi[i]<<" ";
            }
            myflux<<endl;
            myflux.close();
            this->nsave++;
        }
        
        // PRINT
        if (t >= dtprint*nprint){
            cout<<setw(15)<<t<<setw(15)<<T<<setw(15)<<p<<setw(15)<<u<<setw(15)<<rho;
            for (int i=0; i<ns; i++) {
                cout<<setw(15)<<Yi[i];
            }
            cout<<endl;
            this->nprint++;
        }
    }
};





// ###############################################################################
// Main
// ###############################################################################
int main( int argc , char **argv ){
    
    //Initialisation (read input file)
    cout<< "Loading Input file ..."<<endl;
    const char* filename;
    if (argv[1]==NULL) {
        cout <<"WARNING : No input file has been provided... The shockair example will be loaded by default!"<<endl;
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
    
    Mutation::MixtureOptions opts(mixtName);
    opts.setStateModel(stateModel);
    opts.setThermodynamicDatabase(thermoDataBase);
    
    Mutation::Mixture* mix = new Mutation::Mixture(opts);
    int ns = mix->nSpecies();
    double T = 0.0;
    double Yi[ns];
    double Xi[ns];
    for (int i=0; i<ns; i++) {
        Yi[i] = 0.0;
        Xi[i] = 0.0;
    }
    double rho  = 0.0;
    double p    = 0.0;
    double u = 0.0;
    double mdot = 0.0;
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
            ss>>u;
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
        if (line.compare("Time integration options:")==0) {
            getline(in,line);
            int index = 0;
            while (line[0]!='-') {
                getline(in,line);
                stringstream ss(line);
                ss>> timeOptions[index];
                getline(in,line);
                index++;
            }
        }
        if (line.compare("Mesh integration options:")==0) {
            getline(in,line);
            int index = 0;
            while (line[0]!='-') {
                getline(in,line);
                stringstream ss(line);
                ss>> meshOptions[index];
                getline(in,line);
                index++;
            }
        }
    }
    
    double Mm = 0.0;
    if (!massFraction) {
        Mm = mix->mixtureMwMole(Xi);
        for (int i = 0; i<ns; i++) {
            Yi[i] = (mix->speciesMw(i)/Mm)*Xi[i];
        }
    }
    
    Mm = mix->mixtureMwMass(Yi);
    double TP [2];
    if (temperature && density && !pressure) {
        TP[0] = T;
        TP[1] = mix->pressure(TP[0], rho, Yi);
    }
    else if(temperature &&pressure && !density){
        TP[0] = T;
        TP[1] = p;
        rho = p*Mm/(RU*TP[0]);
    }
    else if(pressure && density && !temperature){
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
    
    if (stateModel.compare("ChemNonEq1T")==0) {
        double* rhoi = new double[ns];
        for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
            rhoi[iSpecies] = rho*Yi[iSpecies];
        }
        mix->setState(rhoi,&T,1);
    }
    else{
        cout <<"other state model not yet implemented"<<endl;
    }
    
    if (!velocity) {
        cout << "Not enough initial data , the program will terminated"<<endl;
        exit(0);
    }
    else{
        mdot = u*mix->density();
    }
    
    
    // ------------------------------------------------------------------------------------------------------------
    // Compute jump conditions
    // ------------------------------------------------------------------------------------------------------------
    double Us = u;
    double P1 = TP[1];
    double T1 = TP[0];
    double Rho1 = rho;
    vector_type x_ps(3,0.0);
    rankineHugoniot(Us,P1,T1,Rho1,Yi,mix,x_ps);
    
    // ------------------------------------------------------------------------------------------------------------
    // Initial solution
    // ------------------------------------------------------------------------------------------------------------
    vector_type x (ns+2,0.0);
    
    double* rhoi = new double [ns];
    for (int i = 0; i<ns; i++) {
        x[i] = Yi[i];
        rhoi[i] = Yi[i]*mdot/x_ps[0];
    }
    x[ns] = x_ps[0];
    x[ns+1] = x_ps[2];
    TP[0] = x_ps[2];
    TP[1] = x_ps[1];
    mix->setState(rhoi,&T,1);
    
    
    // ------------------------------------------------------------------------------------------------------------
    double etot = mix->mixtureEnergyMass();
    
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

    //Mesh Integration option
    double x0 = meshOptions[0];
    double xend = meshOptions[1];
    double dx = meshOptions[2];
    double dxsave = meshOptions[3];
    double dxprint = meshOptions[4];
    cout<<"Mesh integration options : "<< "Domain : ["<< x0<<" "<<xend<<"][m] with a dx of "<<dx<<" [m]"<<endl;
    cout<<"The solution is printed every "<<dxprint<<" [m] and saved every "<<dxsave<<" [m]"<<endl;
    cout<<"-----------------------------------------------------------------------------------------------------------"<<endl;
    
    
    std::ofstream myflux;
    myflux.open("output.dat");
    myflux.close();
    
    cout<<"ODE integration ..."<<endl;
    
    cout<<"Position [x] :"<<setw(15)<<"T [K]"<<setw(15)<<"P [Pa]"<<setw(15)<<"U [m/s]"<<setw(15)<<"rho [kg/m^3]";
    for (int i = 0; i<ns; i++) {
        string Yname = "Y_";
        Yname+=mix->speciesName(i);
        cout << setw(15) << Yname;
    }
    cout<<endl;
    
    size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-8 , 1.0e-8 ) ,make_pair( equationsToSolve(mix,mdot) , jacobian(mix,mdot) ) ,x , x0 , xend, dx,intermediateSolution(mdot,mix,dxsave,dxprint));
    
    cout<<"ODE Integration finished in "<<num_of_steps<<" steps"<<endl;

  // ###########################################################################################
  // ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##
  // # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  // ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##
  // ###########################################################################################
  // #                                   Now it's my turn.                                     #
  // ###########################################################################################
  // ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##
  // # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  // ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##
  // ###########################################################################################
  //

cout << endl;
cout << "                                                                    " << endl;
cout << "              _ _ __ __ ________   __         __ \\  |  /           " << endl;
cout << "     _ _ __ __  _______________   /  \\_/\\ /\\_/  \\       -       " << endl;
cout << "              _ _ __ __ ________  \\__/   V   \\__/     __          " << endl; 
cout << "                                                 /   /  \\          " << endl;
cout << "   _      _   ___  ___ ___ _  _                      \\__/          " << endl;
cout << "  | |    /_\\ | _ \\/ __| __| \\| |  LAgrangian               \\    " << endl;
cout << "  | |__ / _ \\|   /\\__ \\ _|| .` |  Reactor for StrEams  \\  \\ \\ " << endl;
cout << "  |____/_/ \\_\\_|_\\|___/___|_|\\_|  in Nonequilibrium     \\  \\  " << endl;
cout << "                                                            \\      " << endl;
cout << "                                                                    " << endl; 
cout << endl;
                                                                                               

  size_t nSpecies = ns;
  size_t ii; // internal index 

  // ============  READING REFERENCE SOLUTION AND SAVING INITIAL CONDITIONS  =================

  ExternalField extField; // creating external field object

  // Reading shocking output file and filling internal variables
  string lineRead;
  ifstream filein ("output.dat");
  if(filein.is_open() == 0)
  {
    cout << "I could not open the file.. Does it exist?\nAborting!\n";
    exit(EXIT_FAILURE);
  }

  getline(filein,lineRead); // read first line
  std::stringstream nowstream(lineRead); // and unpack it in the nowstream variable

  double dummy, xNow, uNow, YiNow, TNow, PNow, rhoNow;     // those are temporary
  std::vector<double> x_vect, T_vect, P_vect;    // those are permanent
  std::vector< std::vector<double> > Yi_vect(nSpecies);    // and also this one

  // -----   Reading first line and saving initial conditions
  // NB: NELLA PARTE FINALE IL VETTORE DA SALVARE AVRA' UNA LUNGHEZZA PREIMPOSTATA TRAMITE L'IMPUT INIZIALE:
  //     INTEGRANDO CON integrate_const INFATTI SO QUANTI OUTPUT AVRO'!!!!!!!!!!!!
  //     MI EVITO IL PUSH BACK!!
  cout << "ATTENTION: this works only after modifing the Shocking part: we need also rho as output!\n";
  nowstream >> xNow >> TNow >> PNow >> rhoNow >> uNow;
  x_vect.push_back(xNow);
  T_vect.push_back(TNow);
  P_vect.push_back(PNow);
  
  // ..and saving initial value for species concentrations
  for(ii = 0; ii < nSpecies; ++ii)
  {
    nowstream >> YiNow;
    Yi_vect[ii].push_back(YiNow);
  }

  // -----   Now saving the externally computed field: I only care about position and velocity
  // The first line has already been extracted, so I save it here.
  extField.x_ext.push_back(xNow);
  extField.u_ext.push_back(uNow);
  extField.rho_ext.push_back(rhoNow);

  // Keep reading now until the EOF
  while(getline(filein,lineRead))
  {
    nowstream.clear();          // clear stream variable nowstream
    nowstream.str(lineRead);    // put newly read line into nowstream

    nowstream >> xNow >> dummy >> dummy >> rhoNow >> uNow; // ...and unpack it into my variables

    extField.x_ext.push_back(xNow);
    extField.u_ext.push_back(uNow);
    extField.rho_ext.push_back(rhoNow);
  }

  filein.close();  // I'm done with this data file

  // ================   CREATING MIXTURE and INITIALIZING SOLVER  ======================

  // Mixture Options Object
  Mutation::MixtureOptions myMixOpts("air5");
  myMixOpts.setStateModel("ChemNonEq1T");
  myMixOpts.setThermodynamicDatabase(thermoDataBase);

  // Mixture Object
  Mutation::Mixture myMix(myMixOpts);
//Mutation::Mixture myJacobiMix(myMixOpts);
// THE FOLLOWING MUST BE MODIFIED ACCORDINGLY TO INITIAL CONDITIONS!!!!!!!!!!!!!!!!!!!
  myMix.addComposition("N:0.8, O:0.2", true); 
//myJacobiMix.addComposition("N:0.8, O:0.2", true);
  // Lagrangian solver, Jacobian and Observer objects
  EquationsToSolve myEqs(&myMix, &extField);
  LagJacobian myJacobi(&myMix, &extField);
  Observer myObs(&myMix); 

  // ==============================   MAIN CYCLE   ===================================
  
  vec_type xState(nSpecies+1);
  double  T_kk, P_kk;
  double* rhoi_kk = new double [nSpecies];  // densities at current step
  double* Xi_kk   = new double [nSpecies];  // species molar fractions
  double* Yi_kk   = new double [nSpecies];  // species mass fractions

  double Yi_sum; // Internal variable to keep sure that the sum of mass fractions gives 1
 
  // Resizing variables to the right dimension (to avoid using push_back)
  P_vect.resize(extField.x_ext.size()-1);
  T_vect.resize(extField.x_ext.size()-1);
  for(ii = 0; ii < nSpecies; ++ii)
    Yi_vect[ii].resize(extField.x_ext.size() - 1);

  // --- BIG CYCLE --- Now re-integrate for every value along the given streamline
  for(int kk = 0; kk < extField.x_ext.size() - 1; ++kk)
  {
    // ---------   Initializing timestep variables 

    // Internal mass fraction (used for readability and since Mutation likes it)
    for(ii = 0; ii < nSpecies; ++ii)
      Yi_kk[ii] = Yi_vect[ii][kk]; 

    // Densities at current step
    for(ii = 0; ii < nSpecies; ++ii)
      rhoi_kk[ii] = Yi_kk[ii] * extField.rho_ext[kk];
  
    // Temperature at current step
    T_kk = T_vect[kk];

    // -----------   Setting mixture at the given densities and temperature
    myMix.setState(rhoi_kk, &T_kk, 1);

    // And then also save pressure
    P_vect[kk]   = myMix.pressure(T_vect[kk], extField.rho_ext[kk], Yi_kk);

    // -----------   Setting current velocities and positions into the solver
    extField.x_now[0]   = extField.x_ext[kk];
    extField.x_now[1]   = extField.x_ext[kk+1];
    extField.u_now[0]   = extField.u_ext[kk];
    extField.u_now[1]   = extField.u_ext[kk+1];
    extField.rho_now[0] = extField.rho_ext[kk];
    extField.rho_now[1] = extField.rho_ext[kk+1];

    // -----------   Setting initial condition for the State variable
    for(ii = 0; ii < nSpecies; ++ii) 
      xState[ii] = Yi_kk[ii];
    
    xState[nSpecies] = T_kk;
    
    // -----------   Integrating from step kk to kk + 1
    typedef runge_kutta4<vec_type> rk4;
    typedef runge_kutta_cash_karp54< vec_type > rkck54;
    typedef controlled_runge_kutta< rkck54 > ctrl_rkck54;

//    integrate_const( ctrl_rkck54(),
//                     myEqs, xState, extField.x_ext[kk], extField.x_ext[kk+1], 0.00000001, myObs );
   // Stiff integration
   double mydx = 0.005*abs(extField.x_ext[kk+1] - extField.x_ext[kk]) + 1e-12; // output step
   integrate_const(make_dense_output< rosenbrock4< double > >(1.0e-13, 1.0e-13),
                   make_pair(myEqs, myJacobi),
                   xState, extField.x_ext[kk], extField.x_ext[kk+1], mydx, myObs);

//   integrate_const( ctrl_rkck54(), myEqs,
//                   xState, extField.x_ext[kk], extField.x_ext[kk+1], 0.00000001, myObs);



    // -----------   Saving new variables at position kk + 1 
    Yi_sum = 0.0; // This ensures the sum of mass fractions gives 1 (numerical round-off?)
    for(ii = 0; ii < nSpecies; ++ii)
      Yi_sum += xState[ii];

    for(ii = 0; ii < nSpecies; ++ii)
      Yi_vect[ii][kk+1] = xState[ii] / Yi_sum;

    T_vect[kk+1] = xState[nSpecies];
    
  }

cout << "BE CAREFUL: THE OUTPUT IS NOT 'UNIQUE'!!! STARTING INTEGRATION POINTS ARE REPEATED\n";
cout << "I COULD PUT A CHECK IN THE OBSERVER MAYBE, SKIPPING SOME POINTS.....\n";

    return 0;
}


