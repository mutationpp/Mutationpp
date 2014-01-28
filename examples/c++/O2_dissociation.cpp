/**
 * @file O2_dissociation.cpp
 *
 * Example program demonstrating the use of the Mutation++ finite-rate 
 * chemistry routine netProductionRates().  This program simulates a constant
 * volume, adiabatic reactor beginning with pure O2 at 4000 K and density of
 * 2.85e-4 kg/m^3.  The program integrates the mass conservation equation for
 * O2 explicitely in time up to 100 seconds until the O/O2 mixture is 
 * essentially at equilibrium.  At each time-step, the temperature is computed
 * using a simple Newton method to solve the total energy equation for T based
 * on a constant density and total energy inside the reactor.  The reactor
 * average temperature, pressure, and mass fractions of O and O2 are output to
 * show the results of the simulation.  The final conditions at 100 seconds are
 * also compared with equilibrium values for the O/O2 mixture at the final
 * temperature and pressure.
 */

// Must include this header file to use the Mutation++ library
#include "mutation++.h"
#include <iostream>
#include <cmath>

using std::setw;
using std::cout;
using std::endl;

using namespace Mutation;
using namespace Mutation::Numerics;

/**
 * This class implements Newton's method for solving the energy equation for
 * temperature given a fixed energy and species mass fractions.
 */
class EnergyEqSolution : public NewtonSolver<double, EnergyEqSolution>
{
public:
    
    EnergyEqSolution(const Mixture& mix, const double etot)
        : m_mix(mix),
          m_etot(etot),
          mp_work(new double [mix.nSpecies()]),
          mp_Y(NULL)
    {
        setEpsilon(1.e-6);
        setJacobianLag(1);
    }
    
    virtual ~EnergyEqSolution() {
        delete [] mp_work;
    }
    
    void setMassFractions(double* const Y) {
        mp_Y = Y;
    }

    void updateFunction(const double& T)
    {
        const int ns = m_mix.nSpecies();
        m_mix.speciesHOverRT(T, mp_work);
        
        m_f = 0.0;
        for (int i = 0; i < ns; ++i)
            m_f += mp_Y[i] * (mp_work[i] - 1.0) / m_mix.speciesMw(i);
        
        m_f = m_f*RU*T - m_etot;
    }
    
    void updateJacobian(const double& T)
    {
        const int ns = m_mix.nSpecies();
        m_mix.speciesCpOverR(T, mp_work);
        
        m_dfdT = 0.0;
        for (int i = 0; i < ns; ++i)
            m_dfdT += mp_Y[i] * (mp_work[i] - 1.0) / m_mix.speciesMw(i);
        
        m_dfdT *= RU;
    }
    
    const double systemSolution() {
        return (m_f / m_dfdT);
    }
    
    const double norm() {
        return std::abs(m_f);
    }

private:
    
    const Mutation::Mixture& m_mix;
    
    double m_f;
    double m_dfdT;
    
    const double  m_etot;
    double* const mp_work;
    double*       mp_Y;
};


/**
 * Performs a simple O2 dissociation simulation.
 */
int main()
{
    Mutation::MixtureOptions opts("O2");
    opts.setStateModel("TPY");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);
    
    // Useful to know the species indices
    const int ns  = mix.nSpecies();
    const int iO2 = mix.speciesIndex("O2");
    const int iO  = mix.speciesIndex("O");
    
    // Initial mass fractions
    double* pY = new double [ns];
    pY[iO2] = 1.0;
    pY[iO]  = 0.0;
    
    // Initial temperature and density
    double TP [2];
    TP[0] = 4000.0;   // K
    double rho = 2.85e-4;  // kg/m^3
    
    // Set the initial mixture state
    TP[1] = mix.pressure(TP[0], rho, pY);
    mix.setState(TP, pY);
    
    // Compute the total energy
    double etot = mix.mixtureEnergyMass();
    
    // Setup the Newton solver for solving the energy equation
    EnergyEqSolution tSolver(mix, etot);
    //tSolver.setWriteConvergenceHistory(true);
    
    // Write the table header
    cout << setw(15) << "Time (s)";
    cout << setw(15) << "T (K)";
    cout << setw(15) << "P (Pa)";
    cout << setw(15) << "rho (kg/m^3)";
    cout << setw(15) << "e (J/kg)";
    cout << setw(15) << "Y_O";
    cout << setw(15) << "Y_O2";
    cout << endl;
    
    // Integrate in time until equilibrium is reached.
    double* pWdot = new double [ns];
    
    double dt   = 0.0;
    double time = 0.0;    
    double tout = 0.0001;
    
    // Write time = 0 conditions
    cout << setw(15) << time;
    cout << setw(15) << TP[0];
    cout << setw(15) << TP[1];
    cout << setw(15) << mix.density();
    cout << setw(15) << mix.mixtureEnergyMass();
    cout << setw(15) << pY[iO];
    cout << setw(15) << pY[iO2];
    cout << endl;
    
    while (time < 100.0) {
        // Get species net production rates due to chemical reactions
        mix.netProductionRates(pWdot);
        
        // Compute a stable time step
        dt = 1.0e-7 * std::abs(rho * pY[iO2] / pWdot[iO2]);
        if (time+dt > tout)
            dt = tout-time;
        
        // Explicitly integrate mass fractions
        time += dt;
        pY[iO2] += pWdot[iO2] * dt / rho;
        pY[iO]  = 1.0 - pY[iO2];
        
        // Compute the temperature at the new state by solving the energy 
        // equation
        //T = newton(mix, pY, etot, T);
        tSolver.setMassFractions(pY);
        tSolver.solve(TP[0]);
        
        // Finally set the mixture state
        TP[1] = mix.pressure(TP[0], rho, pY);
        mix.setState(TP, pY);
        
        // Print results
        if (std::abs(tout-time) < 1.0e-10) {
            cout << setw(15) << time;
            cout << setw(15) << TP[0];
            cout << setw(15) << TP[1];
            cout << setw(15) << mix.density();
            cout << setw(15) << mix.mixtureEnergyMass();
            cout << setw(15) << pY[iO];
            cout << setw(15) << pY[iO2];
            cout << endl;
            
            if (tout > 9.9)          tout += 10.0;
            else if (tout > 0.99)    tout += 1.0;
            else if (tout > 0.099)   tout += 0.1;
            else if (tout > 0.0099)  tout += 0.01;
            else if (tout > 0.00099) tout += 0.001;
            else                     tout += 0.0001;
        }
    }
    
    // Output what the equilibrium values are at these conditions
    /*mix.equilibrate(T, P);
    cout << "Equilibrium values: " << endl;
    cout << setw(30) << TP[0];
    cout << setw(15) << TP[1];
    cout << setw(15) << mix.Y()[iO];
    cout << setw(15) << mix.Y()[iO2];
    cout << endl;
    
    cout << "Total energy given:           " << etot*rho << endl;
    cout << "Total energy at equilibrium:  " << mix.mixtureEnergyMass()*mix.density() << endl;
    
    cout << "Total density given:          " << rho << endl;
    cout << "Total density at equilibrium: " << mix.density() << endl;*/
    
    // Clean up
    delete [] pY;
    delete [] pWdot;
}
