/**
 * @file MultiPhaseEquilSolver.cpp
 *
 * @brief Implementation of the Multiphase equilibrium solver.
 * @see class MultiPhaseEquilSolver
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "MultiPhaseEquilSolver.h"
#include "Thermodynamics.h"
#include "lp.h"

#include <set>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;

//#define VERBOSE
//#define SAVE_PATH
//#define SAVE_IC
//#define SAVE_EIGEN_SCRIPT

#include "Utilities.h"

using std::cout;
using std::endl;
using std::setw;

//using namespace Mutation::Numerics;

namespace Mutation {
    namespace Thermodynamics {
    

// Define static constants for class MultiPhaseEquilSolver
const double MultiPhaseEquilSolver::ms_eps_rel = 0.05;
const double MultiPhaseEquilSolver::ms_eps_abs = 1.0e-12;
const double MultiPhaseEquilSolver::ms_ds_inc  = 4.0;
const double MultiPhaseEquilSolver::ms_ds_dec  = 0.25;
const double MultiPhaseEquilSolver::ms_max_ds  = 1.0;//0.01;

// A change in temperature more than this triggers an update in initial solution
const double MultiPhaseEquilSolver::ms_temp_change_tol = 1.0;

// A change in pressure more than this triggers an update in initial solution
const double MultiPhaseEquilSolver::ms_pres_change_tol = 10000.0;

//==============================================================================

void MultiPhaseEquilSolver::Solution::initialize(int np, int nc, int ns)
{
    // Init sizes
    m_np = np;
    m_nc = nc;
    m_ns = ns;
    m_npr = np;
    m_ncr = nc;
    m_nsr = ns;

    int dsize = 2*ns+np+nc;
    int isize = ns+nc+np+2;

    // Allocate data storage if necessary
    if (dsize > m_dsize) {
        if (mp_ddata != NULL)
            delete [] mp_ddata;
        mp_ddata = new double [dsize];
    }
    m_dsize = dsize;

    if (isize > m_isize) {
        if (mp_idata != NULL)
            delete [] mp_idata;
        mp_idata = new int [isize];
    }
    m_isize = isize;

    // Setup the pointer locations
    mp_g      = mp_ddata;
    mp_y      = mp_g+ns;
    mp_lnNbar = mp_y+ns;
    mp_lambda = mp_lnNbar+np;

    mp_sizes = mp_idata;
    mp_sjr   = mp_sizes+np+2;
    mp_cir   = mp_sjr+ns;

    // Just fill all data with 0
    std::fill(mp_ddata, mp_ddata+m_dsize, 0.0);
    std::fill(mp_idata, mp_idata+m_isize, 0);
}

//==============================================================================

bool MultiPhaseEquilSolver::Solution::setupOrdering(
        int* species_group, bool* zero_constraint)
{
    static std::vector<int> previous_order(m_np+m_nc+4, 0);

    // Count the number of species in each group and order the species such that
    // the groups are contiguous
    for (int i = 0; i < m_np+2; ++i)
        mp_sizes[i] = 0;

    for (int i = 0; i < m_ns; ++i)
        mp_sizes[species_group[i]+1]++;

    for (int i = 1; i < m_np+1; ++i)
        mp_sizes[i] += mp_sizes[i-1];

    for (int i = 0; i < m_ns; ++i)
        mp_sjr[mp_sizes[species_group[i]]++] = i;

    for (int i = m_np+1; i > 0; --i)
        mp_sizes[i] = mp_sizes[i-1];
    mp_sizes[0] = 0;

    m_nsr = mp_sizes[m_np];

    // Update list of reduced constraint indices
    int index = 0;
    for (int i = 0; i < m_nc; ++i)
        if (!zero_constraint[i]) mp_cir[index++] = i;

    m_ncr = index;
    m_npr = m_np;
    int i = m_np-1;
    while (i > 0) {
        if (mp_sizes[i+1] - mp_sizes[i] == 0)
            removePhase(i);
        i--;
    }

    // Check and see if the ordering information has changed since the last call
    // Note this assumes that any change in species order garuntees at least one
    // of the following will change also: npr, ncr, sizes of each group, or
    // ordering of constraints.
    bool order_change =
            (previous_order[0] != m_npr) && (previous_order[1] != m_ncr);
    if (!order_change) {
        for (int i = 0; i < m_np+2; ++i)
            order_change |= (previous_order[i+2] != mp_sizes[i]);
        for (int i = 0; i < m_nc; ++i)
            order_change |= (previous_order[i+4+m_np] != mp_cir[i]);
    }

    // Save the new ordering information
    if (order_change) {
        previous_order[0] = m_npr;
        previous_order[1] = m_ncr;
        for (int i = 0; i < m_np+2; ++i)
            previous_order[i+2] = mp_sizes[i];
        for (int i = 0; i < m_nc; ++i)
            previous_order[i+4+m_np] = mp_cir[i];
    }

    return order_change;
}

//==============================================================================

void MultiPhaseEquilSolver::Solution::setG(
        const double* const p_g0, const double* const p_g, double s)
{
    const double s1 = 1.0 - s;
    for (int j = 0; j < m_ns; ++j)
        mp_g[j] = s1*p_g0[j] + s*p_g[j];
}

//==============================================================================

void MultiPhaseEquilSolver::Solution::setSolution(
    const double* const p_lambda, const double* const p_Nbar, const MatrixXd& B)
{
    for (int i = 0; i < m_ncr; ++i)
        mp_lambda[i] = p_lambda[i];
    for (int m = 0; m < m_npr; ++m)
        mp_lnNbar[m] = std::log(p_Nbar[m]);

    updateY(B);
}

//==============================================================================

void MultiPhaseEquilSolver::Solution::updateY(const MatrixXd& B)
{
    int jk;
    for (int m = 0; m < m_npr; ++m) {
        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j) {
            jk = mp_sjr[j];
            mp_y[j] = mp_lnNbar[m] - mp_g[jk];
            for (int i = 0; i < m_ncr; ++i)
                mp_y[j] += B(jk, mp_cir[i])*mp_lambda[i];
            mp_y[j] = std::exp(0.5*std::min(mp_y[j], 300.0));
        }
    }
}

//==============================================================================

MultiPhaseEquilSolver::Solution&
MultiPhaseEquilSolver::Solution::operator=(const Solution& state)
{
    m_np  = state.m_np;
    m_ns  = state.m_ns;
    m_npr = state.m_npr;
    m_ncr = state.m_ncr;
    m_nsr = state.m_nsr;

    if (state.m_dsize > m_dsize) {
        if (mp_ddata != NULL) delete mp_ddata;
        mp_ddata  = new double [state.m_dsize];
        mp_g      = mp_ddata;
        mp_y      = mp_g+m_ns;
        mp_lnNbar = mp_y+m_ns;
        mp_lambda = mp_lnNbar+m_np;
    }
    m_dsize = state.m_dsize;

    if (state.m_isize > m_isize) {
        if (mp_idata != NULL) delete mp_idata;
        mp_idata = new int [state.m_isize];
        mp_sizes = mp_idata;
        mp_sjr   = mp_sizes+m_np+2;
        mp_cir   = mp_sjr+m_ns;
    }
    m_isize = state.m_isize;

    std::copy(state.mp_ddata, state.mp_ddata+m_dsize, mp_ddata);
    std::copy(state.mp_idata, state.mp_idata+m_isize, mp_idata);

    return *this;
}

//==============================================================================

int MultiPhaseEquilSolver::Solution::removePhase(int phase)
{
//    cout << "removing phase:";
//    for (int j = mp_sizes[phase]; j < mp_sizes[phase+1]; ++j)
//        cout << " " << m_thermo.speciesName(mp_sjr[j]);
//    cout << endl;

    // Check that the phase number is feasible
    assert(phase >= 0);
    assert(phase < m_npr);
    assert(m_npr > 1);

    // If the phase is not the last non-empty phase, we need to shift it
    // to the right to make the non-empty species list contiguous
    const int size = mp_sizes[phase+1] - mp_sizes[phase];
    if (phase != m_npr-1) {
        int temp [size];

        // First copy the phase to be removed into temporary array
        int* ip = temp;
        for (int j = mp_sizes[phase]; j < mp_sizes[phase+1]; ++j)
            *ip++ = mp_sjr[j];

        // Shift all remaining phases to fill the space
        ip = mp_sjr + mp_sizes[phase];
        for (int j = mp_sizes[phase+1]; j < mp_sizes[m_np]; ++j)
            *ip++ = mp_sjr[j];

        // Place the removed phase at the end of the list (before
        // determined species)
        for (int m = phase+1; m < m_np; ++m)
            mp_sizes[m] = mp_sizes[m+1] - size;
        ip = temp;
        for (int j = mp_sizes[m_np-1]; j < mp_sizes[m_np]; ++j)
            mp_sjr[j] = *ip++;

        // Shift the solution to match the new phase ordering
        for (int j = mp_sizes[phase]; j < m_nsr-size; ++j)
            mp_y[j] = mp_y[j+size];
        for (int m = phase; m < m_npr-1; ++m)
            mp_lnNbar[m] = mp_lnNbar[m+1];

        // Set the phase number to the new number to return
        phase = m_np-1;
    }

    // Update reduced species and phase sizes
    m_npr--;
    m_nsr -= size;

    return phase;
}

//==============================================================================

int MultiPhaseEquilSolver::Solution::addPhase(int phase)
{
//	cout << "adding phase:";
//    for (int j = mp_sizes[phase]; j < mp_sizes[phase+1]; ++j)
//        cout << " " << m_thermo.speciesName(mp_sjr[j]);
//    cout << endl;

    assert(phase >= m_npr);
    assert(phase < m_np);

    // If the phase is not the first empty phase, then we need to shift
    // it to the front
    int size = mp_sizes[phase+1] - mp_sizes[phase];
    if (phase > m_npr) {
        int temp [size];

        // First copy the phase to be added into temporary array
        int* p = temp;
        for (int i = mp_sizes[phase]; i < mp_sizes[phase+1]; ++i)
            *p++ = mp_sjr[i];

        // Shift phases that come before to the back
        p = mp_sjr + mp_sizes[phase+1]-1;
        for (int i = mp_sizes[phase]-1; i >= m_nsr; --i)
            *p-- = mp_sjr[i];

        // Place the added phase at the end of the included list
        for (int i = phase+1; i > m_npr; --i)
            mp_sizes[i] = mp_sizes[i-1] + size;
        p = temp;
        for (int i = mp_sizes[m_npr]; i < mp_sizes[m_npr+1]; ++i)
            mp_sjr[i] = *p++;
    }

    // Update reduced species and phase sizes
    m_npr++;
    m_nsr += size;

    return (m_npr-1);
}

//==============================================================================

int MultiPhaseEquilSolver::Solution::checkCondensedPhase(
        const MatrixXd& B)
{
    if (m_np <= m_ncr)
        return -1;

    int    min_m   = -1;
    double min_sum = 0.0;

    for (int m = m_npr; m < m_np; ++m) {
        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j) {
            double sum = mp_g[mp_sjr[j]];
            for (int i = 0; i < m_ncr; ++i)
                sum -= mp_lambda[i]*B(mp_sjr[j], mp_cir[i]);

            if (sum < min_sum) {
                min_m   = m;
                min_sum = sum;
            }
        }
    }

    return min_m;
}

//==============================================================================

void MultiPhaseEquilSolver::Solution::printOrder()
{
    cout << "Species order:" << endl;
    cout << "  Active Phases:" << endl;
    for (int m = 0; m < m_npr; ++m) {
        cout << "    " << m << ":";
        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
            cout << " " << m_thermo.speciesName(mp_sjr[j]);
        cout << endl;
    }
    cout << "  Inactive Phases:" << endl;
    for (int m = m_npr; m < m_np; ++m) {
        cout << "    " << m << ":";
        for (int j = mp_sizes[m]; j < mp_sizes[m+1]; ++j)
            cout << " " << m_thermo.speciesName(mp_sjr[j]);
        cout << endl;
    }
    cout << "  Determined Species:" << endl;
    cout << "   ";
    for (int j = mp_sizes[m_np]; j < m_ns; ++j)
        cout << " " << m_thermo.speciesName(mp_sjr[j]);
    cout << endl;
}

//==============================================================================

void MultiPhaseEquilSolver::Solution::printSolution()
{
    cout << "Solution:" << endl;
    cout << "  lambda = " << endl;
    for (int i = 0; i < m_ncr; ++i)
        cout << "    " << mp_lambda[i] << endl;
    cout << "  Nbar = " << endl;
    for (int m = 0; m < m_npr; ++m)
        cout << "    " << std::exp(mp_lnNbar[m]) << endl;
    cout << "  N = " << endl;
    for (int j = 0; j < m_nsr; ++j)
        cout << "   " << setw(12) << m_thermo.speciesName(mp_sjr[j])
             << ": " << mp_y[j]*mp_y[j] << endl;
}

//==============================================================================

void MultiPhaseEquilSolver::Solution::printG()
{
    cout << "Current G vector" << endl;
    for (int j = 0; j < m_ns; ++j)
        cout << setw(12) << m_thermo.speciesName(j) << ": "
             << mp_g[j] << endl;
}

void MultiPhaseEquilSolver::Solution::unpackMoleFractions(
    double* const p_x, const MoleFracDef mfd)
{
    switch (mfd) {
    case IN_PHASE:
        for (int m = 0; m < npr(); ++m) {
            const double nbar = std::exp(lnNbar()[m]);
            for (int j = sizes()[m]; j < sizes()[m+1]; ++j)
                p_x[sjr()[j]] = y()[j]*y()[j]/nbar;
        }

        for (int j = nsr(); j < m_ns; ++j)
            p_x[sjr()[j]] = 0.0;
        break;

    case GLOBAL:
        for (int j = 0; j < m_ns; ++j)
            p_x[j] = 0.0;
        double sum = 0.0;
        for (int j = 0; j < nsr(); ++j) {
            p_x[sjr()[j]] = y()[j]*y()[j];
            sum += p_x[sjr()[j]];
        }
        for (int j = 0; j < m_ns; ++j)
            p_x[j] /= sum;
        break;
    }
}

//==============================================================================

MultiPhaseEquilSolver::MultiPhaseEquilSolver(
    const Thermodynamics& thermo, const bool pure_condensed) :
        m_thermo(thermo),
        m_pure_condensed(pure_condensed),
        m_solution(thermo),
        m_T(0.0),
        m_P(0.0)
{
    // Sizing information
    m_ns  = m_thermo.nSpecies();
    m_ne  = m_thermo.nElements();
    m_nc  = m_ne;
    
    // System element matrix
    m_B = m_thermo.elementMatrix();
    m_Br = m_B;
    
    // Compute phase information
    mp_phase = new int [m_ns];
    initPhases();

    m_tableau_capacity = std::max((m_nc+2)*(m_ns+2), m_ns+m_nc+m_np*(m_ns+m_nc));
    mp_tableau = new double [m_tableau_capacity];
    mp_ming    = new double [m_ns];
    mp_maxmin  = new double [m_ns];
    mp_g       = new double [m_ns];
    mp_g0      = new double [m_ns];
    mp_c       = new double [m_nc];
    
    std::fill(mp_c, mp_c+m_nc, 0.0);
    
    // Initialize storage for the solution object
    m_solution.initialize(m_np, m_nc, m_ns);
}

//==============================================================================

MultiPhaseEquilSolver::~MultiPhaseEquilSolver()
{
    delete [] mp_phase;
    delete [] mp_tableau;
    delete [] mp_ming;
    delete [] mp_maxmin;
    delete [] mp_g;
    delete [] mp_g0;
    delete [] mp_c;
};

//==============================================================================

void MultiPhaseEquilSolver::addConstraint(const double *const p_A)
{
    // Save constraints
    m_constraints.push_back(Map<const VectorXd>(p_A, m_ns));

    // Update the B matrix
    m_B.resize(m_ns, ++m_nc);
    m_B.block(0,0,m_ns,m_ne) = m_thermo.elementMatrix();

    for (int i = 0; i < m_nc-m_ne; ++i)
        m_B.col(m_ne+i) = m_constraints[i];

    // Resize the tableau used in the simplex solver
    size_t size = (m_nc+2)*(m_ns+2);
    if (m_tableau_capacity < size) {
        delete [] mp_tableau;
        mp_tableau = new double [size];
        m_tableau_capacity = size;
    }

    // Resize the c vector
    delete [] mp_c;
    mp_c = new double [m_nc];
    std::fill(mp_c, mp_c+m_nc, 0.0);

    // The solution should also be reinitialized
    m_solution.initialize(m_np, m_nc, m_ns);
}

//==============================================================================
        
void MultiPhaseEquilSolver::clearConstraints()
{
    m_constraints.clear();
    m_B = m_thermo.elementMatrix();
    m_nc = m_ne;
}

//==============================================================================
        
void MultiPhaseEquilSolver::dNdT(double *const p_dNdt) const
{
    throw NotImplementedError("MultiPhaseEquilSolver::dNdT");
}

//==============================================================================

void MultiPhaseEquilSolver::dXdP(double *const p_dxdp) const
{
    throw NotImplementedError("MultiPhaseEquilSolver::dXdP");
}

//==============================================================================

void MultiPhaseEquilSolver::dSoldg(const double* const p_dg, VectorXd& dx) const
{
    const int ncr = m_solution.ncr();
    const int npr = m_solution.npr();
    const int nsr = m_solution.nsr();
    const int neq = ncr + npr;

    // Compute the Jacobian matrix
    MatrixXd A(neq, neq);
    formSystemMatrix(A);

    // Compute the RHS
    VectorXd rhs = VectorXd::Zero(neq);

    double temp;
    int jk;

    for (int m = 0, j = 0; m < npr; ++m) {
        for (; j < m_solution.sizes()[m+1]; ++j) {
            jk = m_solution.sjr()[j];
            temp = m_solution.y()[j];
            temp *= temp * p_dg[jk];

            for (int i = 0; i < ncr; ++i)
                rhs(i) += temp * m_B(jk, m_solution.cir()[i]);

            rhs(ncr+m) += temp;
        }
    }

    // Get the system solution
    //SVD<double> svd(A);
    //svd.solve(dx, rhs);
    dx = A.selfadjointView<Upper>().ldlt().solve(rhs);
}

//==============================================================================

void MultiPhaseEquilSolver::dXdg(const double* const p_dg, double* const p_dX) const
{
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int npr = m_solution.npr();
    const int neq = ncr + npr;

    // Special case for a single species
    if (nsr == 1) {
        Map<ArrayXd>(p_dX, m_ns).setZero();
        return;
    }

    // Compute the change in the solution variables due to change in dg
    VectorXd dx(neq);
    dSoldg(p_dg, dx);

    double temp;
    int jk;

    // Finally compute the dX/dg for all non zero species
    for (int m = 0, j = 0; m < npr; ++m) {
        double Nbar = std::exp(m_solution.lnNbar()[m]);
        for (; j < m_solution.sizes()[m+1]; ++j) {
            jk = m_solution.sjr()[j];
            p_dX[jk] = -p_dg[jk];
        
            for (int i = 0; i < ncr; ++i)
                p_dX[jk] += m_B(jk,m_solution.cir()[i])*dx(i);
        
            temp = m_solution.y()[j];
            p_dX[jk] *= temp*temp/Nbar;
        }
    }

    // Set the remaining species to zero
    const int* p_jk = m_solution.sjr()+nsr;
    for ( ; p_jk != m_solution.sjr()+m_ns; ++p_jk)
        p_dX[*p_jk] = 0.0;
}

void MultiPhaseEquilSolver::dXdc(int i, double* const p_dX)
{
    const int ncr = m_solution.ncr();
    const int npr = m_solution.npr();
    const int nsr = m_solution.nsr();
    const int neq = ncr + npr;

    // Get the system solution
    VectorXd dx(neq);
    rates_ci(i, dx);

    int jk;
    double temp;

    // Finally compute the dX/dci for all non zero species
    for (int m = 0, j = 0; m < npr; ++m) {
        double Nbar = std::exp(m_solution.lnNbar()[m]);
        for (; j < m_solution.sizes()[m+1]; ++j) {
            jk = m_solution.sjr()[j];
            p_dX[jk] = 0.0;

            for (int i = 0; i < ncr; ++i)
                p_dX[jk] += m_B(jk,m_solution.cir()[i])*dx(i);

            temp = m_solution.y()[j];
            p_dX[jk] *= temp*temp/Nbar;
        }
    }

    // Set the remaining species to zero
    const int* p_jk = m_solution.sjr()+nsr;
    for ( ; p_jk != m_solution.sjr()+m_ns; ++p_jk)
        p_dX[*p_jk] = 0.0;
}

//==============================================================================

void MultiPhaseEquilSolver::dNdg(const double* const p_dg, double* const p_dN) const
{
    const int ncr = m_solution.ncr();
    const int npr = m_solution.npr();
    const int nsr = m_solution.nsr();
    const int neq = ncr + npr;

    // Compute the change in the solution variables due to change in dg
    VectorXd dx(neq);
    dSoldg(p_dg, dx);

    double temp;
    int jk;

    // Finally compute the dN/dg for all non zero species
    for (int m = 0, j = 0; m < npr; ++m) {
        for (; j < m_solution.sizes()[m+1]; ++j) {
            jk = m_solution.sjr()[j];
            p_dN[jk] = dx(ncr+m) - p_dg[jk];

            for (int i = 0; i < ncr; ++i)
                p_dN[jk] += m_B(jk,m_solution.cir()[i])*dx(i);

            temp = m_solution.y()[j];
            p_dN[jk] *= temp*temp;
        }
    }

    // Set the remaining species to zero
    const int* p_jk = m_solution.sjr()+nsr;
    for ( ; p_jk != m_solution.sjr()+m_ns; ++p_jk)
        p_dN[*p_jk] = 0.0;
}

//==============================================================================

void MultiPhaseEquilSolver::initPhases()
{
    std::set<int> phases;
    
    // Assign an integer value to each species corresponding to the phase it
    // belongs to (also keep track of unique phases)
    int npure = 0;
    for (int i = 0; i < m_ns; ++i) {
        mp_phase[i] = (int) m_thermo.species(i).phase();
        if (mp_phase[i] > 0 && m_pure_condensed)
            mp_phase[i] = ++npure;
        phases.insert(mp_phase[i]);
    }
    
    // Number of phases = number of unique phase numbers
    m_np = phases.size();
    
    // Make sure all the phase numbers are in the set {0,...,np-1}
    bool empty;
    for (int p = 0; p < m_np; ++p) {
        empty = true;
        for (int i = 0; i < m_ns; ++i) {
            if (mp_phase[i] == p) {
                empty = false;
                break;
            }
        }
        
        if (empty) {
            for (int i = 0; i < m_ns; ++i)
                if (mp_phase[i] > p) mp_phase[i]--;
            p--;
        }
    }
}

//==============================================================================

std::pair<int, int> MultiPhaseEquilSolver::equilibrate(
    double T, double P, const double* const p_cv, double* const p_sv,
    MoleFracDef mfd)
{

    DEBUG("equilibrate(" << T << " K, " << P << " Pa,")
    for (int i = 0; i < m_nc; ++i)
		DEBUG(m_thermo.elementName(i) << " " << p_cv[i] << (i == m_ne-1 ? ")" : ","))
	DEBUG(endl)
	//exit(1);
    
    // Special case for 1 species
    if (m_ns == 1) {
        DEBUG("only one species..." << endl)
        p_sv[0] = 1.0;
        return std::make_pair(0,0);
    }
    
    // Compute the initial conditions lambda(0), Nbar(0), N(0), and g(0)
    if (!initialConditions(T, P, p_cv)) {
        DEBUG("could not compute the initial conditions!" << endl)
        for (int i = 0; i < m_ns; ++i)
            p_sv[i] = 0;
        return std::make_pair(-1,-1);
    }
#ifdef SAVE_IC
    // Save the initial mole fractions
    std::ofstream of("ic_mf.dat", std::ios_base::app);
    m_solution.unpackMoleFractions(p_sv, mfd);
    of << setw(15) << m_T << setw(15) << m_P;
    for (int i = 0; i < m_ns; ++i)
        of << setw(15) << p_sv[i];
    of << endl;
    of.close();

    // Save the initial solution
    of.open("ic_sol.dat", std::ios_base::app);
    of << setw(15) << m_T << setw(15) << m_P;
    for (int i = 0; i < m_solution.ncr(); ++i)
        of << setw(15) << m_solution.lambda()[i];
    of << endl;
    of.close();
#endif

#ifdef VERBOSE
//    cout << "Gj/RT = " << endl;
//    for (int j = 0; j < m_ns; ++j)
//        cout << mp_g[j] << endl;
//    m_solution.printOrder();
    cout << "Initial conditions: " << endl;
    m_solution.printSolution();
//    m_solution.printG();
#endif

#ifdef SAVE_PATH
    std::ofstream of("path.dat");
    of << setw(12) << "s";
    of << setw(12) << "res";
    of << setw(12) << "G/RT[mol]";
    for (int i = 0; i < m_solution.ncr(); ++i)
        of << setw(12) << "lam_" + m_thermo.elementName(m_solution.cir()[i]);
//    for (int m = 0; m < m_solution.npr(); ++m) {
//        if (m_solution.sizes()[m] == m_solution.sizes()[m+1]-1)
//            of << setw(12) << "N_" + m_thermo.speciesName(m_solution.sjr()[m_solution.sizes()[m]]);
//        else
//            of << setw(11) << "N_" << m;
//    }
    of << endl;
#endif

    double s  = 0.0;
    double ds = ms_max_ds;
    double res, resk = ms_eps_abs;
    Solution last_solution(m_thermo);
    VectorXd dx(m_solution.ncr()+m_solution.npr());
    
    // Integrate the equilibrium solution from s = 0 to s = 1
    m_niters = 0;
    m_nnewts = 0;
    while (s < 1.0) {
        // Save current solution in case we need to take a smaller step
        last_solution = m_solution;
        
        // Compute the rates for lambda and Nbar
        rates(dx);
        #ifdef VERBOSE
        cout << "iter = " << m_niters << ", s = " << s << ", ds = " << ds << endl;
        cout << "dx = " << endl;
        cout << dx << endl;
        #endif

#ifdef SAVE_PATH
        // Compute the Gibbs free energy of the system
        of << setw(12) << s;
        of << setw(12) << resk;

        double G = 0.0;
        for (int i = 0; i < m_solution.ncr(); ++i)
            G += m_solution.lambda()[i] * mp_c[m_solution.cir()[i]];
        of << setw(12) << G;

        for (int i = 0; i < m_solution.ncr(); ++i)
            of << setw(12) << m_solution.lambda()[i];
//        for (int m = 0; m < m_solution.npr(); ++m)
//            of << setw(12) << std::exp(m_solution.lnNbar()[m]);

        of << endl;
#endif

        // Keep trying to obtain a solution at progressively smaller steps
        // until a solution is obtained
        bool solution_obtained = false;
        while (!solution_obtained) {
            // Get trial solution
            m_solution.setG(mp_g0, mp_g, s+ds);
            m_solution.update(dx*ds, m_B);

            // Use newton to reduce the residual
            res = newton();
            solution_obtained =
            		(res < std::max((1.0+ms_eps_rel)*resk, ms_eps_abs));

            // If newton fails to converge then reduce the step size and try
            // again, otherwise we can continue
            if (solution_obtained)
            	resk = res;
            else {
                #ifdef VERBOSE
                cout << "Could not converge to solution, reducing step size" << endl;
                #endif
                ds *= ms_ds_dec;
                m_solution = last_solution;
            }
        }
        
        // We found the solution at s, increase step size and move forward
        #ifdef VERBOSE
        cout << "Step succeded, increasing step size" << endl;
        #endif

        s += ds;
        ds = std::min(ms_max_ds, std::min(ds*ms_ds_inc, 1.0-s));
        m_niters++;

        // If we have finished then check for a phase redistribution
        if (s >= 1.0 && phaseRedistribution()) {
        	s = 0.0;
        	ds = ms_max_ds;
        }
    }
    
    // Check one last newton
    resk = newton();

#ifdef SAVE_PATH
    // Compute the Gibbs free energy of the system
    of << setw(12) << s;
    of << setw(12) << resk;

    double G = 0.0;
    for (int i = 0; i < m_solution.ncr(); ++i)
        G += m_solution.lambda()[i] * mp_c[m_solution.cir()[i]];
    of << setw(12) << G;

    for (int i = 0; i < m_solution.ncr(); ++i)
        of << setw(12) << m_solution.lambda()[i];
//    for (int m = 0; m < m_solution.npr(); ++m)
//        of << setw(12) << std::exp(m_solution.lnNbar()[m]);

    of << endl;
#endif

    if (resk > ms_eps_abs) {
    	cout << "Warning: equilibrium solver finished with residual of "
    	     << resk << "!";
    }

    #ifdef SAVE_EIGEN_SCRIPT
        rates(dx, true);
    #endif

    // Finally, unwrap the solution for the user and return convergence stats
    m_solution.unpackMoleFractions(p_sv, mfd);
    

#ifdef SAVE_PATH
    of.close();
#endif

    return std::make_pair(m_niters, m_nnewts);
}

//==============================================================================

// Simple comparison function for sorting which is used in the next function.
bool sortLargestSpecies(
	const std::pair<int, double>& s1, const std::pair<int, double>& s2)
{
	return (s1.second > s2.second);
}

//==============================================================================

bool MultiPhaseEquilSolver::phaseRedistribution()
{
    const int new_phase = m_solution.checkCondensedPhase(m_B);
    
    // Return false if there is nothing to do
    if (new_phase < 0)
    	return false;

    // Add the phase with the most negative vapor pressure test value to the end
    // of the species list in the solution object
	const int phase = m_solution.addPhase(new_phase);

	// Get the updated sizing information
	const int ncr = m_solution.ncr();
	const int npr = m_solution.npr();
	const int nsr = m_solution.nsr();

	const int* const p_sizes = m_solution.sizes();
	const int* const p_sjr   = m_solution.sjr();
	const int* const p_cir   = m_solution.cir();

	const double* p_y = m_solution.y();

	// Use the tableau as temporary storage
	double* p_N = mp_tableau;
	double* p_Nbar = p_N + nsr;
	double* p_lambda = p_Nbar + npr;

	// Compute the current moles vector for the phases already present
	for (int j = 0; j < p_sizes[npr-1]; ++j)
		p_N[j] = p_y[j]*p_y[j];

	// Add some moles to the new phase which is scaled based on the amount of
	// elements in the new phase
	double Nbar_new = 1.0e-6;
	int size = p_sizes[npr] - p_sizes[npr-1];
	for (int j = p_sizes[npr-1]; j < p_sizes[npr]; ++j)
		p_N[j] = Nbar_new / size;

	// Distribute some moles from ncr largest species into the moles of the new
	// phase. First gather the ncr largest species...
	std::vector< std::pair<int, double> > largest(
		m_solution.ncr(), std::make_pair(-1, 0.0));
	for (int j = 0; j < p_sizes[npr-1]; ++j) {
		if (p_N[j] >= largest.back().second) {
			largest.back() = std::make_pair(j, p_N[j]);
			std::sort(largest.begin(), largest.end(), sortLargestSpecies);
		}
	}

	// Form the system matrix to be solved
	MatrixXd Blarge(ncr, ncr);
	int jk;
	for (int j = 0; j < ncr; ++j) {
		jk = p_sjr[largest[j].first];
		for (int i = 0; i < ncr; ++i)
			Blarge(i,j) = m_B(jk,i);
	}

	// Compute the RHS which is B'*N
	//VectorXd rhs = -asVector(p_N, nsr) * m_solution.reducedMatrix(m_B);
	VectorXd rhs = -Map<VectorXd>(p_N, nsr) * m_solution.reducedMatrix(m_B, m_Br);
	for (int i = 0; i < ncr; ++i)
		rhs(i) += mp_c[p_cir[i]];

	// Solve the linear system for the update in the species moles vector
	//SVD<double> svd(Blarge);
	//RealVector x(ncr);
	//svd.solve(x, rhs);
	VectorXd x = Blarge.jacobiSvd(ComputeThinU | ComputeThinV).solve(rhs);

	// Finally update the species moles vector
	for (int j = 0; j < ncr; ++j)
		p_N[largest[j].first] += x(j);

	// Initialize the solution with residual of zero, based on the redistributed
	// species moles
	initZeroResidualSolution(p_N, p_Nbar, p_lambda);

	// Update the solution
	m_solution.setG(mp_g0, mp_g, 0.0);
	m_solution.setSolution(p_lambda, p_Nbar, m_B);

//        // Adhere to the phase rule (npr <= ncr)
//        if (m_solution.npr() > m_solution.ncr()) {
//            int lowest_phase = 1;
//            int lowest_value = m_solution.lnNbar()[1];
//
//            for (int m = 2; m < m_solution.npr()-1; ++m) {
//                if (m_solution.lnNbar()[m] < lowest_value) {
//                    lowest_value = m_solution.lnNbar()[m];
//                    lowest_phase = m;
//                }
//            }
//
//            m_solution.removePhase(lowest_phase);
//        }

    return true;
}

//==============================================================================

void MultiPhaseEquilSolver::rates(VectorXd& dx, bool save)
{
    // Get some of the solution variables
    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    
    const int* const p_sjr   = m_solution.sjr();
    const int* const p_sizes = m_solution.sizes();
    
    //const double* const p_y = m_solution.y();
    Map<const VectorXd> y(m_solution.y(), nsr);

    // Compute a least squares factorization of H
    static MatrixXd H; H = y.asDiagonal()*m_solution.reducedMatrix(m_B, m_Br);
    static JacobiSVD<MatrixXd> svd; svd.compute(H, ComputeThinU | ComputeThinV);

    // Use tableau for temporary storage
    Map<VectorXd> ydg(mp_tableau, nsr);
    Map<VectorXd> dlamg(ydg.data()+ydg.size(), ncr);
    Map<MatrixXd> dlamy(dlamg.data()+dlamg.size(), ncr, npr);
    Map<MatrixXd> P(dlamy.data()+dlamy.size(), nsr, npr);

    // Compute dlamg
    for (int j = 0, jk = p_sjr[0]; j < nsr; jk = p_sjr[++j])
        ydg[j] = y[j] * (mp_g[jk] - mp_g0[jk]);
    dlamg = svd.solve(ydg);
    
    // Compute dlambda_m for each phase m
//    for (int m = 0; m < npr; ++m) {
//        rhs.setZero();
//        rhs.segment(p_sizes[m], p_sizes[m+1]-p_sizes[m]) =
//            y.segment(p_sizes[m], p_sizes[m+1]-p_sizes[m]);
//        Map<VectorXd>(p_dlamy+m*ncr, ncr) = svd.solve(rhs);
//    }
    P.setZero();
    for (int m = 0; m < npr; ++m)
        P.block(p_sizes[m], m, p_sizes[m+1]-p_sizes[m], 1) =
            y.segment(p_sizes[m], p_sizes[m+1]-p_sizes[m]);
    dlamy = svd.solve(P);

    // Compute the linear system to be solved for the d(lnNbar)/ds variables
    MatrixXd A = MatrixXd::Zero(npr, npr);

    double sum;
    for (int p = 0; p < npr; ++p) {
        for (int m = p; m < npr; ++m) {
            for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j) {
                sum = 0.0;
                for (int i = 0; i < ncr; ++i)
                    sum += H(j,i)*dlamy(i,p);
                A(m,p) += y[j]*sum;
            }
        }
    }
    
    VectorXd b(npr);
    int j, n;
    for (int m = 0; m < npr; j = ++m) {
        j = p_sizes[m]; n = p_sizes[m+1]-j;
        b(m) = y.segment(j,n).dot((H*dlamg - ydg).segment(j,n));
    }

    // Solve for d(lnNbar)/ds
    dx.tail(npr) = A.selfadjointView<Upper>().ldlt().solve(b);

    // Finally compute the d(lambda)/ds vector
    dx.head(ncr) = dlamg - dlamy*dx.tail(npr);
}

void MultiPhaseEquilSolver::rates_ci(int k, VectorXd& dx)
{
    // Get some of the solution variables
    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();

    const int* const p_sjr   = m_solution.sjr();
    const int* const p_sizes = m_solution.sizes();

    VectorXd y = Map<const ArrayXd>(m_solution.y(), nsr).max(1.e-6);

    // Compute a least squares factorization of H
    static MatrixXd H; H = y.asDiagonal()*m_solution.reducedMatrix(m_B, m_Br);
    static JacobiSVD<MatrixXd> svd; svd.compute(H, ComputeThinU | ComputeThinV);

    // Use tableau for temporary storage
    Map<VectorXd> phi(mp_tableau, nsr);
    Map<VectorXd> dlamg(phi.data()+phi.size(), ncr);
    Map<MatrixXd> dlamy(dlamg.data()+dlamg.size(), ncr, npr);
    Map<MatrixXd> P(dlamy.data()+dlamy.size(), nsr, npr);

    // Compute dlamg
    phi.setZero();
    for (int j = 0; j < nsr; ++j) {
        if (m_Br.row(j).array().abs().sum() == std::abs(m_Br(j,k))) {
            phi(j) = 1.0/(y(j)*m_Br(j,k));
            break;
        }
    }
    cout << "phi = \n" << phi << endl;
    dlamg = svd.solve(phi);
    cout << "dlamg = \n" << dlamg << endl;

    // Compute dlambda_m for each phase m
//    for (int m = 0; m < npr; ++m) {
//        rhs.setZero();
//        rhs.segment(p_sizes[m], p_sizes[m+1]-p_sizes[m]) =
//            y.segment(p_sizes[m], p_sizes[m+1]-p_sizes[m]);
//        Map<VectorXd>(p_dlamy+m*ncr, ncr) = svd.solve(rhs);
//    }
    P.setZero();
    for (int m = 0; m < npr; ++m)
        P.block(p_sizes[m], m, p_sizes[m+1]-p_sizes[m], 1) =
            y.segment(p_sizes[m], p_sizes[m+1]-p_sizes[m]);
    dlamy = svd.solve(P);

    // Compute the linear system to be solved for the d(lnNbar)/ds variables
    MatrixXd A = MatrixXd::Zero(npr, npr);

    double sum;
    for (int p = 0; p < npr; ++p) {
        for (int m = p; m < npr; ++m) {
            for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j) {
                sum = 0.0;
                for (int i = 0; i < ncr; ++i)
                    sum += H(j,i)*dlamy(i,p);
                A(m,p) += y[j]*sum;
            }
        }
    }

    VectorXd b(npr);
    int j, n;
    for (int m = 0; m < npr; j = ++m) {
        j = p_sizes[m]; n = p_sizes[m+1]-j;
        b(m) = y.segment(j,n).dot((H*dlamg).segment(j,n));
    }

    // Solve for d(lnNbar)/ds
    dx.tail(npr) = A.selfadjointView<Upper>().ldlt().solve(b);

    // Finally compute the d(lambda)/ds vector
    dx.head(ncr) = dlamg - dlamy*dx.tail(npr);
}

//==============================================================================

double MultiPhaseEquilSolver::newton()
{
    using namespace Eigen;
    #ifdef VERBOSE
    cout << "newton:" << endl;
    #endif
    const int    max_iters = 5;
    const double phase_tol = std::log(1.0e-6);
    
    int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();
    const double* const p_lnNbar = m_solution.lnNbar();
    
    // First compute the residual
    static VectorXd r; r.resize(ncr+npr);
    computeResidual(r);
    
    double res = r.norm();
    if (res > 1.0)
        return res;

    static MatrixXd A;  A.resize(ncr+npr,ncr+npr);
    static VectorXd dx; dx.resize(ncr+npr);
    
    int iter = 0;
    while (res > ms_eps_abs && iter < max_iters) {
        #ifdef VERBOSE
        cout << "newton, iter = " << iter << ", " << "res = " << res << endl;
        m_solution.printSolution();
        #endif
        // First check if a phase needs to be removed (always include gas phase)
        int m = 1;
        while(m < m_solution.npr()) {
            if (m_solution.lnNbar()[m] < phase_tol)
                m_solution.removePhase(m);
            else
                m++;
        }
        
        // Update residual if we need to
        if (npr != m_solution.npr()) {
            npr = m_solution.npr();
            r.resize(ncr+npr);
            computeResidual(r);
            res = r.norm();

            A.resize(ncr+npr,ncr+npr);
            //dx.resize(ncr+npr);
        }
        
        // Compute the system jacobian
        formSystemMatrix(A);
        #ifdef VERBOSE
        cout << "jacobian matrix = " << endl;
        cout << A << endl;
        cout << "residual vector = " << endl;
        cout << r << endl;
        #endif
        
        // Solve the linear system (if it is singular then don't bother)
        static LDLT<MatrixXd, Upper> ldlt;
        ldlt.compute(A);
        dx = ldlt.solve(-r);
        
        #ifdef VERBOSE
        cout << "dx = " << endl;
        cout << dx << endl;
        #endif
        
        // Update the solution
        m_solution.update(dx, m_B);
        
        // Update residual
        computeResidual(r);
        double new_res = r.norm();
        #ifdef VERBOSE
        cout << "new res = " << new_res << endl;
        #endif
        iter++;
        if (new_res > res)
            break;
        else
            res = new_res;
    }
    
    m_nnewts += iter;
    return res;
}

//==============================================================================

bool MultiPhaseEquilSolver::checkForDeterminedSpecies()
{
    // First determine which species must be zero (either due to constraints or
    // to temperature limits on condensed species)
    int  species_group[m_ns];
    bool zero_constraint[m_nc];
    
    // Group each species
    for (int i = 0; i < m_ns; ++i) {
        species_group[i] = mp_phase[i];
        if (species_group[i] > 0 && !m_thermo.speciesThermoValidAtT(i, m_T))
            species_group[i] = m_np;
    }
    
    // Check for constraints which force species to be zero
    for (int i = 0; i < m_nc; ++i) {
        zero_constraint[i] = false;
        if (mp_c[i] == 0.0) {
            // Determine first non-zero sign
            bool pos;
            int j;
            for (j = 0; j < m_ns; ++j) {
                if (m_B(j,i) != 0.0) {
                    pos = m_B(j,i) > 0.0;
                    break;
                }
            }
            
            // Compare all the other signs to the first non-zero sign
            for ( ; j < m_ns; ++j)
                if ((m_B(j,i) != 0.0) && (pos != (m_B(j,i) > 0.0)))
                    break;
            
            if (j == m_ns) {
                // Keep track of the species which will be zeroed
                for (j = 0; j < m_ns; ++j) {
                    if (m_B(j,i) != 0.0)
                        species_group[j] = m_np;
                }
                
                // Keep track of the constraint that is zero
                zero_constraint[i] = true;
            }
        }
    }

    m_solution.setupOrdering(&species_group[0], zero_constraint);
    return true;
}

//==============================================================================

bool MultiPhaseEquilSolver::initialConditions(
        const double T, const double P, const double* const p_c)
{
    DEBUG("entering initialConditions()" << endl)
    const double alpha = 1.0e-3;
    
    // Determine which parameters have changed since the last equilibrate call
    bool composition_change = false;
    for (int i = 0; i < m_nc; ++i)
        composition_change |= (p_c[i] != mp_c[i]);
    bool temperature_change = (std::abs(T - m_T) > ms_temp_change_tol);
    bool pressure_change    = (std::abs(P - m_P) > ms_pres_change_tol);

    // Initialize input variables
    m_T  = T;
    m_P  = P;
    std::copy(p_c, p_c+m_nc, mp_c);

    // Initialize the Gibbs energy vector (regardless of tolerance)
    m_thermo.speciesGOverRT(m_T, m_P, mp_g);
    DEBUG("Species G: " << endl)
    for (int i = 0; i < m_ns; ++i)
        DEBUG(setw(20) << m_thermo.speciesName(i) << mp_g[i] << endl)
    DEBUG(endl)

    // Setup species ordering and remove "determined" species from consideration
    // (depends on temperature and constraint changes regardless of tolerance)
    bool order_change = checkForDeterminedSpecies();

    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();

    // Check to see if we can reuse the previous solution
    if (!(composition_change || temperature_change || pressure_change || m_np != 1)) {//order_change)) {
        // Can use the previous solution (just need to get g(0) which is the
        // g* from the previous solution, which is the current g because s = 1)
        std::copy(m_solution.g(), m_solution.g()+m_ns, mp_g0);

        // Nothing left to do
        return true;
    }

    // Use the tableau as temporary storage
    double* p_N = mp_tableau;
    double* p_Nbar = p_N + nsr;
    double* p_lambda = p_Nbar + npr;
    int j, jk;

    // A composition change or order change triggers a complete reinitialization
    if (composition_change || order_change)
        if (!updateMaxMinSolution()) return false;

    // Otherwise only update the MinG solution
    if (!updateMinGSolution(mp_g)) return false;

    // Form estimate of species moles from min-g and max-min vectors
    for (j = 0; j < nsr; ++j)
        p_N[j] = mp_ming[j]*(1.0-alpha) + mp_maxmin[j]*alpha;

    // Compute the initial solution which makes sure the residual is zero
    initZeroResidualSolution(p_N, p_Nbar, p_lambda);
    
    m_solution.setG(mp_g0, mp_g, 0.0);
    m_solution.setSolution(p_lambda, p_Nbar, m_B);
    
    // Make sure that we initially satisfy the phase rule
    while (m_solution.npr() > m_solution.ncr()) {
        double min_value = 0.0;
        int min_index = 1;
        for (int j = p_sizes[1]; j < p_sizes[2]; ++j)
            min_value += m_solution.y()[j];
        for (int m = 2; m < m_solution.npr(); ++m) {
            double value = 0.0;
            for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j)
                value += m_solution.y()[j];
            if (value < min_value) {
                min_value = value;
                min_index = m;
            }
        }
        m_solution.removePhase(min_index);
    }

    return true;
}

//==============================================================================

void MultiPhaseEquilSolver::initZeroResidualSolution(
	double* const p_N, double* const p_Nbar, double* const p_lambda)
{
	int npr = m_solution.npr();
	int nsr = m_solution.nsr();
	int ncr = m_solution.ncr();

	const int* const p_sizes = m_solution.sizes();
	const int* const p_sjr   = m_solution.sjr();
	const int* const p_cir   = m_solution.cir();

	// Compute initial phase moles
	for (int m = 0, j = 0; m < npr; ++m) {
		p_Nbar[m] = 0.0;
		for ( ; j < p_sizes[m+1]; ++j)
			p_Nbar[m] += p_N[j];
	}

	// Compute RHS of lambda least-squares problem (temporarily store in g0)
	for (int m = 0, j = 0; m < npr; ++m) {
		const double nbar = p_Nbar[m];
		for ( ; j < p_sizes[m+1]; ++j) {
			p_N[j] = std::log(p_N[j] / nbar);
			mp_g0[j] = p_N[j] + mp_g[p_sjr[j]];
		}
	}

	// Compute SVD of Br
	Map<VectorXd>(p_lambda, ncr) =
	    m_solution.reducedMatrix(m_B, m_Br).jacobiSvd(
	        ComputeThinU | ComputeThinV).solve(Map<VectorXd>(mp_g0, nsr));

	// Now compute the g0 which satisfies the constraints
	int jk;
	for (int j = 0; j < nsr; ++j) {
		jk = p_sjr[j];
		mp_g0[jk] = -p_N[j];
		for (int i = 0; i < ncr; ++i)
			mp_g0[jk] += m_B(jk,p_cir[i])*p_lambda[i];
	}
}

//==============================================================================

bool MultiPhaseEquilSolver::updateMinGSolution(const double* const p_g)
{
    // Use the current solution's ordering
    const int nsr = m_solution.nsr();
    const int ncr = m_solution.ncr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    
    // Build the tableau
    // 0 -f'
    double* p = mp_tableau;
    *p++ = 0;
    for (int i = 0; i < nsr; ++i)
        *p++ = -p_g[p_sjr[i]];
    
    // c -B'
    int ik;
    for (int i = 0; i < ncr; ++i) {
        ik = p_cir[i];
        *p++ = mp_c[ik];
        for (int j = 0; j < nsr; ++j)
            *p++ = -m_B(p_sjr[j], ik);
    }
    
    // 0 0
    for (int i = 0; i < nsr+1; ++i)
        *p++ = 0.0;
    
    DEBUG("Tableau for Min-G composition:" << endl)
    for (int i = 0; i < ncr+2; ++i) {
        for (int j = i*(nsr+1); j < (i+1)*(nsr+1); ++j)
            DEBUG(mp_tableau[j] << " ")
        DEBUG(endl)
    }
    DEBUG(endl)

    // Use the simplex algorithm to get the min-g solution
    int izrov [nsr];
    int iposv [ncr];
    int ret = Numerics::simplex(mp_tableau, ncr, nsr, 0, 0, izrov, iposv, 1.0E-9);
    
    // Error check
    if (ret != 0) {
        cout << "Error in computing the min-g solution in equilibrium solver!" << endl;
        if (ret < 0)
            cout << "--> no solution exists for the given problem" << endl;
        else
            cout << "--> solution is unbounded" << endl;
        return false;
    }
    
    // Unravel the solution
    for (int i = 0; i < nsr; ++i)
        mp_ming[i] = 0.0;
    for (int i = 0; i < ncr; ++i) {
        if (iposv[i] < nsr)
            mp_ming[iposv[i]] = mp_tableau[(nsr+1)*(i+1)];
        else {
            cout << "Linearly dependent in min-g!" << endl;
            return false;
        }
    }
    
    DEBUG("Successfully computed Min-G solution." << endl)

    return true;

//    cout << "min-g solution:" << endl;
//    for (int i = 0; i < nsr; ++i)
//        cout << m_thermo.speciesName(p_sjr[i]) << " " << mp_ming[i] << endl;
}

bool MultiPhaseEquilSolver::updateMaxMinSolution()
{    
    // Use the current solution's ordering
    const int nsr = m_solution.nsr();
    const int ncr = m_solution.ncr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    
    // Build the tableau
    // 0 0' 1
    double* p = mp_tableau;
    for (int i = 0; i < nsr+1; ++i)
        *p++ = 0.0;
    *p++ = 1.0;
    
    // c -B' -sum(B')
    double sum;
    int ik;
    for (int i = 0; i < ncr; ++i) {
        ik = p_cir[i];
        *p++ = mp_c[ik];
        sum = 0.0;
        for (int j = 0; j < nsr; ++j) {
            *p = -m_B(p_sjr[j], ik);
            sum += *p++;
        }
        *p++ = sum;
    }
    
    // 0 0 0
    for (int i = 0; i < nsr+2; ++i)
        *p++ = 0.0;

    DEBUG("Tableau for Max-Min composition:" << endl)
    for (int i = 0; i < ncr+2; ++i) {
        for (int j = i*(nsr+2); j < (i+1)*(nsr+2); ++j)
            DEBUG(mp_tableau[j] << " ")
        DEBUG(endl)
    }
    DEBUG(endl)

    // Use the simplex algorithm to get the max-min solution
    int izrov [nsr+1];
    int iposv [ncr];
    int ret = Numerics::simplex(mp_tableau, ncr, nsr+1, 0, 0, izrov, iposv, 1.0E-9);

    // Error check
    if (ret != 0) {
        cout << "Error in computing the max-min solution in equilibrium solver!" << endl;
        if (ret < 0)
            cout << "--> no solution exists for the given problem" << endl;
        else
            cout << "--> solution is unbounded" << endl;
        return false;
    }
    
    // Unravel the solution
    for (int i = 0; i < nsr; ++i)
        mp_maxmin[i] = *mp_tableau;
    for (int i = 0; i < ncr; ++i) {
        if (iposv[i] < nsr)
            mp_maxmin[iposv[i]] += mp_tableau[(nsr+2)*(i+1)];
    }
    
    DEBUG("Successfully computed Max-Min solution." << endl)

    return true;

//    cout << "max-min solution:" << endl;
//    for (int i = 0; i < nsr; ++i)
//        cout << m_thermo.speciesName(p_sjr[i]) << " " << mp_maxmin[i] << endl;
}

//==============================================================================

void MultiPhaseEquilSolver::formSystemMatrix(Eigen::MatrixXd& A) const
{
    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();
    const double* const p_y = m_solution.y();
    const double* const p_lnNbar = m_solution.lnNbar();
    
    A = MatrixXd::Zero(ncr+npr, ncr+npr);
    
    double temp, Nj;
    int mk, jk;
    
    for (int m = 0; m < npr; ++m) {
        mk = ncr+m;
        
        for (int j = p_sizes[m]; j < p_sizes[m+1]; ++j) {
            jk = p_sjr[j];
            Nj = p_y[j]*p_y[j];
            for (int i = 0; i < ncr; ++i) {
                temp = Nj * m_B(jk, p_cir[i]);
                for (int k = i; k < ncr; ++k)
                    A(i,k) += temp * m_B(jk, p_cir[k]);
                A(i,mk) += temp;
            }
            A(mk,mk) += Nj;
        }
        
        A(mk,mk) -= std::exp(p_lnNbar[m]);
    }
}

//==============================================================================

void MultiPhaseEquilSolver::computeResidual(Eigen::VectorXd& r) const
{
    const int npr = m_solution.npr();
    const int ncr = m_solution.ncr();
    const int nsr = m_solution.nsr();
    const int* const p_sjr = m_solution.sjr();
    const int* const p_cir = m_solution.cir();
    const int* const p_sizes = m_solution.sizes();
    const double* const p_y = m_solution.y();
    const double* const p_lnNbar = m_solution.lnNbar();
    
    for (int i = 0; i < ncr; ++i)
        r(i) = -mp_c[p_cir[i]];
    
    int jk;
    for (int j = 0; j < nsr; ++j) {
        jk = p_sjr[j];
        for (int i = 0; i < ncr; ++i)
            r(i) += m_B(jk, p_cir[i]) * p_y[j] * p_y[j];
    }
    
    for (int m = 0; m < npr; ++m) {
        jk = ncr + m;
        r(jk) = -std::exp(p_lnNbar[m]);
        for (int j = p_sizes[m]; j < p_sizes[m+1]; j++)
            r(jk) += p_y[j] * p_y[j];
    }
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
