/**
 * @file mppshock.cpp
 *
 * @brief Utility which provides the post shock equilibrium state
 * for a given mixture.
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

#include "mutation++.h"

#include <string>

using namespace std;
using namespace Eigen;
using namespace Mutation;

class State
{
public:
    explicit State(const double& p, const double& u, const double& T) :
        m_P(p), m_V(u), m_T(T) {}
    State (const State&) = delete;
    State (State&&) = default;
    State& operator=( const State &) = default;

    double getP() const { return m_P; }
    double getV() const { return m_V; }
    double getT() const { return m_T; }

private:
    const double m_P;
    const double m_V;
    const double m_T;
};

// This class computes the 1T Rankine-Hugoniot relations.
// It will be abstracted to treat additional state models.
class RankineHugoniot
{
public:
    explicit RankineHugoniot(Mixture& mix) : m_mix(mix) {}

    // void setPreShockState(const State& pre_shock)
    State computeRH(const State& pre_shock)
    {
        const double P_pre = pre_shock.getP();
        const double V_pre = pre_shock.getV();
        const double T_pre = pre_shock.getT();

        m_mix.setState(&P_pre, &T_pre, set_state_PT);

        // Thermodynamic state based on the pre-shock conditions
        const double gamma = m_mix.mixtureEquilibriumGamma();
        const double gp1 = gamma + 1.;
        const double gm1 = gamma - 1.;
        const double a_speed = m_mix.equilibriumSoundSpeed();
        const double Ms = V_pre/a_speed;
        const double Ms2 = Ms*Ms;

        //  Rankine-Hugoniot relations
        const double P_post = P_pre*(1. + 2.*gamma/gp1*(Ms2 - 1.));
        const double V_post = V_pre - a_speed*(2./gp1*(Ms - 1./Ms));
        const double T_post = T_pre*(1. + (2.*gm1/(gp1*gp1)*((gamma*Ms2 + 1.)/Ms2)*(Ms2 - 1.)));

        return move(State(P_post, V_post, T_post));
    }

    void setPreShockState(
        const double& P_pre, const double& V_pre, const double T_pre)
        // const ArrayXd& v_ye) // Add me in the state
    {
        m_mix.setState(&P_pre, &T_pre, set_state_PT);

        // Thermodynamic state based on the pre-shock conditions
        const double gamma = m_mix.mixtureEquilibriumGamma();
        const double gp1 = gamma + 1.;
        const double gm1 = gamma - 1.;
        const double a_speed = m_mix.equilibriumSoundSpeed();
        const double Ms = V_pre/a_speed;
        const double Ms2 = Ms*Ms;

        //  Rankine-Hugoniot relations
        m_P_post = P_pre*(1. + 2.*gamma/gp1*(Ms2 - 1.));
        m_V_post = V_pre - a_speed*(2./gp1*(Ms - 1./Ms));
        m_T_post = T_pre*(1. + (2.*gm1/(gp1*gp1)*((gamma*Ms2 + 1.)/Ms2)*(Ms2 - 1.)));
    }

    double postShockP() const { return m_P_post; }
    double postShockV() const { return m_V_post; }
    double postShockT() const { return m_T_post; }

private:
    Mixture& m_mix;
    const int set_state_PT = 1;

    double m_P_post = 1.;
    double m_V_post = 1.;
    double m_T_post = 1.;
};

class ShockOptions {
public:
    explicit ShockOptions(const vector<string>& args)
    {
        // Error Checking
        // for ( int i = 1; i < n_args; i++) v_args.push_back(args[i]);

        // if (v_args.size() == 0)
    // if (optionExists(argc, argv, "-h") || optionExists(argc, argv, "--help"))
    //     printHelpMessage(argv[0]);
        printHelpMessage();

    // Get the temperature
    // if (optionExists(argc, argv, "-T")) {
    //     if (!parseRange(
    //         getOption(argc, argv, "-T"), opts.T1, opts.T2, opts.dT)) {
    //         cout << "Bad format for temperature range!" << endl;
    //         printHelpMessage(argv[0]);
    //     }

        m_P_pre = 1000.;
        m_V_pre = 10000.;
        m_T_pre = 300.;

        m_mix_name = "air_5";

    }


    double preShockP() const { return m_P_pre; }
    double preShockV() const { return m_V_pre; }
    double preShockT() const { return m_T_pre; }
    const string& mixtureName() const { return m_mix_name; }
    int postShockEquilMethod() const { return m_equil_method; }

private:
//     bool optionExists(int argc, char** argv, const std::string& option)
// {
//     return (std::find(argv, argv+argc, option) != argv+argc);
// }

    void printHelpMessage()
    {
        cout << "mppshock is work in progress.\n";
        cout << "Algorithms (Bisection/Newton) are implemented.\n";
        cout << "Parser and refactoring will be finalized soon.\n";
        cout << "mppshock -T 300 -P 1000 -V 10000 -m 0 air_5\n";
    }

private:
    double m_P_pre;
    double m_V_pre;
    double m_T_pre;

    string m_mix_name;
    const int m_equil_method = 0;
};

enum EquilibriumMethod { bisection = 0, newton = 1 };

class PostShockEquilibrium {
public:
    virtual ~PostShockEquilibrium() = default;
    virtual State computePostShockEquilibrium(const State& pre_shock) = 0;
};

class PostShockEquilibriumBisection : public PostShockEquilibrium {
public:
    explicit PostShockEquilibriumBisection(Mixture& mix) : m_mix(mix) {}

    State computePostShockEquilibrium(const State& pre_shock)
    {
        const double p_pre = pre_shock.getP();
        const double u_pre = pre_shock.getV();
        const double T_pre = pre_shock.getT();

        m_mix.setState(&p_pre, &T_pre, set_state_PT);
        const double rho_pre = m_mix.density();
        const double h_pre = m_mix.mixtureHMass();
        const double mdot = rho_pre*u_pre;

        // Setting up solver
        double r = 0.1; // Initial density ratio
        double p_eq, T_eq;

        // Solving LTE based on bisection method
        double err_r = 1.;
        while (err_r > m_tol) {
            p_eq = p_pre + mdot*u_pre*(1. - r);
            double h_eq = h_pre + 0.5*u_pre*u_pre*(1 - r*r);

            double Tl = T_pre;
            double Tu = 60000.;

            T_eq = 0.5*(Tl + Tu);

            double err_T = 1.;
            while (err_T > m_tol) {
                // Lower
                m_mix.setState(&p_eq, &Tl, set_state_PT);
                const double fl = -h_eq + m_mix.mixtureHMass();

                // Upper
                m_mix.setState(&p_eq, &Tu, set_state_PT);
                const double fu = -h_eq + m_mix.mixtureHMass();

                // Solution
                m_mix.setState(&p_eq, &T_eq, set_state_PT);
                const double f  = -h_eq + m_mix.mixtureHMass();

                double T_ref = T_eq;
                if ((f*fu) < 0.){
                    Tl = T_eq;
                    T_eq  = .5*(Tl + Tu);
                } else if ((f*fl) < 0.) {
                    Tu = T_eq;
                    T_eq  = .5*(Tl + Tu);
                }
                err_T = abs(T_eq-T_ref)/T_ref;
            }

            // Density update
            const double rho_eq = m_mix.density();

            double r_old = r;
            r = rho_pre / rho_eq;
            err_r = abs(r - r_old)/r_old;
        }

        m_mix.setState(&p_eq, &T_eq, set_state_PT);
        p_eq = p_pre + mdot*u_pre*(1. - r);
        const double rho_eq = m_mix.density();
        const double u_eq = mdot / rho_eq;

        return State(p_eq, u_eq, T_eq);
    }

private:
    Mixture& m_mix;

    const int set_state_PT = 1;
    const double m_tol = 1.e-12;
    const double m_T_high = 90000.;
};

class PostShockEquilibriumNewton : public PostShockEquilibrium {
public:
    explicit PostShockEquilibriumNewton(Mixture& mix) : m_mix(mix) {}

    State computePostShockEquilibrium(const State& pre_shock)
    {
        const double p_pre = pre_shock.getP();
        const double u_pre = pre_shock.getV();
        const double T_pre = pre_shock.getT();

        m_mix.setState(&p_pre, &T_pre, set_state_PT);
        const double rho_pre = m_mix.density();
        const double h_pre = m_mix.mixtureHMass();

        // Conserved quanities
        const double mdot = rho_pre*u_pre;
        const double momentum = mdot*u_pre+p_pre;
        const double E = h_pre + .5*u_pre*u_pre;

        // Initial Guess
        double ratio = .99;
        double T_eq = T_pre;
        double p_eq = p_pre;
        double u_eq = u_pre;

        // Outer loop for Mass/Momentum
        double dR = 1.;
        while (dR > m_tol) {
            p_eq = p_pre + mdot*u_pre*(1.-ratio);
            double rho_eq = rho_pre/ratio;
            double RHS  = h_pre + .5*u_pre*u_pre*(1. - ratio*ratio);

            // Inner loop for Energy
            double dT = 1.;
            while (abs(dT) > m_tol) {
                m_mix.setState(&p_eq, &T_eq, set_state_PT);

                const double h_eq = m_mix.mixtureHMass();
                const double cp_eq = m_mix.mixtureFrozenCpMass();

                dT = (h_eq-RHS)/cp_eq;
                T_eq = T_eq - dT;
            }

            m_mix.setState(&p_eq, &T_eq, set_state_PT);
            rho_eq = m_mix.density();

            // Density ratio update and residual
            const double old_ratio = ratio;
            ratio = rho_pre/rho_eq;
            dR = abs(ratio - old_ratio)/ratio;
        }
        u_eq = u_pre*ratio;

        return move(State(p_eq, u_eq, T_eq));
    }

private:
    Mixture& m_mix;

    const int set_state_PT = 1;
    const double m_tol = 1.e-10;
};

class PrintResults {
public:
    PrintResults(
        const State& free_stream,
        const State& RH,
        const State& equil) :
            m_fs(free_stream), m_RH(RH), m_eq(equil){ }

    void output() noexcept {
        cout.precision(m_precision);

        cout << endl;
        cout
            << tab << setw(width) << " "
            << setw(width) << "Pressure (Pa)"
            << setw(width) << "Velocity (m/s)"
            << setw(width) << "Temperature (K)" << endl;;
        cout
            << tab << setw(width) << left << "Free Stream:" << fixed
            << setw(width) << right << m_fs.getP()
            << setw(width) << right << m_fs.getV()
            << setw(width) << right << m_fs.getT() << endl;
        cout
            << tab << setw(width) << left << "Rankine Hugoniot:" << fixed
            << setw(width) << right << m_RH.getP()
            << setw(width) << right << m_RH.getV()
            << setw(width) << right << m_RH.getT() << endl;
        cout
            << tab << setw(width) << left << "Equilibrated:" << fixed
            << setw(width) << right << m_eq.getP()
            << setw(width) << right << m_eq.getV()
            << setw(width) << right << m_eq.getT() << endl;
        cout << endl;
        }

private:
    const State& m_fs;
    const State& m_RH;
    const State& m_eq;

    const int m_precision = 1;
    const string tab = "    ";
    const int width = 20;
};

int main(int argc, char** argv)
{
    vector<string> args(argv, argv + argc);
    ShockOptions mppshock_opts(args);

    MixtureOptions mix_opts(mppshock_opts.mixtureName());

    const string state_model = "Equil";
    mix_opts.setStateModel(state_model);
    Mixture mix(mix_opts);

    State pre_shock(
        mppshock_opts.preShockP(),
        mppshock_opts.preShockV(),
        mppshock_opts.preShockT());

    RankineHugoniot rh(mix);
    State post_shock_RH = rh.computeRH(pre_shock);

    const int method = newton;
    unique_ptr<PostShockEquilibrium> p_method = nullptr;
    if (method == bisection) {
        p_method = unique_ptr<PostShockEquilibrium> {
            new PostShockEquilibriumBisection(mix) };
    } else if (method == newton) {
        p_method = unique_ptr<PostShockEquilibrium> {
            new PostShockEquilibriumNewton(mix) };
    }

    State post_shock_equil = p_method->computePostShockEquilibrium(pre_shock);

    PrintResults print_results(pre_shock, post_shock_RH, post_shock_equil);
    print_results.output();

    return 0;
}


