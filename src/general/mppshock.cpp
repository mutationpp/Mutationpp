/**
 * @file mppshock.cpp
 *
 * @brief Utility which provides the post shock equilibrium state
 * for a given mixture.
 */

/*
 * Copyright 2020 von Karman Institute for Fluid Dynamics (VKI)
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

// Class storing the current gas state in pressure, temperature and velocity.
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
    explicit RankineHugoniot(Mutation::Mixture& mix) : m_mix(mix) {}

    State computeRH(const State& pre_shock) noexcept
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

        return std::move(State(P_post, V_post, T_post));
    }

    double postShockP() const { return m_P_post; }
    double postShockV() const { return m_V_post; }
    double postShockT() const { return m_T_post; }

private:
    Mutation::Mixture& m_mix;
    const int set_state_PT = 1;

    double m_P_post = 1.;
    double m_V_post = 1.;
    double m_T_post = 1.;
};

enum EquilibriumMethod { bisection = 0, newton = 1 };

// Class responsible for parsing and storing the input options.
class ShockOptions {
public:
    explicit ShockOptions(const std::vector<std::string>& args)
    {
        if ( args.size() < n_args_min ) printHelpMessage();

        // Printing help
        if (std::find(args.begin(), args.end(), "-h") != args.end()) printHelpMessage();
        if (std::find(args.begin(), args.end(), "--help") != args.end()) printHelpMessage();

        // Parsing options
        const auto p_begin = args.begin();
        const auto p_end = args.end();

        const auto opts_p = find(p_begin, p_end, "-P");
        if (opts_p != p_end && opts_p+1 != p_end)
            try { m_P_pre  = stod(args[opts_p-p_begin+1]);
            } catch(std::exception& error) {
                std::cerr << "Error in option -P. The input value is not a double!\n";
            }

        const auto opts_V = find(p_begin, p_end, "-V");
        if (opts_V != p_end && opts_V+1 != p_end)
            try { m_V_pre  = stod(args[opts_V-p_begin+1]);
            } catch(std::exception& error) {
                std::cerr << "Error in option -V. The input value is not a double!\n";
            }

        const auto opts_T = find(args.begin(), args.end(), "-T");
        if (opts_T != p_end && opts_T+1 != p_end)
            try { m_T_pre  = stod(args[opts_T-p_begin+1]);
            } catch(std::exception& error) {
                std::cerr << "Error in option -T. The input value is not a double!\n";
            }

        const auto opts_m = find(args.begin(), args.end(), "-m");
        if (opts_m != p_end && opts_m+1 != p_end)
            try { m_equil_method  = stoi(args[opts_m-p_begin+1]);
            } catch(std::exception& error) {
                std::cerr << "Error in option -m. The input value is not an int!\n";
            }

        if (m_equil_method > 1) {
            std::cerr << "Unknown option " << m_equil_method
                 << " for equilibrium method.\n";
            printHelpMessage();
        }

        // Parsing mixture
        m_mix_name = args.back();
    }

    double preShockP() const { return m_P_pre; }
    double preShockV() const { return m_V_pre; }
    double preShockT() const { return m_T_pre; }
    const std::string& mixtureName() const { return m_mix_name; }
    int postShockEquilMethod() const { return m_equil_method; }

private:
    void printHelpMessage() {
        std::cout << "Computes the Rankine-Hugoniot and post-shock equilibrium \n"
             << "given the pre-shock state for a give mixture.\n";
        std::cout << "Usage:   mppshock [OPTIONS] mixture\n";
        std::cout << tab << "-h, --help          prints this help message\n";
        std::cout << tab << "-P                  Free stream pressure range in Pa (default: 1 atm)\n";
        std::cout << tab << "-T                  Free stream temperature range in K (default: 300 K)\n";
        std::cout << tab << "-V                  Free stream velocity range in m/s (default: 10000 m/s)\n";
        std::cout << tab << "-m                  Solution algorithm; 0: Bisection - 1: Newton-Raphson (default)\n";
        std::cout << "Example: mppshock -P 1000 -V 10000 -T 300 -m 0 air_5\n";
        exit(0);
    }

private:
    double m_P_pre = Mutation::ONEATM;
    double m_V_pre = 10000.;
    double m_T_pre = 300.;
    int m_equil_method = newton;

    std::string m_mix_name;

    const std::string tab = "    ";
    const size_t n_args_min = 2;
};

// Abstract class responsible for the computation of the post shock equilibrium.
// The solid classes implement different algorithms for the equilibrium calculation.
class PostShockEquilibrium {
public:
    // Virtual destructor
    virtual ~PostShockEquilibrium() = default;

    // This function is the core of the class. It takes as input the pre-shock state
    // and returns the relaxed equilibrated post-shock state.
    virtual State computePostShockEquilibrium(const State& pre_shock) = 0;
};

// Computes the post shock equilibrium based on a bisection algorithm. This
// method should be used only for validation and in stiff cases where the Newton
// solver fails to converge.
class PostShockEquilibriumBisection : public PostShockEquilibrium {
public:
    explicit PostShockEquilibriumBisection(Mutation::Mixture& mix) : m_mix(mix) {}

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
            double Tu = m_T_high;

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
    Mutation::Mixture& m_mix;

    const int set_state_PT = 1;
    const double m_tol = 1.e-12;
    const double m_T_high = 100000.;
};

// Computing the post shock equilibrium based on a Newton-Raphson loop. Even
// though this methodology is less robust, it is very efficient.
class PostShockEquilibriumNewton : public PostShockEquilibrium {
public:
    explicit PostShockEquilibriumNewton(Mutation::Mixture& mix) : m_mix(mix) {}

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

        return std::move(State(p_eq, u_eq, T_eq));
    }

private:
    Mutation::Mixture& m_mix;

    const int set_state_PT = 1;
    const double m_tol = 1.e-10;
};

// Helper class that prints
class PrintResults {
public:
    explicit PrintResults(
        const State& free_stream,
        const State& RH,
        const State& equil) :
            m_fs(free_stream), m_RH(RH), m_eq(equil) {}

    void output() noexcept {
        std::cout.precision(m_precision);

        std::cout << std::endl;
        std::cout
            << tab << std::setw(width) << " "
            << std::setw(width) << "Pressure (Pa)"
            << std::setw(width) << "Velocity (m/s)"
            << std::setw(width) << "Temperature (K)" << std::endl;;
        std::cout
            << tab << std::setw(width) << std::left
            << "Free Stream:" << std::fixed
            << std::setw(width) << std::right << m_fs.getP()
            << std::setw(width) << std::right << m_fs.getV()
            << std::setw(width) << std::right << m_fs.getT() << std::endl;
        std::cout
            << tab << std::setw(width) << std::left
            << "Rankine-Hugoniot:" << std::fixed
            << std::setw(width) << std::right << m_RH.getP()
            << std::setw(width) << std::right << m_RH.getV()
            << std::setw(width) << std::right << m_RH.getT() << std::endl;
        std::cout
            << tab << std::setw(width) << std::left
            << "Equilibrated:" << std::fixed
            << std::setw(width) << std::right << m_eq.getP()
            << std::setw(width) << std::right << m_eq.getV()
            << std::setw(width) << std::right << m_eq.getT() << std::endl;
        std::cout << std::endl;
    }

private:
    const State& m_fs;
    const State& m_RH;
    const State& m_eq;

    const int m_precision = 1;
    const std::string tab = "    ";
    const int width = 20;
};

/**
 * @page mppshock Mutation++ Post Shock Equilibrium Solver (mppshock)
 *
 * __Usage:__
 *
 * mppshock [OPTIONS] mixture
 *
 * Compute equilibrium properties for mixture over a set of temperatures and
 * pressures using the Mutation++ library.  Use
 *
 *     mppshock -h
 *
 * for a full list of options.
 */

int main(int argc, char** argv)
{
    std::vector<std::string> args(argv, argv + argc);
    ShockOptions mppshock_opts(args);

    Mutation::MixtureOptions mix_opts(mppshock_opts.mixtureName());

    const std::string state_model = "Equil";
    mix_opts.setStateModel(state_model);
    Mutation::Mixture mix(mix_opts);

    State pre_shock(
        mppshock_opts.preShockP(),
        mppshock_opts.preShockV(),
        mppshock_opts.preShockT());

    RankineHugoniot rh(mix);
    State post_shock_RH = rh.computeRH(pre_shock);

    const int method = mppshock_opts.postShockEquilMethod();
    std::unique_ptr<PostShockEquilibrium> p_method = nullptr;
    if (method == bisection) {
        p_method = std::unique_ptr<PostShockEquilibrium> {
            new PostShockEquilibriumBisection(mix) };
    } else if (method == newton) {
        p_method = std::unique_ptr<PostShockEquilibrium> {
            new PostShockEquilibriumNewton(mix) };
    }

    State post_shock_equil = p_method->computePostShockEquilibrium(pre_shock);

    PrintResults print_results(pre_shock, post_shock_RH, post_shock_equil);
    print_results.output();

    return 0;
}


