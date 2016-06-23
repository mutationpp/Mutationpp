#include "mutation++.h"

int main()
{
    // Initialization
    Mutation::MixtureOptions opts("air5_frc_catalysis_cu");
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    const int ns   = mix.nSpecies();

    Eigen::VectorXd v_rhoi(ns);
    Eigen::VectorXd v_y(ns);
    Eigen::VectorXd v_x(ns);
    double p_static;
    double wall_temperature;

    Eigen::VectorXd v_Mw(ns);
    double mixture_Mw;
    double rho;

    Eigen::VectorXd v_wall_prod_rates(ns);

    const int iN   = mix.speciesIndex("N");
    const int iO   = mix.speciesIndex("O");
    const int iN2  = mix.speciesIndex("N2");
    const int iO2  = mix.speciesIndex("O2");
    const int iNO  = mix.speciesIndex("NO");

    // ID 38
    p_static = 5000.;         // Pa
    wall_temperature = 350.0; // K

    v_y[iN]  = 2.29420728e-06;
    v_y[iO]  = 5.98294058e-02;
    v_y[iN2] = 6.91336995e-01;
    v_y[iO2] = 2.21596802e-01;
    v_y[iNO] = 2.72264675e-02;
    double yNOp     = 8.03534124e-06;
    double yem      = 1.47019541e-10;
    v_y[iNO] += yNOp + yem;
    
    // Yi to Rhoi
    mix.convert<Mutation::Thermodynamics::Y_TO_X>(v_y.data(), v_x.data());
    mixture_Mw = v_Mw.dot(v_x);
    rho =  p_static / (Mutation::RU / mixture_Mw * wall_temperature );

    v_rhoi = v_y * rho;

    //
    int set_state_with_rhoi_T = 1;

    mix.setState(v_rhoi.data(), &wall_temperature, set_state_with_rhoi_T);
    mix.setWallState(v_rhoi.data(), &wall_temperature, set_state_with_rhoi_T);

    mix.surfaceProductionRates(v_wall_prod_rates.data());
    std::cout << v_wall_prod_rates << std::endl;

}
