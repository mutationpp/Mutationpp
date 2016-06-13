#include "SurfaceBalanceSolver.h"

namespace Mutation {
    namespace GasSurfaceInteraction {

//======================================================================================

        class SurfaceBalanceSolverEquilibrium : public SurfaceBalanceSolver
{

public:
    SurfaceBalanceSolverEquilibrium(ARGS args)
        :m_thermo(args.s_thermo),
         m_surf_props(args.s_surf_props),
         m_wall_state(args.s_wall_state),
         mp_diff_vel_calc(NULL),
         mp_mass_blowing_rate(NULL),
         equil_object (NULL),
         m_ns(args.s_thermo.nGas()),
         v_rhoi(args.s_thermo.nGas()),
         set_state_with_rhoi_T(1),
         v_elementalMolePyro(args.s_thermo.nElements()),
         v_elementalMoleChar(args.s_thermo.nElements()),
         v_elementalMoleEdge(args.s_thermo.nElements()),
         v_elementalMassPyro(args.s_thermo.nElements()),
         v_elementalMassChar(args.s_thermo.nElements()),
         v_elementalMassEdge(args.s_thermo.nElements()),
         v_elementalMassWall(args.s_thermo.nElements()),
         v_elementalMoleWall(args.s_thermo.nElements()),
         surface_constraint(m_surf_props.surfaceConstraint()),
         Bprime_pyro(m_surf_props.BprimePyro())
    { 

        Mutation::Utilities::IO::XmlElement::const_iterator iter_prod_terms = args.s_node_prod_terms.begin();

//@toremove @BD        std::basic_ostream(*iter_prod_terms);

        DataWallProductionTerms l_data_wall_production_terms = {m_thermo, 
                                                                args.s_transport,
                                                                args.s_gsi_mechanism, 
                                                                *iter_prod_terms, 
                                                                m_surf_props, 
                                                                m_wall_state};
        for(; iter_prod_terms != args.s_node_prod_terms.end(); ++iter_prod_terms){
            addSurfaceProductionTerm(Mutation::Utilities::Config::Factory<WallProductionTerms>::create(
                  iter_prod_terms->tag(), l_data_wall_production_terms));
        }

        // DiffusionVelocityCalculator
        mp_diff_vel_calc = new DiffusionVelocityCalculator(m_thermo, args.s_transport); //@BD Remove? or set it to Null

        // Setting 
        m_surf_props.getSurfaceCompositions(v_elementalMolePyro, v_elementalMoleChar);

        // Thermodynamic Object for equilibrium
        equil_object = new Mutation::Thermodynamics::Thermodynamics(m_surf_props.surfaceSpecies(), "NASA-9", "ChemNonEq1T"); // check if ChemNonEq works for equilbrate

        GetCondesedElementIndex();

        if(surface_constraint == true){
        addConstraint();}

        // Mass Blowing Rate?

    }

    //======================================================================================

    ~SurfaceBalanceSolverEquilibrium()
    {
        Solid_index.clear();
        v_speciesMoleChar.clear();
        if (mp_diff_vel_calc != NULL){ delete mp_diff_vel_calc; }
        if (mp_mass_blowing_rate != NULL){ delete mp_mass_blowing_rate; }
        if (equil_object != NULL){ delete equil_object; }
        if (surface_constraint != false){ m_thermo.clearEquilibriumContraints(); }
    }

    //======================================================================================

    Eigen::VectorXd& computeGSIProductionRate()
    {
        std::cerr << "Wall treated in chemical equilibrium" << std::endl;
        exit(1);
    }

    //======================================================================================

    void setDiffusionModel(
        const Eigen::VectorXd& lv_mole_frac_edge, 
        const double& l_dx)
    {
        mp_diff_vel_calc->setDiffusionModel( lv_mole_frac_edge, l_dx );
    }

    //======================================================================================

    void solveSurfaceBalance()
    {
        m_Twall = m_wall_state.getWallT()(0);
        m_Pwall = m_wall_state.getWallP();

        computeEquilibriumSurface( );

        m_wall_state.setWallState(v_rhoi.data(), &m_Twall, set_state_with_rhoi_T);
    }
    //======================================================================================

     void getBprimeCondensedSpecies(std::vector<std::string>& s_condensed_species)
     {
         for (int i = equil_object->nGas(); i < equil_object->nSpecies(); ++i ){
             s_condensed_species.push_back(equil_object->speciesName(i));}
      }

   //======================================================================================

    void getBprimeParameters(double & Bprime_char, std::vector<double>& CondensedMoleFrac)
    {
            CondensedMoleFrac = v_speciesMoleChar;
            Bprime_char = ComputeBprimeChar();
            v_speciesMoleChar.clear();
    }

   //======================================================================================

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;
    Mutation::Thermodynamics::Thermodynamics* equil_object;

    SurfaceProperties& m_surf_props;
    WallState& m_wall_state;

    DiffusionVelocityCalculator* mp_diff_vel_calc;
    MassBlowingRate* mp_mass_blowing_rate;

    std::vector<WallProductionTerms*> v_surf_prod;

    // VARIABLES FOR SOLVER
    const size_t m_ns;
    double m_Twall;
    double m_Pwall;
    Eigen::VectorXd v_rhoi;
    Eigen::VectorXd v_elementalMoleWall;
    Eigen::VectorXd v_elementalMassWall;
    Eigen::VectorXd v_elementalMassEdge;
    Eigen::VectorXd v_elementalMassPyro;
    Eigen::VectorXd v_elementalMassChar;
    Eigen::VectorXd v_elementalMoleEdge;
    Eigen::VectorXd v_elementalMolePyro;
    Eigen::VectorXd v_elementalMoleChar;

    std::vector<double> v_speciesMoleChar;

    std::vector<std::string> Solid_index;
    const double  Bprime_pyro;
    bool surface_constraint;

    const int set_state_with_rhoi_T;

    //======================================================================================

    void computeMoleFracfromPartialDens( Eigen::VectorXd& lv_rhoi, Eigen::VectorXd& lv_xi ){

        m_thermo.setState( &lv_rhoi(0), &m_Twall, set_state_with_rhoi_T );
        for ( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){
            lv_xi(i_ns) = m_thermo.X()[i_ns];
        }

    }

    //======================================================================================
    /**
     * @todo Replace this function with the default in Thermodynamics
     */
    void computePartialDensfromMoleFrac( Eigen::VectorXd& lv_xi, Eigen::VectorXd& lv_rhoi )
    {
        for( int i_ns = 0 ; i_ns < m_ns ; i_ns++ ){ // @BD Remove the loop
            lv_rhoi(i_ns) = lv_xi(i_ns) * m_Pwall * m_thermo.speciesMw(i_ns) / ( m_Twall * Mutation::RU );
        }

    }

    //======================================================================================

    void computeMassFracfromMoleFrac(Eigen::VectorXd& lv_xi, Eigen::VectorXd& lv_yi)
    {
        m_thermo.convert<Mutation::Thermodynamics::X_TO_Y>(lv_xi.data(), lv_yi.data());
    }

    //======================================================================================

    void computeMoleFracfromMassFrac(Eigen::VectorXd& lv_yi, Eigen::VectorXd& lv_xi)
    {
        m_thermo.convert<Mutation::Thermodynamics::Y_TO_X>(lv_yi.data(), lv_xi.data());
    }

    //======================================================================================

    void computeElementalMassFracFromMassFrac(Eigen::VectorXd& lv_yi, Eigen::VectorXd& lv_ye)
    {
        m_thermo.convert<Mutation::Thermodynamics::Y_TO_YE>(lv_yi.data(), lv_ye.data());
    }
    //======================================================================================

    void computeElementalMoleFracFromElementalMassFrac(Eigen::VectorXd& lv_ye, Eigen::VectorXd& lv_xe)
    {
       m_thermo.convert<Mutation::Thermodynamics::YE_TO_XE>(lv_ye.data(), lv_xe.data());
    }
    //======================================================================================

    void computeElementalMassFracFromElementalMoleFrac(Eigen::VectorXd& lv_xe, Eigen::VectorXd& lv_ye)
    {
       m_thermo.convert<Mutation::Thermodynamics::XE_TO_YE>(lv_xe.data(), lv_ye.data());
    }
    //======================================================================================

    void addConstraint()
    {
        Eigen::VectorXd v_constraintVector(equil_object->nSpecies());
        size_t m_nConstraints = 0;


        int m_nCondensedElements = Solid_index.size();

        if(Solid_index.size()>2){
            m_nConstraints = m_nCondensedElements-2;
        }


        for (int j = 0; j<m_nConstraints; ++j){
            double alpha = v_elementalMoleChar(equil_object->elementIndex(Solid_index[j]))
                    /v_elementalMoleChar(equil_object->elementIndex(Solid_index[j+1]));

            for(int i=0; i< equil_object->nGas(); ++i){
                v_constraintVector(i)= 0;
            }

            for (int i = equil_object->nGas(); i <equil_object->nSpecies(); ++i){
                v_constraintVector(i) = equil_object->elementMatrix()(i,equil_object->elementIndex(Solid_index[j]))
                        -alpha*equil_object->elementMatrix()(i,equil_object->elementIndex(Solid_index[j+1]));
            }

            equil_object->addEquilibriumConstraint(&v_constraintVector(0));
        }
    }
    //======================================================================================

    void computeEquilibriumSurface()
    {
        Eigen::VectorXd v_speciesMoleWall (equil_object->nSpecies());
        Eigen::VectorXd v_speciesMassWall (equil_object->nSpecies());

        double Bprime_char_ini = 10000.0;

        for(int i = 0; i< m_thermo.nElements(); ++i){
            v_elementalMoleEdge(i) = m_thermo.getDefaultComposition(i);}

        computeElementalMassFracFromElementalMoleFrac(v_elementalMoleEdge, v_elementalMassEdge);
        computeElementalMassFracFromElementalMoleFrac(v_elementalMolePyro, v_elementalMassPyro);
        computeElementalMassFracFromElementalMoleFrac(v_elementalMoleChar, v_elementalMassChar);

        double sum = 0.0E0;

        for(int i = 0; i<equil_object->nElements(); ++i){
            v_elementalMassWall(i) = (v_elementalMassEdge(i)+Bprime_pyro*v_elementalMassPyro(i)+Bprime_char_ini*v_elementalMassChar(i))
                                    /(1.+Bprime_pyro+Bprime_char_ini);
            sum += v_elementalMassWall(i);}

        for (int i = 0; i < equil_object->nElements(); ++i){
            v_elementalMassWall(i) /= sum;
        }

        computeElementalMoleFracFromElementalMassFrac(v_elementalMassWall,v_elementalMoleWall);


        equil_object->equilibriumComposition(m_Twall, m_Pwall, &v_elementalMoleWall(0), &v_speciesMoleWall(0),Mutation::Thermodynamics::GLOBAL);

        for (int i = equil_object->nGas(); i < equil_object->nSpecies(); ++i ){
            v_speciesMoleChar.push_back(v_speciesMoleWall(i));}


       sum = 0.0E0;

        for (int i = 0; i < equil_object->nGas(); ++i){
            sum += v_speciesMoleWall(i);
         }


        for (int i = 0; i < equil_object->nGas(); ++i){
            v_speciesMoleWall(i) /= sum;
        }

        computePartialDensfromMoleFrac(v_speciesMoleWall, v_rhoi);
        computeMassFracfromMoleFrac(v_speciesMoleWall, v_speciesMassWall);
        computeElementalMassFracFromMassFrac(v_speciesMassWall,v_elementalMassWall);
        computeElementalMoleFracFromElementalMassFrac(v_elementalMassWall, v_elementalMoleWall);

    }

    //======================================================================================

    double ComputeBprimeChar(){

        int ic = equil_object->elementIndex(Solid_index[0]);
        double Bprime_char;

        Bprime_char = (v_elementalMassEdge(ic) + Bprime_pyro*v_elementalMassPyro(ic) -
                       v_elementalMassWall(ic)*(1.0 + Bprime_pyro)) / (v_elementalMassWall(ic) - v_elementalMassChar(ic));
        Bprime_char = std::max(Bprime_char, 0.0);

        return Bprime_char;

    }

    //======================================================================================

    void GetCondesedElementIndex()
    {
        for(int i = 0; i<equil_object->nElements(); ++i){
            if (v_elementalMoleChar(i)>0){
                Solid_index.push_back(equil_object->elementName(i));
            }
        }
    }
    //======================================================================================

    void errorWallStateNotSet() const {

        if ( m_wall_state.isWallStateSet() == 0 ){
            std::cerr << "The wall state must have been set!" << std::endl; /** @todo better error */
            exit(1);
        }

    }

    //======================================================================================

    void addSurfaceProductionTerm( WallProductionTerms* lp_wall_prod_term ){
        v_surf_prod.push_back( lp_wall_prod_term );
    }

    //======================================================================================

}; // class SurfaceBalanceSolverEquilibrium

Mutation::Utilities::Config::ObjectProvider<SurfaceBalanceSolverEquilibrium, SurfaceBalanceSolver> surface_balance_solver_equilibrium("equilibrium");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
