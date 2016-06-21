#include "Solver.h"
#include "dfdt.h"
#include "jacobian.h"

typedef boost::numeric::ublas::vector<double> vector_type;

//=============================================================================================================

Solver::Solver(size_t l_n_args, std::string l_file_name)
        : m_n_args(l_n_args),
          s_file_name(l_file_name),
          p_mix(NULL),
          p_setup_props(NULL),
          p_setup_problem(NULL),
          p_problem(NULL),
          p_mesh(NULL),
          p_shock_relations(NULL),
          p_data_post(NULL),
          p_data_pre(NULL),
          p_output(NULL)
{ }

//=============================================================================================================

void Solver::initialize(){

    errorInsufficientNumberOfArguments();

    // Parse Input File
    p_setup_props = new SetupProperties(s_file_name);

    //Initializing Mutation++
    setMutationpp();

    // Initializing Abstract Factory Method
    p_setup_problem = SetupProblem::createFactory(p_setup_props->getProblemType(), p_setup_props->getStateModel());

    // Data Pre-Post Shock
    p_data_pre = p_setup_problem->getDataPreShock(*p_mix, p_setup_props->getFreeStreamConditions());
    p_data_post = p_setup_problem->getDataPostShock(*p_mix);

    // Mesh
    p_mesh = new Mesh(p_setup_props->getMesh());

    // PostShock
    p_shock_relations = p_setup_problem->getShockRelations(*p_mix);

    // Problem
    p_problem = p_setup_problem->getProblem(*p_mix, *p_data_post);

    // Setup Output
    p_output = new BasicOutput(*p_data_post);

    // Deallocate Some Memory
    p_setup_props->cleanMeshInfo();
    p_setup_props->cleanFreeStreamConditionsInfo();

    // Apply Shock Relations
    p_shock_relations->applyShockRelations(*p_data_pre, *p_data_post);

}

//=============================================================================================================

void Solver::solve(){

    // Call solver to solve the problem
    vector_type x = p_data_post->getInitialState();
    size_t num_of_steps = integrate_const(make_dense_output<rosenbrock4<double> >(1.0e-8, 1.0e-8),
                                          std::make_pair(dfdt(*p_problem), jacobian(*p_problem)),
                                          x,
                                          p_mesh->getXinitial(), p_mesh->getXfinal(), p_mesh->getdX(),
                                          *p_output);

}

//=============================================================================================================

Solver::~Solver(){
    if (p_mix != NULL) delete p_mix;
    if (p_setup_props != NULL) delete p_setup_props;
    if (p_setup_problem != NULL) delete p_setup_problem;
    if (p_problem != NULL) delete p_problem;
    if (p_shock_relations != NULL) delete p_shock_relations;
    if (p_data_pre != NULL) delete p_data_pre;
    if (p_data_post != NULL) delete p_data_post;
    if (p_output != NULL) delete p_output;
}

//=============================================================================================================

inline void Solver::setMutationpp(){
    Mutation::MixtureOptions opts(p_setup_props->getMixture());
    opts.setStateModel(p_setup_props->getStateModel());
    opts.setThermodynamicDatabase(p_setup_props->getThermoDB());
    p_mix = new Mutation::Mixture(opts);
}

//=============================================================================================================

inline void Solver::errorInsufficientNumberOfArguments(){
    if ( m_n_args < 2 ) {
        std::cout << "Error: The name of the input file for the solver should be provided!" << std::endl;
        exit(1);
    }
}

