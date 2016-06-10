// This file contains the Lagrangian Solver headers




// =================================================================================
//    current solution
// I pass this to both the EquationsToSolve and Jacobian objects, so that I don't 
// have to fix their internal parameters. They just link here.
class ExternalField
{
  public:
    std::vector<double> x_ext, u_ext, rho_ext;

    double x_now[2]; // <--- Questi li aggiorno passo dopo passo
    double u_now[2];
    double rho_now[2];
};




// =================================================================================
// function implementing ODEs in Lagrangian form

void LagrProblem( const vec_type &x , vec_type &dxdt , double t, \
                  Mutation::Mixture* mix, ExternalField* extField )
{
// This function creates the ODE system in the given state "x".
// The input mixture basically stores the info on how many nuclei for each species are 
// present and is set to the current state "x" via "setState()".
// Values for thermodynamics quantities and chemical rates are then extracted and used 
// to create the system (compute derivatives).
// Info on the outer field is given via the extField pointer, that stores many things, 
// among which the positions, velocities and densities of the two streamline points 
// among which we are integrating.
//
// Function arguments:
// x        -> state of the system: nSpecies mass fractions + 1 Temperature
// dxdt     -> derivative 
// t        -> current step, abscissa
// mix      -> pointer to mixture from which picking parameters such as chemical rates
// extField -> pointer to the externally computed flow field
// Du2Dt    -> lagrangian derivative of the squared velocity
  size_t ii;                       // internal counter
  size_t ns = mix->nSpecies();     // ..shorthand

  double x_1     = extField->x_now[1];  // some more shorthands..
  double x_0     = extField->x_now[0];
  double u_1     = extField->u_now[1];
  double u_0     = extField->u_now[0];
  double rho_1   = extField->rho_now[1];
  double rho_0   = extField->rho_now[0];

  double u_now   = u_0;   // should interpolate... but for now keep it simple..
  double rho_now = rho_0; // idem, keep it simple for now: oder 0 interpolation
// ------------   qualcosa tipo u_now = u_ext[0] + (u_ext[1] - u_ext[0]) / (x_ext[1]-x_ext[0]) * x
// -------------   MA METTENDO UN CHECK CHE I DUE PUNTI NON SIANO COINCIDENTI, senno' divido per 0!!

  // Some arrays that store per-species values
  double * Wdot = new double[ns];
  double * hi   = new double[ns];
  double * cpi  = new double[ns];
  double * Mwi  = new double[ns];
  double * rhoi = new double[ns];

  // Fill some vector from the state x, needed to set the state
  for (ii = 0; ii < ns; ++ii)
  {
    rhoi[ii] = x[ii] * rho_now;
  }
  double T_now = x[ns];

  // .. and set the state!
  mix->setState(rhoi, &T_now, 1);
 
  // Extract values from the mixture
  mix->netProductionRates(Wdot);          // Extract chemical rates from mixture object
  mix->getEnthalpiesMass(hi);             // Extract enthalpies per unit mass
  mix->getCpsMass(cpi);                   // Extract Cp_i per unit mass

  // ---------  Now compute the system  ---------
  // As many mass balance equations as many species in the mixture
  for (ii = 0; ii < ns; ++ii)
  {
    dxdt[ii] = Wdot[ii] / rho_now;
  }

  // Now the temperature equation
  double den      = 0.; // denominator
  double num_part = 0; // numerator (part)

  for (ii = 0; ii < ns; ++ii) 
  {
    den      += x[ii]*cpi[ii];
    num_part += hi[ii]*Wdot[ii]/rho_now;
  }

  // Numerically approximate Du^2 / Dt
  double Du2Dt = u_now*( pow(u_1,2) - pow(u_0,2) )/(x_1 - x_0);

  dxdt[ns] = - (num_part + 0.5*Du2Dt)/den; 

  // Making it a derivative w.r.t. x
  for (ii = 0; ii < ns+1; ++ii)  // note that i = 0:ns! All the states, also the last
    dxdt[ii] = dxdt[ii]/u_now;

  // Freeding space!
  delete [] Wdot;
  delete [] hi;
  delete [] cpi; 
  delete [] Mwi;
  delete [] rhoi;

}



// ==================================================================================
//    Functor called by the Boost integration function

struct EquationsToSolve
{
  // Pointers to Mixture and external field objects
  Mutation::Mixture* mix;
  ExternalField* extField;

  // Constructor, initializing the mixture and external field
  EquationsToSolve( Mutation::Mixture* p_mix, ExternalField* p_extField)
  {
    this->mix      = p_mix;
    this->extField = p_extField;
  }

  // Overloading operator, so that it works as a Boost-friendly function 
  void operator()( const vec_type &x, vec_type &dxdt, double t )
  {
    // Calls the equations builder, also passing him the pointer to mixture and
    // to the external field
    LagrProblem( x, dxdt, t, mix, extField );
  }
};



// ==================================================================================
// Jacobian struct
struct LagJacobian
{
  // Pointers to Mixture and external field objects
  Mutation::Mixture* mix;
  ExternalField* extField;

  // Constructor, initializing the mixture and external field
  LagJacobian( Mutation::Mixture* p_mix, ExternalField* p_extField )
  {
    this->mix      = p_mix;
    this->extField = p_extField;
  }

  // Overloading operator: when called by Boost it must compose the Jacobian
  void operator()( const vec_type &x, mat_type &J, double t, vec_type &dfdt )
  {
    // THIS WILL COMPUTE A FINITE DIFFERENCE JACOBIAN
// DANGER!!! OCCHIO!!! what happens when one molar fraction have reached unity? Does it crash??
// Maybe I should introduce a check and.. if it is reached switch to backward finite differences..
    
  // In this overloading, the state "x" is used to evaluate the reference derivative "dxdt".
  // "x" is then perturbed in one direction and the perturbated dxdt_pert is obtained.
  // This allows to compute the derivative of "dxdt" with respect to the perturbation.
  // Base system in latex notation: d \vec{x} / dt = \vec{F}(x)
  // so the Jacobian is (index notation): J(i,j) = d F(i) / d x(j)
  // First of all I perturbate the j-th component of x and compute the resulting F, then I save 
  // at each J(i,j) the value d F(i) / d x(j).
  //
  // The dimension of the Jacobian is (nSpecies+1)*(nSpecies+1) since the state is made of the 
  // species plus one temperature
    size_t Jdim = mix->nSpecies() + 1; // dimension of the Jacobian 
    double pert = 1e-8;                // entity of perturbation
    vec_type x_pert(Jdim);
    vec_type dxdt(Jdim);
    vec_type dxdt_pert(Jdim);

    size_t ii,jj;

    // Computing the reference solution
    LagrProblem( x, dxdt, t, mix, extField );

    // Initializing the perturbation with the same value of the reference state "x"
    for( ii = 0; ii < Jdim; ++ii )
    {
      x_pert[ii] = x[ii];
    }

    // Now for every component of x, I have to compute J(i,j)
    for( jj = 0; jj < Jdim; ++jj )
    {
      x_pert[jj] *= 1.0 + pert;                       // perturbating the i-th component
      LagrProblem( x_pert, dxdt_pert, t, mix, extField ); // computing the perturbated solution
      x_pert[jj]  = x[jj];                            // restoring the previous value
     
      // And now saving all the components of that COLUMN
      for( ii = 0; ii < Jdim; ++ii )
      { 
        J(ii,jj) = ( dxdt_pert[ii] - dxdt[ii] )/(x[jj]*pert + 1.0e-16);
      }
    }
  }
};


// ==================================================================================
// Observer, prints time and state when called (during integration)
struct Observer
{
   Mutation::Mixture* mix;

   // Constructor, assigning a mixture object, from which to know number of species & co
   Observer(Mutation::Mixture* gasMix)
   {
     this->mix = gasMix;
   }

   void operator()( const vec_type &x, const double t )
   {
      // print current abscissa
      std::cout << t << "  ";

      // print mass fractions
      for (int ii = 0; ii < mix->nSpecies(); ++ii)
         std::cout << x[ii] << "  ";

      // and finally print temperature, that is the last field
      std::cout << x[mix->nSpecies()] << std::endl;
   }
};




////    // =======================================================================
////    //    Lagrangian Solver Class
////    
////    class LagSolver
////    {
////      public:
////    
////        size_t nSpecies;
////        // Solution is saved in the following variables
////        std::vector<double> x_vect, u_vect, T_vect, P_vect;    
////        std::vector< std::vector<double> > Yi_vect;
////        std::vector<double> x_ext, u_ext;                        // External field imported
////    
////        // Constructor, initializes nSpecies from input and Yi_vect from nSpecies
////    
////    };
