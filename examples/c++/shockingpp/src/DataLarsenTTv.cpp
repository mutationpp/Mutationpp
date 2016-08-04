#include "DataLarsenTTv.h"

DataLarsenTTv::DataLarsenTTv(Mutation::Mixture& l_mix)
                        : m_mix(l_mix),
                          n_sp(l_mix.nSpecies()),
                          n_meq(0), // An 1D solver can have 1 component of velocity
                          n_eneq(l_mix.nEnergyEqns()),
                          n_eq(n_sp+n_eneq),
                          pos_V(n_sp),  // NOT NEEDED, just initialized
                          pos_T(n_sp),  // NOT NEEDED, just initialized
                          m_P(0.0),
                          m_V(0.0),     // NOT NEEDED, just initialized
                          m_rho(0.0),
                          m_mdot(0.0),  // NOT NEEDED, just initialized
                          v_rhoi(n_sp, 0.0),
                          v_xi(n_sp, 0.0),
                          v_yi(n_sp, 0.0),
                          v_T(n_eneq, 0.0),
                          v_X(n_eq, 0.0)
{
  Y0_ext.resize(m_mix.nSpecies()); // Right dimension to Y0
}

// -----------------------------------------------------------------------------

DataLarsenTTv::~DataLarsenTTv() {}

// -----------------------------------------------------------------------------

vector_type DataLarsenTTv::getInitialState()
{
  vector_type x(n_eq); // Species + 1 temperature
  
  // chemical species
  for(size_t ii = 0; ii < m_mix.nSpecies(); ++ii) {
    x[ii] = Y0_ext[ii];
  }
  // translational temperature
  x[m_mix.nSpecies()] = T_ext.at(0);
  // vibrational temperature
  x[m_mix.nSpecies()+1] = Tv_ext.at(0);

  return x;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::setCurrentStep(int pos)
{
  currentStep = pos;

  return; 
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepValues(double x_now[], double u_now[], double rho_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function
  x_now[0]   = x_ext.at(currentStep);
  x_now[1]   = x_ext.at(currentStep+1);
  u_now[0]   = u_ext.at(currentStep);
  u_now[1]   = u_ext.at(currentStep+1);
  rho_now[0] = rho_ext.at(currentStep);
  rho_now[1] = rho_ext.at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

int DataLarsenTTv::getStrNele()
{
  // Returns the number of elements of the imported streamline
  return x_ext.size();
}
// -----------------------------------------------------------------------------

void DataLarsenTTv::readDataFile()
{
  std::string lineRead;
  std::ifstream filein("output.dat");
  if(filein.is_open() == 0)
  {
    std::cout << "I could not open the file.. Does it exist?\nAborting!\n";
    exit(EXIT_FAILURE);
  }

  getline(filein,lineRead); // read first line
  std::stringstream nowstream(lineRead); // and unpack it in the nowstream variable

  double xNow, TNow, TvNow, PNow, rhoNow, uNow, YiNow, dummy;     // those are temporary

  // -----   Reading first line and saving initial conditions
  nowstream >> xNow >> TNow >> TvNow >> PNow >> rhoNow >> uNow;

  // ..and saving initial value for species mass fractions
  for(size_t ii = 0; ii < m_mix.nSpecies(); ++ii)
  {
    nowstream >> YiNow;
    Y0_ext.at(ii) = YiNow;
  }

  // -----   Now saving the externally computed field: I only care about position and velocity
  // The first line has already been extracted, so I save it here.
  x_ext.push_back(xNow);
  T_ext.push_back(TNow);
  Tv_ext.push_back(TvNow);
  u_ext.push_back(uNow);
  rho_ext.push_back(rhoNow);

  // Keep reading now until the EOF
  while(getline(filein,lineRead))
  {
    nowstream.clear();          // clear stream variable nowstream
    nowstream.str(lineRead);    // put newly read line into nowstream

    nowstream >> xNow >> TNow >> TvNow >> dummy >> rhoNow >> uNow; // ...and unpack it into my variables

    x_ext.push_back(xNow);
    T_ext.push_back(TNow);
    Tv_ext.push_back(TvNow);
    u_ext.push_back(uNow);
    rho_ext.push_back(rhoNow);
  }

  filein.close();  // I'm done with this data file

  return;
}

