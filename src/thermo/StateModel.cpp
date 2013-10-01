#include "StateModel.h"
#include "Utilities.h"
#include <algorithm>

namespace Mutation {
    namespace Thermodynamics {

void StateModel::setStateTPX(
    const double* const T, const double* const P, const double* const X)
{
    std::copy(T, T+nT(), mp_T);
    std::copy(P, P+nP(), mp_P);
    std::copy(X, X+m_ns, mp_X);
}

void StateModel::setStateTPX(
    const double T, const double P, const double* const X)
{
    std::fill(mp_T, mp_T+nT(), T);
    std::fill(mp_P, mp_P+nP(), P);
    std::copy(X, X+m_ns, mp_X);
}

void StateModel::init() 
{
    // Set temperature indices
    m_it = m_itr = m_itv = m_itel = m_ite = 0;
    m_nt = 1;

    // Allocate storage
    mp_T = new double [nT()];
    mp_P = new double [nP()];
    mp_X = new double [m_ns];
    
    // Initialize everything to 0
    std::fill(mp_T, mp_T+nT(), 0.0);
    std::fill(mp_P, mp_P+nP(), 0.0);
    std::fill(mp_X, mp_X+m_ns, 0.0);
}

// Register the basic state model
Utilities::Config::ObjectProvider<StateModel, StateModel> tpModel("T");


/**
 * Represents a mixture using two temperatures: \f$T_1 = T = T_e = T_{el}\f$ and 
 * \f$T_2 = T_v = T_r\f$.  
 */
class TTvStateModel : public StateModel
{
public:
    
    /**
     * Constructor.
     */
    TTvStateModel(int ns)
        : StateModel(ns)
    {
        // Set the temperature indices
        m_it  = m_itr  =         0;
        m_itv = m_itel = m_ite = 1;
        m_nt = 2;
    }
    
}; // class TTvStateModel

// Register the T,Tv model
Utilities::Config::ObjectProvider<TTvStateModel, StateModel> ttvModel("T,Tv");

    } // namespace Thermodynamics
} // namespace Mutation


