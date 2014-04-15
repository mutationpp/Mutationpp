#include "StateModel.h"

namespace Mutation {
    namespace Thermodynamics {

/**
 * StateModel which has only one variable set {T, P} and sets the mixture at the
 * equilibrium composition for the given T and P.
 */
class EquilTPStateModel : public StateModel
{
public:
    /**
     * Constructor.
     */
    EquilTPStateModel(const Thermodynamics& thermo)
        : StateModel(thermo)
    { }

    /**
     * Sets the mixture state by computing the equilibrium composition at the
     * given temperature and pressure using the default elemental composition.
     *
     * @param p_T - pointer to temperature value (K)
     * @param p_P - pointer to pressure value (Pa)
     * @param vars - should be ignored
     */
    void setState(
        const double* const p_T, const double* const p_P, const int vars = 0)
    {
        if (vars != 0) {
            cout << "Only T,P variable set implemented in EquilTP StateModel!"
                 << endl;
            exit(1);
        }

        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_T[0];
        m_P = p_P[0];
        m_thermo.equilibriumComposition(m_T, m_P, mp_X);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    EquilTPStateModel, StateModel> equil_tp("EquilTP");

    } // namespace Thermodynamics
} // namespace Mutation
