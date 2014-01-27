#include "StateModel.h"
#include "Utilities.h"
#include <algorithm>

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

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
     */
    virtual void setState(const double* const p_T, const double* const p_P)
    {
        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_T[0];
        m_P = p_P[0];
        m_thermo.equilibriumComposition(m_T, m_P, mp_X);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    EquilTPStateModel, StateModel> equil_tp("EquilTP");

//==============================================================================

class EquilTPXeStateModel : public StateModel
{
public:
    /**
     * Constructor.
     */
    EquilTPXeStateModel(const Thermodynamics& thermo)
        : StateModel(thermo)
    { }

    /**
     * Sets the mixture state by computing the equilibrium composition at the
     * given temperature and pressure using the default elemental composition.
     *
     * @param p_T - pointer to temperature value (K)
     * @param p_P - pointer to pressure value (Pa)
     */
    virtual void setState(const double* const p_TP, const double* const p_Xe)
    {
        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_TP[0];
        m_P = p_TP[1];
        m_thermo.equilibriumComposition(m_T, m_P, p_Xe, mp_X);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    EquilTPXeStateModel, StateModel> equil_tpxe("EquilTPXe");

//==============================================================================

class TPXStateModel : public StateModel
{
public:
    /**
     * Constructor.
     */
    TPXStateModel(const Thermodynamics& thermo)
        : StateModel(thermo)
    { }
    
    /**
     * Sets the mixture state from mixture temperature, pressure, and species
     * mole fractions.
     *
     * @param p_TP - pointer to {T (K), P (Pa)} array
     * @param p_X  - pointer to species mole fraction array
     */
    virtual void setState(const double* const p_TP, const double* const p_X)
    {
        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_TP[0];
        m_P = p_TP[1];
        std::copy(p_X, p_X+m_thermo.nSpecies(), mp_X);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    TPXStateModel, StateModel> tpx("TPX");

//==============================================================================

class TPYStateModel : public StateModel
{
public:
    /**
     * Constructor.
     */
    TPYStateModel(const Thermodynamics& thermo)
        : StateModel(thermo)
    { }
    
    /**
     * Sets the mixture state from mixture temperature, pressure, and species
     * mass fractions.
     *
     * @param p_TP - pointer to {T (K), P (Pa)} array
     * @param p_Y  - pointer to species mass fraction array
     */
    virtual void setState(const double* const p_TP, const double* const p_Y)
    {
        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_TP[0];
        m_P = p_TP[1];
        m_thermo.convert<Y_TO_X>(p_Y, mp_X);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    TPYStateModel, StateModel> tpy("TPY");

//==============================================================================

class TRhoiStateModel : public StateModel
{
public:
    /**
     * Constructor.
     */
    TRhoiStateModel(const Thermodynamics& thermo)
        : StateModel(thermo)
    { }
    
    /**
     * Sets the mixture state from mixture temperature and species densities.
     *
     * @param p_T    - pointer to mixture temperature
     * @param p_rhoi - pointer to species densities array
     */
    virtual void setState(const double* const p_T, const double* const p_rhoi)
    {
//        cout << "T,Rhoi : setState(" << p_T[0];
//        for (int i = 0; i < m_thermo.nSpecies(); ++i)
//            cout << ", " << p_rhoi[i];
//        cout << ")" << endl;
        m_T = m_Tr = m_Tv = m_Tel = m_Te = p_T[0];
        
        m_P = 0.0;
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
            m_P += p_rhoi[i];
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
            mp_X[i] = p_rhoi[i] / m_P;
        m_P = m_thermo.pressure(m_T, m_P, mp_X);
        
        m_thermo.convert<Y_TO_X>(mp_X, mp_X);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    TRhoiStateModel, StateModel> trhoi("TRhoi");

//==============================================================================

class TTvRhoiStateModel : public StateModel
{
public:
    /**
     * Constructor.
     */
    TTvRhoiStateModel(const Thermodynamics& thermo)
        : StateModel(thermo)
    { }
    
    /**
     * Sets the mixture state from mixture temperature, pressure, and species
     * mass fractions.
     *
     * @param p_T    - pointer to {T (K), Tv (K)} array
     * @param p_rhoi - pointer to species densities array
     */
    virtual void setState(const double* const p_T, const double* const p_rhoi)
    {
        m_T = m_Tr = p_T[0];
        m_Tv = m_Tel = m_Te = p_T[1];
        
        m_P = 0.0;
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
            m_P += p_rhoi[i];
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
            mp_X[i] = p_rhoi[i] / m_P;
        m_P = m_thermo.pressure(m_T, m_P, mp_X);
        
        m_thermo.convert<Y_TO_X>(mp_X, mp_X);
    }

    virtual Mutation::Transfer::TransferModel* createTransferModel(
        Thermodynamics& thermo,
        Mutation::Transport::Transport& transport,
        Mutation::Kinetics::Kinetics& kinetics)
    {
        return new Transfer::TransferModelVT(thermo);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    TTvRhoiStateModel, StateModel> ttvrhoi("TTvRhoi");

//==============================================================================

class TTERhoiStateModel : public StateModel
{
public:
    /**
     * Constructor.
     */
    TTERhoiStateModel(const Thermodynamics& thermo)
        : StateModel(thermo)
    { }
    
    /**
     * Sets the mixture state from mixture temperature, pressure, and species
     * mass fractions.
     *
     * @param p_T    - pointer to {T (K), Tv (K)} array
     * @param p_rhoi - pointer to species densities array
     */
    virtual void setState(const double* const p_T, const double* const p_rhoi)
    {

        m_T = m_Tr = p_T[0];
        m_Tv = m_Tel = m_Te = p_T[1];
        
        m_P = 0.0;
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
            m_P += p_rhoi[i];
        for (int i = 0; i < m_thermo.nSpecies(); ++i)
            mp_X[i] = p_rhoi[i] / m_P;
        m_P = m_thermo.pressure(m_T, m_P, mp_X);
        
        m_thermo.convert<Y_TO_X>(mp_X, mp_X);



    }

    virtual Mutation::Transfer::TransferModel* createTransferModel(
        Thermodynamics& thermo,
        Mutation::Transport::Transport& transport,
        Mutation::Kinetics::Kinetics& kinetics)
    {
        return new Transfer::TransferModelTE(thermo,transport);
    }
};

// Register the state model
Utilities::Config::ObjectProvider<
    TTERhoiStateModel, StateModel> tterhoi("TTERhoi");


    } // namespace Thermodynamics
} // namespace Mutation


