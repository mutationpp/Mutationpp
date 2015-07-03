#ifndef JACOBIANMANAGERGSI_H
#define JACOBIANMANAGERGSI_H

#include "GSIStoichiometryManager.h"
#include "JacobianManagerGSI.h"

namespace Mutation{
    namespace gsi{

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

class JacobianManagerGSI{
public:
    virtual ~JacobianManagerGSI() = 0;

    virtual double getJacobian(const int& i_row, const int& i_column) = 0;

};

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

class JacPhysisorption : public JacobianManagerGSI{
public:
    JacPhysisorption(const std::vector<double>& v_reactants);
    ~JacPhysisorption(){ }

    double getJacobian(const int& i_row, const int& i_column);
private:

const std::vector<double>& v_Jac_reac; // Possibly here have a stoichiometry manager also...



};

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

class JacThermalDesorption : public JacobianManagerGSI{
public:
    JacThermalDesorption();
    ~JacThermalDesorption(){ }

    double getJacobian(const int& i_row, const int& i_column);
private:
};

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

class JacChemisorption : public JacobianManagerGSI{
public:
    JacChemisorption();
    ~JacChemisorption(){ }

    double getJacobian(const int& i_row, const int& i_column);
private:
};

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

class JacERRecombination : public JacobianManagerGSI{
public:
    JacERRecombination();
    ~JacERRecombination(){ }

    double getJacobian(const int& i_row, const int& i_column);
private:
};

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

//class Jac : public JacobianManagerGSI{
//
//};
//class Jac : public JacobianManagerGSI{
//
//};

    } // namespace gsi
} // namespace Mutation


#endif // JACOBIANMANAGERGSI_H
