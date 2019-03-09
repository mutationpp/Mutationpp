<a id="top"></a>

_This is an incomplete file brought over from the original Wiki and needs to be heavily edited and updated._

# State Models

## Introduction

## Generalized Governing Equations

### Mass Conservation
A general set of mass conservation equations can be written as
\f[
\frac{\partial \rho_i}{\partial t} + \nabla\cdot(\rho_i\vec{v}) +
    \nabla\cdot\vec{J}_i = \dot{\omega}_i \quad \forall\; i\in\mathscr{M},
\f]
where the index \f$i\;\in\mathscr{M}\f$ refers to the \f$i^\text{th}\f$ mass
component being modeled.

\f$\vec{J}_i\f$ is the diffusion mass flux of component
\f$i\f$.  When chemical nonequilibrium is assumed, components may be individual
species and/or energy states.  In this case \f$J_i = \rho_i \vec{V}_i \;\forall\;i\in\mathscr{S}\f$
and \f$\mathscr{M} = \mathscr{S}\f$ is the set of species/energy state indices being
considered.  For chemical equilibrium, \f$\vec{J}_i = \sum_{k\in\mathscr{S}}
\nu_{ik} M_{w,i}/M_{w,k} \rho_k \vec{V}_k\f$ if elemental diffusion is allowed, or \f$\vec{J}_i = 0\f$
otherwise.

The production rates of component \f$i\f$, \f$\dot{\omega}_i\f$ represent
any process which may produce or destroy a mass component such as elementary
chemical reactions.  In general, both mass diffusion fluxes and production rates
are linearly dependent following \f$\sum_{i\in\mathscr{M}}\vec{J}_i = \vec{0}\f$ and
\f$\sum_{i\in\mathscr{M}}\dot{\omega}_i = 0\f$.

### Momentum Conservation
Momentum conservation is written as
\f[
\frac{\partial\rho\vec{v}}{\partial t} + \nabla\cdot(\rho\vec{v}\otimes\vec{v})
 + \nabla\cdot(p\matrix{I}) - \nabla\cdot\matrix{\tau} = \vec{0},
\f]
where \f$\matrix{\tau}\f$ is the second order shear stress tensor
\f[
\matrix{\tau} = \eta \left[ \nabla \vec{v} + (\nabla \vec{v})^T -
    \frac{2}{3} (\nabla\cdot \vec{v}) \matrix{I} \right].
\f]
The total density is defined simply as \f$\rho = \sum_{i\in\mathscr{M}} \rho_i\f$.
Pressure is given by Dalton's Law as \f$p = \sum_{k\in\mathscr{S}} p_k =
\sum_{k\in\mathscr{H}} n_k k_B T_h + n_e k_B T_e\f$, where the heavy particle
translation temperature \f$T_h\f$ may or may not equal the electron temperature
\f$T_e\f$ depending on the model in use.

### Total Energy Conservation
The internal energy of the system may be considered split into \f$n^\mathscr{E}\f$ different
energy modes which follow Boltzmann distributions according to an individual temperature
prescribed to each mode. Total energy conservation is then written accordingly as
\f[
\frac{\partial \rho E}{\partial t} + \nabla\cdot(\rho\vec{v}H) =
    \nabla\cdot(\matrix{\tau}\cdot\vec{v}) -\nabla\cdot\vec{q},
\f]
where the total energy and enthalpy densities are the summation of each internal
energy mode and the bulk kinetic energy density.
\f[
\rho E = \sum_{m\in\mathscr{E}} \rho e^m + \frac{1}{2}\rho\vec{v}\cdot\vec{v}
\f]
\f[
\rho H = \sum_{m\in\mathscr{E}} \rho h^m + \frac{1}{2}\rho\vec{v}\cdot\vec{v}
\f]
\f$\mathscr{E} = \{1,\dots,n^\mathscr{E}\}\f$ is the set of energy mode indices.
Likewise, the total heat flux vector is the summation of heat fluxes due to each
energy mode.
\f[
\vec{q} = \sum_{m\in\mathscr{E}} \vec{q}^m
\f]
\f[
\vec{q}^m = \sum_{k\in\mathscr{S}}\rho_k h_k^m\vec{V}_k-\lambda^m\nabla T^m
\f]

### Internal Energy Conservation
\f$n^\mathscr{E}-1\f$ additional energy conservation equations are necessary to
close the system.  They may be written as
\f[
\frac{\partial\rho e^m}{\partial t} + \nabla\cdot(\rho\vec{v}e^m) =
    -\nabla\cdot\vec{q}^m + \Omega^m - \delta_{m\mathscr{I_e}}\; p_e\nabla\cdot\vec{v}
    \quad \forall\; m \in \{2,\dots,n^\mathscr{E}\}
\f]
where \f$\delta\f$ is the Kronecker delta function and \f$\mathscr{I_e}\f$ is
the index of the energy mode which includes the free electron energy contribution.
\f$\Omega^m\f$ is the energy transfered into the energy mode \f$m\f$ (and out of
all other modes).

## StateModel Interface

Symbol          | Description           | Function
:--------------:|:----------------------|:-------------------------
\f$n^\mathscr{E}\f$ | Number of energy modes | (not available)
\f$n^\mathscr{M}\f$ | Number of mass components | (not available)
\f$\vec{J}_i\f$ | Mass diffusion fluxes | (not available)
\f$\dot{\omega}_i\f$ | Mass production rates | (not available)
\f$T^m\f$       | Temperature for each mode  | Mutation::Thermodynamics::StateModel::getTemperatures()
\f$T_h\f$ | Heavy particle translation temperature | (not available)
\f$T_e\f$ | Electron translation temperature | (not available)
\f$\lambda^m\f$ | Thermal conductivity of each energy mode | (not available)
\f$e^m\f$ | Energy per mass of each mode | (not available)
\f$h^m\f$ | Enthalpy per mass of each mode | (not available)
\f$\Omega^m\f$ | Energy transfer source terms | (not available)
\f$\mathscr{I_e}\f$ | Index of free electron energy mode | (not available)

