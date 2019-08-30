/*!
 *  \file VariableOrderReac.h
 *    \brief Standard kernel for coupling multiple species via a reaction based mechanism
 *    \details This file creates a standard MOOSE kernel for the coupling species in a
 *                simulation based on the following reaction scheme...
 *
 *                sum(i, v_i*C_i) <-- --> n*q + sum(j, v_j*C_j)
 *
 *                In this reaction scheme, the i-th species may interact to form n number of
 *                species (q). In return, some other gas species may be produced (C_j) as
 *                by-products from the reaction. This expression then formulates the following
 *                rate equation for the coupled species (q)...
 *
 *                dq/dt = n*k_f*product(i, C_i^v_i) - n*k_r*(q)^n*product(j, C_j^v_j)
 *
 *                Parameters are as follows...
 *
 *                v_i,j = stoichiometric coefficients for reactants/products
 *                C_i,j = concentrations of reactants/products (mol/L or mol/kg)
 *                n = number of product species (q) produced
 *                q = concentration of coupled species produced (mol/kg)
 *                k_f,r = rate constant for the forward/reverse reaction
 *
 *  \author Alexander Wiechert
 *    \date 08/29/2019
 *    \copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *             by Idaho National Laboratory and Oak Ridge National Laboratory
 *             engineers and scientists. Portions Copyright (c) 2018, all
 *             rights reserved.
 *
 *             Alexander Wiechert does not claim any ownership or copyright to the
 *             MOOSE framework in which these kernels are constructed, only
 *             the kernels themselves. The MOOSE framework copyright is held
 *             by the Battelle Energy Alliance, LLC (c) 2010, all rights reserved.
 */

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#pragma once

#include "TimeDerivative.h"

///VariableOrderReac class object forward declaration
class VariableOrderReac;

template<>
InputParameters validParams<VariableOrderReac>();

/// VariableOrderReac class object inherits from TimeDerivative object
/** This class object inherits from the TimeDerivative object in the MOOSE framework.
 All public and protected members of this class are required function overrides.
 The kernel interfaces with the non-linear variables for the concentrations of various
 species. */

class VariableOrderReac : public TimeDerivative
{
public:
    /// Required constructor for objects in MOOSE
    VariableOrderReac(const InputParameters & parameters);
    
protected:
    /// Function to compute the rate function for the reaction
    Real computeRateFunction();
    
    /// Function to compute the diagonal Jacobi for the rate function
    Real computeRateFunctionJacobi();
    
    /// Function to compute the off-diagonal Jacobi for the rate function
    Real computeRateFunctionSpeciesOffDiagJacobi(int i);
    
    /// Required residual function for standard kernels in MOOSE
    /** This function returns a residual contribution for this object.*/
    virtual Real computeQpResidual();
    
    /// Required Jacobian function for standard kernels in MOOSE
    /** This function returns a Jacobian contribution for this object. The Jacobian being
     computed is the associated diagonal element in the overall Jacobian matrix for the
     system and is used in preconditioning of the linear sub-problem. */
    virtual Real computeQpJacobian();
    
    /// Not Required, but aids in the preconditioning step
    /** This function returns the off diagonal Jacobian contribution for this object. By
     returning a non-zero value we will hopefully improve the convergence rate for the
     cross coupling of the variables. */
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);
    
    Real _forward;                                       ///< Rate constant for the forward reaction (per hour)
    Real _reverse;                                       ///< Rate constant for the reverse reaction (per hour)
    Real _m_stoich;                                      ///< Stoichiometry of the primary species
    Real _m_order;                                       ///< Reaction Order of the primary species
    
    std::vector<Real> _stoich;                           ///< Vector of species stoichiometry: (-) = reactant, (+) = product
    std::vector<Real> _order;                            ///< Reaction order of Species
    
    std::vector<const VariableValue *> _coupled_species;     ///< Pointer list to the coupled species (mol/L for gases, mol/kg for all else)
    std::vector<unsigned int> _coupled_species_vars;         ///< Indices for the species in the system

    unsigned int _coupled_var_i;                        ///< Variable index for the coupled species
    
private:
    
};
