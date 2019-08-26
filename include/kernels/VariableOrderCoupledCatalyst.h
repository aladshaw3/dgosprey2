/*!
 *  \file VariableOrderCoupledCatalyst.h
 *	\brief Standard kernel for coupling multiple gas and adsorbed species together via a reaction based mechanism
 *	\details This file creates a standard MOOSE kernel for the coupling multiple gas, adsorption, and catalytic species in a
 *				simulation based on the following reaction scheme...
 *
 *				sum(i, v_i*C_i) + m*C + n*L  <-- --> x*D + m*C + sum(j, v_j*C_j)
 *
 *				In this reaction scheme, the i-th species and m number of a catalytic species (C) may interact with
 *              m number of available surface sites (L) to form x number of deactivated sites (D). The catalytic speices
 *              is returned along with some other gas species that may be produced (C_j) as by-products from the reaction.
 *              This expression then formulates the following rate equation for the deactivated site (D)...
 *
 *				dD/dt = x*(C)^a*[k_f*(L)^b*product(i, C_i^z_i) - k_r*(D)^f*product(j, C_j^z_j)]
 *
 *				Parameters are as follows...
 *
 *				v_i,j = stoichiometric coefficients for gas reactants/products
 *              z_i,j = reaction oder with respect to the gas reactants/products
 *				C_i,j = concentrations of gas reactants/products (mol/L)
 *				m = number of catalysts used in the reaction
 *				n = number of adsorption sites (L) reduced
 *              x = number of reduced adsortion sites (D) produced
 *				L = concentration of available sites (mol/kg)
 *				D = concentration of reduced species produced (mol/kg)
 *				k_f,r = rate constant for the forward/reverse reaction
 *              a = reaction order with repect to the catalyst
 *              b = reaction order with respect to available adsorption sites
 *              f = reaction order with respect to the reduced species
 *
 *				Rate expression must be coupled with all involved adsorbed species and uses a site-balance to
 *				account for the loss of adsorption sites during multi-species adsorption
 *
 *				SiteBalance: (L) = Lmax - (D) - sum(i, m_i/n_i*q_i)
 *
 *				where Lmax is the maximum capacity for adsorption (or the maximum available sites), m_i is the
 *				number of sites the m-th adsorbate occupies, n_i is the number of adsorbed species produced from
 *				the reaction consuming m_i sites, and q_i is the adsorbed concentration of the i-th species.
 *
 *  \author Alexander Wiechert
 *	\date 8/23/2018
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2018, all
 *             rights reserved.
 *
 *			   Alexander Wiechert does not claim any ownership or copyright to the
 *			   MOOSE framework in which these kernels are constructed, only
 *			   the kernels themselves. The MOOSE framework copyright is held
 *			   by the Battelle Energy Alliance, LLC (c) 2010, all rights reserved.
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

#include "CoupledCatalyst.h"

/// CoupledCatalyst class object forward declarationss
class VariableOrderCoupledCatalyst;

template<>
InputParameters validParams<VariableOrderCoupledCatalyst>();

/// CoupledCatalyst class object inherits from CoupledConstChemisorption object
/** This class object inherits from the Kernel object in the MOOSE framework.
    All public and protected members of this class are required function overrides.
    The kernel interfaces with the non-linear variables for gas concentrations and
    adsorbed concentrations. */
class VariableOrderCoupledCatalyst : public CoupledCatalyst
{
public:
    /// Required constructor for objects in MOOSE
    VariableOrderCoupledCatalyst(const InputParameters & parameters);

protected:
    /// Function to compute the rate function for the reaction
    Real computeRateFunction();

    /// Function to compute the diagonal Jacobi for the rate function
    Real computeRateFunctionJacobi();

    /// Function to compute the gas concentration off-diagonal Jacobi for the rate function
    Real computeRateFunctionGasOffDiagJacobi(int i);

    /// Function to compute the catalyst off-diagonal Jacobi for the rate function
    Real computeRateFunctionCatalystOffDiagJacobi(int i);

    /// Function to compute the adsorption concentration off-diagonal Jacobi for the rate function
    Real computeRateFunctionAdsOffDiagJacobi(int i);

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

    std::vector<Real> _c_order;                      ///< Catalyst reaction order
    std::vector<Real> _g_order;                      ///< Gas reaction order
    std::vector<Real> _a_order;                      ///< Adsorbed/Deactivated reaction order
    std::vector<Real> _s_order;                      ///< Adsorbent reaction order

private:

};
