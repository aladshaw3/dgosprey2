/*!
 *  \file CoupledConstChemisorption.h
 *	\brief Standard kernel for coupling multiple gas and adsorbed species together via a reaction based mechanism
 *	\details This file creates a standard MOOSE kernel for the coupling multiple gas and adsorption species in a
 *				simulation based on the following reaction scheme...
 *
 *				sum(i, v_i*C_i) + m*L  <-- --> n*q + sum(j, v_j*C_j)
 *
 *				In this reaction scheme, the i-th species may interact with m number of available surface sites (L)
 *				to form n number of adsorbed species (q). In return, some other gas species may be produced (C_j) as
 *				by-products from the reaction. This expression then formulates the following rate equation for the
 *				adsorbed species (q)...
 *
 *				dq/dt = n*k_f*(L)^m*product(i, C_i^v_i) - n*k_r*(q)^n*product(j, C_j^v_j)
 *
 *				Parameters are as follows...
 *
 *				v_i,j = stoichiometric coefficients for gas reactants/products
 *				C_i,j = concentrations of gas reactants/products (mol/L)
 *				m = number of adsorption sites used in the reaction
 *				n = number of adsorption product species (q) produced
 *				L = concentration of available sites (mol/kg)
 *				q = concentration of adsorbed species produced (mol/kg)
 *				k_f,r = rate constant for the forward/reverse reaction
 *
 *				Rate expression must be coupled with all involved adsorbed species and uses a site-balance to 
 *				account for the loss of adsorption sites during multi-species adsorption
 *
 *				SiteBalance: (L) = Lmax - sum(i, m_i/n_i*q_i)
 *
 *				where Lmax is the maximum capacity for adsorption (or the maximum available sites), m_i is the 
 *				number of sites the m-th adsorbate occupies, n_i is the number of adsorbed species produced from
 *				the reaction consuming m_i sites, and q_i is the adsorbed concentration of the i-th species.
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 02/05/2018
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

#include "TimeDerivative.h"

//#ifndef CoupledConstChemisorption_h
//#define CoupledConstChemisorption_h

/// CoupledConstChemisorption class object forward declarationss
class CoupledConstChemisorption;

template<>
InputParameters validParams<CoupledConstChemisorption>();

/// CoupledConstChemisorption class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces with the non-linear variables for gas concentrations and
	adsorbed concentrations. */
class CoupledConstChemisorption : public TimeDerivative
{
public:
	/// Required constructor for objects in MOOSE
	CoupledConstChemisorption(const InputParameters & parameters);
	
protected:
	/// Function to compute the rate function for the reaction
	Real computeRateFunction();
	
	/// Function to compute the diagonal Jacobi for the rate function
	Real computeRateFunctionJacobi();
	
	/// Function to compute the gas concentration off-diagonal Jacobi for the rate function
	Real computeRateFunctionGasOffDiagJacobi(int i);
	
	/// Function to compute the adsorption concentration off-diagonal Jacobi for the rate function
	Real computeRateFunctionAdsOffDiagJacobi(int i);
	
	/// Function to compute the site balance
	Real computeSiteBalance();
	
	/// Function to compute the off-diagonal Jacobi for the site balance
	Real computeSiteBalanceOffDiagJacobi(int i);
	
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
	
	Real _maxcap;										///< Maximum number of adsorption sites (mol/kg)
	Real _forward;										///< Rate constant for the forward reaction (per hour)
	Real _reverse;										///< Rate constant for the reverse reaction (per hour)
	
	std::vector<Real> _gas_stoich;						///< Vector of gas species stoichiometry: (-) = reactant, (+) = product
	std::vector<Real> _ads_sites;						///< Number of sites each adsorbent uses (m_i)
	std::vector<Real> _ads_stoich;						///< Number of adsorbed species each reaction forms (n_i)
	
	std::vector<const VariableValue *> _coupled_gas;	///< Pointer list to the coupled gases (mol/L)
	std::vector<unsigned int> _coupled_gas_vars;		///< Indices for the gas species in the system
	std::vector<const VariableValue *> _coupled_ads;	///< Pointer list to the coupled adsorbed concentrations (mol/kg)
	std::vector<unsigned int> _coupled_ads_vars;		///< Indices for the adsorbed species in the system
	int _ads_index;										///< Index for the primary adsorption species
	unsigned int _coupled_var_i;						///< Variable index for the coupled species (same species as _ads_index)
	
private:
	
};

//#endif /* CoupledConstChemisorption_h */
