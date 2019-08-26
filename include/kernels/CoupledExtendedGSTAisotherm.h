/*!
 *  \file CoupledExtendedGSTAisotherm.h
 *	\brief Standard kernel for coupling a vector non-linear variables via an extended gsta isotherm
 *	\details This file creates a standard MOOSE kernel for the coupling of a vector non-linear variables
 *			together via an extended gsta forcing function...
 *
 *			                                    SUM(n, Kn_i * (coupled_variable_i/Co)^n)
 *					variable = qmax_i * ----------------------------------------------------
 *                                       1 + SUM(j, SUM(n, Kn_j * (coupled_variable_j/Co)^n)
 *
 *						where Co = 100.0 / (8.3144621 * _coupled_temp)
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 05/10/2018
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2018, all
 *             rights reserved.
 *
 *			   Austin Ladshaw does not claim any ownership or copyright to the
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

#include "Kernel.h"

/// CoupledExtendedGSTAisotherm class object forward declarationss
class CoupledExtendedGSTAisotherm;

template<>
InputParameters validParams<CoupledExtendedGSTAisotherm>();

/// CoupledExtendedGSTAisotherm class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the set of non-linear variables to couple an extended GSTA
	forcing function between given objects. */
class CoupledExtendedGSTAisotherm : public Kernel
{
public:
	/// Required constructor for objects in MOOSE
	CoupledExtendedGSTAisotherm(const InputParameters & parameters);
	
protected:
	/// Function to compute the Extended GSTA Equilibrium value
	Real computeExtGSTAEquilibrium();
	
	/// Function to compute the Jacobi for the main coupled concentration
	Real computeExtGSTAConcJacobi();
	
	/// Function to compute the off-diagonal Jacobi for the other coupled concentrations
	Real computeExtGSTA_ConcOffJacobi(int i);
	
	/// Function to compute the off-diagonal Jacobi for the other coupled concentrations
	Real computeExtGSTA_TempOffJacobi();
	
	/// Function to compute the top of the extended GSTA function
	Real computeTopExtGSTA();
	
	/// Function to compute the bottom of the extended GSTA function
	Real computeBottomExtGSTA();
	
	/// Function to compute derivative of top with respect to given concentration
	Real computeTopDerivativeConc(int i);
	
	/// Function to compute derivative of bottom with respect to given concentration
	Real computeBottomDerivativeConc(int i);
	
	/// Function to compute derivative of top with respect to temperature
	Real computeTopDerivativeTemp();
	
	/// Function to compute derivative of bottom with respect to temperature
	Real computeBottomDerivativeTemp();
	
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
	
	Real _maxcap;										///< Maximum Capacity for the primary adsorbed species
	std::vector<int> _num_sites;						///< Number of adsorption sites for each variable
	std::vector< std::vector<Real> > _gstaparams;		///< GSTA Coefficients for the all coupled variables
	std::vector< Real > _param_1;						///< GSTA parameter 1 for all variables
	std::vector< Real > _param_2;						///< GSTA parameter 2 for all variables
	std::vector< Real > _param_3;						///< GSTA parameter 3 for all variables
	std::vector< Real > _param_4;						///< GSTA parameter 4 for all variables
	std::vector< Real > _param_5;						///< GSTA parameter 5 for all variables
	std::vector< Real > _param_6;						///< GSTA parameter 6 for all variables
	std::vector<const VariableValue *> _coupled;		///< Pointer list to the coupled gases
	std::vector<unsigned int> _coupled_vars;			///< Indices for the gas species in the system
	const VariableValue & _coupled_i;					///< Primary Coupled variable index
	const unsigned int _coupled_var_i;					///< Variable identification for the primary coupled variable
	int _main_index;									///< Index for primary species
	const VariableValue & _coupled_temp;				///< Coupled gas temperature variable
	const unsigned int _coupled_var_temp;				///< Variable identification for the coupled temperature variable
	
private:
	
};
