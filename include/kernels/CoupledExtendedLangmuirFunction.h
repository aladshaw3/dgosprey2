/*!
 *  \file CoupledExtendedLangmuirFunction.h
 *	\brief Standard kernel for coupling a vector non-linear variables via an extended langmuir function
 *	\details This file creates a standard MOOSE kernel for the coupling of a vector non-linear variables
 *			together via an extended langmuir forcing function, 
 *			i.e., variable = b_i * K_i * coupled_variable_i / 1 + sum(j, K_j * coupled_variable_j).
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 06/21/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
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

/// CoupledExtendedLangmuirFunction class object forward declarationss
class CoupledExtendedLangmuirFunction;

template<>
InputParameters validParams<CoupledExtendedLangmuirFunction>();

/// CoupledExtendedLangmuirFunction class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the set of non-linear variables to couple an extended Langmuir 
	forcing function between given objects. */
class CoupledExtendedLangmuirFunction : public TimeDerivative
{
public:
	/// Required constructor for objects in MOOSE
	CoupledExtendedLangmuirFunction(const InputParameters & parameters);
	
protected:
	/// Function to compute the Extended Langmuir Equilibrium value
	Real computeExtLangmuirEquilibrium();
	
	/// Function to compute the Jacobi for the main coupled concentration
	Real computeExtLangmuirConcJacobi();
	
	/// Function to compute the off-diagonal Jacobi for the other coupled concentrations
	Real computeExtLangmuirOffJacobi(int i);
	
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
	std::vector<Real> _langmuircoef;					///< Langmuir Coefficients for the coupled variables
	std::vector<const VariableValue *> _coupled;		///< Pointer list to the coupled gases
	std::vector<unsigned int> _coupled_vars;			///< Indices for the gas species in the system
	const VariableValue & _coupled_i;					///< Primary Coupled variable
	const unsigned int _coupled_var_i;					///< Variable identification for the primary coupled variable
	int _lang_index;									///< Index for primary langmuir coefficient
	
private:
	
};
