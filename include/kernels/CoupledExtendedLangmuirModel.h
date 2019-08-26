/*!
 *  \file CoupledExtendedLangmuirModel.h
 *	\brief Standard kernel for coupling a vector non-linear variables via an extended langmuir function
 *	\details This file creates a standard MOOSE kernel for the coupling of a vector non-linear variables
 *			together via an extended langmuir forcing function,
 *			i.e., variable = b_i * K_i * coupled_variable_i / 1 + sum(j, K_j * coupled_variable_j).
 *			In this kernel, the langmuir coefficients {K_i} are calculated as a function of temperature
 *			using the van't Hoff expression: ln(K_i) = -dH_i/(R*T) + dS_i/R
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 10/09/2017
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

#include "CoupledExtendedLangmuirFunction.h"
#include "flock.h"

/// CoupledExtendedLangmuirModel class object forward declarationss
class CoupledExtendedLangmuirModel;

template<>
InputParameters validParams<CoupledExtendedLangmuirModel>();

/// CoupledExtendedLangmuirFunction class object inherits from CoupledExtendedLangmuirFunction object
/** This class object inherits from the CoupledExtendedLangmuirFunction object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the set of non-linear variables to couple an extended Langmuir
	forcing function between given objects. */
class CoupledExtendedLangmuirModel : public CoupledExtendedLangmuirFunction
{
public:
	/// Required constructor for objects in MOOSE
	CoupledExtendedLangmuirModel(const InputParameters & parameters);
	
protected:
	/// Function to compute all langmuir coefficients from temperature
	void computeAllLangmuirCoeffs();
	
	/// Function to compute the Jacobi for the coupled temperature
	Real computeExtLangmuirTempJacobi();
	
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
	
	std::vector<Real> _enthalpies;						///< Vector of enthalpies for all langmuir coefficients (J/mol)
	std::vector<Real> _entropies;						///< Vector of entropies for all langmuir coefficients (J/K/mol)
	const VariableValue & _coupled_temp;				///< Coupled variable for temperature
	const unsigned int _coupled_var_temp;				///< Index for the coupled temperature variable
	
private:
	
};
