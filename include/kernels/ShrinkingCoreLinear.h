/*!
 *  \file ShrinkingCoreLinear.h
 *	\brief Standard kernel for a standard costant shrinking core
 *	\details This file creates a standard MOOSE kernel for constant shrinking core mechanism that
 *			can be added to the non-linear residuals. It contains all the same parameters as the more
 *			generic base class, but couples with another non-linear variable via a linear relationship.
 *			The linear coefficient for that relationship is added as an additional parameter to be
 *			set by the user.
 *
 *  \author Alexander Wiechert, Austin Ladshaw
 *	\date 12/19/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2015, all
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

#include "ShrinkingCoreConstant.h"

/// ShrinkingCoreLinear class object forward declarations
class ShrinkingCoreLinear;

template<>
InputParameters validParams<ShrinkingCoreLinear>();

/// ShrinkingCoreLinear class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel has all the protected members from Kernel, but also
	includes coefficients for the shrinking core model.
 */
class ShrinkingCoreLinear : public ShrinkingCoreConstant
{
public:
	/// Required constructor for objects in MOOSE
	ShrinkingCoreLinear(const InputParameters & parameters);
	
protected:
	/// Function to compute equilibria adsorption as a linear function
	void computeLinearIsotherm ();
	/// Function to compute the concentration derivative of qe
	Real computeConcentrationOffDiag ();
	/// Function to compute rate function with respect to coupled variable
	Real computeShrinkingCoreRateConcOffDiagJac ();
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
	
private:
	Real _partitioncoef; ///< Coefficient for a linear isotherm
	const VariableValue & _coupled_conc;	///< Coupled gas concentration variable
	const unsigned int _coupled_var_conc;	///< Variable identification for the coupled concentration variable
	
};
