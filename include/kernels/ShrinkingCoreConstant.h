/*!
 *  \file ShrinkingCoreConstant.h
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

#include "TimeDerivative.h"
#include <float.h>

/// ShrinkingCoreConstant class object forward declarations
class ShrinkingCoreConstant;

template<>
InputParameters validParams<ShrinkingCoreConstant>();

/// ShrinkingCoreConstant class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel has all the protected members from Kernel, but also
	includes coefficients for the shrinking core model.
 */
class ShrinkingCoreConstant : public TimeDerivative
{
public:
	/// Required constructor for objects in MOOSE
	ShrinkingCoreConstant(const InputParameters & parameters);
	
protected:
	///Function to compute the fraction of adsorption
	Real computeAdsorbedFraction ();
	///Function to compute the Jacobian of the fraction of adsorption
	Real computeAdsorbedFractionJacobian ();
	///Function to compute the current rate of the shrinking core model
	Real computeShrinkingCoreRate ();
	///Function to compute Jacobian of rate function
	Real computeShrinkingCoreRateJacobian();
	/// Required residual function for standard kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	/// Required Jacobian function for standard kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
	
	Real _maxcap;			///< Coefficient for the max capacity
	Real _tau1;				///< Coefficient for the time to reach equilibrium if gas film mass transfer is controlling
	Real _tau2;				///< Coefficient for the time to reach equilibrium if gas pore diffusion is controlling
	Real _tau3;				///< Coefficient for the time to reach equilibrium if gas-solid reaction is controlling
	
private:
	
};
