/*!
 *  \file LinearDrivingForce.h
 *	\brief Standard kernel for a generic linear driving force mechanism
 *	\details This file creates a standard MOOSE kernel for a linear driving force type of mechanism that
 *			can be added to the non-linear residuals. It contains a boolean argument to determine whether
 *			the driving force is gaining or losing, a coefficient for the rate of the driving force, and
 *			a driving value to where the non-linear coupled variable is heading towards.
 *
 *  \author Austin Ladshaw
 *	\date 11/20/2015
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2015, all
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

/// LinearDrivingForce class object forward declarations
class LinearDrivingForce;

template<>
InputParameters validParams<LinearDrivingForce>();

/// LinearDrivingForce class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel has several protected members including: a boolean for gaining or 
	losing mechanisms, a coefficient for the rate or strength of the driving force,
	a driving value to where the coupled non-linear variable is driving toward, and
	the coupled non-linear variable. 
 
	\note To create a specific linear driving force kernel, inherit from this class
	and use other non-linear variables or material properties to change the protected
	member values to reflect the physics for your problem. */
class LinearDrivingForce : public Kernel
{
public:
	/// Required constructor for objects in MOOSE
	LinearDrivingForce(const InputParameters & parameters);

protected:
	/// Required residual function for standard kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	/// Required Jacobian function for standard kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
	
	bool _gaining;			///< Boolean to mark whether the driving force is gaining or losing (True = gaining)
	Real _coef;				///< Coefficient for the strength or rate of the driving force
	Real _driving_value;	///< Value the coupled variable is driving towards

private:

};
