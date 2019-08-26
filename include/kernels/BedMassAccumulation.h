/*!
 *  \file BedMassAccumulation.h
 *	\brief Time Derivative kernel for the accumulation of mass of a species in a fixed-bed column
 *	\details This file creates a time derivative kernel to be used in the mass transport equations for adsorption
 *			in a fixed-bed column. It combines the retardation coefficient from a material property with the standard
 *			time derivative kernel object in MOOSE.
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

#include "TimeDerivative.h"

/// BedMassAccumulation class object forward declarations
class BedMassAccumulation;

template<>
InputParameters validParams<BedMassAccumulation>();

/// BedMassAccumulation class object inherits from TimeDerivative object
/** This class object inherits from the TimeDerivative object.
	All public and protected members of this class are required function overrides.
	The kernel interfaces with the mass retardation coefficient calculated in a materials property
	file and calls the standard TimeDerivative functions while appending the retardation coefficient
	to those values. */
class BedMassAccumulation : public TimeDerivative
{
public:
	/// Required constructor for objects in MOOSE
	BedMassAccumulation(const InputParameters & parameters);
  
protected:
	/// Required residual function for standard kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	/// Required Jacobian function for standard kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
  
private:
	const MaterialProperty<Real> & _retardation;	///< Reference to the mass retardation coefficient material property
};
