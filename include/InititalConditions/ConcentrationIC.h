/*!
 *  \file ConcentrationIC.h
 *	\brief Initial Condition kernel for initial concentration of a species in a fixed-bed column
 *	\details This file creates an initial condition for the concentration of a species in the bed. The initial condition for
 *			concentration is assumed a constant value at all points in the bed. However, this can be modified later to include
 *			varying initial conditions.
 *
 *	\note If you want to have spatially varying initial conditions, you will need to modify the virtual value function of
 *		this kernel. Otherwise, it is assumed that the non-linear variable is initially constant at all points in the domain.
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

#include "InitialCondition.h"

/// ConcentrationIC class object forward declarations
class ConcentrationIC;

template<> InputParameters validParams<ConcentrationIC>();

/// ConcentrationIC class object inherits from InitialCondition object
/** This class object inherits from the InitialCondition object in the MOOSE framework.
	All public and protected members of this class are required function overrides. The object
	will establish the initial conditions for a species' concentration as constant throughout the
	domain.
 
	\note You can have the non-linear variable vary spatially in the domain by inheriting from
	and or modifying this file to do so. */
class ConcentrationIC : public InitialCondition
{
public:
	/// Required constructor for objects in MOOSE
	ConcentrationIC(const InputParameters & parameters);
	/// Required function override for setting the value of the non-linear variable at a given point
	/** This function passes a point p as an argument. The return value will be the value of the non-linear
		variable at that point. That information is used to establish the spatially varying initial condition
		for the given non-linear variable. */
	virtual Real value(const Point & p);
  
private:
  Real _y_IC;				///< Initial molefraction of the species in the gas phase
  Real _PT_IC;				///< Initial total pressure in the column (kPa)
  Real _T_IC;				///< Initial temperature in the column (K)
};
