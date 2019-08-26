/*!
 *  \file TotalPressureIC.h
 *	\brief Initial Condition kernel for initial temperature in a fixed-bed column
 *	\details This file creates an initial condition for the temperature in the bed. The initial condition for temperature
 *			is assumed a constant value at all points in the bed. However, this can be modified later to include spatially
 *			varying initial conditions for temperature.
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

/// TotalPressureIC class object forward declarations
class TotalPressureIC;

template<> InputParameters validParams<TotalPressureIC>();

/// TotalPressureIC class object inherits from InitialCondition object
/** This class object inherits from the InitialCondition object in the MOOSE framework.
	All public and protected members of this class are required function overrides. The object
	will establish the initial conditions for total pressure as constant throughout the
	domain.
 
	\note You can have the non-linear variable vary spatially in the domain by inheriting from
	and or modifying this file to do so. */
class TotalPressureIC : public InitialCondition
{
public:
	/// Required constructor for objects in MOOSE
	TotalPressureIC(const InputParameters & parameters);
	/// Required function override for setting the value of the non-linear variable at a given point
	/** This function passes a point p as an argument. The return value will be the value of the non-linear
		variable at that point. That information is used to establish the spatially varying initial condition
		for the given non-linear variable. */
	virtual Real value(const Point & p);
  
private:
	Real _PT_IC;			///< Initial condition for the total pressure in the column (kPa)
};
