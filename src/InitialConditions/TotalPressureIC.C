/*!
 *  \file TotalPressureIC.h
 *	\brief Initial Condition kernel for initial temperature in a fixed-bed column
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

#include "TotalPressureIC.h"
registerMooseObject("dgospreyApp", TotalPressureIC);

template<> InputParameters validParams<TotalPressureIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("initial_pressure","The value of the initial pressure in the column (kPa).");
  return params;
}

TotalPressureIC::TotalPressureIC(const InputParameters & parameters) :
InitialCondition(parameters),
_PT_IC(getParam<Real>("initial_pressure"))
{
  
}

Real TotalPressureIC::value(const Point & p)
{
  return _PT_IC;
}
