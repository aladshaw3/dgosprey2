/*!
 *  \file WallHeatAccumulation.h
 *	\brief Time Derivative kernel for the accumulation of heat in a walls of the column
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

#include "WallHeatAccumulation.h"
registerMooseObject("dgospreyApp", WallHeatAccumulation);

template<>
InputParameters validParams<WallHeatAccumulation>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}


WallHeatAccumulation::WallHeatAccumulation(const InputParameters & parameters)
:TimeDerivative(parameters),
_wall_density(getMaterialProperty<Real>("wall_density")),
_wall_heat_capacity(getMaterialProperty<Real>("wall_heat_capacity"))
{
  
}

Real WallHeatAccumulation::computeQpResidual()
{
  return _wall_density[_qp] * _wall_heat_capacity[_qp] * TimeDerivative::computeQpResidual();
}

Real WallHeatAccumulation::computeQpJacobian()
{
  return _wall_density[_qp] * _wall_heat_capacity[_qp] * TimeDerivative::computeQpJacobian();
}
