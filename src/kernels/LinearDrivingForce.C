/*!
 *  \file LinearDrivingForce.h
 *	\brief Standard kernel for a generic linear driving force mechanism
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

#include "LinearDrivingForce.h"
registerMooseObject("dgospreyApp", LinearDrivingForce);

template<>
InputParameters validParams<LinearDrivingForce>()
{
	InputParameters params = validParams<Kernel>();
	params.addParam<bool>("gaining",true,"True if driving force is a gaining term and false if driving force is a loss term");
	params.addParam<Real>("ldf_coef",1.0,"Coefficient multiplied by driving force");
	params.addParam<Real>("driving_value",1.0,"Value of driving force for the conserved quantity");
  return params;
}


LinearDrivingForce::LinearDrivingForce(const InputParameters & parameters)
  :Kernel(parameters),
   _gaining(getParam<bool>("gaining")),
   _coef(getParam<Real>("ldf_coef")),
   _driving_value(getParam<Real>("driving_value"))
{
	if (_gaining == true)
		_coef = -_coef;
}

Real LinearDrivingForce::computeQpResidual()
{
	return _test[_i][_qp] * _coef * (_driving_value - _u[_qp]);
}

Real LinearDrivingForce::computeQpJacobian()
{
  return -_phi[_j][_qp]*_test[_i][_qp]*_coef;
}

