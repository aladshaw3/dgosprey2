/*!
 *  \file GAdvection.h
 *	\brief Kernel for use with the corresponding DGAdvection object
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

#include "GAdvection.h"
registerMooseObject("dgospreyApp", GAdvection);

template<>
InputParameters validParams<GAdvection>()
{
	InputParameters params = validParams<Kernel>();
	params.addParam<Real>("vx",0, "x-component of velocity vector");
	params.addParam<Real>("vy",0,"y-component of velocity vector");
	params.addParam<Real>("vz",0,"z-component of velocity vector");
	return params;
}

GAdvection::GAdvection(const InputParameters & parameters) :
Kernel(parameters),
_vx(getParam<Real>("vx")),
_vy(getParam<Real>("vy")),
_vz(getParam<Real>("vz"))

{
	_velocity(0)=_vx;
	_velocity(1)=_vy;
	_velocity(2)=_vz;
}

Real GAdvection::computeQpResidual()
{
	return -_u[_qp]*(_velocity*_grad_test[_i][_qp]);
}

Real GAdvection::computeQpJacobian()
{
	return -_phi[_j][_qp]*(_velocity*_grad_test[_i][_qp]);
}
