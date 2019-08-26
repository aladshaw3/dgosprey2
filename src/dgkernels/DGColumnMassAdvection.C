/*!
 *  \file DGColumnMassAdvection.h
 *	\brief Discontinous Galerkin kernel for mass advection in a fixed-bed column
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

#include "DGColumnMassAdvection.h"
registerMooseObject("dgospreyApp", DGColumnMassAdvection);

template<>
InputParameters validParams<DGColumnMassAdvection>()
{
	InputParameters params = validParams<DGAdvection>();
	return params;
}

DGColumnMassAdvection::DGColumnMassAdvection(const InputParameters & parameters) :
DGAdvection(parameters),
_vel(getMaterialProperty<Real>("velocity"))
{

}

Real DGColumnMassAdvection::computeQpResidual(Moose::DGResidualType type)
{
	
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp];
	_velocity(2)=0.0;
	
	return DGAdvection::computeQpResidual(type);
}

Real DGColumnMassAdvection::computeQpJacobian(Moose::DGJacobianType type)
{
	
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp];
	_velocity(2)=0.0;
	
	return DGAdvection::computeQpJacobian(type);
}
