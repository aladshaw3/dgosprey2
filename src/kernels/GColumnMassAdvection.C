/*!
 *  \file GColumnMassAdvection.h
 *	\brief Kernel for use with the corresponding DGColumnMassAdvection object
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

#include "GColumnMassAdvection.h"
registerMooseObject("dgospreyApp", GColumnMassAdvection);

template<>
InputParameters validParams<GColumnMassAdvection>()
{
	InputParameters params = validParams<GAdvection>();
	return params;
}

GColumnMassAdvection::GColumnMassAdvection(const InputParameters & parameters) :
GAdvection(parameters),
_vel(getMaterialProperty<Real>("velocity"))

{

}

Real GColumnMassAdvection::computeQpResidual()
{
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp];
	_velocity(2)=0.0;
	
	return GAdvection::computeQpResidual();
}

Real GColumnMassAdvection::computeQpJacobian()
{
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp];
	_velocity(2)=0.0;
	
	return GAdvection::computeQpJacobian();
}
