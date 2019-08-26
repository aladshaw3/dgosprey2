/*!
 *  \file CoupledLinearLDF.h
 *	\brief Standard kernel for a coupled linear driving force mechanism
 *	\details This file creates a standard MOOSE kernel for a coupled linear driving force mechanism that
 *			can be added to the non-linear residuals. It contains all the same parameters as the more
 *			generic base class, but couples with another non-linear variable via a linear relationship.
 *			The linear coefficient for that relationship is added as an additional parameter to be
 *			set by the user.
 *
 *  \author Austin Ladshaw
 *	\date 03/30/2017
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

#include "CoupledLinearLDF.h"
registerMooseObject("dgospreyApp", CoupledLinearLDF);

template<>
InputParameters validParams<CoupledLinearLDF>()
{
	InputParameters params = validParams<LinearDrivingForce>();
	params.addParam<Real>("linear_coef",1.0,"Coefficient multiplied by the coupled variable");
	params.addRequiredCoupledVar("coupled","Name of the variable being coupled");
	return params;
}


CoupledLinearLDF::CoupledLinearLDF(const InputParameters & parameters)
:LinearDrivingForce(parameters),
_coupled_coef(getParam<Real>("linear_coef")),
_coupled_u(coupledValue("coupled"))
{
}

Real CoupledLinearLDF::computeQpResidual()
{
	_driving_value = _coupled_coef*_coupled_u[_qp];
	return LinearDrivingForce::computeQpResidual();
}

Real CoupledLinearLDF::computeQpJacobian()
{
	return LinearDrivingForce::computeQpJacobian();
}

