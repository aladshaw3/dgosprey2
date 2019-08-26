/*!
 *  \file CoupledLinearForcingFunction.h
 *	\brief Standard kernel for coupling non-linear variables via a linear function
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via a linear forcing function, i.e., variable = K * coupled_variable .
 *
 *  \author Austin Ladshaw
 *	\date 03/30/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
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

#include "CoupledLinearForcingFunction.h"
registerMooseObject("dgospreyApp", CoupledLinearForcingFunction);

template<>
InputParameters validParams<CoupledLinearForcingFunction>()
{
	InputParameters params = validParams<Kernel>();
	params.addParam<bool>("gaining",true,"If coupled function is positive, then gaining = true");
	params.addParam<Real>("coeff",1.0,"Coefficient for the linear function");
	params.addRequiredCoupledVar("coupled","Name of the variable being coupled");
	return params;
}

CoupledLinearForcingFunction::CoupledLinearForcingFunction(const InputParameters & parameters)
: Kernel(parameters),
_gaining(getParam<bool>("gaining")),
_coef(getParam<Real>("coeff")),
_coupled_u(coupledValue("coupled")),
_coupled_var(coupled("coupled"))
{
	if (_gaining == false)
		_coef = -_coef;
}

Real CoupledLinearForcingFunction::computeQpResidual()
{
	return _u[_qp]*_test[_i][_qp] - _coef*_coupled_u[_qp]*_test[_i][_qp];
}

Real CoupledLinearForcingFunction::computeQpJacobian()
{
	return _phi[_j][_qp]*_test[_i][_qp];
}

Real CoupledLinearForcingFunction::computeQpOffDiagJacobian(unsigned int jvar)
{
    
	if (jvar == _coupled_var)
		return -_coef*_phi[_j][_qp]*_test[_i][_qp];
    
	return 0.0;
}


