/*!
 *  \file CoupledExtendedLangmuirFunction.h
 *	\brief Standard kernel for coupling a vector non-linear variables via an extended langmuir function
 *	\details This file creates a standard MOOSE kernel for the coupling of a vector non-linear variables
 *			together via an extended langmuir forcing function,
 *			i.e., variable = b_i * K_i * coupled_variable_i / 1 + sum(j, K_j * coupled_variable_j).
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 06/21/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
 *             rights reserved.
 *
 *			   Alexander Wiechert does not claim any ownership or copyright to the
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

#include "CoupledExtendedLangmuirFunction.h"
registerMooseObject("dgospreyApp", CoupledExtendedLangmuirFunction);

template<>
InputParameters validParams<CoupledExtendedLangmuirFunction>()
{
	InputParameters params = validParams<TimeDerivative>();
	params.addParam<Real>("max_capacity",0.0,"Maximum Capacity for Langmuir Function (mol/kg)");
	params.addParam< std::vector<Real> >("langmuir_coeff","Coefficient for the langmuir function (L/mol)");
	params.addRequiredCoupledVar("coupled_list","List of names of the variables being coupled");
	params.addRequiredCoupledVar("main_coupled","Name of the primary variable being coupled");
	return params;
}

CoupledExtendedLangmuirFunction::CoupledExtendedLangmuirFunction(const InputParameters & parameters)
: TimeDerivative(parameters),
_maxcap(getParam<Real>("max_capacity")),
_langmuircoef(getParam<std::vector<Real> >("langmuir_coeff")),
_coupled_i(coupledValue("main_coupled")),
_coupled_var_i(coupled("main_coupled"))
{
	unsigned int n = coupledComponents("coupled_list");
	_coupled_vars.resize(n);
	_coupled.resize(n);
	
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		_coupled_vars[i] = coupled("coupled_list",i);
		_coupled[i] = &coupledValue("coupled_list",i);
		if (_coupled_vars[i] == _coupled_var_i)
			_lang_index = i;
	}

}

Real CoupledExtendedLangmuirFunction::computeExtLangmuirEquilibrium()
{
	double sum = 0.0;
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		sum = sum + _langmuircoef[i] * (*_coupled[i])[_qp];
	}
	return _maxcap*((_langmuircoef[_lang_index]*_coupled_i[_qp])/(1.0+sum));
}

Real CoupledExtendedLangmuirFunction::computeExtLangmuirConcJacobi()
{
	double sum = 0.0;
	for (unsigned int j = 0; j<_coupled.size(); ++j)
	{
		sum = sum + _langmuircoef[j] * (*_coupled[j])[_qp];
	}
	double numerator = _langmuircoef[_lang_index]*_phi[_j][_qp] * (1.0 + sum - _langmuircoef[_lang_index]*_coupled_i[_qp]);
	double denom = (1.0+sum)*(1.0+sum);
	return _maxcap*numerator/denom;
}

Real CoupledExtendedLangmuirFunction::computeExtLangmuirOffJacobi(int i)
{
	double sum = 0.0;
	for (unsigned int j = 0; j<_coupled.size(); ++j)
	{
		sum = sum + _langmuircoef[j] * (*_coupled[j])[_qp];
	}
	double numerator = _langmuircoef[_lang_index]*_coupled_i[_qp]*_phi[_j][_qp]*_langmuircoef[i];
	double denom = (1.0+sum)*(1.0+sum);
	return _maxcap*numerator/denom;

}

Real CoupledExtendedLangmuirFunction::computeQpResidual()
{
	return _u[_qp]*_test[_i][_qp]-_test[_i][_qp]*computeExtLangmuirEquilibrium();
}

Real CoupledExtendedLangmuirFunction::computeQpJacobian()
{
	return _phi[_j][_qp]*_test[_i][_qp];
}

Real CoupledExtendedLangmuirFunction::computeQpOffDiagJacobian(unsigned int jvar)
{
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		if (jvar == _coupled_vars[i] && jvar != _coupled_var_i)
		{
			return _test[_i][_qp]*computeExtLangmuirOffJacobi(i);
		}
	}
	
	if (jvar == _coupled_var_i)
	{
		return -_test[_i][_qp]*computeExtLangmuirConcJacobi();
	}
	
	return 0.0;
}


