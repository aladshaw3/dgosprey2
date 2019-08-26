/*!
 *  \file CoupledExtendedGSTAisotherm.h
 *	\brief Standard kernel for coupling a vector non-linear variables via an extended gsta isotherm
 *	\details This file creates a standard MOOSE kernel for the coupling of a vector non-linear variables
 *			together via an extended gsta forcing function...
 *
 *			                                    SUM(n, Kn_i * (coupled_variable_i/Co)^n)
 *					variable = qmax_i * ----------------------------------------------------
 *                                       1 + SUM(j, SUM(n, Kn_j * (coupled_variable_j/Co)^n)
 *
 *						where Co = 100.0 / (8.3144621 * _coupled_temp)
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 05/10/2018
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2018, all
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

#include "CoupledExtendedGSTAisotherm.h"
registerMooseObject("dgospreyApp", CoupledExtendedGSTAisotherm);

template<>
InputParameters validParams<CoupledExtendedGSTAisotherm>()
{
	InputParameters params = validParams<Kernel>();
	params.addParam<Real>("max_capacity",0.0,"Maximum Capacity for Langmuir Function (mol/kg)");
	params.addParam< std::vector<int> >("num_sites","Number of adsorption sites for each species");
	params.addParam< std::vector<Real> >("gsta_param_1","1st GSTA param for each variable");
	params.addParam< std::vector<Real> >("gsta_param_2","2nd GSTA param for each variable");
	params.addParam< std::vector<Real> >("gsta_param_3","3rd GSTA param for each variable");
	params.addParam< std::vector<Real> >("gsta_param_4","4th GSTA param for each variable");
	params.addParam< std::vector<Real> >("gsta_param_5","5th GSTA param for each variable");
	params.addParam< std::vector<Real> >("gsta_param_6","6th GSTA param for each variable");
	params.addRequiredCoupledVar("coupled_list","List of names of the variables being coupled");
	params.addRequiredCoupledVar("main_coupled","Name of the primary variable being coupled");
	params.addRequiredCoupledVar("coupled_temp","Name of the temperature variable to couple with");
	return params;
}

CoupledExtendedGSTAisotherm::CoupledExtendedGSTAisotherm(const InputParameters & parameters)
: Kernel(parameters),
_maxcap(getParam<Real>("max_capacity")),
_num_sites(getParam<std::vector<int> >("num_sites")),
_param_1(getParam<std::vector<Real> >("gsta_param_1")),
_param_2(getParam<std::vector<Real> >("gsta_param_2")),
_param_3(getParam<std::vector<Real> >("gsta_param_3")),
_param_4(getParam<std::vector<Real> >("gsta_param_4")),
_param_5(getParam<std::vector<Real> >("gsta_param_5")),
_param_6(getParam<std::vector<Real> >("gsta_param_6")),
_coupled_i(coupledValue("main_coupled")),
_coupled_var_i(coupled("main_coupled")),
_coupled_temp(coupledValue("coupled_temp")),
_coupled_var_temp(coupled("coupled_temp"))
{
	
	if (_maxcap < 0.0)
		_maxcap = 0.0;
	
	unsigned int n = coupledComponents("coupled_list");
	_coupled_vars.resize(n);
	_coupled.resize(n);
	_gstaparams.resize(n);
	_num_sites.resize(n);
	_param_1.resize(n);
	_param_2.resize(n);
	_param_3.resize(n);
	_param_4.resize(n);
	_param_5.resize(n);
	_param_6.resize(n);
	
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		_coupled_vars[i] = coupled("coupled_list",i);
		_coupled[i] = &coupledValue("coupled_list",i);
		if (_coupled_vars[i] == _coupled_var_i)
			_main_index = i;
		
		if (_num_sites[i] <= 0.0)
		{
			_num_sites[i] = 1.0;
		}
	}
	
	if (_num_sites[_main_index] <= 0)
		_maxcap = 0.0;
	
	//Initialize gstaparams
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		_gstaparams[i].resize(_num_sites[i]);
		
		for (int n=0; n<_num_sites[i]; n++)
		{
			if (n == 0)
			{
				_gstaparams[i][n] = _param_1[i];
			}
			else if (n == 1)
			{
				_gstaparams[i][n] = _param_2[i];
			}
			else if (n == 2)
			{
				_gstaparams[i][n] = _param_3[i];
			}
			else if (n == 3)
			{
				_gstaparams[i][n] = _param_4[i];
			}
			else if (n == 4)
			{
				_gstaparams[i][n] = _param_5[i];
			}
			else if (n == 5)
			{
				_gstaparams[i][n] = _param_6[i];
			}
			else
			{
				_gstaparams[i][n] = 0.0;
			}
		}
	}
}

Real CoupledExtendedGSTAisotherm::computeExtGSTAEquilibrium()
{
	double top = computeTopExtGSTA();
	double bot = computeBottomExtGSTA();
	return (_maxcap/_num_sites[_main_index]) * (top/bot);
}

Real CoupledExtendedGSTAisotherm::computeExtGSTAConcJacobi()
{
	double top = computeTopDerivativeConc(_main_index)*computeBottomExtGSTA() - computeTopExtGSTA()*computeBottomDerivativeConc(_main_index);
	double bot = computeBottomExtGSTA()*computeBottomExtGSTA();
	return (_maxcap/_num_sites[_main_index]) * _phi[_j][_qp]*(top/bot);
}

Real CoupledExtendedGSTAisotherm::computeExtGSTA_ConcOffJacobi(int i)
{
	double top = computeTopDerivativeConc(i)*computeBottomExtGSTA() - computeTopExtGSTA()*computeBottomDerivativeConc(i);
	double bot = computeBottomExtGSTA()*computeBottomExtGSTA();
	return (_maxcap/_num_sites[_main_index]) * _phi[_j][_qp]*(top/bot);
}

Real CoupledExtendedGSTAisotherm::computeExtGSTA_TempOffJacobi()
{
	double top = computeTopDerivativeTemp()*computeBottomExtGSTA() - computeTopExtGSTA()*computeBottomDerivativeTemp();
	double bot = computeBottomExtGSTA()*computeBottomExtGSTA();
	return (_maxcap/_num_sites[_main_index]) * _phi[_j][_qp]*(top/bot);
}

Real CoupledExtendedGSTAisotherm::computeTopExtGSTA()
{
	double sum = 0.0;
	double Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int n=0; n<_num_sites[_main_index]; n++)
	{
		sum += (double)(n+1) * _gstaparams[_main_index][n] * std::pow((_coupled_i[_qp]/Co),(double)(n+1));
	}
	
	return sum;
}

Real CoupledExtendedGSTAisotherm::computeBottomExtGSTA()
{
	double sum = 1.0;
	double Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int i=0; i<_coupled.size(); ++i)
	{
		for (int n=0; n<_num_sites[i]; n++)
		{
			sum += _gstaparams[i][n] * std::pow(((*_coupled[i])[_qp]/Co),(double)(n+1));
		}
	}
	
	return sum;
}

Real CoupledExtendedGSTAisotherm::computeTopDerivativeConc(int i)
{
	double sum = 0.0;
	double Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	if (i == _main_index)
	{
		for (int n=0; n<_num_sites[_main_index]; n++)
		{
			sum += (double)(n+1) * (double)(n+1) * _gstaparams[_main_index][n] * std::pow((_coupled_i[_qp]/Co),(double)(n)) / Co;
		}
	}
	else
	{
		return 0.0;
	}
	
	return sum;
}

Real CoupledExtendedGSTAisotherm::computeBottomDerivativeConc(int i)
{
	double sum = 0.0;
	double Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int n=0; n<_num_sites[i]; n++)
	{
		sum += (double)(n+1) * _gstaparams[i][n] * std::pow(((*_coupled[i])[_qp]/Co),(double)(n)) / Co;
	}
	
	return sum;
}

Real CoupledExtendedGSTAisotherm::computeTopDerivativeTemp()
{
	double sum = 0.0;
	double Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int n=0; n<_num_sites[_main_index]; n++)
	{
		sum += (double)(n+1) * (double)(n+1) * _gstaparams[_main_index][n] * std::pow((_coupled_i[_qp]/Co),(double)(n)) * (_coupled_i[_qp]*8.3144621/100.0);
	}
	
	return sum;
}

Real CoupledExtendedGSTAisotherm::computeBottomDerivativeTemp()
{
	double sum = 0.0;
	double Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int i=0; i<_coupled.size(); ++i)
	{
		for (int n=0; n<_num_sites[i]; n++)
		{
			sum += (double)(n+1) * _gstaparams[i][n] * std::pow(((*_coupled[i])[_qp]/Co),(double)(n)) * ((*_coupled[i])[_qp]*8.3144621/100.0);
		}
	}
	
	return sum;
}

Real CoupledExtendedGSTAisotherm::computeQpResidual()
{
	return _u[_qp]*_test[_i][_qp]-_test[_i][_qp]*computeExtGSTAEquilibrium();
}

Real CoupledExtendedGSTAisotherm::computeQpJacobian()
{
	return _phi[_j][_qp]*_test[_i][_qp];
}

Real CoupledExtendedGSTAisotherm::computeQpOffDiagJacobian(unsigned int jvar)
{
	// Off-diagonal element for all other coupled variables
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		if (jvar == _coupled_vars[i] && jvar != _coupled_var_i)
		{
			return -_test[_i][_qp]*computeExtGSTA_ConcOffJacobi(i);
		}
	}
	
	// Off-diagonal element for this coupled variable
	if (jvar == _coupled_var_i)
	{
		return -_test[_i][_qp]*computeExtGSTAConcJacobi();
	}
	
	// Off-diagonal element for coupled temperature
	if (jvar == _coupled_var_temp)
	{
		return -_test[_i][_qp]*computeExtGSTA_TempOffJacobi();
	}
	
	return 0.0;
}


