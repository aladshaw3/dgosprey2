/*!
 *  \file CoupledExtendedGSTAmodel.h
 *	\brief Standard kernel for coupling a vector non-linear variables via an extended gsta model
 *	\details This file creates a standard MOOSE kernel for the coupling of a vector non-linear variables
 *			together via an extended gsta forcing function and computes gsta parameters as a function
 *			of temperature by coupling with the Thermodynamic Properties Material property.
 *
 *			                                    SUM(n, Kn_i * (coupled_variable_i/Co)^n)
 *					variable = qmax_i * ----------------------------------------------------
 *                                       1 + SUM(j, SUM(n, Kn_j * (coupled_variable_j/Co)^n)
 *
 *						where Co = 100.0 / (8.3144621 * _coupled_temp)
 *						and  Kn_i = exp( -dHn_i/RT + dSn_i/R )
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

#include "CoupledExtendedGSTAmodel.h"
registerMooseObject("dgospreyApp", CoupledExtendedGSTAmodel);

template<>
InputParameters validParams<CoupledExtendedGSTAmodel>()
{
	InputParameters params = validParams<CoupledExtendedGSTAisotherm>();
	return params;
}

CoupledExtendedGSTAmodel::CoupledExtendedGSTAmodel(const InputParameters & parameters)
: CoupledExtendedGSTAisotherm(parameters),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data"))
{

}

void CoupledExtendedGSTAmodel::computeGSTAparams()
{
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		_gstaparams[i].resize(_magpie_dat[_qp].gsta_dat[i].m);
		_num_sites[i] = _magpie_dat[_qp].gsta_dat[i].m;
	}
	_maxcap = _magpie_dat[_qp].gsta_dat[_main_index].qmax;
	
	
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		for (int n = 0; n<_num_sites[i]; n++)
		{
			_gstaparams[i][n] = std::exp( lnKo(_magpie_dat[_qp].gsta_dat[i].dHo[n], _magpie_dat[_qp].gsta_dat[i].dSo[n], _coupled_temp[_qp]) );
		}
	}
}

Real CoupledExtendedGSTAmodel::computeExtGSTA_TempOffJacobi()
{
	double top = computeTopDerivativeTemp()*CoupledExtendedGSTAisotherm::computeBottomExtGSTA() - CoupledExtendedGSTAisotherm::computeTopExtGSTA()*computeBottomDerivativeTemp();
	double bot = CoupledExtendedGSTAisotherm::computeBottomExtGSTA()*CoupledExtendedGSTAisotherm::computeBottomExtGSTA();
	return (_maxcap/_num_sites[_main_index]) * _phi[_j][_qp]*(top/bot);
}

Real CoupledExtendedGSTAmodel::computeTopDerivativeTemp()
{
	double sum = 0.0;
	double Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int n=0; n<_num_sites[_main_index]; n++)
	{
		sum += (double)(n+1) * (double)(n+1) * _gstaparams[_main_index][n] * std::pow((_coupled_i[_qp]/Co),(double)(n)) * (_coupled_i[_qp]*8.3144621/100.0);
		sum += (double)(n+1) * _gstaparams[_main_index][n] * std::pow((_coupled_i[_qp]/Co),(double)(n+1)) * _magpie_dat[_qp].gsta_dat[_main_index].dHo[n] / 8.3144621 / std::pow(_coupled_temp[_qp], 2.0);
	}
	
	return sum;
}

Real CoupledExtendedGSTAmodel::computeBottomDerivativeTemp()
{
	double sum = 0.0;
	double Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int i=0; i<_coupled.size(); ++i)
	{
		for (int n=0; n<_num_sites[i]; n++)
		{
			sum += (double)(n+1) * _gstaparams[i][n] * std::pow(((*_coupled[i])[_qp]/Co),(double)(n)) * ((*_coupled[i])[_qp]*8.3144621/100.0);
			sum += _gstaparams[i][n] * std::pow(((*_coupled[i])[_qp]/Co),(double)(n+1)) * _magpie_dat[_qp].gsta_dat[i].dHo[n] / 8.3144621 / std::pow(_coupled_temp[_qp], 2.0);
		}
	}
	
	return sum;
}

Real CoupledExtendedGSTAmodel::computeQpResidual()
{
	computeGSTAparams();
	return CoupledExtendedGSTAisotherm::computeQpResidual();
}

Real CoupledExtendedGSTAmodel::computeQpJacobian()
{
	return CoupledExtendedGSTAisotherm::computeQpJacobian();
}

Real CoupledExtendedGSTAmodel::computeQpOffDiagJacobian(unsigned int jvar)
{
	computeGSTAparams();
	
	// Off-diagonal element for all other coupled variables
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		if (jvar == _coupled_vars[i] && jvar != _coupled_var_i)
		{
			return -_test[_i][_qp]*CoupledExtendedGSTAisotherm::computeExtGSTA_ConcOffJacobi(i);
		}
	}
	
	// Off-diagonal element for this coupled variable
	if (jvar == _coupled_var_i)
	{
		return -_test[_i][_qp]*CoupledExtendedGSTAisotherm::computeExtGSTAConcJacobi();
	}
	
	// Off-diagonal element for coupled temperature
	if (jvar == _coupled_var_temp)
	{
		return -_test[_i][_qp]*computeExtGSTA_TempOffJacobi();
	}
	
	
	return 0.0;
}
