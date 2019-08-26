/*!
 *  \file CoupledGSTAisotherm.h
 *	\brief Standard kernel for coupling non-linear variables via the GSTA isotherm model
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via the GSTA isotherm model.
 *
 *			GSTA isotherm: q = (q_max / m) * SUM(n*Kno*(p/Po)^n)/(1+SUM(Kno*(p/Po)^n))
 *			where q is amount adsorbed, q_max is the maximum capacity, m is the number of adsorption sites
 *			and Kno are the dimensionless equilibrium parameters. Also, p is partial pressure of gas and Po
 *			is taken as a reference state pressure (100 kPa).
 *
 *	\note	For the use of this kernel, our coupled variable with be a gas concentration in mol/L (C), therefore,
 *			we need to use ideal gas law to rewrite the GSTA model in terms of C as opposed to p. Thus, we are also
 *			forced to couple with kernel with column temperature.
 *
 *			Ideal Gas Law: p = C*R*T
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 08/24/2017
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

#include "CoupledGSTAisotherm.h"
registerMooseObject("dgospreyApp", CoupledGSTAisotherm);

template<>
InputParameters validParams<CoupledGSTAisotherm>()
{
	InputParameters params = validParams<TimeDerivative>();
	params.addParam<Real>("max_capacity",0.0,"Maximum Capacity for GSTA isotherm (mol/kg)");
	params.addParam<Real>("num_sites",1.0,"Number of adsorption sites for GSTA isotherm");
	params.addParam< std::vector<Real> >("gsta_params","Isotherm parameters for the GSTA isotherm");
	params.addRequiredCoupledVar("coupled_gas","Name of the gas species to couple with");
	params.addRequiredCoupledVar("coupled_temp","Name of the temperature variable to couple with");
	return params;
}

CoupledGSTAisotherm::CoupledGSTAisotherm(const InputParameters & parameters)
: TimeDerivative(parameters),
_maxcap(getParam<Real>("max_capacity")),
_numsites(getParam<Real>("num_sites")),
_gstaparam(getParam<std::vector<Real> >("gsta_params")),
_coupled_u(coupledValue("coupled_gas")),
_coupled_var_u(coupled("coupled_gas")),
_coupled_temp(coupledValue("coupled_temp")),
_coupled_var_temp(coupled("coupled_temp"))
{
	if (_numsites <= 0.0)
	{
		_numsites = 1.0;
		_maxcap = 0.0;
	}
	
	if (_maxcap < 0.0)
		_maxcap = 0.0;
}

Real CoupledGSTAisotherm::computeGSTAequilibrium()
{
	double top = 0.0, bot = 1.0, Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int n = 0; n<(int)_numsites; n++)
	{
		top = top + ( (double)(n+1) * _gstaparam[n] * std::pow((_coupled_u[_qp]/Co),(double)(n+1)) );
		bot = bot + ( _gstaparam[n] * std::pow((_coupled_u[_qp]/Co),(double)(n+1)) );
	}
	
	return (_maxcap/_numsites)*(top/bot);
}

Real CoupledGSTAisotherm::computeGSTAconcDerivative()
{
	double a = 0.0, b = 1.0, c = 0.0, d = 0.0, Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int n = 0; n<(int)_numsites; n++)
	{
		d = d + ( (double)(n+1) * _gstaparam[n] * std::pow((_coupled_u[_qp]/Co),(double)(n+1)) );
		b = b + ( _gstaparam[n] * std::pow((_coupled_u[_qp]/Co),(double)(n+1)) );
		
		a = a + ( (double)(n+1) * (double)(n+1) * _gstaparam[n] * std::pow((1.0/Co),(double)(n+1)) * std::pow((_coupled_u[_qp]),(double)(n)) );
		c = c + ( (double)(n+1) * _gstaparam[n] * std::pow((1.0/Co),(double)(n+1)) * std::pow((_coupled_u[_qp]),(double)(n)) );
	}

	return (_maxcap/_numsites)*_phi[_j][_qp]*( ((a*b) - (c*d)) / (b*b) );
}

Real CoupledGSTAisotherm::computeQpResidual()
{
	return _u[_qp]*_test[_i][_qp]-_test[_i][_qp]*computeGSTAequilibrium();
}

Real CoupledGSTAisotherm::computeQpJacobian()
{
	return _phi[_j][_qp]*_test[_i][_qp];
}

Real CoupledGSTAisotherm::computeQpOffDiagJacobian(unsigned int jvar)
{
	// Off-diagonal element for coupled gas
	if (jvar == _coupled_var_u)
	{
		return -_test[_i][_qp]*computeGSTAconcDerivative();
	}
	
	// Off-diagonal element for coupled temperature
	if (jvar == _coupled_var_temp)
	{
		double a = 0.0, b = 1.0, c = 0.0, d = 0.0, Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
		
		for (int n = 0; n<(int)_numsites; n++)
		{
			d = d + ( (double)(n+1) * _gstaparam[n] * std::pow((_coupled_u[_qp]/Co),(double)(n+1)) );
			b = b + ( _gstaparam[n] * std::pow((_coupled_u[_qp]/Co),(double)(n+1)) );
			
			a = a + ( (double)(n+1) * (double)(n+1) * _gstaparam[n] * std::pow((_coupled_u[_qp]*8.3144621/100.0),(double)(n+1)) * std::pow((_coupled_temp[_qp]),(double)(n)) );
			c = c + ( (double)(n+1) * _gstaparam[n] * std::pow((_coupled_u[_qp]*8.3144621/100.0),(double)(n+1)) * std::pow((_coupled_temp[_qp]),(double)(n)) );
		}
		
		return -_test[_i][_qp]*(_maxcap/_numsites)*_phi[_j][_qp]*( ((a*b) - (c*d)) / (b*b) );
	}
	
	
	return 0.0;
}


