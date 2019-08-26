/*!
 *  \file CoupledGSTAmodel.h
 *	\brief Standard kernel for coupling non-linear variables via the GSTA model
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via the GSTA model, which includes estimating model parameters from material
 *			properties files. The primary difference between this kernel and the CoupledGSTAisotherm
 *			kernel is the coupling of the material properties.
 *
 *			This kernel extends the CoupledGSTAisotherm kernel by calculating the model parameters from
 *			information in the ThermodynamicProperties material. The Kno parameters (described below) are
 *			to be estimated from the site enthalpies (dHno) and entropies (dSno) using the van't Hoff
 *			expression (shown below). Thus, this model is inherently a function of temperature and will
 *			require a different form of coupling with the temperature parameter.
 *
 *			van't Hoff: ln(Kno) = -dHno/(R*T) + dSno/R
 *			where R is the gas law constant and T is the column temperature.
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
 *	\date 08/28/2017
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

#include "CoupledGSTAmodel.h"
registerMooseObject("dgospreyApp", CoupledGSTAmodel);

template<>
InputParameters validParams<CoupledGSTAmodel>()
{
	InputParameters params = validParams<CoupledGSTAisotherm>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

CoupledGSTAmodel::CoupledGSTAmodel(const InputParameters & parameters)
: CoupledGSTAisotherm(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data"))
{
	
}

void CoupledGSTAmodel::computeGSTAparams()
{
	_gstaparam.resize(_magpie_dat[_qp].gsta_dat[_index].m);
	_numsites = _magpie_dat[_qp].gsta_dat[_index].m;
	_maxcap = _magpie_dat[_qp].gsta_dat[_index].qmax;
	
	for (int n = 0; n<(int)_numsites; n++)
	{
		_gstaparam[n] = std::exp( lnKo(_magpie_dat[_qp].gsta_dat[_index].dHo[n], _magpie_dat[_qp].gsta_dat[_index].dSo[n], _coupled_temp[_qp]) );
	}
}

Real CoupledGSTAmodel::computeGSTAtempDerivative()
{
	double a = 0.0, b = 1.0, c = 0.0, d = 0.0, Co = 100.0 / (8.3144621 * _coupled_temp[_qp]);
	
	for (int n = 0; n<(int)_numsites; n++)
	{
		d = d + ( (double)(n+1) * _gstaparam[n] * std::pow((_coupled_u[_qp]/Co),(double)(n+1)) );
		b = b + ( _gstaparam[n] * std::pow((_coupled_u[_qp]/Co),(double)(n+1)) );
		
		a = a + ( (double)(n+1) * _gstaparam[n] * std::pow((_coupled_u[_qp]*8.3144621/100.0),(double)(n+1)) * (( (double)(n+1) * std::pow((_coupled_temp[_qp]),(double)(n)) ) + (_magpie_dat[_qp].gsta_dat[_index].dHo[n]/8.3144621 * std::pow((_coupled_temp[_qp]),(double)(n-1)) )) );
		c = c + ( _gstaparam[n] * std::pow((_coupled_u[_qp]*8.3144621/100.0),(double)(n+1)) * (( (double)(n+1) * std::pow((_coupled_temp[_qp]),(double)(n)) ) + (_magpie_dat[_qp].gsta_dat[_index].dHo[n]/8.3144621 * std::pow((_coupled_temp[_qp]),(double)(n-1)) )) );
	}
	
	return (_maxcap/_numsites)*_phi[_j][_qp]*( ((a*b) - (c*d)) / (b*b) );
}

Real CoupledGSTAmodel::computeQpResidual()
{
	computeGSTAparams();
	return CoupledGSTAisotherm::computeQpResidual();
}

Real CoupledGSTAmodel::computeQpJacobian()
{
	return CoupledGSTAisotherm::computeQpJacobian();
}

Real CoupledGSTAmodel::computeQpOffDiagJacobian(unsigned int jvar)
{
	computeGSTAparams();
	
	// Off-diagonal element for coupled gas
	if (jvar == _coupled_var_u)
	{
		return CoupledGSTAisotherm::computeQpOffDiagJacobian(jvar);
	}
	
	// Off-diagonal element for coupled temperature
	if (jvar == _coupled_var_temp)
	{
		return -_test[_i][_qp]*computeGSTAtempDerivative();
	}
	
	
	return 0.0;
}


