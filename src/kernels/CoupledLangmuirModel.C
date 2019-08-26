/*!
 *  \file CoupledLangmuirModel.h
 *	\brief Standard kernel for coupling non-linear variables via a langmuir model
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via a langmuir forcing function, i.e., variable = b * K * coupled_variable / 1 + K * coupled_variable.
 *			In this kernel, the Langmuir parameter (K) is a function of the non-linear variable temperature through
 *			the van't Hoff expression: ln(K) = -dH/(R*T) + dS/R
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 10/09/2017
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

#include "CoupledLangmuirModel.h"
registerMooseObject("dgospreyApp", CoupledLangmuirModel);

template<>
InputParameters validParams<CoupledLangmuirModel>()
{
	InputParameters params = validParams<CoupledLangmuirForcingFunction>();
	params.addParam<Real>("enthalpy",0.0,"Langmuir Reaction Enthalpy [dH] (J/mol)");
	params.addParam<Real>("entropy",0.0,"Langmuir Reaction Entropy [dS] (J/K/mol)");
	params.addRequiredCoupledVar("coupled_temp","Name of the temperature variable being coupled");
	return params;
}

CoupledLangmuirModel::CoupledLangmuirModel(const InputParameters & parameters)
: CoupledLangmuirForcingFunction(parameters),
_enthalpy(getParam<Real>("enthalpy")),
_entropy(getParam<Real>("entropy")),
_coupled_temp(coupledValue("coupled_temp")),
_coupled_var_temp(coupled("coupled_temp"))
{
	
}

void CoupledLangmuirModel::computeLangmuirCoefficient()
{
	_langmuircoef = std::exp( lnKo(_enthalpy, _entropy, _coupled_temp[_qp]) );
}

Real CoupledLangmuirModel::computeLangTempDerivative()
{
	return _maxcap*((_langmuircoef*_phi[_j][_qp]*_coupled_u[_qp]*(_enthalpy/Rstd/_coupled_temp[_qp]/_coupled_temp[_qp]))/((1.0+_langmuircoef*_coupled_u[_qp])*(1.0+_langmuircoef*_coupled_u[_qp])));
}

Real CoupledLangmuirModel::computeQpResidual()
{
	computeLangmuirCoefficient();
	return CoupledLangmuirForcingFunction::computeQpResidual();
}

Real CoupledLangmuirModel::computeQpJacobian()
{
	return _phi[_j][_qp]*_test[_i][_qp];
}

Real CoupledLangmuirModel::computeQpOffDiagJacobian(unsigned int jvar)
{
	computeLangmuirCoefficient();
	
	//Off-diagonal element for coupled concentration
	if (jvar == _coupled_var)
		return -_test[_i][_qp]*CoupledLangmuirForcingFunction::computeLangConcDerivative();
	
	// Off-diagonal element for coupled temperature
	if (jvar == _coupled_var_temp)
		return -_test[_i][_qp]*computeLangTempDerivative();
	
	return 0.0;
}

