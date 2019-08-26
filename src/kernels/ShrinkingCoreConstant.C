/*!
 *  \file ShrinkingCoreConstant.h
 *	\brief Standard kernel for a standard costant shrinking core
 *	\details This file creates a standard MOOSE kernel for constant shrinking core mechanism that
 *			can be added to the non-linear residuals. It contains all the same parameters as the more
 *			generic base class, but couples with another non-linear variable via a linear relationship.
 *			The linear coefficient for that relationship is added as an additional parameter to be
 *			set by the user.
 *
 *  \author Alexander Wiechert, Austin Ladshaw
 *	\date 12/19/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2015, all
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

#include "ShrinkingCoreConstant.h"
registerMooseObject("dgospreyApp", ShrinkingCoreConstant);

template<>
InputParameters validParams<ShrinkingCoreConstant>()
{
	InputParameters params = validParams<TimeDerivative>();
	params.addParam<Real>("max_capacity",1.0,"Equilibrium adsorption capacity (mol/kg)");
	params.addParam<Real>("filmtrans_coef",1.0,"Time to reach equilibrium if gas film mass transfer is controlling (h)");
	params.addParam<Real>("porediff_coef",1.0,"Time to reach equilibrium if gas pore diffusion is controlling (h)");
	params.addParam<Real>("reaction_coef",1.0,"Time to reach equilibrium if gas-solid reaction is controlling (h)");
	return params;
}


ShrinkingCoreConstant::ShrinkingCoreConstant(const InputParameters & parameters)
:TimeDerivative(parameters),
_maxcap(getParam<Real>("max_capacity")),
_tau1(getParam<Real>("filmtrans_coef")),
_tau2(getParam<Real>("porediff_coef")),
_tau3(getParam<Real>("reaction_coef"))
{
}

Real ShrinkingCoreConstant::computeAdsorbedFraction ()
{
	if (_maxcap <= 0.0)
		return 1.0;
	return 1.0 - _u[_qp]/_maxcap;
}

Real ShrinkingCoreConstant::computeAdsorbedFractionJacobian ()
{
	if (_maxcap <= 0.0)
		return 0.0;
	return -1.0/_maxcap;
}

Real ShrinkingCoreConstant::computeShrinkingCoreRate ()
{
	double X = computeAdsorbedFraction ();
	if (X <= 0.0)
		return 0.0;
	return (3.0*_maxcap*X)/(3.0*(_tau1-2.0*_tau2)*X + 6.0*_tau2*std::pow(X,2.0/3.0) + _tau3*std::pow(X,1.0/3.0));
}

Real ShrinkingCoreConstant::computeShrinkingCoreRateJacobian ()
{
	double g, gp, h, hp;
	double X = computeAdsorbedFraction ();
	if (X <= 0.0)
		return 0.0;
	g = 3.0*_maxcap*X;
	h = 3.0*(_tau1-2.0*_tau2)*X + 6.0*_tau2*std::pow(X,2.0/3.0) + _tau3*std::pow(X,1.0/3.0);
	
	double Xp = computeAdsorbedFractionJacobian ();
	gp = 3.0 * _maxcap * Xp;
	hp = 3.0*(_tau1-2.0*_tau2)*Xp + 4.0*_tau2*std::pow(X,-1.0/3.0)*Xp + (1.0/3.0)*_tau3*std::pow(X,-2.0/3.0)*Xp;
	
	return (gp*h - hp*g)/(h*h);
}

Real ShrinkingCoreConstant::computeQpResidual()
{

	return _test[_i][_qp]*_u_dot[_qp] - computeShrinkingCoreRate()*_test[_i][_qp];
}

Real ShrinkingCoreConstant::computeQpJacobian()
{
	return _test[_i][_qp] * _phi[_j][_qp] * _du_dot_du[_qp] - _test[_i][_qp]*_phi[_j][_qp]*computeShrinkingCoreRateJacobian();
}

