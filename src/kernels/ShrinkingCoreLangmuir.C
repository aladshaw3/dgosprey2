/*!
 *  \file ShrinkingCoreLangmuir.C
 *	\brief Standard kernel for a standard costant shrinking core
 *	\details This file creates a standard MOOSE kernel for constant shrinking core mechanism that
 *			can be added to the non-linear residuals. It contains all the same parameters as the more
 *			generic base class, but couples with another non-linear variable via a langmuir function.
 *			The langmuir coefficient for that relationship is added as an additional parameter to be
 *			set by the user.
 *
 *  \author Alexander Wiechert, Austin Ladshaw
 *	\date 1/11/2018
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

#include "ShrinkingCoreLangmuir.h"
registerMooseObject("dgospreyApp", ShrinkingCoreLangmuir);

template<>
InputParameters validParams<ShrinkingCoreLangmuir>()
{
    InputParameters params = validParams<ShrinkingCoreConstant>();
    params.addParam<Real>("langmuir_coef",1.0,"Equilibrium adsorption constant (L/kg)");
    params.addParam<Real>("total_capacity",1.0,"Maximum adsorption capacity (mol/kg)");
    params.addRequiredCoupledVar("coupled_conc","Concentration to couple adsorption variable");
    return params;
}


ShrinkingCoreLangmuir::ShrinkingCoreLangmuir(const InputParameters & parameters)
:ShrinkingCoreConstant(parameters),
_langmuircoef(getParam<Real>("langmuir_coef")),
_absmaxcap(getParam<Real>("total_capacity")),
_coupled_conc(coupledValue("coupled_conc")),
_coupled_var_conc(coupled("coupled_conc"))
{
}

void ShrinkingCoreLangmuir::computeLangmuirIsotherm ()
{
    _maxcap = _absmaxcap * (_langmuircoef * _coupled_conc[_qp]) / (1.0 + _langmuircoef * _coupled_conc[_qp]);
}

Real ShrinkingCoreLangmuir::computeConcentrationOffDiag ()
{
    if (_coupled_conc[_qp] <= 0.0)
        return 0.0;
    return _u[_qp]/(_absmaxcap * _langmuircoef * _coupled_conc[_qp] * _coupled_conc[_qp]);
}


Real ShrinkingCoreLangmuir::computeShrinkingCoreRateConcOffDiagJac ()
{
    computeLangmuirIsotherm ();
    double g, gp, h, hp;
    double X = computeAdsorbedFraction ();
    if (X <= 0.0)
        return 0.0;
    g = 3.0*_maxcap*X;
    h = 3.0*(_tau1-2.0*_tau2)*X + 6.0*_tau2*std::pow(X,2.0/3.0) + _tau3*std::pow(X,1.0/3.0);
    
    double Xp = computeConcentrationOffDiag ();
    gp = 3.0 * _maxcap * Xp - 3.0 * X * _absmaxcap * _langmuircoef / ((1.0 + _langmuircoef * _coupled_conc[_qp]) * (1.0 + _langmuircoef * _coupled_conc[_qp]));
    hp = 3.0*(_tau1-2.0*_tau2)*Xp + 4.0*_tau2*std::pow(X,-1.0/3.0)*Xp + (1.0/3.0)*_tau3*std::pow(X,-2.0/3.0)*Xp;
    
    return (gp*h - hp*g)/(h*h);
}

Real ShrinkingCoreLangmuir::computeQpResidual()
{
    computeLangmuirIsotherm ();
    return ShrinkingCoreConstant::computeQpResidual ();
}

Real ShrinkingCoreLangmuir::computeQpJacobian()
{
    computeLangmuirIsotherm ();
    return ShrinkingCoreConstant::computeQpJacobian ();
}

Real ShrinkingCoreLangmuir::computeQpOffDiagJacobian(unsigned int jvar)
{
    computeLangmuirIsotherm ();
    
    //Off-diagonal for concentration
    if (jvar == _coupled_var_conc)
    {
        return -_test[_i][_qp]*_phi[_j][_qp]*computeShrinkingCoreRateConcOffDiagJac();
    }
    
    return 0.0;
}


