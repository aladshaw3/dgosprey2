/*!
 *  \file CoupledLangmuirForcingFunction.h
 *	\brief Standard kernel for coupling non-linear variables via a langmuir function
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via a langmuir forcing function, i.e., variable = b * K * coupled_variable / 1 + K * coupled_variable.
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 06/07/2017
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

#include "CoupledLangmuirForcingFunction.h"
registerMooseObject("dgospreyApp", CoupledLangmuirForcingFunction);

template<>
InputParameters validParams<CoupledLangmuirForcingFunction>()
{
    InputParameters params = validParams<TimeDerivative>();
    params.addParam<Real>("max_capacity",0.0,"Maximum Capacity for Langmuir Function (mol/kg)");
    params.addParam<Real>("langmuir_coeff",1.0,"Coefficient for the langmuir function (L/mol)");
    params.addRequiredCoupledVar("coupled","Name of the variable being coupled");
    return params;
}

CoupledLangmuirForcingFunction::CoupledLangmuirForcingFunction(const InputParameters & parameters)
: TimeDerivative(parameters),
_maxcap(getParam<Real>("max_capacity")),
_langmuircoef(getParam<Real>("langmuir_coeff")),
_coupled_u(coupledValue("coupled")),
_coupled_var(coupled("coupled"))
{
    
}

Real CoupledLangmuirForcingFunction::computeLangmuirEquilibrium()
{
	return _maxcap*((_langmuircoef*_coupled_u[_qp])/(1.0+_langmuircoef*_coupled_u[_qp]));
}

Real CoupledLangmuirForcingFunction::computeLangConcDerivative()
{
	return _maxcap*((_langmuircoef*_phi[_j][_qp])/((1.0+_langmuircoef*_coupled_u[_qp])*(1.0+_langmuircoef*_coupled_u[_qp])));
}

Real CoupledLangmuirForcingFunction::computeQpResidual()
{
    return _u[_qp]*_test[_i][_qp]-_test[_i][_qp]*computeLangmuirEquilibrium();
}

Real CoupledLangmuirForcingFunction::computeQpJacobian()
{
    return _phi[_j][_qp]*_test[_i][_qp];
}

Real CoupledLangmuirForcingFunction::computeQpOffDiagJacobian(unsigned int jvar)
{
    if (jvar == _coupled_var)
        return -_test[_i][_qp]*computeLangConcDerivative();
    
    return 0.0;
}


