/*!
 *  \file HeatofAdsorption.h
 *	\brief Kernel for calculation of heat of adsorption based on Thermodynamic Properties and Solid Concentration
 *	\details This file creates a kernel for calculating the heat of adsorption as a non-linear forcing function of
 *				the thermodynamic properties and the solid concentration. It inherits from the CoupledLinearForcingFunction
 *				kernel and overrides its computeQpResidual and computeQpJacobian functions. Input parameters include
 *				the species index and the coupled solid concentration.
 *
 *				Mathematical Description: heat_i = solid_i * isostericheat_i     where isostericheat_i is a parameter
 *				calculated from the thermodynamic properties.
 *
 *  \author Austin Ladshaw
 *	\date 04/28/2017
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

#include "HeatofAdsorption.h"
registerMooseObject("dgospreyApp", HeatofAdsorption);

template<>
InputParameters validParams<HeatofAdsorption>()
{
	InputParameters params = validParams<CoupledLinearForcingFunction>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

HeatofAdsorption::HeatofAdsorption(const InputParameters & parameters)
: CoupledLinearForcingFunction(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data")),
_ads_heat(getMaterialProperty<std::vector<Real> >("adsorption_heat"))
{
	if (_gaining == false)
		_coef = -_coef;
}

Real HeatofAdsorption::computeQpResidual()
{
	_coef = _ads_heat[_qp][_index];
	return CoupledLinearForcingFunction::computeQpResidual();
}

Real HeatofAdsorption::computeQpJacobian()
{
	return CoupledLinearForcingFunction::computeQpJacobian();
}

Real HeatofAdsorption::computeQpOffDiagJacobian(unsigned int jvar)
{
	_coef = _ads_heat[_qp][_index];
	//return CoupledLinearForcingFunction::computeQpOffDiagJacobian(jvar);
	return 0.0;
}

