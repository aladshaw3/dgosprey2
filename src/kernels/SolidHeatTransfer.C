/*!
 *  \file SolidHeatTransfer.h
 *	\brief Standard kernel for the transfer of heat via adsorption
 *	\details This file creates a standard MOOSE kernel for the transfer of heat between the
 *			bulk gas of the fixed-bed and the adsorbent material in the column. The
 *			heat transfer is based on the amount of material in the bed and the heat of adsorption variables.
 *
 *  \author Austin Ladshaw
 *	\date 04/28/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2016, all
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

#include "SolidHeatTransfer.h"
registerMooseObject("dgospreyApp", SolidHeatTransfer);

template<>
InputParameters validParams<SolidHeatTransfer>()
{
	InputParameters params = validParams<CoupledCoeffTimeDerivative>();
	return params;
}

SolidHeatTransfer::SolidHeatTransfer(const InputParameters & parameters)
: CoupledCoeffTimeDerivative(parameters),
_porosity(getMaterialProperty<Real>("porosity")),
_pellet_density(getMaterialProperty<Real>("pellet_density"))
{
	
}

Real SolidHeatTransfer::computeQpResidual()
{
	_time_coef = -(1.0-_porosity[_qp])*_pellet_density[_qp];
	return CoupledCoeffTimeDerivative::computeQpResidual();
}

Real SolidHeatTransfer::computeQpJacobian()
{
	return CoupledCoeffTimeDerivative::computeQpJacobian();
}

Real SolidHeatTransfer::computeQpOffDiagJacobian(unsigned int jvar)
{
	_time_coef = -(1.0-_porosity[_qp])*_pellet_density[_qp];
	return CoupledCoeffTimeDerivative::computeQpOffDiagJacobian(jvar);
}

