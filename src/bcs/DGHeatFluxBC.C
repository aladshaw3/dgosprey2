/*!
 *  \file DGHeatFluxBC.h
 *	\brief Boundary Condition kernel for the heat flux in and out of the ends of the fixed-bed column
 *  \author Austin Ladshaw
 *	\date 11/20/2015
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2015, all
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

#include "DGHeatFluxBC.h"
registerMooseObject("dgospreyApp", DGHeatFluxBC);

template<>
InputParameters validParams<DGHeatFluxBC>()
{
	InputParameters params = validParams<DGFluxBC>();
	params.addParam<Real>("input_temperature","Inlet temperature of column (K)");
	
	return params;
}

DGHeatFluxBC::DGHeatFluxBC(const InputParameters & parameters) :
DGFluxBC(parameters),
_input_temperature(getParam<Real>("input_temperature")),
_vel(getMaterialProperty<Real>("velocity")),
_gas_density(getMaterialProperty<Real>("gas_density")),
_gas_heat_capacity(getMaterialProperty<Real>("gas_heat_capacity"))
{
	
}

Real
DGHeatFluxBC::computeQpResidual()
{
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp]*_gas_heat_capacity[_qp]*_gas_density[_qp];
	_velocity(2)=0.0;
	
	_u_input = _input_temperature;
	
	return DGFluxBC::computeQpResidual();
}

Real
DGHeatFluxBC::computeQpJacobian()
{
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp]*_gas_heat_capacity[_qp]*_gas_density[_qp];
	_velocity(2)=0.0;
	
	_u_input = _input_temperature;
	
	return DGFluxBC::computeQpJacobian();
}

