/*!
 *  \file DGMassFluxBC.h
 *	\brief Boundary Condition kernel for the mass flux in and out of the ends of the fixed-bed column
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

#include "DGMassFluxBC.h"
registerMooseObject("dgospreyApp", DGMassFluxBC);

template<>
InputParameters validParams<DGMassFluxBC>()
{
	InputParameters params = validParams<DGFluxBC>();
	
	params.addParam<Real>("input_temperature","Inlet temperature of column (K)");
	params.addParam<Real>("input_pressure","Inlet pressure of column (kPa)");
	params.addParam<Real>("input_molefraction","Inlet mole fraction of the gas");
	params.addParam<unsigned int>("index", 0, "The index of the coupling variable. Must be given in same order of appearance as in the FlowProperties Material block. Indexing starts from 0. 0 is default value.");
	
	return params;
}

DGMassFluxBC::DGMassFluxBC(const InputParameters & parameters) :
DGFluxBC(parameters),
_input_temperature(getParam<Real>("input_temperature")),
_input_pressure(getParam<Real>("input_pressure")),
_input_molefraction(getParam<Real>("input_molefraction")),
_vel(getMaterialProperty<Real>("velocity")),
_index(getParam<unsigned int>("index"))
{

}

Real
DGMassFluxBC::computeQpResidual()
{
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp];
	_velocity(2)=0.0;
	
	_u_input = (_input_pressure * _input_molefraction) / (8.3144621 * _input_temperature);
		
	return DGFluxBC::computeQpResidual();
}

Real
DGMassFluxBC::computeQpJacobian()
{
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp];
	_velocity(2)=0.0;
	
	_u_input = (_input_pressure * _input_molefraction) / (8.3144621 * _input_temperature);
	
	return DGFluxBC::computeQpJacobian();
}

