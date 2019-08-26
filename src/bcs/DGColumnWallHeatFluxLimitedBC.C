/*!
 *  \file DGColumnWallHeatFluxLimitedBC.h
 *	\brief Boundary Condition kernel for a dirichlet-like boundary condition of heat on the column wall
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

#include "DGColumnWallHeatFluxLimitedBC.h"
registerMooseObject("dgospreyApp", DGColumnWallHeatFluxLimitedBC);

template<>
InputParameters validParams<DGColumnWallHeatFluxLimitedBC>()
{
	InputParameters params = validParams<DGFluxLimitedBC>();
	params.addCoupledVar("wall_temp", "Name of the wall temperature variable to couple with");
	
	return params;
}

DGColumnWallHeatFluxLimitedBC::DGColumnWallHeatFluxLimitedBC(const InputParameters & parameters) :
DGFluxLimitedBC(parameters),
_wall_temp(coupledValue("wall_temp")),
_bed_wall_transfer_coeff(getMaterialProperty<Real>("bed_wall_transfer_coeff")),
_conductivity(getMaterialProperty<Real>("conductivity"))
{
	
}

Real
DGColumnWallHeatFluxLimitedBC::computeQpResidual()
{
	_velocity(0)=_bed_wall_transfer_coeff[_qp];
	_velocity(1)=0.0;
	_velocity(2)=0.0;
	
	//Output
	if ((_velocity)*_normals[_qp] > 0.0)
	{
		_velocity(0)=0.0;
	}
	//Input
	else
	{
		_velocity(0)=_bed_wall_transfer_coeff[_qp];
	}
	
	_Diffusion(0,0) =  _conductivity[_qp];
	_Diffusion(0,1) = 0.0;
	_Diffusion(0,2) = 0.0;
	
	_Diffusion(1,0) = 0.0;
	_Diffusion(1,1) = _conductivity[_qp];
	_Diffusion(1,2) = 0.0;
	
	_Diffusion(2,0) = 0.0;
	_Diffusion(2,1) = 0.0;
	_Diffusion(2,2) = 0.0;
	
	_u_input = _wall_temp[_qp];
	
	return DGFluxLimitedBC::computeQpResidual();
}

Real
DGColumnWallHeatFluxLimitedBC::computeQpJacobian()
{
	_velocity(0)=_bed_wall_transfer_coeff[_qp];
	_velocity(1)=0.0;
	_velocity(2)=0.0;
	
	//Output
	if ((_velocity)*_normals[_qp] > 0.0)
	{
		_velocity(0)=0.0;
	}
	//Input
	else
	{
		_velocity(0)=_bed_wall_transfer_coeff[_qp];
	}
	
	_Diffusion(0,0) =  _conductivity[_qp];
	_Diffusion(0,1) = 0.0;
	_Diffusion(0,2) = 0.0;
	
	_Diffusion(1,0) = 0.0;
	_Diffusion(1,1) = _conductivity[_qp];
	_Diffusion(1,2) = 0.0;
	
	_Diffusion(2,0) = 0.0;
	_Diffusion(2,1) = 0.0;
	_Diffusion(2,2) = 0.0;
	
	_u_input = _wall_temp[_qp];
	
	return DGFluxLimitedBC::computeQpJacobian();
}
