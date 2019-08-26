/*!
 *  \file GColumnHeatDispersion.h
 *	\brief Kernel for use with the corresponding DGColumnHeatDispersion object
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

#include "GColumnHeatDispersion.h"
registerMooseObject("dgospreyApp", GColumnHeatDispersion);

template<>
InputParameters validParams<GColumnHeatDispersion>()
{
	InputParameters params = validParams<GAnisotropicDiffusion>();
	return params;
}

GColumnHeatDispersion::GColumnHeatDispersion(const InputParameters & parameters) :
GAnisotropicDiffusion(parameters),
_conductivity(getMaterialProperty<Real>("conductivity"))
{
	
}

Real
GColumnHeatDispersion::computeQpResidual()
{
	_Diffusion(0,0) =  _conductivity[_qp];
	_Diffusion(0,1) = 0.0;
	_Diffusion(0,2) = 0.0;
	
	_Diffusion(1,0) = 0.0;
	_Diffusion(1,1) = _conductivity[_qp];
	_Diffusion(1,2) = 0.0;
	
	_Diffusion(2,0) = 0.0;
	_Diffusion(2,1) = 0.0;
	_Diffusion(2,2) = 0.0;
	
	return GAnisotropicDiffusion::computeQpResidual();
}

Real
GColumnHeatDispersion::computeQpJacobian()
{
	_Diffusion(0,0) =  _conductivity[_qp];
	_Diffusion(0,1) = 0.0;
	_Diffusion(0,2) = 0.0;
	
	_Diffusion(1,0) = 0.0;
	_Diffusion(1,1) = _conductivity[_qp];
	_Diffusion(1,2) = 0.0;
	
	_Diffusion(2,0) = 0.0;
	_Diffusion(2,1) = 0.0;
	_Diffusion(2,2) = 0.0;
	
	return GAnisotropicDiffusion::computeQpJacobian();
}


