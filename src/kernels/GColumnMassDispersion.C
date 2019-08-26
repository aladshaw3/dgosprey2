/*!
 *  \file GColumnMassDispersion.h
 *	\brief Kernel for use with the corresponding DGColumnMassDispersion object
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

#include "GColumnMassDispersion.h"
registerMooseObject("dgospreyApp", GColumnMassDispersion);

template<>
InputParameters validParams<GColumnMassDispersion>()
{
	InputParameters params = validParams<GAnisotropicDiffusion>();
	params.addParam<unsigned int>("index", 0, "The index of the coupling variable. Must be given in same order of appearance as in the FlowProperties Material block. Indexing starts from 0. 0 is default value.");
	return params;
}

GColumnMassDispersion::GColumnMassDispersion(const InputParameters & parameters) :
GAnisotropicDiffusion(parameters),
_index(getParam<unsigned int>("index")),
_dispersion(getMaterialProperty<std::vector<Real> >("dispersion")),
_molecular_diffusion(getMaterialProperty<std::vector<Real> >("molecular_diffusion"))
{

}

Real
GColumnMassDispersion::computeQpResidual()
{
	_Diffusion(0,0) =  _molecular_diffusion[_qp][_index];
	_Diffusion(0,1) = 0.0;
	_Diffusion(0,2) = 0.0;
	
	_Diffusion(1,0) = 0.0;
	_Diffusion(1,1) = _dispersion[_qp][_index];
	_Diffusion(1,2) = 0.0;
	
	_Diffusion(2,0) = 0.0;
	_Diffusion(2,1) = 0.0;
	_Diffusion(2,2) = 0.0;
	
	return GAnisotropicDiffusion::computeQpResidual();
}

Real
GColumnMassDispersion::computeQpJacobian()
{
	_Diffusion(0,0) =  _molecular_diffusion[_qp][_index];
	_Diffusion(0,1) = 0.0;
	_Diffusion(0,2) = 0.0;
	
	_Diffusion(1,0) = 0.0;
	_Diffusion(1,1) = _dispersion[_qp][_index];
	_Diffusion(1,2) = 0.0;
	
	_Diffusion(2,0) = 0.0;
	_Diffusion(2,1) = 0.0;
	_Diffusion(2,2) = 0.0;
	
	return GAnisotropicDiffusion::computeQpJacobian();
}

