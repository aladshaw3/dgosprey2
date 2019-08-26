/*!
 *  \file DGAdvection.h
 *	\brief Discontinous Galerkin kernel for advection
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

#include "DGAdvection.h"
registerMooseObject("dgospreyApp", DGAdvection);

template<>
InputParameters validParams<DGAdvection>()
{
	InputParameters params = validParams<DGKernel>();
	params.addParam<Real>("vx",0,"x-component of velocity vector");
	params.addParam<Real>("vy",0,"y-component of velocity vector");
	params.addParam<Real>("vz",0,"z-component of velocity vector");
	return params;
}

DGAdvection::DGAdvection(const InputParameters & parameters) :
DGKernel(parameters),
_vx(getParam<Real>("vx")),
_vy(getParam<Real>("vy")),
_vz(getParam<Real>("vz"))
{
	_velocity(0)=_vx;
	_velocity(1)=_vy;
	_velocity(2)=_vz;
}

Real DGAdvection::computeQpResidual(Moose::DGResidualType type)
{
	Real r = 0;
	
	switch (type)
	{
		case Moose::Element:
			if ( (_velocity * _normals[_qp]) >= 0.0)
				r += (_velocity * _normals[_qp]) * _u[_qp] * _test[_i][_qp];
			else
				r += (_velocity * _normals[_qp]) * _u_neighbor[_qp] * _test[_i][_qp];
			break;
			
		case Moose::Neighbor:
			if ( (_velocity * _normals[_qp]) >= 0.0)
				r += -(_velocity * _normals[_qp]) * _u[_qp] * _test_neighbor[_i][_qp];
			else
				r += -(_velocity * _normals[_qp]) * _u_neighbor[_qp] * _test_neighbor[_i][_qp];
			break;
	}
	return r;
}

Real DGAdvection::computeQpJacobian(Moose::DGJacobianType type)
{
	Real r = 0;
	
	switch (type)
	{
			
		case Moose::ElementElement:
			if ( (_velocity * _normals[_qp]) >= 0.0)
				r += (_velocity * _normals[_qp]) * _phi[_j][_qp] * _test[_i][_qp];
			else
				r += 0.0;
			break;
			
		case Moose::ElementNeighbor:
			if ( (_velocity * _normals[_qp]) >= 0.0)
				r += 0.0;
			else
				r += (_velocity * _normals[_qp]) * _phi_neighbor[_j][_qp] * _test[_i][_qp];
			break;
			
		case Moose::NeighborElement:
			if ( (_velocity * _normals[_qp]) >= 0.0)
				r += -(_velocity * _normals[_qp]) * _phi[_j][_qp] * _test_neighbor[_i][_qp];
			else
				r += 0.0;
			break;
			
		case Moose::NeighborNeighbor:
			if ( (_velocity * _normals[_qp]) >= 0.0)
				r += 0.0;
			else
				r += -(_velocity * _normals[_qp]) * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
			break;
	}
	return r;
}
