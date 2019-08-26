/*!
 *  \file WallAmbientHeatTransfer.h
 *	\brief Standard kernel for the transfer of heat from the column wall to the ambient air
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

#include "WallAmbientHeatTransfer.h"
registerMooseObject("dgospreyApp", WallAmbientHeatTransfer);

template<>
InputParameters validParams<WallAmbientHeatTransfer>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("coupled", "Name of the ambient temperature variable to be coupled.");
  return params;
}


WallAmbientHeatTransfer::WallAmbientHeatTransfer(const InputParameters & parameters)
:Kernel(parameters),
_wall_exterior_transfer_coeff(getMaterialProperty<Real>("wall_exterior_transfer_coeff")),
_inner_dia(getMaterialProperty<Real>("inner_dia")),
_outer_dia(getMaterialProperty<Real>("outer_dia")),
_ambient_temp(coupledValue("coupled"))
{
}

Real WallAmbientHeatTransfer::computeQpResidual()
{
  //Note: if _inner_dia >= _outer_dia, this will lead to an error
  Real _coef = (4.0 * _wall_exterior_transfer_coeff[_qp] * _outer_dia[_qp]) / ( (_outer_dia[_qp]*_outer_dia[_qp]) - (_inner_dia[_qp]*_inner_dia[_qp]) );
  return _test[_i][_qp] * _coef * (_u[_qp] - _ambient_temp[_qp]);
}

Real WallAmbientHeatTransfer::computeQpJacobian()
{
  Real _coef = (4.0 * _wall_exterior_transfer_coeff[_qp] * _outer_dia[_qp]) / ( (_outer_dia[_qp]*_outer_dia[_qp]) - (_inner_dia[_qp]*_inner_dia[_qp]) );
  
  return _test[_i][_qp] * _coef * _phi[_j][_qp];
}

