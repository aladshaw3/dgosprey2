/*!
 *  \file BedProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with the fixed-bed
 *	\warning THIS KERNEL IS INCOMPLETE! ONLY USED FOR DATA STORAGE FOR VARIOUS BED PARAMETERS!
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

#include "BedProperties.h"
registerMooseObject("dgospreyApp", BedProperties);

template<>
// input parameters are the parameters that are constant and not calculated from other parameters
InputParameters validParams<BedProperties>()
{
	InputParameters params = validParams<Material>();
  
	params.addParam<Real>("length", "Bed length (cm)");
	params.addParam<Real>("inner_diameter", "Wall inner diameter (cm)");
	params.addParam<Real>("outer_diameter","Wall outer diameter (cm)");
	params.addParam<Real>("bulk_porosity","Bed bulk porosity");
	params.addParam<Real>("wall_density","Wall material density (g/cm^3)");
	params.addParam<Real>("wall_heat_capacity","Wall heat capacity (J/g/K)");
	params.addParam<Real>("wall_heat_trans_coef","Bed-Wall Heat Transfer Coefficient (J/hr/cm^2/K)");
	params.addParam<Real>("extern_heat_trans_coef","Wall-External Heat Transfer Coefficient (J/hr/cm^2/K)");
  
	return params;
}

BedProperties::BedProperties(const InputParameters & parameters)
  :Material(parameters),
   
   _length(getParam<Real>("length")),
   _din(getParam<Real>("inner_diameter")),
   _dout(getParam<Real>("outer_diameter")),
   _eb(getParam<Real>("bulk_porosity")),
   _rhow(getParam<Real>("wall_density")),
   _hw(getParam<Real>("wall_heat_capacity")),
   _Uw(getParam<Real>("wall_heat_trans_coef")),
   _Ua(getParam<Real>("extern_heat_trans_coef")),
   _inner_dia(declareProperty<Real>("inner_dia")),
   _outer_dia(declareProperty<Real>("outer_dia")),
   _bed_length(declareProperty<Real>("bed_length")),
   _porosity(declareProperty<Real>("porosity")),
   _wall_density(declareProperty<Real>("wall_density")),
   _wall_heat_capacity(declareProperty<Real>("wall_heat_capacity")),
   _bed_wall_transfer_coeff(declareProperty<Real>("bed_wall_transfer_coeff")),
   _wall_exterior_transfer_coeff(declareProperty<Real>("wall_exterior_transfer_coeff"))
   
{

}

void
BedProperties::computeQpProperties()
{
	//For constant bed properties...
	_inner_dia[_qp] = _din;
	_outer_dia[_qp] = _dout;
	_bed_length[_qp] = _length;
	_porosity[_qp] = _eb;
	_wall_density[_qp] = _rhow;
	_wall_heat_capacity[_qp] = _hw;
	_bed_wall_transfer_coeff[_qp] = _Uw;
	_wall_exterior_transfer_coeff[_qp] = _Ua;
	//Note: some of these may vary with space or temperature or concentration, but for now we assume constant 
}
