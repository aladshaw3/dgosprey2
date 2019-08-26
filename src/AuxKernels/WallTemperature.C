/*!
 *  \file WallTemperature.h
 *	\brief Auxillary kernel to calculate wall temperature from ambient and column temperatures
 *	\details This file is responsible for calculating the wall temperature in a fixed-bed adsorber
 *			given the column temperature and ambient temperature or other temperature. It is linked
 *			with material properties to account for parameters of wall density, wall heat capacity,
 *			heat transfer coefficients, and column dimensions.
 *
 *			The calculation is carried out with a semi-implicit, first order linear approximation at
 *			the start and end of each time step.
 *
 *  \author Austin Ladshaw
 *	\date 06/05/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2015, all
 *             rights reserved.
 *
 *			   Austin Ladshaw does not claim any owernship or copyright to the
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

#include "WallTemperature.h"
registerMooseObject("dgospreyApp", WallTemperature);

template<>
InputParameters validParams<WallTemperature>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addCoupledVar("column_temp","Coupled variable for column temperature");
	params.addCoupledVar("ambient_temp", "Coupled variable for ambient temperature");
	return params;
}

WallTemperature::WallTemperature(const InputParameters & parameters) :
AuxKernel(parameters),
_column_temp(coupledValue("column_temp")),
_ambient_temp(coupledValue("ambient_temp")),
_hbw(getMaterialProperty< Real >("bed_wall_transfer_coeff")),
_haw(getMaterialProperty< Real >("wall_exterior_transfer_coeff")),
_din(getMaterialProperty< Real >("inner_dia")),
_dout(getMaterialProperty< Real >("outer_dia")),
_rhow(getMaterialProperty< Real >("wall_density")),
_cw(getMaterialProperty< Real >("wall_heat_capacity"))
{
	//Forces specific execution behavior of the auxkernel
	//_exec_flags.clear();
	//_exec_flags.push_back(EXEC_INITIAL);
	//_exec_flags.push_back(EXEC_TIMESTEP_END);
}

Real WallTemperature::computeValue()
{
	Real value = 0.0;
	
	double Kbw, Kaw;
	
	Kbw = 4.0*_hbw[_qp]*_din[_qp]/((_dout[_qp]*_dout[_qp])-(_din[_qp]*_din[_qp]));
	Kaw = 4.0*_haw[_qp]*_dout[_qp]/((_dout[_qp]*_dout[_qp])-(_din[_qp]*_din[_qp]));
	
	value = ((_rhow[_qp]*_cw[_qp]*_u_old[_qp])+(_dt*Kbw*_column_temp[_qp])+(_dt*Kaw*_ambient_temp[_qp]))/((_rhow[_qp]*_cw[_qp])+(_dt*Kbw)+(_dt*Kaw));
	
	return value;
}
