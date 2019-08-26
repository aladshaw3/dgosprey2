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

#pragma once

#include "AuxKernel.h"


/// Wall Temperature class object forward declaration
class WallTemperature;

template<>
InputParameters validParams<WallTemperature>();

/// Wall Temperature class inherits from AuxKernel
/** This class object creates an AuxKernel for use in the MOOSE framework. The AuxKernel will
	calculate the wall temperature (in K) based on the non-linear variables of column temperature
	and ambient temperature, as well as other parameters from material properties. */
class WallTemperature : public AuxKernel
{
public:
	/// Standard MOOSE public constructor
	WallTemperature(const InputParameters & parameters);
	
protected:
	/// Required MOOSE function override
	/** This is the function that is called by the MOOSE framework when a calculation of the total
		system pressure is needed. You are required to override this function for any inherited
		AuxKernel. */
	virtual Real computeValue();
	
private:
	const VariableValue & _column_temp;					///< Reference to the temperature non-linear variable
	const VariableValue & _ambient_temp;				///< Reference to the temperature non-linear variable
	const MaterialProperty<Real> & _hbw;				///< Reference to the bed-wall transfer coefficient material property
	const MaterialProperty<Real> & _haw;				///< Reference to the ambient-wall transfer coefficient material property
	const MaterialProperty<Real> & _din;				///< Reference to the inner diameter of the column
	const MaterialProperty<Real> & _dout;				///< Reference to the outer diameter of the column
	const MaterialProperty<Real> & _rhow;				///< Reference to the wall density material property
	const MaterialProperty<Real> & _cw;					///< Reference to the heat capacity coefficient material property
	
};
