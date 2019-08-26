/*!
 *  \file TotalColumnPressure.h
 *	\brief Auxillary kernel to calculate total column pressure based on temperature and concentrations
 *	\details This file is responsible for calculating the total column pressure in a fixed-bed adsorber
 *			given the temperature and the concentrations of all species in the gas phase. The gas phase
 *			is assumed to behave ideally and ideal gas law is employed to estimate the pressure.
 *
 *  \author Austin Ladshaw
 *	\date 11/20/2015
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


/// Total Column Pressure class object forward declaration
class TotalColumnPressure;

template<>
InputParameters validParams<TotalColumnPressure>();

/// Total Column Pressure class inherits from AuxKernel
/** This class object creates an AuxKernel for use in the MOOSE framework. The AuxKernel will
	calculate the total column pressure (in kPa) based on the non-linear variables of temperature
	and concentration of each species in the gas phase. Total pressure is calculated based on the
	ideal gas law. */
class TotalColumnPressure : public AuxKernel
{
public:
	/// Standard MOOSE public constructor
	TotalColumnPressure(const InputParameters & parameters);
  
protected:
	/// Required MOOSE function override
	/** This is the function that is called by the MOOSE framework when a calculation of the total
		system pressure is needed. You are required to override this function for any inherited 
		AuxKernel. */
	virtual Real computeValue();
	
private:
	const VariableValue & _temperature;					///< Reference to the temperature non-linear variable
	std::vector<unsigned int> _index;				///< Indices of the gaseous species coupled to the object
	std::vector<const VariableValue *> _gas_conc;			///< Pointer list for the non-linear concentration variables
  
};
