/*!
 *  \file DGOSPREY_TimeStepper.h
 *	\brief Standard executioner for choosing the time step for DGOSPREY simulations
 *	\details This file creates a standard MOOSE executioner for determining time step sizes based
 *			on flow rates, particle sizes, dispersion coefficients, etc. It inherits from the 
 *			SolutionTimeAdaptiveDT MOOSE object and alters the time step choosen based on system
 *			parameters for a given simulation case.
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 06/22/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
 *             rights reserved.
 *
 *			   Alexander Wiechert does not claim any ownership or copyright to the
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

#include "SolutionTimeAdaptiveDT.h"

class DGOSPREY_TimeStepper;

template <>
InputParameters validParams<DGOSPREY_TimeStepper>();

/**
 *
 */
class DGOSPREY_TimeStepper : public SolutionTimeAdaptiveDT
{
public:
	DGOSPREY_TimeStepper(const InputParameters & parameters);
	virtual ~DGOSPREY_TimeStepper();
	
	virtual void step() override;
	
	virtual void rejectStep() override;
	
protected:
	virtual Real computeInitialDT() override;
	virtual Real computeDT() override;
	
private:
	Real _flow_rate;						///< Inlet flow rate for the fixed-bed column (cm^3/hr)
	Real _pellet_dia;						///< Nominal diameter of the adsorbent pellets in the system (cm)
	Real _din;								///< Column inner diameter (cm)
	
};
