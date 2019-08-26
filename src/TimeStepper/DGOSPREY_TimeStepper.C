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

#include "DGOSPREY_TimeStepper.h"
#include "FEProblem.h"
#include "Transient.h"

#include <chrono>

registerMooseObject("dgospreyApp", DGOSPREY_TimeStepper);

template <>
InputParameters
validParams<DGOSPREY_TimeStepper>()
{
	InputParameters params = validParams<SolutionTimeAdaptiveDT>();
	params.addParam<Real>("flow_rate","Volumetric flow rate in the system (cm^3/hr)");
	params.addParam<Real>("pellet_diameter","Diameter of the adsorbent pellet (cm)");
	params.addParam<Real>("inner_diameter", "Wall inner diameter (cm)");
	
	return params;
}

DGOSPREY_TimeStepper::DGOSPREY_TimeStepper(const InputParameters & parameters)
: SolutionTimeAdaptiveDT(parameters),
_flow_rate(getParam<Real>("flow_rate")),
_pellet_dia(getParam<Real>("pellet_diameter")),
_din(getParam<Real>("inner_diameter"))
{

}

DGOSPREY_TimeStepper::~DGOSPREY_TimeStepper()
{
	
}

void
DGOSPREY_TimeStepper::step()
{
	SolutionTimeAdaptiveDT::step();
}

Real
DGOSPREY_TimeStepper::computeInitialDT()
{
	double timestep = getParam<Real>("dt");
	double alpha = 10000.0;
	double vel = (_flow_rate / ((M_PI/4.0) * _din * _din));
	if (timestep < (alpha/(360.0+(100.0*_pellet_dia*vel))))
		timestep = (alpha/(360.0+(100.0*_pellet_dia*vel)));
	return timestep;
}

Real
DGOSPREY_TimeStepper::computeDT()
{
	return SolutionTimeAdaptiveDT::computeDT();
}

void
DGOSPREY_TimeStepper::rejectStep()
{
	SolutionTimeAdaptiveDT::rejectStep();
}

