/*!
 *  \file ConcentrationIC.h
 *	\brief Initial Condition kernel for initial concentration of a species in a fixed-bed column
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

#include "ConcentrationIC.h"
registerMooseObject("dgospreyApp", ConcentrationIC);

template<> InputParameters validParams<ConcentrationIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("initial_mole_frac","The initial mole fraction of the gas species in the column. Note: All mole fractions must sum to one, otherwise errors will occur that will not be caught!!!");
  params.addRequiredParam<Real>("initial_press","Initial condition for total pressure");
  params.addRequiredParam<Real>("initial_temp","Initial condition for temperature");
  return params;
}

ConcentrationIC::ConcentrationIC(const InputParameters & parameters) :
InitialCondition(parameters),
_y_IC(getParam<Real>("initial_mole_frac")),
_PT_IC(getParam<Real>("initial_press")),
_T_IC(getParam<Real>("initial_temp"))
{
  
}

Real ConcentrationIC::value(const Point & p)
{
  //Note: we will not check to ensure the mole fraction is valid. This may cause error.
  return (_PT_IC * _y_IC) / (8.3144621 * _T_IC);
}
