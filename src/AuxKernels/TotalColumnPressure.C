/*!
 *  \file TotalColumnPressure.h
 *	\brief Auxillary kernel to calculate total column pressure based on temperature and concentrations
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

#include "TotalColumnPressure.h"
registerMooseObject("dgospreyApp", TotalColumnPressure);

template<>
InputParameters validParams<TotalColumnPressure>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("temperature","Coupled variable for temperature");
  params.addCoupledVar("coupled_gases", "Gas concentrations variables being coupled");
  return params;
}

TotalColumnPressure::TotalColumnPressure(const InputParameters & parameters) :
AuxKernel(parameters),
_temperature(coupledValue("temperature"))
{
	//Forces specific execution behavior of the auxkernel
	//_exec_flags.clear();
	//_exec_flags.push_back(EXEC_INITIAL);
	//_exec_flags.push_back(EXEC_TIMESTEP_END);
	
	unsigned int n = coupledComponents("coupled_gases");
	_index.resize(n);
	_gas_conc.resize(n);
	
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		_index[i] = coupled("coupled_gases",i); 
		_gas_conc[i] = &coupledValue("coupled_gases",i);
	}

}

Real
TotalColumnPressure::computeValue()
{
	Real _PT = 0.0;
  
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		_PT = _PT + ( (*_gas_conc[i])[_qp] * 8.3144621 * _temperature[_qp] );
	}
	return _PT;
}
