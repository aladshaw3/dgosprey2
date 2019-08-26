/*!
 *  \file DGHeatFluxLimitedBC.h
 *	\brief Boundary Condition kernel to mimic a dirichlet boundary condition at the column inlet
 *	\details This file creates a dirichlet-like boundary condition kernel for the column temperature
 *			at the inlet of the system. The outlet boundary condition would remain unchanged from the
 *			standard DG form of the boundaries. Only the inlet BC is affected by this file. See
 *			FluxLimitedBC.h for more details.
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

#pragma once

#include "DGFluxLimitedBC.h"

/// DGHeatFluxLimitedBC class object forward declaration
class DGHeatFluxLimitedBC;

template<>
InputParameters validParams<DGHeatFluxLimitedBC>();

/// DGHeatFluxLimitedBC class object inherits from DGFluxLimitedBC object
/** This class object inherits from the DGFluxLimitedBC object (see DGFluxLimitedBC.h for more details).
	All public and protected members of this class are required function overrides. The object
	will take in the given variable and material properties to override some objects declared
	for the generic DGFluxLimitedBC to fit this particular boundary condition. Then, it just calls the
	appropriate DGFluxLimitedBC functions. */
class DGHeatFluxLimitedBC : public DGFluxLimitedBC
{
public:
	/// Required constructor for BC objects in MOOSE
	DGHeatFluxLimitedBC(const InputParameters & parameters);
	
protected:
	/// Required function override for BC objects in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	/// Required function override for BC objects in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
	
private:
	/// Value of the column temperature at the inlet of the system
	Real _input_temperature;
	
	const MaterialProperty<Real> & _vel;				///< Reference to the velocity material property
	const MaterialProperty<Real> & _gas_density;		///< Reference to the gas density material property
	const MaterialProperty<Real> & _gas_heat_capacity;	///< Reference to the gas heat capacity material property
	const MaterialProperty<Real> & _conductivity;		///< Reference to the thermal conductivity material property
	
};
