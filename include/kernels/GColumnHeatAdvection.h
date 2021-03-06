/*!
 *  \file GColumnHeatAdvection.h
 *	\brief Kernel for use with the corresponding DGColumnHeatAdvection object
 *	\details This file creates a standard MOOSE kernel that is to be used in conjunction with DGColumnHeatAdvection
 *			for the discontinous Galerkin formulation of the heat advection physics for a fixed-bed adsorber. It couples
 *			with material properties to override the velocity parameter in the inherited GAdection kernel, then simply
 *			calls the corresponding methods of the base class.
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

#include "GAdvection.h"

/// GColumnHeatAdvection class object forward declarations
class GColumnHeatAdvection;

template<>
InputParameters validParams<GColumnHeatAdvection>();

/// GColumnHeatAdvection class object inherits from GAdvection object
/** This class object inherits from the generic GAdvection kernel for use with 
	the corresponding DGColumnHeatAdvection kernel to complete the physical 
	description of DG methods in MOOSE. It is coupled with the material properties
	of linear velocity, gas density, and gas heat capacity and uses those parameters
	to override the components of the velocity vector of the more generic GAdvection 
	class. */
class GColumnHeatAdvection : public GAdvection
{
public:
	/// Required constructor for objects in MOOSE
	GColumnHeatAdvection(const InputParameters & parameters);
	
protected:
	/// Required residual function for standard kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	/// Required Jacobian function for standard kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
	
private:
	const MaterialProperty<Real> & _vel;				///< Reference to the linear velocity material property
	const MaterialProperty<Real> & _gas_density;		///< Reference to the gas density material property
	const MaterialProperty<Real> & _gas_heat_capacity;	///< Reference to the gas heat capacity material property

};
