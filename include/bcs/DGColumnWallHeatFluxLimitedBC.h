/*!
 *  \file DGColumnWallHeatFluxLimitedBC.h
 *	\brief Boundary Condition kernel for a dirichlet-like boundary condition of heat on the column wall
 *	\details This file creates a boundary condition that mimics a Dirichlet boundary condition for the
 *			heat of the column at the wall. A true Dirichlet boundary condition does not exist in DG
 *			methods. However, this will create a weak form that imposes a constraint that the solution
 *			must be of a certain value at the boundary.
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

/// DGColumnWallHeatFluxLimitedBC class object forward declaration
class DGColumnWallHeatFluxLimitedBC;

template<>
InputParameters validParams<DGColumnWallHeatFluxLimitedBC>();

/// DGColumnWallHeatFluxLimitedBC class object inherits from DGFluxLimitedBC object
/** This class object inherits from the DGFluxLimitedBC object (see DGFluxLimitedBC.h for more details).
	All public and protected members of this class are required function overrides. The object
	will take in the given variable and material properties to override some objects declared
	for the generic DGFluxBC to fit this particular boundary condition. Then, it just calls the
	appropriate DGFluxLimitedBC functions. */
class DGColumnWallHeatFluxLimitedBC : public DGFluxLimitedBC
{
public:
	/// Required constructor for BC objects in MOOSE
	DGColumnWallHeatFluxLimitedBC(const InputParameters & parameters);
	
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
	
  	const VariableValue & _wall_temp;								///< Reference to the coupled variable for wall temperature of the column
	const MaterialProperty<Real> & _bed_wall_transfer_coeff;///< Reference to the bed-wall transfer coefficient material property
	const MaterialProperty<Real> & _conductivity;			///< Reference to the thermal conductivity material property
	
};
