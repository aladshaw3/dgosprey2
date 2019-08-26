/*!
 *  \file BedWallHeatTransfer.h
 *	\brief Standard kernel for the transfer of heat from the fixed-bed to the column wall
 *	\details This file creates a standard MOOSE kernel for the transfer of energy as heat between the 
 *			bulk gas temperature of the fixed-bed and the temperature of the walls of the column. The 
 *			heat transfer is based on the thickness of the wall and a bed-wall heat transfer coefficient. 
 *			It is coupled to the heat of the gas in the column and is a primary kernel used in determining
 *			the heat of the wall.
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

#include "Kernel.h"

/// BedWallHeatTransfer class object forward declarations
class BedWallHeatTransfer;

template<>
InputParameters validParams<BedWallHeatTransfer>();

/// BedWallHeatTransfer class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the material properties for the size of the column, as well
	as the heat transfer coefficient for the exchange of energy from the gas to the wall,
	in order to form a residuals and Jacobians for the wall temperature variable. */
class BedWallHeatTransfer : public Kernel
{
public:
	/// Required constructor for objects in MOOSE
	BedWallHeatTransfer(const InputParameters & parameters);
  
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
	const MaterialProperty<Real> & _bed_wall_transfer_coeff;	///< Reference to the bed-wall heat transfer material property
	const MaterialProperty<Real> & _inner_dia;					///< Reference to the wall inner diameter material property
	const MaterialProperty<Real> & _outer_dia;					///< Reference to the wall outer diameter material property
  
	const VariableValue & _column_temp;								///< Reference to the gas temperature coupled non-linear variable 
  
};
