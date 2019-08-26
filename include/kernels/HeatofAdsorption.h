/*!
 *  \file HeatofAdsorption.h
 *	\brief Kernel for calculation of heat of adsorption based on Thermodynamic Properties and Solid Concentration
 *	\details This file creates a kernel for calculating the heat of adsorption as a non-linear forcing function of
 *				the thermodynamic properties and the solid concentration. It inherits from the CoupledLinearForcingFunction
 *				kernel and overrides its computeQpResidual and computeQpJacobian functions. Input parameters include
 *				the species index and the coupled solid concentration.
 *
 *				Mathematical Description: heat_i = solid_i * isostericheat_i     where isostericheat_i is a parameter
 *				calculated from the thermodynamic properties.
 *
 *  \author Austin Ladshaw
 *	\date 04/28/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
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

#include "CoupledLinearForcingFunction.h"
#include "flock.h"
#include "DataStruct_StoreLoad.h"

/// HeatofAdsorption class object forward declarationss
class HeatofAdsorption;

template<>
InputParameters validParams<HeatofAdsorption>();

/// HeatofAdsorption class object inherits from CoupledLinearForcingFunction object
/** This class object inherits from the CoupledLinearForcingFunction object in DGOSPREY.
	All public and protected members of this class are required function overrides.
	The kernel interfaces between the solid concentration non-linear variable and the
	adsorption heat non-linear variable and computes adsorption heat as a linear function
	of the solid concentration. */
class HeatofAdsorption : public CoupledLinearForcingFunction
{
public:
	/// Required constructor for objects in MOOSE
	HeatofAdsorption(const InputParameters & parameters);
	
protected:
	/// Required residual function for standard kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	
	/// Required Jacobian function for standard kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
	
	/// Not Required, but aids in the preconditioning step
	/** This function returns the off diagonal Jacobian contribution for this object. By
		returning a non-zero value we will hopefully improve the convergence rate for the
		cross coupling of the variables. */
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);
	
private:
	unsigned int _index;									///< Index of the gaseous species to calculate adsorbtion heat for
	const MaterialProperty< MAGPIE_DATA > & _magpie_dat;	///< Material Property holding the MAGPIE data structure
	const MaterialProperty<std::vector<Real> > & _ads_heat;	///< Reference to the heat of adsorption coefficient
	
};
