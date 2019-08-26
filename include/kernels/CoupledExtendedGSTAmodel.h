/*!
 *  \file CoupledExtendedGSTAmodel.h
 *	\brief Standard kernel for coupling a vector non-linear variables via an extended gsta model
 *	\details This file creates a standard MOOSE kernel for the coupling of a vector non-linear variables
 *			together via an extended gsta forcing function and computes gsta parameters as a function
 *			of temperature by coupling with the Thermodynamic Properties Material property.
 *
 *			                                    SUM(n, Kn_i * (coupled_variable_i/Co)^n)
 *					variable = qmax_i * ----------------------------------------------------
 *                                       1 + SUM(j, SUM(n, Kn_j * (coupled_variable_j/Co)^n)
 *
 *						where Co = 100.0 / (8.3144621 * _coupled_temp)
 *						and  Kn_i = exp( -dHn_i/RT + dSn_i/R )
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 05/10/2018
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2018, all
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

#include "CoupledExtendedGSTAisotherm.h"
#include "flock.h"
#include "DataStruct_StoreLoad.h"

/// CoupledExtendedGSTAmodel class object forward declarationss
class CoupledExtendedGSTAmodel;

template<>
InputParameters validParams<CoupledExtendedGSTAmodel>();

/// CoupledExtendedGSTAisotherm class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the set of non-linear variables to couple an extended GSTA
	forcing function between given objects. */
class CoupledExtendedGSTAmodel : public CoupledExtendedGSTAisotherm
{
public:
	/// Required constructor for objects in MOOSE
	CoupledExtendedGSTAmodel(const InputParameters & parameters);
	
protected:
	/// Function to compute the off-diagonal Jacobi for the other coupled concentrations
	Real computeExtGSTA_TempOffJacobi();
	
	/// Function to compute derivative of top with respect to temperature
	Real computeTopDerivativeTemp();
	
	/// Function to compute derivative of bottom with respect to temperature
	Real computeBottomDerivativeTemp();
	
	/// Function to compute the isotherm parameters based on temperature
	void computeGSTAparams();
	
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
	
	const MaterialProperty< MAGPIE_DATA > & _magpie_dat;	///< Material Property holding the MAGPIE data structure
	
private:
	
};
