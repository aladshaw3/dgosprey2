/*!
 *  \file CoupledGSTAmodel.h
 *	\brief Standard kernel for coupling non-linear variables via the GSTA model
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via the GSTA model, which includes estimating model parameters from material
 *			properties files. The primary difference between this kernel and the CoupledGSTAisotherm
 *			kernel is the coupling of the material properties. 
 *
 *			This kernel extends the CoupledGSTAisotherm kernel by calculating the model parameters from
 *			information in the ThermodynamicProperties material. The Kno parameters (described below) are
 *			to be estimated from the site enthalpies (dHno) and entropies (dSno) using the van't Hoff 
 *			expression (shown below). Thus, this model is inherently a function of temperature and will
 *			require a different form of coupling with the temperature parameter. 
 *
 *			van't Hoff: ln(Kno) = -dHno/(R*T) + dSno/R
 *			where R is the gas law constant and T is the column temperature.
 *
 *			GSTA isotherm: q = (q_max / m) * SUM(n*Kno*(p/Po)^n)/(1+SUM(Kno*(p/Po)^n))
 *			where q is amount adsorbed, q_max is the maximum capacity, m is the number of adsorption sites
 *			and Kno are the dimensionless equilibrium parameters. Also, p is partial pressure of gas and Po
 *			is taken as a reference state pressure (100 kPa).
 *
 *	\note	For the use of this kernel, our coupled variable with be a gas concentration in mol/L (C), therefore,
 *			we need to use ideal gas law to rewrite the GSTA model in terms of C as opposed to p. Thus, we are also
 *			forced to couple with kernel with column temperature.
 *
 *			Ideal Gas Law: p = C*R*T
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 08/28/2017
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

#include "CoupledGSTAisotherm.h"
#include "flock.h"
#include "DataStruct_StoreLoad.h"

/// CoupledGSTAmodel class object forward declarationss
class CoupledGSTAmodel;

template<>
InputParameters validParams<CoupledGSTAmodel>();

/// CoupledGSTAmodel class object inherits from CoupledGSTAisotherm object
/** This class object inherits from the CoupledGSTAisotherm object in the DGOSPREY framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the two non-linear variables to couple the GSTA isotherm
	model between given objects. Parameters of the GSTA model are determined through the
	parameters in the ThermodynamicProperties file. */
class CoupledGSTAmodel : public CoupledGSTAisotherm
{
public:
	/// Required constructor for objects in MOOSE
	CoupledGSTAmodel(const InputParameters & parameters);
	
protected:
	/// Function to compute the isotherm parameters based on temperature
	void computeGSTAparams();
	
	/// Function to compute the isotherm derivative with respect to temperature
	Real computeGSTAtempDerivative();
	
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
	
	unsigned int _index;									///< Index of the gaseous species to calculate adsorbtion for
	const MaterialProperty< MAGPIE_DATA > & _magpie_dat;	///< Material Property holding the MAGPIE data structure
	
private:
	
};
