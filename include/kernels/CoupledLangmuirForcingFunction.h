/*!
 *  \file CoupledLangmuirForcingFunction.h
 *	\brief Standard kernel for coupling non-linear variables via a langmuir function
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via a langmuir forcing function, i.e., variable = b * K * coupled_variable / 1 + K * coupled_variable.
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 06/07/2017
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

#include "TimeDerivative.h"

/// CoupledLangmuirForcingFunction class object forward declarationss
class CoupledLangmuirForcingFunction;

template<>
InputParameters validParams<CoupledLangmuirForcingFunction>();

/// CoupledLangmuirForcingFunction class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the two non-linear variables to couple a Langmuir forcing
	function between given objects. */
class CoupledLangmuirForcingFunction : public TimeDerivative
{
public:
    /// Required constructor for objects in MOOSE
    CoupledLangmuirForcingFunction(const InputParameters & parameters);
    
protected:
	/// Function to compute the Langmuir equilibrium
	Real computeLangmuirEquilibrium();
	
	/// Function to compute the Langmuir derivative with respect to concentration
	Real computeLangConcDerivative();
	
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
    
    Real _maxcap;							///< Manimum Capacity
    Real _langmuircoef;						///< Langmuir Coefficient for the coupled function
    const VariableValue & _coupled_u;		///< Coupled variable
    const unsigned int _coupled_var;		///< Variable identification for the coupled variable
    
private:
    
};
