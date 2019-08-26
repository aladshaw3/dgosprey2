/*!
 *  \file DGMassFluxBC.h
 *	\brief Boundary Condition kernel for the mass flux in and out of the ends of the fixed-bed column
 *	\details This file creates a boundary condition kernel for the mass flux across the boundary of
 *			the ends of the column in the fixed-bed adsorber. It inherits from the DGFluxBC, which
 *			acts as a generic flux BC module. This kernel is coupled to the column temperature variable,
 *			as well as the material properties for thermal conductivity, gas density and heat capacity,
 *			and the velocity in the domain.
 *
 *			This type of boundary condition for DG kernels applies the true flux boundary condition.
 *			Alternatively, you can use the "FluxLimitedBC" to impose a Dirichlet boundary condition
 *			on the system. Although, in true finite volumes or DG methods, there is no Dirichlet
 *			boundary conditions, because the solutions are based on fluxes into and out of cells in
 *			a domain.
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

#include "DGFluxBC.h"

/// DGMassFluxBC class object forward declaration
class DGMassFluxBC;

template<>
InputParameters validParams<DGMassFluxBC>();

/// DGMassFluxBC class object inherits from DGFluxBC object
/** This class object inherits from the DGFluxBC object (see DGFluxBC.h for more details).
	All public and protected members of this class are required function overrides. The object
	will take in the given variable and material properties to override some objects declared
	for the generic DGFluxBC to fit this particular boundary condition. Then, it just calls the
	appropriate DGFluxBC functions. */
class DGMassFluxBC : public DGFluxBC
{
public:
	/// Required constructor for BC objects in MOOSE
	DGMassFluxBC(const InputParameters & parameters);
	
protected:
	/// Required function override for BC objects in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	/// Required function override for BC objects in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
    
    /// Value of the column temperature at the inlet of the system
    Real _input_temperature;
    /// Value of the column total pressure at the inlet of the system
    Real _input_pressure;
    /// Value of the molefraction of the specific species at the inlet of the system
    Real _input_molefraction;
    
    const MaterialProperty<Real> & _vel;								///< Reference to the velocity material property
    unsigned int _index;												///< Index of the species of interest at the boundary
	
private:
	
};
