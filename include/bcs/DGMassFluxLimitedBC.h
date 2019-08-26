/*!
 *  \file DGMassFluxLimitedBC.h
 *	\brief Boundary Condition kernel to mimic a dirichlet boundary condition at the column inlet
 *	\details This file creates a dirichlet-like boundary condition kernel for the gas species concentration
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

/// DGMassFluxLimitedBC class object forward declaration
class DGMassFluxLimitedBC;

template<>
InputParameters validParams<DGMassFluxLimitedBC>();

/// DGMassFluxLimitedBC class object inherits from DGFluxLimitedBC object
/** This class object inherits from the DGFluxLimitedBC object (see DGFluxLimitedBC.h for more details).
	All public and protected members of this class are required function overrides. The object
	will take in the given variable and material properties to override some objects declared
	for the generic DGFluxLimitedBC to fit this particular boundary condition. Then, it just calls the
	appropriate DGFluxLimitedBC functions. */
class DGMassFluxLimitedBC : public DGFluxLimitedBC
{
public:
	/// Required constructor for BC objects in MOOSE
	DGMassFluxLimitedBC(const InputParameters & parameters);
	
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
	/// Value of the column total pressure at the inlet of the system
	Real _input_pressure;
	/// Value of the molefraction of the specific species at the inlet of the system
	Real _input_molefraction;
	
	const MaterialProperty<Real> & _vel;								///< Reference to the velocity material property
	unsigned int _index;												///< Index of the species of interest at the boundary
	const MaterialProperty<std::vector<Real> > & _dispersion;			///< Reference to the dispersion coefficient material property
	const MaterialProperty<std::vector<Real> > & _molecular_diffusion;	///< Reference to the molecular diffusion material property
	
};
