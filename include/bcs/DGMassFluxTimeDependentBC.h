/*!
 *  \file DGMassFluxTimeDependentBC.h
 *	\brief Boundary Condition kernel for time dependent mass flux in and out of the ends of the fixed-bed column
 *	\details This file creates a boundary condition kernel for the time dependent mass flux across the boundary of
 *			the ends of the column in the fixed-bed adsorber. It inherits from the DGMassFluxBC, which
 *			is coupled to the column temperature variable, as well as the material properties for
 *          thermal conductivity, gas density and heat capacity, and the velocity in the domain.
 *
 *			This type of boundary condition for DG kernels applies the true flux boundary condition.
 *			Alternatively, you can use the "FluxLimitedBC" to impose a Dirichlet boundary condition
 *			on the system. Although, in true finite volumes or DG methods, there is no Dirichlet
 *			boundary conditions, because the solutions are based on fluxes into and out of cells in
 *			a domain.
 *
 *  \author Alexander Wiechert
 *	\date 7/2/2018
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2018, all
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

#include "DGMassFluxBC.h"

/// DGMassFluxTimeDependentBC class object forward declaration
class DGMassFluxTimeDependentBC;

template<>
InputParameters validParams<DGMassFluxTimeDependentBC>();

/// DGMassFluxTimeDependentBC class object inherits from DGMassFluxBC object
/** This class object inherits from the DGMassFluxBC object (see DGMassFluxBC.h for more details).
    All public and protected members of this class are required function overrides. The object
    will take in the given variable and material properties to override some objects declared
    for the generic DGMassFluxBC to fit this particular boundary condition. Then, it just calls the
    appropriate DGMassFluxBC functions. */
class DGMassFluxTimeDependentBC : public DGMassFluxBC
{
public:
    /// Required constructor for BC objects in MOOSE
    DGMassFluxTimeDependentBC(const InputParameters & parameters);

protected:
    /// Required function override for BC objects in MOOSE
    /** This function returns a residual contribution for this object.*/
    virtual Real computeQpResidual();
    /// Required function override for BC objects in MOOSE
    /** This function returns a Jacobian contribution for this object. The Jacobian being
    computed is the associated diagonal element in the overall Jacobian matrix for the
    system and is used in preconditioning of the linear sub-problem. */
    virtual Real computeQpJacobian();

    /// Value of the start time for input into the system
    Real _s_time;
    /// Value of the end time for input into the system
    Real _e_time;

private:

};
