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


#include "DGMassFluxTimeDependentBC.h"
registerMooseObject("dgospreyApp", DGMassFluxTimeDependentBC);

template<>
InputParameters validParams<DGMassFluxTimeDependentBC>()
{
    InputParameters params = validParams<DGMassFluxBC>();

    params.addParam<Real>("start_time","Starting time for the Boundary Condition (hr)");
    params.addParam<Real>("end_time","Ending time for the Boundary Condition (hr)");

    return params;
}

DGMassFluxTimeDependentBC::DGMassFluxTimeDependentBC(const InputParameters & parameters) :
DGMassFluxBC(parameters),
_s_time(getParam<Real>("start_time")),
_e_time(getParam<Real>("end_time"))
{

}

Real DGMassFluxTimeDependentBC::computeQpResidual()
{
    if (_t >= _s_time && _t <= _e_time)
    {
        return DGMassFluxBC::computeQpResidual();
    }
    else
    {
        //_input_molefraction = 0.0;
        //_u_input = (_input_pressure * _input_molefraction) / (8.3144621 * _input_temperature);
        //return DGFluxBC::computeQpResidual();
        _u_input = 0.0;
        return DGFluxBC::computeQpResidual();
    }
}

Real
DGMassFluxTimeDependentBC::computeQpJacobian()
{
    return DGMassFluxBC::computeQpJacobian();
}

