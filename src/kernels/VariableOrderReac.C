/*!
 *  \file VariableOrderReac.h
 *    \brief Standard kernel for coupling multiple species via a reaction based mechanism
 *    \details This file creates a standard MOOSE kernel for the coupling species in a
 *                simulation based on the following reaction scheme...
 *
 *                sum(i, v_i*C_i) <-- --> n*q + sum(j, v_j*C_j)
 *
 *                In this reaction scheme, the i-th species may interact to form n number of
 *                species (q). In return, some other gas species may be produced (C_j) as
 *                by-products from the reaction. This expression then formulates the following
 *                rate equation for the coupled species (q)...
 *
 *                dq/dt = n*k_f*product(i, C_i^v_i) - n*k_r*(q)^n*product(j, C_j^v_j)
 *
 *                Parameters are as follows...
 *
 *                v_i,j = stoichiometric coefficients for reactants/products
 *                C_i,j = concentrations of reactants/products (mol/L or mol/kg)
 *                n = number of product species (q) produced
 *                q = concentration of coupled species produced (mol/kg)
 *                k_f,r = rate constant for the forward/reverse reaction
 *
 *  \author Alexander Wiechert
 *    \date 08/29/2019
 *    \copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *             by Idaho National Laboratory and Oak Ridge National Laboratory
 *             engineers and scientists. Portions Copyright (c) 2018, all
 *             rights reserved.
 *
 *             Alexander Wiechert does not claim any ownership or copyright to the
 *             MOOSE framework in which these kernels are constructed, only
 *             the kernels themselves. The MOOSE framework copyright is held
 *             by the Battelle Energy Alliance, LLC (c) 2010, all rights reserved.
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

#include "VariableOrderReac.h"
registerMooseObject("dgospreyApp", VariableOrderReac);

template<>
InputParameters validParams<VariableOrderReac>()
{
    InputParameters params = validParams<TimeDerivative>();
    params.addParam<Real>("forward_rate",0.0,"Forward reaction rate constant (per hour)");
    params.addParam<Real>("reverse_rate",0.0,"Reverse reaction rate constant (per hour)");
    params.addParam<Real>("main_stoichiometry",0.0,"Stoichiometry of the main product variable of this reaction");
    params.addParam<Real>("main_order",0.0,"Order of the main product variable of this reaction");
    
    params.addParam< std::vector<Real> >("stoichiometry","Species stoichiometric coefficients: (-) sign = reactants and (+) sign for products");
    params.addParam< std::vector<Real> >("order","Order of species in reaction");
    
    params.addRequiredCoupledVar("coupled_species","List of names of variables being coupled for the rate equation");
    params.addRequiredCoupledVar("main_variable","Name of the non-linear variable that this kernel acts on");
    
    return params;
}

VariableOrderReac::VariableOrderReac(const InputParameters & parameters)
: TimeDerivative(parameters),
_forward(getParam<Real>("forward_rate")),
_reverse(getParam<Real>("reverse_rate")),
_m_stoich(getParam<Real>("main_stoichiometry")),
_m_order(getParam<Real>("main_order")),

_stoich(getParam<std::vector<Real> >("stoichiometry")),
_order(getParam<std::vector<Real> >("order")),

_coupled_var_i(coupled("main_variable"))
{
    unsigned int species_n = coupledComponents("coupled_species");
    _coupled_species_vars.resize(species_n);
    _coupled_species.resize(species_n);
    
    for (unsigned int i = 0; i<_coupled_species.size(); ++i)
    {
        _coupled_species_vars[i] = coupled("coupled_species",i);
        _coupled_species[i] = &coupledValue("coupled_species",i);
    }
    
    if (_coupled_species.size() != _stoich.size())
        Moose::out << "ERROR!!! Vectors for coupled species and species stoichiometry do not match in size!\n\n";
    if (_coupled_species.size() != _order.size())
        Moose::out << "ERROR!!! Vectors for coupled species and species reaction order do not match in size!\n\n";
}

Real VariableOrderReac::computeRateFunction()
{
    double react = 1.0, prod = 1.0;
    for (unsigned int i=0; i<_coupled_species.size(); i++)
    {
        if (_stoich[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_species[i])[_qp],_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_species[i])[_qp],_order[i]);
        }
    }
    return _m_stoich*(_forward*react - _reverse*pow(_u[_qp],_m_order)*prod);
}

Real VariableOrderReac::computeRateFunctionJacobi()
{
    double react = 1.0, prod = 1.0;
    for (unsigned int i=0; i<_coupled_species.size(); i++)
    {
        if (_stoich[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_species[i])[_qp],_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_species[i])[_qp],_order[i]);
        }
    }
    return -_m_stoich*_m_order*_reverse*pow(_u[_qp],_m_order-1)*prod*_phi[_j][_qp];
}

Real VariableOrderReac::computeRateFunctionSpeciesOffDiagJacobi(int i)
{
    //Is index i a reactant or product?
    //If _species_stoich[i] > 0, prod; else, react;
    double react = 1.0, prod = 1.0;
    for (unsigned int j=0; j<_coupled_species.size(); j++)
    {
        if (j != i)
        {
            if (_stoich[i] >= 0.0)
            {
                prod = prod * pow((*_coupled_species[i])[_qp],_order[i]);
            }
            else
            {
                react = react * pow((*_coupled_species[i])[_qp],_order[i]);
            }
        }
    }
    
    if (_stoich[i] >= 0.0)
    {
        return -_m_stoich*_reverse*pow(_u[_qp],_m_order)*prod*_order[i]*pow((*_coupled_species[i])[_qp],_order[i]-1)*_phi[_j][_qp];
    }
    else
    {
        return _m_stoich*_forward*react*_order[i]*pow((*_coupled_species[i])[_qp],_order[i]-1)*_phi[_j][_qp];
    }
}

Real VariableOrderReac::computeQpResidual()
{
    return _test[_i][_qp]*_u_dot[_qp] - _test[_i][_qp]*computeRateFunction();
}

Real VariableOrderReac::computeQpJacobian()
{
    return _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp] - _test[_i][_qp]*computeRateFunctionJacobi();
}

Real VariableOrderReac::computeQpOffDiagJacobian(unsigned int jvar)
{
    //Off-diagonals for concentrations
    for (unsigned int i = 0; i<_coupled_species.size(); ++i)
    {
        if (jvar == _coupled_species_vars[i] && jvar != _coupled_var_i)
        {
            return -_test[_i][_qp]*computeRateFunctionSpeciesOffDiagJacobi(i);
        }
    }
    return 0.0;
}



