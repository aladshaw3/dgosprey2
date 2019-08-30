/*!
 *  \file VariableOrderCoupledSolid.C
 *  \brief Standard kernel for coupling multiple gas, solid and adsorbed species
 *   together via a reaction based mechanism
 *
 *  \details This file creates a standard MOOSE kernel for the coupling multiple gas,
 *   adsorption, solid and catalytic species in a simulation based on the following
 *   reaction scheme...
 *
 *                sum(i, v_i*C_i) + m*C + n*L  <-- --> x*D + m*C + sum(j, v_j*C_j)
 *
 *                In this reaction scheme, the i-th species and m number of a catalytic
 *                species (C) may interact with m number of available surface sites (L)
 *                to form x number of deactivated sites (D). The catalytic speices is
 *                returned along with some other gas species that may be produced (C_j)
 *                as by-products from the reaction.
 *
 *                This expression then formulates the following rate equation for the
 *                deactivated site (D)...
 *
 *        dD/dt = x*(C)^a*[k_f*(L)^b*product(i, C_i^z_i) - k_r*(D)^f*product(j, C_j^z_j)]
 *
 *                Parameters are as follows...
 *
 *                v_i,j = stoichiometric coefficients for the reactants/products
 *                z_i,j = reaction oder with respect to the reactants/products
 *                C_i,j = concentrations of reactants/products (mol/L or mol/kg)
 *                m = number of catalysts used in the reaction
 *                n = number of adsorption sites (L) reduced
 *                x = number of reduced adsortion sites (D) produced
 *                L = concentration of available sites (mol/kg)
 *                D = concentration of reduced species produced (mol/kg)
 *                k_f,r = rate constant for the forward/reverse reaction
 *                a = reaction order with repect to the catalyst
 *                b = reaction order with respect to available adsorption sites
 *                f = reaction order with respect to the reduced species
 *
 *                Rate expression must be coupled with all involved adsorbed species and
 *                uses a site-balance to account for the loss of adsorption sites during
 *                multi-species adsorption
 *
 *                SiteBalance: (L) = Lmax - (D) - sum(i, m_i/n_i*q_i)
 *
 *                where Lmax is the maximum capacity for adsorption (or the maximum available
 *                sites), m_i is the number of sites the m-th adsorbate occupies, n_i is the
 *                number of adsorbed species produced from the reaction consuming m_i sites,
 *                and q_i is the adsorbed concentration of the i-th species.
 *
 *  \author Alexander Wiechert
 *   \date 8/30/2019
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

#include "VariableOrderCoupledSolid.h"
registerMooseObject("dgospreyApp", VariableOrderCoupledSolid);

template<>
InputParameters validParams<VariableOrderCoupledSolid>()
{
    InputParameters params = validParams<VariableOrderCoupledCatalyst>();
    
    params.addParam< std::vector<Real> >("solids_order","Reaction order for each solid species");
    params.addParam< std::vector<Real> >("solids_stoichiometry","Stoichiometry of each solid species");
    
    params.addRequiredCoupledVar("coupled_solids","List of names of the solid variables being coupled");
    
    return params;
}

VariableOrderCoupledSolid::VariableOrderCoupledSolid(const InputParameters & parameters)
: VariableOrderCoupledCatalyst(parameters),
_sol_order(getParam<std::vector<Real> >("solids_order")),
_sol_stoic(getParam<std::vector<Real> >("solids_stoichiometry"))
{
    unsigned int solid_n = coupledComponents("coupled_solids");
    _coupled_solid_vars.resize(solid_n);
    _coupled_solid.resize(solid_n);
    
    for (unsigned int i = 0; i<_coupled_solid.size(); ++i)
    {
        _coupled_solid_vars[i] = coupled("coupled_solids",i);
        _coupled_solid[i] = &coupledValue("coupled_solids",i);
        
    }
    
    if (_coupled_solid.size() != _sol_stoic.size())
        Moose::out << "ERROR!!! Vectors for coupled solids and solids stoichiometry do not match in size!\n\n";
    if (_coupled_solid.size() != _sol_order.size())
        Moose::out << "ERROR!!! Vectors for coupled solids and solids reaction order do not match in size!\n\n";
}

Real VariableOrderCoupledSolid::computeRateFunction()
{
    
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int i=0; i<_coupled_gas.size(); i++)
    {
        if (_gas_stoich[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_gas[i])[_qp],_g_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_gas[i])[_qp],_g_order[i]);
        }
    }
    
    for (unsigned int i=0; i<_coupled_solid.size(); i++)
    {
        if (_sol_stoic[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_c_order[j]);
    }
    
    return _ads_stoich[_ads_index]*cat*( (_forward*pow(CoupledConstChemisorption::computeSiteBalance(),_s_order[_ads_index])*react) - (_reverse*pow(_u[_qp],_a_order[_ads_index])*prod) );
}

Real VariableOrderCoupledSolid::computeRateFunctionJacobi()
{
    
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int i=0; i<_coupled_gas.size(); i++)
    {
        if (_gas_stoich[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_gas[i])[_qp],_g_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_gas[i])[_qp],_g_order[i]);
        }
    }
    
    for (unsigned int i=0; i<_coupled_solid.size(); i++)
    {
        if (_sol_stoic[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_c_order[j]);
    }
    
    return _ads_stoich[_ads_index]*cat*(_s_order[_ads_index]*_forward*pow(CoupledConstChemisorption::computeSiteBalance(),_s_order[_ads_index]-1)*react*CoupledConstChemisorption::computeSiteBalanceOffDiagJacobi(_ads_index) - _a_order[_ads_index]*_reverse*pow(_u[_qp],_a_order[_ads_index]-1)*prod*_phi[_j][_qp]);
    
}

Real VariableOrderCoupledSolid::computeRateFunctionGasOffDiagJacobi(int i)
{
    //Is index i a reactant or product?
    //If _gas_stoich[i] > 0, prod; else, react;
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int j=0; j<_coupled_gas.size(); j++)
    {
        if (j != i)
        {
            if (_gas_stoich[j] >= 0.0)
            {
                prod = prod * pow((*_coupled_gas[j])[_qp],_g_order[j]);
            }
            else
            {
                react = react * pow((*_coupled_gas[j])[_qp],_g_order[j]);
            }
        }
    }
    
    for (unsigned int i=0; i<_coupled_solid.size(); i++)
    {
        if (_sol_stoic[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_c_order[j]);
    }
    
    if (_gas_stoich[i] >= 0.0)
    {
        return -_ads_stoich[_ads_index]*cat*_reverse*pow(_u[_qp],_a_order[_ads_index])*prod*_g_order[i]*pow((*_coupled_gas[i])[_qp],_g_order[i]-1)*_phi[_j][_qp];
    }
    else
    {
        return _ads_stoich[_ads_index]*cat*_forward*pow(CoupledConstChemisorption::computeSiteBalance(),_s_order[_ads_index])*react*_g_order[i]*pow((*_coupled_gas[i])[_qp],_g_order[i]-1)*_phi[_j][_qp];
    }
}

Real VariableOrderCoupledSolid::computeRateFunctionCatalystOffDiagJacobi(int i)
{
    
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int j=0; j<_coupled_gas.size(); j++)
    {
        if (_gas_stoich[j] >= 0.0)
        {
            prod = prod * pow((*_coupled_gas[j])[_qp],_g_order[j]);
        }
        else
        {
            react = react * pow((*_coupled_gas[j])[_qp],_g_order[j]);
        }
    }
    
    for (unsigned int i=0; i<_coupled_solid.size(); i++)
    {
        if (_sol_stoic[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        if (j != i)
        {
            cat = cat * pow((*_coupled_cat[j])[_qp],_c_order[j]);
        }
    }
    
    return _ads_stoich[_ads_index]*_c_order[i]*pow((*_coupled_cat[i])[_qp],(_c_order[i]-1))*cat*((_forward*pow(CoupledConstChemisorption::computeSiteBalance(),_s_order[_ads_index])*react) - (_reverse*pow(_u[_qp],_a_order[_ads_index])*prod));
    
}

Real VariableOrderCoupledSolid::computeRateFunctionSolidsOffDiagJacobi(int i)
{
    //Is index i a reactant or product?
    //If _sol_stoic[i] > 0, prod; else, react;
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int j=0; j<_coupled_solid.size(); j++)
    {
        if (j != i)
        {
            if (_sol_stoic[i] >= 0.0)
            {
                prod = prod * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
            }
            else
            {
                react = react * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
            }
        }
    }
    
    for (unsigned int j=0; j<_coupled_gas.size(); j++)
    {
        if (_gas_stoich[j] >= 0.0)
        {
            prod = prod * pow((*_coupled_gas[j])[_qp],_g_order[j]);
        }
        else
        {
            react = react * pow((*_coupled_gas[j])[_qp],_g_order[j]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_c_order[j]);
    }
    
    if (_sol_stoic[i] >= 0.0)
    {
        return -_ads_stoich[_ads_index]*cat*_reverse*pow(_u[_qp],_a_order[_ads_index])*prod*_sol_order[i]*pow((*_coupled_solid[i])[_qp],_sol_order[i]-1)*_phi[_j][_qp];
    }
    else
    {
        return _ads_stoich[_ads_index]*cat*_forward*pow(CoupledConstChemisorption::computeSiteBalance(),_s_order[_ads_index])*react*_sol_order[i]*pow((*_coupled_solid[i])[_qp],_sol_order[i]-1)*_phi[_j][_qp];
    }
}

Real VariableOrderCoupledSolid::computeRateFunctionAdsOffDiagJacobi(int i)
{
    
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int i=0; i<_coupled_gas.size(); i++)
    {
        if (_gas_stoich[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_gas[i])[_qp],_g_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_gas[i])[_qp],_g_order[i]);
        }
    }
    
    for (unsigned int i=0; i<_coupled_solid.size(); i++)
    {
        if (_sol_stoic[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_solid[i])[_qp],_sol_order[i]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_c_order[j]);
    }
    
    return _ads_stoich[_ads_index]*cat*_s_order[_ads_index]*_forward*pow(CoupledConstChemisorption::computeSiteBalance(),_s_order[_ads_index]-1)*react*CoupledConstChemisorption::computeSiteBalanceOffDiagJacobi(i);
}

Real VariableOrderCoupledSolid::computeQpResidual()
{
    return _test[_i][_qp]*_u_dot[_qp] - _test[_i][_qp]*computeRateFunction();
}

Real VariableOrderCoupledSolid::computeQpJacobian()
{
    return _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp] - _test[_i][_qp]*computeRateFunctionJacobi();
}

Real VariableOrderCoupledSolid::computeQpOffDiagJacobian(unsigned int jvar)
{
    //Off-diagonals for gas concentrations
    for (unsigned int i = 0; i<_coupled_gas.size(); ++i)
    {
        if (jvar == _coupled_gas_vars[i] && jvar != _coupled_var_i)
        {
            return -_test[_i][_qp]*computeRateFunctionGasOffDiagJacobi(i);
        }
    }
    //Off-diagonals for catalyst concentrations
    for (unsigned int i = 0; i<_coupled_cat.size(); ++i)
    {
        if (jvar == _coupled_cat_vars[i] && jvar != _coupled_var_i)
        {
            return -_test[_i][_qp]*computeRateFunctionCatalystOffDiagJacobi(i);
        }
    }
    //Off-diagonals for solid concentrations
    for (unsigned int i = 0; i<_coupled_solid.size(); ++i)
    {
        if (jvar == _coupled_solid_vars[i] && jvar != _coupled_var_i)
        {
            return -_test[_i][_qp]*computeRateFunctionSolidsOffDiagJacobi(i);
        }
    }
    //Off-diagonals for adsorption
    for (unsigned int i = 0; i<_coupled_ads.size(); ++i)
    {
        if (jvar == _coupled_ads_vars[i] && jvar != _coupled_var_i)
        {
            return -_test[_i][_qp]*computeRateFunctionAdsOffDiagJacobi(i);
        }
    }
    return 0.0;
}
