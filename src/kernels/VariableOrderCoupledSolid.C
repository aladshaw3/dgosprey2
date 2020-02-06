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
    InputParameters params = validParams<TimeDerivative>();
    
    params.addRequiredCoupledVar("main_variable","Name of the non-linear variable that this kernel acts on");
    params.addParam<Real>("max_capacity",0.0,"Maximum Capacity for the Rate Function (mol/kg)");
    params.addParam<Real>("forward_rate",0.0,"Forward reaction rate constant (per hour)");
    params.addParam<Real>("reverse_rate",0.0,"Reverse reaction rate constant (per hour)");
    params.addParam<Real>("site_order",0.0,"Reaction order for the adsorption site");
    params.addParam<Real>("main_order",0.0,"Reaction order for the main variable site");
    params.addParam<Real>("main_stoichiometry",0.0,"Stoichiometry of the main variable");
    
    params.addParam< std::vector<Real> >("species_stoichiometry","Stoichiometric coefficients for the reaction species: (-) sign = reactants and (+) sign for products");
    params.addParam< std::vector<Real> >("species_order","Reaction order for the reaction species: (-) sign = reactants and (+) sign for products");
    params.addParam< std::vector<Real> >("catalyst_order","Reaction order for catalyst in the reaction: (-) sign = reactants and (+) sign for products");
    params.addParam< std::vector<Real> >("adsorbed_sites","Adsorption site coefficients for the number of sites used in each reaction");
    params.addParam< std::vector<Real> >("adsorbed_stoichiometry","Adsorption stoichiometric coefficients for the number of adsorbed species produced in each reaction");
    params.addParam< std::vector<Real> >("adsorbed_order","Reaction order for the adsorbed species in each reaction");
    
    params.addRequiredCoupledVar("coupled_species","List of names of the variables being coupled");
    params.addRequiredCoupledVar("coupled_adsorption","List of names of the adsorbed variables being coupled in this reaction");
    params.addRequiredCoupledVar("coupled_all_adsorption","List of names of  all adsorbed variables that occupy adsorption sites");
    params.addRequiredCoupledVar("coupled_catalysts","List of names of the catalytic variables being coupled");
    
    return params;
}

VariableOrderCoupledSolid::VariableOrderCoupledSolid(const InputParameters & parameters)
: TimeDerivative(parameters),
_maxcap(getParam<Real>("max_capacity")),
_forward(getParam<Real>("forward_rate")),
_reverse(getParam<Real>("reverse_rate")),
_so(getParam<Real>("site_order")),
_main_order(getParam<Real>("main_order")),
_main_stoich(getParam<Real>("main_stoichiometry")),
_spec_stoich(getParam<std::vector<Real> >("species_stoichiometry")),
_spec_order(getParam<std::vector<Real> >("species_order")),
_cat_order(getParam<std::vector<Real> >("catalyst_order")),
_ads_sites(getParam<std::vector<Real> >("adsorbed_sites")),
_ads_stoich(getParam<std::vector<Real> >("adsorbed_stoichiometry")),
_ads_order(getParam<std::vector<Real> >("adsorbed_order")),
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
    
    if (_coupled_species.size() != _spec_stoich.size())
        Moose::out << "ERROR!!! Vectors for coupled solids and solids stoichiometry do not match in size!\n\n";
    if (_coupled_species.size() != _spec_order.size())
        Moose::out << "ERROR!!! Vectors for coupled solids and solids reaction order do not match in size!\n\n";
    
    unsigned int cat_n = coupledComponents("coupled_catalysts");
    _coupled_cat_vars.resize(cat_n);
    _coupled_cat.resize(cat_n);
    
    for (unsigned int i = 0; i<_coupled_cat.size(); ++i)
    {
        _coupled_cat_vars[i] = coupled("coupled_catalysts",i);
        _coupled_cat[i] = &coupledValue("coupled_catalysts",i);
    }
    
    if (_coupled_cat.size() != _cat_order.size())
        Moose::out << "ERROR!!! Vectors for coupled catalysts and catalyst stoichiometry do not match in size!\n\n";
    
    unsigned int ads_n = coupledComponents("coupled_adsorption");
    _coupled_ads_vars.resize(ads_n);
    _coupled_ads.resize(ads_n);
    
    for (unsigned int i = 0; i<_coupled_ads.size(); ++i)
    {
        _coupled_ads_vars[i] = coupled("coupled_adsorption",i);
        _coupled_ads[i] = &coupledValue("coupled_adsorption",i);
    }
    
    if (_coupled_ads.size() != _ads_stoich.size())
        Moose::out << "ERROR!!! Vectors for coupled adsorption and adsorption stoichiometry do not match in size!\n\n";
    if (_coupled_ads.size() != _ads_order.size())
        Moose::out << "ERROR!!! Vectors for coupled adsorption and adsorption reaction order do not match in size!\n\n";
 
    unsigned int all_ads_n = coupledComponents("coupled_all_adsorption");
    _coupled_all_ads_vars.resize(all_ads_n);
    _coupled_all_ads.resize(all_ads_n);
    
    for (unsigned int i = 0; i<_coupled_all_ads.size(); ++i)
    {
        _coupled_all_ads_vars[i] = coupled("coupled_all_adsorption",i);
        _coupled_all_ads[i] = &coupledValue("coupled_all_adsorption",i);
        
        if (_coupled_all_ads_vars[i] == _coupled_var_i)
            _ads_index = i;
    }
    
    if (_coupled_all_ads.size() != _ads_sites.size())
        Moose::out << "ERROR!!! Vectors for all coupled adsorption and adsorption sites do not match in size!\n\n";
}

Real VariableOrderCoupledSolid::computeRateFunction()
{
    
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int i=0; i<_coupled_species.size(); i++)
    {
        if (_spec_stoich[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_species[i])[_qp],_spec_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_species[i])[_qp],_spec_order[i]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_cat_order[j]);
    }
    
    for (unsigned int k=0; k<_coupled_ads.size(); k++)
    {
            if (_ads_stoich[k] >= 0.0)
            {
                prod = prod * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
            else
            {
                react = react * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
    }
    
    return _main_stoich*cat*(react*_forward*pow(computeSiteBalance(),_so) - _reverse*prod*pow(_u[_qp],_main_order));
}

Real VariableOrderCoupledSolid::computeRateFunctionJacobi()
{
    
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int i=0; i<_coupled_species.size(); i++)
    {
        if (_spec_stoich[i] >= 0.0)
        {
            prod = prod * pow((*_coupled_species[i])[_qp],_spec_order[i]);
        }
        else
        {
            react = react * pow((*_coupled_species[i])[_qp],_spec_order[i]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_cat_order[j]);
    }
    
    for (unsigned int k=0; k<_coupled_ads.size(); k++)
    {
            if (_ads_stoich[k] >= 0.0)
            {
                prod = prod * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
            else
            {
                react = react * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
    }
    
    return _main_stoich*cat*(_so*_forward*pow(computeSiteBalance(),_so-1) *react*computeSiteBalanceOffDiagJacobi(_ads_index) - _main_order*_reverse*pow(_u[_qp],_main_order-1)*prod*_phi[_j][_qp]);
    
}

Real VariableOrderCoupledSolid::computeRateFunctionSpeciesOffDiagJacobi(int i)
{
    //Is index i a reactant or product?
    //If _spec_stoich[i] > 0, prod; else, react;
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int j=0; j<_coupled_species.size(); j++)
    {
        if (j != i)
        {
            if (_spec_stoich[j] >= 0.0)
            {
                prod = prod * pow((*_coupled_species[j])[_qp],_spec_order[j]);
            }
            else
            {
                react = react * pow((*_coupled_species[j])[_qp],_spec_order[j]);
            }
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_cat_order[j]);
    }
    
    for (unsigned int k=0; k<_coupled_ads.size(); k++)
    {
            if (_ads_stoich[k] >= 0.0)
            {
                prod = prod * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
            else
            {
                react = react * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
    }
    
    if (_spec_stoich[i] >= 0.0)
    {
        return -_main_stoich*cat*_reverse*pow(_u[_qp],_main_order)*prod*_spec_order[i]*pow((*_coupled_species[i])[_qp],_spec_order[i]-1)*_phi[_j][_qp];
    }
    else
    {
        return _main_stoich*cat*_forward*pow(computeSiteBalance(),_so)*react*_spec_order[i]*pow((*_coupled_species[i])[_qp],_spec_order[i]-1)*_phi[_j][_qp];
    }
}

Real VariableOrderCoupledSolid::computeRateFunctionCatalystOffDiagJacobi(int i)
{
    
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int j=0; j<_coupled_species.size(); j++)
    {
        if (_spec_stoich[j] >= 0.0)
        {
            prod = prod * pow((*_coupled_species[j])[_qp],_spec_order[j]);
        }
        else
        {
            react = react * pow((*_coupled_species[j])[_qp],_spec_order[j]);
        }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        if (j != i)
        {
            cat = cat * pow((*_coupled_cat[j])[_qp],_cat_order[j]);
        }
    }
    
    for (unsigned int k=0; i<_coupled_ads.size(); k++)
    {
            if (_ads_stoich[k] >= 0.0)
            {
                prod = prod * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
            else
            {
                react = react * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
    }
    
    return _main_stoich*_cat_order[i]*pow((*_coupled_cat[i])[_qp],(_cat_order[i]-1))*cat*((_forward*pow(computeSiteBalance(),_so)*react) - (_reverse*pow(_u[_qp],_main_order)*prod))*_phi[_j][_qp];
    
}

Real VariableOrderCoupledSolid::computeRateFunctionAdsOffDiagJacobi(int i)
{
    
    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int j=0; j<_coupled_species.size(); j++)
    {
            if (_spec_stoich[j] >= 0.0)
            {
                prod = prod * pow((*_coupled_species[j])[_qp],_spec_order[j]);
            }
            else
            {
                react = react * pow((*_coupled_species[j])[_qp],_spec_order[j]);
            }
    }
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_cat_order[j]);
    }
    
    for (unsigned int k=0; i<_coupled_ads.size(); k++)
    {
        if (k != i)
        {
            if (_ads_stoich[k] >= 0.0)
            {
                prod = prod * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
            else
            {
                react = react * pow((*_coupled_ads[k])[_qp],_ads_order[k]);
            }
        }
    }
    
    if (_ads_stoich[i] >= 0.0)
    {
        return  _main_stoich*cat*(react*_forward*_so*(pow(computeSiteBalance(),_so-1))*computeSiteBalanceOffDiagJacobi(i) - _reverse*prod*_ads_order[i]*pow((*_coupled_ads[i])[_qp],_ads_order[i]-1)*pow(_u[_qp],_main_order)*_phi[_j][_qp]);
    }
    
    else
    {
        return _main_stoich*cat*react*_forward*(_so*(pow(computeSiteBalance(),_so-1))*computeSiteBalanceOffDiagJacobi(i)*pow((*_coupled_ads[i])[_qp],_ads_order[i])+_ads_order[i]*pow((*_coupled_ads[i])[_qp],_ads_order[i]-1)*computeSiteBalance()*_phi[_j][_qp]);
    }

}

Real VariableOrderCoupledSolid::computeSiteBalance()
{
    double sum = 0.0;
    for (unsigned int i=0; i<_coupled_all_ads.size(); i++)
    {
        sum += _ads_sites[i]*(*_coupled_all_ads[i])[_qp];
    }
    return _maxcap - sum;
}

Real VariableOrderCoupledSolid::computeSiteBalanceOffDiagJacobi(int i)
{
    return -_ads_sites[i]*_phi[_j][_qp];
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
    //Off-diagonals for species concentrations
    for (unsigned int i = 0; i<_coupled_species.size(); ++i)
    {
        if (jvar == _coupled_species_vars[i] && jvar != _coupled_var_i)
        {
            return -_test[_i][_qp]*computeRateFunctionSpeciesOffDiagJacobi(i);
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
