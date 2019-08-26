/*!
 *  \file CoupledCatalyst.h
 *	\brief Standard kernel for coupling multiple gas, adsorbed, and catalytic species together via a reaction based mechanism
 *	\details This file creates a standard MOOSE kernel for the coupling multiple gas, adsorption, and catalytic species in a
 *				simulation based on the following reaction scheme...
 *
 *				sum(i, v_i*C_i) + m*C + n*L  <-- --> x*D + m*C + sum(j, v_j*C_j)
 *
 *				In this reaction scheme, the i-th species and m number of a catalytic species (C) may interact with
 *              m number of available surface sites (L) to form x number of deactivated sites (D). The catalytic speices
 *              is returned along with some other gas species that may be produced (C_j) as by-products from the reaction.
 *              This expression then formulates the following rate equation for the deactivated site (D)...
 *
 *				dD/dt = x*(C)^m*[k_f*(L)^n*product(i, C_i^v_i) - k_r*(D)^x*product(j, C_j^v_j)]
 *
 *				Parameters are as follows...
 *
 *				v_i,j = stoichiometric coefficients for gas reactants/products
 *				C_i,j = concentrations of gas reactants/products (mol/L)
 *				m = number of catalysts used in the reaction
 *				n = number of adsorption sites (L) reduced
 *              x = number of reduced adsortion sites (D) produced
 *				L = concentration of available sites (mol/kg)
 *				D = concentration of reduced species produced (mol/kg)
 *				k_f,r = rate constant for the forward/reverse reaction
 *
 *				Rate expression must be coupled with all involved adsorbed species and uses a site-balance to
 *				account for the loss of adsorption sites during multi-species adsorption
 *
 *				SiteBalance: (L) = Lmax - (D) - sum(i, m_i/n_i*q_i)
 *
 *				where Lmax is the maximum capacity for adsorption (or the maximum available sites), m_i is the
 *				number of sites the m-th adsorbate occupies, n_i is the number of adsorbed species produced from
 *				the reaction consuming m_i sites, and q_i is the adsorbed concentration of the i-th species.
 *
 *  \author Alexander Wiechert
 *	\date 6/28/2018
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

#include "CoupledCatalyst.h"
registerMooseObject("dgospreyApp", CoupledCatalyst);

template<>
InputParameters validParams<CoupledCatalyst>()
{
    InputParameters params = validParams<CoupledConstChemisorption>();

    params.addParam< std::vector<Real> >("catalyst_stoichiometry","Stoichiometric coefficients for each catalyst species");

    params.addRequiredCoupledVar("coupled_catalysts","List of names of the catalytic variables being coupled");

    return params;
}

CoupledCatalyst::CoupledCatalyst(const InputParameters & parameters)
: CoupledConstChemisorption(parameters),
_cat_stoich(getParam<std::vector<Real> >("catalyst_stoichiometry"))

{
    unsigned int cat_n = coupledComponents("coupled_catalysts");
    _coupled_cat_vars.resize(cat_n);
    _coupled_cat.resize(cat_n);

    for (unsigned int i = 0; i<_coupled_cat.size(); ++i)
    {
        _coupled_cat_vars[i] = coupled("coupled_catalysts",i);
        _coupled_cat[i] = &coupledValue("coupled_catalysts",i);

    }

    if (_coupled_cat.size() != _cat_stoich.size())
    Moose::out << "ERROR!!! Vectors for coupled catalysts and catalyst stoichiometry do not match in size!\n\n";
}

Real CoupledCatalyst::computeRateFunction()
{
    double no_cat = 1.0, cat = 1.0;
    no_cat = CoupledConstChemisorption::computeRateFunction();
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_cat_stoich[j]);
    }

    return no_cat * cat;
}

Real CoupledCatalyst::computeRateFunctionJacobi()
{
    double no_cat = 1.0, cat = 1.0;
    no_cat = CoupledConstChemisorption::computeRateFunctionJacobi();
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_cat_stoich[j]);
    }
    
    return no_cat * cat;
}

Real CoupledCatalyst::computeRateFunctionGasOffDiagJacobi(int i)
{
    double no_cat = 1.0, cat = 1.0;
    no_cat = CoupledConstChemisorption::computeRateFunctionGasOffDiagJacobi(i);
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_cat_stoich[j]);
    }
    
    return no_cat * cat;
}

Real CoupledCatalyst::computeRateFunctionCatalystOffDiagJacobi(int i)
{

    double react = 1.0, prod = 1.0, cat = 1.0;
    for (unsigned int j=0; j<_coupled_gas.size(); j++)
    {
        if (_gas_stoich[j] >= 0.0)
        {
            prod = prod * pow((*_coupled_gas[j])[_qp],_gas_stoich[j]);
        }
        else
        {
            react = react * pow((*_coupled_gas[j])[_qp],fabs(_gas_stoich[j]));
        }
    }

    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        if (j != i)
        {
            cat = cat * pow((*_coupled_cat[j])[_qp],_cat_stoich[j]);
        }
    }

    return _ads_stoich[_ads_index]*_cat_stoich[i]*pow((*_coupled_cat[i])[_qp],(_cat_stoich[i]-1))*cat*((_forward*pow(CoupledConstChemisorption::computeSiteBalance(),_ads_sites[_ads_index])*react) - (_reverse*pow(_u[_qp],_ads_stoich[_ads_index])*prod));

}

Real CoupledCatalyst::computeRateFunctionAdsOffDiagJacobi(int i)
{
    double no_cat = 1.0, cat = 1.0;
    no_cat = CoupledConstChemisorption::computeRateFunctionAdsOffDiagJacobi(i);
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_cat_stoich[j]);
    }
    
    return no_cat * cat;
}

Real CoupledCatalyst::computeQpResidual()
{
    return _test[_i][_qp]*_u_dot[_qp] - _test[_i][_qp]*computeRateFunction();
}

Real CoupledCatalyst::computeQpJacobian()
{
    return _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp] - _test[_i][_qp]*computeRateFunctionJacobi();
}

Real CoupledCatalyst::computeQpOffDiagJacobian(unsigned int jvar)
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


