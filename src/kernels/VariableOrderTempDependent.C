/*!
 *  \file VariableOrderTempDependent.C
 *	\brief Standard kernel for coupling multiple gas and adsorbed species together via a reaction based mechanism
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
 *				dD/dt = x*(C)^a*[k_f*(L)^b*product(i, C_i^z_i) - k_r*(D)^f*product(j, C_j^z_j)]
 *
 *				Parameters are as follows...
 *
 *				v_i,j = stoichiometric coefficients for gas reactants/products
 *              z_i,j = reaction oder with respect to the gas reactants/products
 *				C_i,j = concentrations of gas reactants/products (mol/L)
 *				m = number of catalysts used in the reaction
 *				n = number of adsorption sites (L) reduced
 *              x = number of reduced adsortion sites (D) produced
 *				L = concentration of available sites (mol/kg)
 *				D = concentration of reduced species produced (mol/kg)
 *				k_f,r = rate constant for the forward/reverse reaction
 *              a = reaction order with repect to the catalyst
 *              b = reaction order with respect to available adsorption sites
 *              f = reaction order with respect to the reduced species
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
 *              The forward and reverse reaction rate constants are temperature dependent and set throught the
 *              Arrhenius equaiton ...
 *
 *              k = A*e^(-Ea/RT)
 *
 *              where k is the reaction rate constant, A is the pre-eponential factor, Ea is the activation energy
 *              (J/mol), R is the gas constant (J/K*mol), and T is the temperature (K).
 *
 *  \author Alexander Wiechert
 *	\date 8/30/2018
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

#include "VariableOrderTempDependent.h"
registerMooseObject("dgospreyApp", VariableOrderTempDependent);

template<>
InputParameters validParams<VariableOrderTempDependent>()
{
	InputParameters params = validParams<VariableOrderCoupledCatalyst>();

	params.addParam<Real>("forward_prefactor",0.0,"Pre-exponential factor for forward reaction");
	params.addParam<Real>("reverse_prefactor",0.0,"Pre-exponential factor for reverse reaction");
	params.addParam<Real>("forward_activation_energy",0.0,"Activation Energy for the forward reaciton (J/mol)");
	params.addParam<Real>("reverse_activation_energy",0.0,"Activation Energy for the reverse reaciton (J/mol)");
    params.addRequiredCoupledVar("coupled_temp","Name of the temperature variable being coupled");

	return params;
}

VariableOrderTempDependent::VariableOrderTempDependent(const InputParameters & parameters)
: VariableOrderCoupledCatalyst(parameters),
_for_a(getParam<Real>("forward_prefactor")),
_rev_a(getParam<Real>("reverse_prefactor")),
_for_activation(getParam<Real>("forward_activation_energy")),
_rev_activation(getParam<Real>("reverse_activation_energy")),
_coupled_temp(coupledValue("coupled_temp")),
_coupled_var_temp(coupled("coupled_temp"))
{

}

Real VariableOrderTempDependent::computeForwardRateConstant()
{
    return _for_a*pow(2.7182818284,-_for_activation/(8.3144598*_coupled_temp[_qp]));
}

Real VariableOrderTempDependent::computeReverseRateConstant()
{
    return _rev_a*pow(2.7182818284,-_rev_activation/(8.3144598*_coupled_temp[_qp]));
}

Real VariableOrderTempDependent::computeForwardTempDerivative()
{
    return (_for_activation/(8.3144598*pow(_coupled_temp[_qp],2)))*_for_a*pow(2.7182818284,-_for_activation/(8.3144598*_coupled_temp[_qp]));
}

Real VariableOrderTempDependent::computeReverseTempDerivative()
{
    return (_rev_activation/(8.3144598*pow(_coupled_temp[_qp],2)))*_rev_a*pow(2.7182818284,-_rev_activation/(8.3144598*_coupled_temp[_qp]));
}

Real VariableOrderTempDependent::computeRateFunction()
{
    _forward = computeForwardRateConstant();
    _reverse = computeReverseRateConstant();
    return VariableOrderCoupledCatalyst::computeRateFunction();
}

Real VariableOrderTempDependent::computeRateFunctionJacobi()
{
    _forward = computeForwardRateConstant();
    _reverse = computeReverseRateConstant();
    return VariableOrderCoupledCatalyst::computeRateFunctionJacobi();
}

Real VariableOrderTempDependent::computeRateFunctionGasOffDiagJacobi(int i)
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
    
    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_c_order[j]);
    }

	if (_gas_stoich[i] >= 0.0)
	{
        _reverse = computeReverseRateConstant();
		return -_ads_stoich[_ads_index]*cat*_reverse*pow(_u[_qp],_a_order[_ads_index])*prod*_g_order[i]*pow((*_coupled_gas[i])[_qp],_g_order[i]-1)*_phi[_j][_qp];
	}
	else
	{
        _forward = computeForwardRateConstant();
        return _ads_stoich[_ads_index]*cat*_forward*pow(CoupledConstChemisorption::computeSiteBalance(),_s_order[_ads_index])*react*_g_order[i]*pow((*_coupled_gas[i])[_qp],_g_order[i]-1)*_phi[_j][_qp];
	}
}

Real VariableOrderTempDependent::computeRateFunctionCatalystOffDiagJacobi(int i)
{
    _forward = computeForwardRateConstant();
    _reverse = computeReverseRateConstant();
    return VariableOrderCoupledCatalyst::computeRateFunctionCatalystOffDiagJacobi(i);
}

Real VariableOrderTempDependent::computeRateFunctionAdsOffDiagJacobi(int i)
{
    _forward = computeForwardRateConstant();
    return VariableOrderCoupledCatalyst::computeRateFunctionAdsOffDiagJacobi(i);
}

Real VariableOrderTempDependent::computeRateFunctionTempOffDiagJacobi()
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

    for (unsigned int j=0; j<_coupled_cat.size(); j++)
    {
        cat = cat * pow((*_coupled_cat[j])[_qp],_c_order[j]);
    }

return _ads_stoich[_ads_index]*cat*( (computeForwardTempDerivative()*pow(CoupledConstChemisorption::computeSiteBalance(),_s_order[_ads_index])*react) - (computeReverseTempDerivative()*pow(_u[_qp],_a_order[_ads_index])*prod) );
}

Real VariableOrderTempDependent::computeQpResidual()
{
    return _test[_i][_qp]*_u_dot[_qp] - _test[_i][_qp]*computeRateFunction();
}

Real VariableOrderTempDependent::computeQpJacobian()
{
    return _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp] - _test[_i][_qp]*computeRateFunctionJacobi();
}

Real VariableOrderTempDependent::computeQpOffDiagJacobian(unsigned int jvar)
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
    //Off-diagonals for temperature
    if (jvar == _coupled_var_temp)
    {
        return -_test[_i][_qp]*computeRateFunctionTempOffDiagJacobi();
    }
    
    return 0.0;
}

