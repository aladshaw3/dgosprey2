/*!
 *  \file CoupledConstChemisorption.h
 *	\brief Standard kernel for coupling multiple gas and adsorbed species together via a reaction based mechanism
 *	\details This file creates a standard MOOSE kernel for the coupling multiple gas and adsorption species in a
 *				simulation based on the following reaction scheme...
 *
 *				sum(i, v_i*C_i) + m*L  <-- --> n*q + sum(j, v_j*C_j)
 *
 *				In this reaction scheme, the i-th species may interact with m number of available surface sites (L)
 *				to form n number of adsorbed species (q). In return, some other gas species may be produced (C_j) as
 *				by-products from the reaction. This expression then formulates the following rate equation for the
 *				adsorbed species (q)...
 *
 *				dq/dt = n*k_f*(L)^m*product(i, C_i^v_i) - n*k_r*(q)^n*product(j, C_j^v_j)
 *
 *				Parameters are as follows...
 *
 *				v_i,j = stoichiometric coefficients for gas reactants/products
 *				C_i,j = concentrations of gas reactants/products (mol/L)
 *				m = number of adsorption sites used in the reaction
 *				n = number of adsorption product species (q) produced
 *				L = concentration of available sites (mol/kg)
 *				q = concentration of adsorbed species produced (mol/kg)
 *				k_f,r = rate constant for the forward/reverse reaction
 *
 *				Rate expression must be coupled with all involved adsorbed species and uses a site-balance to
 *				account for the loss of adsorption sites during multi-species adsorption
 *
 *				SiteBalance: (L) = Lmax - sum(i, m_i/n_i*q_i)
 *
 *				where Lmax is the maximum capacity for adsorption (or the maximum available sites), m_i is the
 *				number of sites the m-th adsorbate occupies, n_i is the number of adsorbed species produced from
 *				the reaction consuming m_i sites, and q_i is the adsorbed concentration of the i-th species.
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 02/05/2018
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

#include "CoupledConstChemisorption.h"
registerMooseObject("dgospreyApp", CoupledConstChemisorption);

template<>
InputParameters validParams<CoupledConstChemisorption>()
{
	InputParameters params = validParams<TimeDerivative>();
	params.addParam<Real>("max_capacity",0.0,"Maximum Capacity for the Rate Function (mol/kg)");
	params.addParam<Real>("forward_rate",0.0,"Forward reaction rate constant (per hour)");
	params.addParam<Real>("reverse_rate",0.0,"Reverse reaction rate constant (per hour)");
	
	params.addParam< std::vector<Real> >("gases_stoichiometry","Stoichiometric coefficients for the gas species: (-) sign = reactants and (+) sign for products");
	params.addParam< std::vector<Real> >("adsorbed_sites","Adsorption site coefficients for the number of sites used in each reaction");
	params.addParam< std::vector<Real> >("adsorbed_stoichiometry","Adsorption stoichiometric coefficients for the number of adsorbed species produced in each reaction");
	
	params.addRequiredCoupledVar("coupled_gases","List of names of the gas variables being coupled");
	params.addRequiredCoupledVar("coupled_adsorption","List of names of the adsorbed variables being coupled");
	params.addRequiredCoupledVar("main_variable","Name of the non-linear variable that this kernel acts on");
	
	return params;
}

CoupledConstChemisorption::CoupledConstChemisorption(const InputParameters & parameters)
: TimeDerivative(parameters),
_maxcap(getParam<Real>("max_capacity")),
_forward(getParam<Real>("forward_rate")),
_reverse(getParam<Real>("reverse_rate")),

_gas_stoich(getParam<std::vector<Real> >("gases_stoichiometry")),
_ads_sites(getParam<std::vector<Real> >("adsorbed_sites")),
_ads_stoich(getParam<std::vector<Real> >("adsorbed_stoichiometry")),
_coupled_var_i(coupled("main_variable"))
{
	unsigned int gas_n = coupledComponents("coupled_gases");
	_coupled_gas_vars.resize(gas_n);
	_coupled_gas.resize(gas_n);
	
	for (unsigned int i = 0; i<_coupled_gas.size(); ++i)
	{
		_coupled_gas_vars[i] = coupled("coupled_gases",i);
		_coupled_gas[i] = &coupledValue("coupled_gases",i);
	}
	
	unsigned int ads_n = coupledComponents("coupled_adsorption");
	_coupled_ads_vars.resize(ads_n);
	_coupled_ads.resize(ads_n);
	
	for (unsigned int i = 0; i<_coupled_ads.size(); ++i)
	{
		_coupled_ads_vars[i] = coupled("coupled_adsorption",i);
		_coupled_ads[i] = &coupledValue("coupled_adsorption",i);
		
		if (_coupled_ads_vars[i] == _coupled_var_i)
			_ads_index = i;
	}
	
	if (_coupled_gas.size() != _gas_stoich.size())
		Moose::out << "ERROR!!! Vectors for coupled gases and gas stoichiometry do not match in size!\n\n";
	
	if (_coupled_ads.size() != _ads_sites.size())
		Moose::out << "ERROR!!! Vectors for coupled adsorption and adsorbed sites do not match in size!\n\n";
	
	if (_coupled_ads.size() != _ads_stoich.size())
		Moose::out << "ERROR!!! Vectors for coupled adsorption and adsorption stoichiometry do not match in size!\n\n";
}

Real CoupledConstChemisorption::computeRateFunction()
{
	double react = 1.0, prod = 1.0;
	for (unsigned int i=0; i<_coupled_gas.size(); i++)
	{
		if (_gas_stoich[i] >= 0.0)
		{
			prod = prod * pow((*_coupled_gas[i])[_qp],_gas_stoich[i]);
		}
		else
		{
			react = react * pow((*_coupled_gas[i])[_qp],fabs(_gas_stoich[i]));
		}
	}
	return _ads_stoich[_ads_index]*( (_forward*pow(computeSiteBalance(),_ads_sites[_ads_index])*react) - (_reverse*pow(_u[_qp],_ads_stoich[_ads_index])*prod) );
}

Real CoupledConstChemisorption::computeRateFunctionJacobi()
{
	double react = 1.0, prod = 1.0;
	for (unsigned int i=0; i<_coupled_gas.size(); i++)
	{
		if (_gas_stoich[i] >= 0.0)
		{
			prod = prod * pow((*_coupled_gas[i])[_qp],_gas_stoich[i]);
		}
		else
		{
			react = react * pow((*_coupled_gas[i])[_qp],fabs(_gas_stoich[i]));
		}
	}
	return _ads_stoich[_ads_index]*_ads_sites[_ads_index]*_forward*pow(computeSiteBalance(),_ads_sites[_ads_index]-1)*react*computeSiteBalanceOffDiagJacobi(_ads_index) - _ads_stoich[_ads_index]*_ads_stoich[_ads_index]*_reverse*pow(_u[_qp],_ads_stoich[_ads_index]-1)*prod*_phi[_j][_qp];
}

Real CoupledConstChemisorption::computeRateFunctionGasOffDiagJacobi(int i)
{
	//Is index i a reactant or product?
	//If _gas_stoich[i] > 0, prod; else, react;
	double react = 1.0, prod = 1.0;
	for (unsigned int j=0; j<_coupled_gas.size(); j++)
	{
		if (j != i)
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
	}
	if (_gas_stoich[i] >= 0.0)
	{
		return -_ads_stoich[_ads_index]*_reverse*pow(_u[_qp],_ads_stoich[_ads_index])*prod*_gas_stoich[i]*pow((*_coupled_gas[i])[_qp],_gas_stoich[i]-1)*_phi[_j][_qp];
	}
	else
	{
		return _ads_stoich[_ads_index]*_forward*pow(computeSiteBalance(),_ads_sites[_ads_index])*react*fabs(_gas_stoich[i])*pow((*_coupled_gas[i])[_qp],fabs(_gas_stoich[i])-1)*_phi[_j][_qp];
	}
}

Real CoupledConstChemisorption::computeRateFunctionAdsOffDiagJacobi(int i)
{
	double react = 1.0, prod = 1.0;
	for (unsigned int i=0; i<_coupled_gas.size(); i++)
	{
		if (_gas_stoich[i] >= 0.0)
		{
			prod = prod * pow((*_coupled_gas[i])[_qp],_gas_stoich[i]);
		}
		else
		{
			react = react * pow((*_coupled_gas[i])[_qp],fabs(_gas_stoich[i]));
		}
	}
	return _ads_stoich[_ads_index]*_ads_sites[_ads_index]*_forward*pow(computeSiteBalance(),_ads_sites[_ads_index]-1)*react*computeSiteBalanceOffDiagJacobi(i);
}

Real CoupledConstChemisorption::computeSiteBalance()
{
	double sum = 0.0;
	for (unsigned int i=0; i<_coupled_ads.size(); i++)
	{
		sum += (_ads_sites[i]/_ads_stoich[i])*(*_coupled_ads[i])[_qp];
	}
	return _maxcap - sum;
}

Real CoupledConstChemisorption::computeSiteBalanceOffDiagJacobi(int i)
{
	return -(_ads_sites[i]/_ads_stoich[i])*_phi[_j][_qp];
}

Real CoupledConstChemisorption::computeQpResidual()
{
	return _test[_i][_qp]*_u_dot[_qp] - _test[_i][_qp]*computeRateFunction();
}

Real CoupledConstChemisorption::computeQpJacobian()
{
	return _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp] - _test[_i][_qp]*computeRateFunctionJacobi();
}

Real CoupledConstChemisorption::computeQpOffDiagJacobian(unsigned int jvar)
{
	//Off-diagonals for gas concentrations
	for (unsigned int i = 0; i<_coupled_gas.size(); ++i)
	{
		if (jvar == _coupled_gas_vars[i] && jvar != _coupled_var_i)
		{
			return -_test[_i][_qp]*computeRateFunctionGasOffDiagJacobi(i);
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


