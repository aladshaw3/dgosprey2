/*!
 *  \file ThermodynamicProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with MAGPIE simulations
 *  \author Austin Ladshaw
 *	\date 04/28/2017
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

#include "ThermodynamicProperties.h"
registerMooseObject("dgospreyApp", ThermodynamicProperties);

template<>
// input parameters are the parameters that are constant and not calculated from other parameters
InputParameters validParams<ThermodynamicProperties>()
{
	InputParameters params = validParams<Material>();
	
	params.addCoupledVar("temperature","Coupled variable for temperature");
	params.addCoupledVar("total_pressure","Coupled variable for total pressure");
	params.addCoupledVar("coupled_gases", "Gas concentrations variables being coupled");
	
	params.addParam<std::vector<int> >("number_sites","The number of adsorption sites for each species in GSTA isotherm");
	params.addParam<std::vector<Real> >("maximum_capacity","The adsorption capacity for each species in the GSTA isotherm (mol/kg)");
	params.addParam<std::vector<Real> >("molar_volume","The van der Waals molar volume of each species (mol/cm^3)");
	params.addParam<std::vector<Real> >("enthalpy_site_1","The molar enthalpy of the species on site 1");
	params.addParam<std::vector<Real> >("enthalpy_site_2","The molar enthalpy of the species on site 2");
	params.addParam<std::vector<Real> >("enthalpy_site_3","The molar enthalpy of the species on site 3");
	params.addParam<std::vector<Real> >("enthalpy_site_4","The molar enthalpy of the species on site 4");
	params.addParam<std::vector<Real> >("enthalpy_site_5","The molar enthalpy of the species on site 5");
	params.addParam<std::vector<Real> >("enthalpy_site_6","The molar enthalpy of the species on site 6");
	
	params.addParam<std::vector<Real> >("entropy_site_1","The molar entropy of the species on site 1");
	params.addParam<std::vector<Real> >("entropy_site_2","The molar entropy of the species on site 2");
	params.addParam<std::vector<Real> >("entropy_site_3","The molar entropy of the species on site 3");
	params.addParam<std::vector<Real> >("entropy_site_4","The molar entropy of the species on site 4");
	params.addParam<std::vector<Real> >("entropy_site_5","The molar entropy of the species on site 5");
	params.addParam<std::vector<Real> >("entropy_site_6","The molar entropy of the species on site 6");
	
	
	return params;
}

ThermodynamicProperties::ThermodynamicProperties(const InputParameters & parameters)
:Material(parameters),

_temperature(coupledValue("temperature")),
_temperature_old(coupledValueOld("temperature")),
_total_pressure(coupledValue("total_pressure")),
_num_sites(getParam<std::vector<int> >("number_sites")),
_max_capacity(getParam<std::vector<Real> >("maximum_capacity")),
_molar_volume(getParam<std::vector<Real> >("molar_volume")),
_enthalpy_1(getParam<std::vector<Real> >("enthalpy_site_1")),
_enthalpy_2(getParam<std::vector<Real> >("enthalpy_site_2")),
_enthalpy_3(getParam<std::vector<Real> >("enthalpy_site_3")),
_enthalpy_4(getParam<std::vector<Real> >("enthalpy_site_4")),
_enthalpy_5(getParam<std::vector<Real> >("enthalpy_site_5")),
_enthalpy_6(getParam<std::vector<Real> >("enthalpy_site_6")),

_entropy_1(getParam<std::vector<Real> >("entropy_site_1")),
_entropy_2(getParam<std::vector<Real> >("entropy_site_2")),
_entropy_3(getParam<std::vector<Real> >("entropy_site_3")),
_entropy_4(getParam<std::vector<Real> >("entropy_site_4")),
_entropy_5(getParam<std::vector<Real> >("entropy_site_5")),
_entropy_6(getParam<std::vector<Real> >("entropy_site_6")),

_magpie_dat(declareProperty< MAGPIE_DATA >("magpie_data")),
_magpie_dat_old(getMaterialPropertyOld< MAGPIE_DATA >("magpie_data")),
_ads_heat(declareProperty<std::vector<Real> >("adsorption_heat"))
{
	unsigned int n = coupledComponents("coupled_gases");
	_index.resize(n);
	_gas_conc.resize(n);
	_gas_conc_old.resize(n);
	
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		_index[i] = coupled("coupled_gases",i);
		_gas_conc[i] = &coupledValue("coupled_gases",i);
		_gas_conc_old[i] = &coupledValueOld("coupled_gases",i);
	}
}

void ThermodynamicProperties::initQpStatefulProperties()
{
	_magpie_dat[_qp].sys_dat.N = _gas_conc.size();
	_magpie_dat[_qp].gpast_dat.resize(_magpie_dat[_qp].sys_dat.N);
	_magpie_dat[_qp].gsta_dat.resize(_magpie_dat[_qp].sys_dat.N);
	_magpie_dat[_qp].mspd_dat.resize(_magpie_dat[_qp].sys_dat.N);
	_ads_heat[_qp].resize(_gas_conc.size());
	
	for (int i=0; i<_magpie_dat[_qp].sys_dat.N; i++)
	{
		_magpie_dat[_qp].mspd_dat[i].eta.resize(_magpie_dat[_qp].sys_dat.N);
		_magpie_dat[_qp].gpast_dat[i].gama_inf.resize(_magpie_dat[_qp].sys_dat.N);
		_magpie_dat[_qp].gpast_dat[i].po.resize(_magpie_dat[_qp].sys_dat.N);
		_magpie_dat[_qp].gsta_dat[i].qmax = _max_capacity[i];
		
		//This species is adsorbable
		if (_magpie_dat[_qp].gsta_dat[i].qmax > 0.0)
		{
			_magpie_dat[_qp].mspd_dat[i].v = _molar_volume[i];
			_magpie_dat[_qp].gsta_dat[i].m = _num_sites[i];
			if (_magpie_dat[_qp].gsta_dat[i].m < 1)
				_magpie_dat[_qp].gsta_dat[i].m = 1;
			if (_magpie_dat[_qp].gsta_dat[i].m > 6)
				_magpie_dat[_qp].gsta_dat[i].m = 6;
			_magpie_dat[_qp].gsta_dat[i].dHo.resize(_magpie_dat[_qp].gsta_dat[i].m);
			_magpie_dat[_qp].gsta_dat[i].dSo.resize(_magpie_dat[_qp].gsta_dat[i].m);
			for (int n=0; n<_magpie_dat[_qp].gsta_dat[i].m; n++)
			{
				if (n == 0)
				{
					_magpie_dat[_qp].gsta_dat[i].dHo[n] = _enthalpy_1[i];
					_magpie_dat[_qp].gsta_dat[i].dSo[n] = _entropy_1[i];
				}
				else if (n == 1)
				{
					_magpie_dat[_qp].gsta_dat[i].dHo[n] = _enthalpy_2[i];
					_magpie_dat[_qp].gsta_dat[i].dSo[n] = _entropy_2[i];
				}
				else if (n == 2)
				{
					_magpie_dat[_qp].gsta_dat[i].dHo[n] = _enthalpy_3[i];
					_magpie_dat[_qp].gsta_dat[i].dSo[n] = _entropy_3[i];
				}
				else if (n == 3)
				{
					_magpie_dat[_qp].gsta_dat[i].dHo[n] = _enthalpy_4[i];
					_magpie_dat[_qp].gsta_dat[i].dSo[n] = _entropy_4[i];
				}
				else if (n == 4)
				{
					_magpie_dat[_qp].gsta_dat[i].dHo[n] = _enthalpy_5[i];
					_magpie_dat[_qp].gsta_dat[i].dSo[n] = _entropy_5[i];
				}
				else if (n == 5)
				{
					_magpie_dat[_qp].gsta_dat[i].dHo[n] = _enthalpy_6[i];
					_magpie_dat[_qp].gsta_dat[i].dSo[n] = _entropy_6[i];
				}
				else
				{
					_magpie_dat[_qp].gsta_dat[i].dHo[n] = 0.0;
					_magpie_dat[_qp].gsta_dat[i].dSo[n] = 0.0;
				}
			}
			
			_ads_heat[_qp][i] = -_enthalpy_1[i];
		}
		//This species will not adsorb
		else
		{
			_magpie_dat[_qp].gsta_dat[i].qmax = 0.0;
			_magpie_dat[_qp].mspd_dat[i].v = 0;
			_magpie_dat[_qp].gsta_dat[i].m = 1;
			_magpie_dat[_qp].gsta_dat[i].dHo.resize(_magpie_dat[_qp].gsta_dat[i].m);
			_magpie_dat[_qp].gsta_dat[i].dSo.resize(_magpie_dat[_qp].gsta_dat[i].m);
			for (int n=0; n<_magpie_dat[_qp].gsta_dat[i].m; n++)
			{
				_magpie_dat[_qp].gsta_dat[i].dHo[n] = 0.0;
				_magpie_dat[_qp].gsta_dat[i].dSo[n] = 0.0;
			}
			_ads_heat[_qp][i] = 0.0;
		}
	}
}

void ThermodynamicProperties::computeQpProperties()
{
	_magpie_dat[_qp].sys_dat.total_eval = 0;
	_magpie_dat[_qp].sys_dat.avg_norm = 0;
	_magpie_dat[_qp].sys_dat.max_norm = 0;
	_magpie_dat[_qp].sys_dat.Recover = false;
	_magpie_dat[_qp].sys_dat.Carrier = false;
	_magpie_dat[_qp].sys_dat.Ideal = false;
	_magpie_dat[_qp].sys_dat.Output = false;
	
	_magpie_dat[_qp].sys_dat.PT = _total_pressure[_qp];
	_magpie_dat[_qp].sys_dat.T = _temperature_old[_qp];
	
	double tempPT = 0.0;
	
	//Loop over all gas species
	for (int i=0; i<_magpie_dat[_qp].sys_dat.N; i++)
	{
		double pi = 0.0;
		pi = (*_gas_conc_old[i])[_qp] * 8.3144621 * _temperature_old[_qp];
		if (pi < 0.0)
			pi = 0.0;
		tempPT = pi + tempPT;
	}
	_magpie_dat[_qp].sys_dat.PT = tempPT;
	
	for (int i=0; i<_magpie_dat[_qp].sys_dat.N; i++)
	{
		double pi = 0.0;
		pi = (*_gas_conc_old[i])[_qp] * 8.3144621 * _temperature_old[_qp];
		if (pi < 0.0)
			pi = 0.0;
		
		if (_magpie_dat[_qp].gsta_dat[i].qmax == 0.0)
		{
			_magpie_dat[_qp].sys_dat.Carrier = true;
			_magpie_dat[_qp].gpast_dat[i].y = 0.0;
		}
		else
		{
			_magpie_dat[_qp].gpast_dat[i].y = pi / tempPT;
		}
		
		if (_magpie_dat[_qp].gpast_dat[i].y < 0.0)
		{
			_magpie_dat[_qp].sys_dat.Carrier = true;
			_magpie_dat[_qp].gpast_dat[i].y = 0.0;
		}
		
		_ads_heat[_qp][i] = -_enthalpy_1[i];
	}
	
}

