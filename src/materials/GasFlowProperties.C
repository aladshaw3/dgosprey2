/*!
 *  \file GasFlowProperties.h
 *	\brief Material Properties kernel that will setup and calculate gas flow properties based on physical characteristics
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

#include "GasFlowProperties.h"
registerMooseObject("dgospreyApp", GasFlowProperties);

template<>
// input parameters are the parameters that are constant and not calculated from other parameters
InputParameters validParams<GasFlowProperties>()
{
	InputParameters params = validParams<Material>();
	
	params.addParam<std::vector<Real> >("molecular_weight","Molecular wieghts of each gas phase component (g/mol)");
	params.addParam<std::vector<Real> >("comp_heat_capacity","Heat capacity of each gas phase component (J/g/K)");
	params.addParam<std::vector<Real> >("comp_ref_viscosity","Reference viscosities for each gas phase component (g/cm/s)");
	params.addParam<std::vector<Real> >("comp_ref_temp","Reference temperature for each gas phase component (K)");
	params.addParam<std::vector<Real> >("comp_Sutherland_const","Sutherland's constant for each gas phase component (K)");
	params.addParam<Real>("flow_rate","Volumetric flow rate in the system (cm^3/hr)");
	params.addCoupledVar("temperature","Coupled variable for temperature");
	params.addCoupledVar("total_pressure","Coupled variable for total pressure");
	params.addCoupledVar("coupled_gases", "Gas concentration variables being coupled");
	
	return params;
}

GasFlowProperties::GasFlowProperties(const InputParameters & parameters)
:Material(parameters),

_molecular_weight(getParam<std::vector<Real> >("molecular_weight")),
_comp_heat_capacity(getParam<std::vector<Real> >("comp_heat_capacity")),
_comp_ref_viscosity(getParam<std::vector<Real> >("comp_ref_viscosity")),
_comp_ref_temp(getParam<std::vector<Real> >("comp_ref_temp")),
_comp_Sutherland_const(getParam<std::vector<Real> >("comp_Sutherland_const")),
_flow_rate(getParam<Real>("flow_rate")),

_velocity(declareProperty<Real>("velocity")),
_gas_density(declareProperty<Real>("gas_density")),
_gas_viscosity(declareProperty<Real>("gas_viscosity")),
_gas_heat_capacity(declareProperty<Real>("gas_heat_capacity")),
_gas_molecular_wieght(declareProperty<Real>("gas_molecular_wieght")),
_inner_dia(getMaterialProperty<Real>("inner_dia")),
_porosity(getMaterialProperty<Real>("porosity")),
_pellet_density(getMaterialProperty<Real>("pellet_density")),
_pellet_heat_capacity(getMaterialProperty<Real>("pellet_heat_capacity")),
_pellet_diameter(getMaterialProperty<Real>("pellet_diameter")),
_binder_porosity(getMaterialProperty<Real>("binder_porosity")),
_pore_size(getMaterialProperty<Real>("pore_size")),

_heat_retardation(declareProperty<Real>("heat_retardation")),
_conductivity(declareProperty<Real>("conductivity")),
_molecular_diffusion(declareProperty<std::vector<Real> >("molecular_diffusion")),
_dispersion(declareProperty<std::vector<Real> >("dispersion")),
_retardation(declareProperty<Real>("retardation")),
_mixed_gas(declareProperty< MIXED_GAS >("mixed_gas_data")),
_mixed_gas_old(getMaterialPropertyOld< MIXED_GAS >("mixed_gas_data")),
_film_transfer(declareProperty<std::vector<Real> >("film_transfer")),
_pore_diffusion(declareProperty<std::vector<Real> >("pore_diffusion")),
_temperature(coupledValue("temperature")),
_temperature_old(coupledValueOld("temperature")),
_total_pressure(coupledValue("total_pressure"))

{
	unsigned int n = coupledComponents("coupled_gases");
	_index.resize(n);
	_gas_conc.resize(n);
	_gas_conc_old.resize(n);
	
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		_index[i] = coupled("coupled_gases",i); //may only be useful for compute Jacobian Off Diagonal (see ~/projects/truck/moose/modules/chemical_reactions/ .../ CoupledConvectionReactionSub.C)
		_gas_conc[i] = &coupledValue("coupled_gases",i);
		_gas_conc_old[i] = &coupledValueOld("coupled_gases",i);
	}
	
	/*
	 Note: When using _gas_conc[i], it must be appropriately dereferenced as follows...
	 (*_gas_conc[i])[_qp] = concentration of species i at node _qp
	 
	 Opposite for materials
	 _molecular_diffusion[_qp][i] = Dm[_qp][i]
	 */
}

void GasFlowProperties::initQpStatefulProperties()
{
	int success = 0;
	int num_species = (int)_gas_conc.size();
	_mixed_gas[_qp].CheckMolefractions = false;
	success = initialize_data(num_species,&_mixed_gas[_qp]);
	
	_molecular_diffusion[_qp].resize(_gas_conc.size());
	_dispersion[_qp].resize(_gas_conc.size());
	_film_transfer[_qp].resize(_gas_conc.size());
	_pore_diffusion[_qp].resize(_gas_conc.size());
	_molefrac.resize(_gas_conc.size());
}

void GasFlowProperties::computeQpProperties()
{
	//Check to see if MIXED_GAS has been initialized
	int success = 0;
	
	
	_mixed_gas[_qp].CheckMolefractions = false;
	_velocity[_qp] = _porosity[_qp] * (_flow_rate / (_porosity[_qp] * (M_PI/4.0) * _inner_dia[_qp] * _inner_dia[_qp]));
	_retardation[_qp] = _porosity[_qp];
	
	_gas_molecular_wieght[_qp] = 0.0;
	_gas_viscosity[_qp] = 0.0;
	_gas_heat_capacity[_qp] = 0.0;
	_heat_retardation[_qp] = 0.0;
	
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		
		Real _yi = ( ((*_gas_conc_old[i])[_qp] * 8.3144621 * _temperature_old[_qp])/_total_pressure[_qp] );
		
		if (_yi > 1.0) _yi = 1.0;
		if (_yi < 0.0) _yi = 0.0;
		
		_molefrac[i] = _yi;
		_mixed_gas[_qp].species_dat[i].molecular_weight = _molecular_weight[i];
		_mixed_gas[_qp].species_dat[i].Sutherland_Temp = _comp_ref_temp[i];
		_mixed_gas[_qp].species_dat[i].Sutherland_Const = _comp_Sutherland_const[i];
		_mixed_gas[_qp].species_dat[i].Sutherland_Viscosity = _comp_ref_viscosity[i];
		_mixed_gas[_qp].species_dat[i].specific_heat= _comp_heat_capacity[i];
		
		Real _sum_yj_over_Dij = 0.0, _sum_yi_over_Dij_prime = 0.0;
		
		Real _rho_i = (_total_pressure[_qp] * _molecular_weight[i]) / (8.3144621 * _temperature_old[_qp]);
		
		Real _mu_i = _comp_ref_viscosity[i] * ( ( _comp_ref_temp[i] + _comp_Sutherland_const[i]) / ( _temperature_old[_qp] + _comp_Sutherland_const[i]) ) * std::pow(_temperature_old[_qp]/_comp_ref_temp[i],(3.0/2.0));
		
		Real _phi = 0.873143 + (0.000072375 * _temperature_old[_qp]);
		
		for (unsigned int j = 0; j<_gas_conc.size(); j++)
		{
			if ( j != i)
			{
				Real _yj = ( ((*_gas_conc_old[j])[_qp] * 8.3144621 * _temperature_old[_qp])/_total_pressure[_qp] );
				
				if (_yj > 1.0) _yj = 1.0;
				if (_yj < 0.0) _yj = 0.0;
				
				Real _rho_j = (_total_pressure[_qp] * _molecular_weight[j]) / (8.3144621 * _temperature_old[_qp]);
		  
				Real _mu_j = _comp_ref_viscosity[j] * ( ( _comp_ref_temp[j] + _comp_Sutherland_const[j]) / ( _temperature_old[_qp] + _comp_Sutherland_const[j]) ) * std::pow(_temperature_old[_qp]/_comp_ref_temp[j],(3.0/2.0));
		  
				Real _Dij = ( (4.0/std::pow(2.0,0.5))*std::pow((1.0/_molecular_weight[i])+(1.0/_molecular_weight[j]),0.5) );
				_Dij = _Dij / std::pow(std::pow((_rho_i*_rho_i)/(1.385*1.385*_mu_i*_mu_i*_molecular_weight[i]),0.25)+std::pow((_rho_j*_rho_j)/(1.385*1.385*_mu_j*_mu_j*_molecular_weight[j]),0.25),2.0);
		  
				Real _Dij_prime = (_total_pressure[_qp]*_Dij) / 100.0;
				
				_sum_yj_over_Dij = _sum_yj_over_Dij + (_yj/_Dij);
				_sum_yi_over_Dij_prime = _sum_yi_over_Dij_prime + (_yi/_Dij_prime);
				
			}
			
		} //jth Loop
		
		if (_sum_yj_over_Dij == 0.0 || _yi >= 1.0)
		{
			_molecular_diffusion[_qp][i] = ( (4.0/std::pow(2.0,0.5))*std::pow((1.0/_molecular_weight[i])+(1.0/_molecular_weight[i]),0.5) );
			_molecular_diffusion[_qp][i] = _molecular_diffusion[_qp][i] / std::pow(std::pow((_rho_i*_rho_i)/(1.385*1.385*_mu_i*_mu_i*_molecular_weight[i]),0.25)+std::pow((_rho_i*_rho_i)/(1.385*1.385*_mu_i*_mu_i*_molecular_weight[i]),0.25),2.0);
		}
		else
		{
			_molecular_diffusion[_qp][i] = (1.0 - _yi) / _sum_yj_over_Dij;
		}
		
		if (_yi != 0.0)
			_gas_viscosity[_qp] = _gas_viscosity[_qp] + (_mu_i / (1.0 + ( (113.65*_phi*_mu_i*_temperature_old[_qp])/(_yi*_molecular_weight[i]) ) * _sum_yi_over_Dij_prime) );
		else
			_gas_viscosity[_qp] = _gas_viscosity[_qp] + 0.0;
		
		_gas_molecular_wieght[_qp] = _gas_molecular_wieght[_qp] + (_yi * _molecular_weight[i]);
		_gas_heat_capacity[_qp] = _gas_heat_capacity[_qp] + (_yi * _comp_heat_capacity[i]);
		
	} //ith Loop
	
	//Set variables for EGRET and calculate properties
	success = set_variables(_total_pressure[_qp],_temperature_old[_qp],fabs(_velocity[_qp]/3600.0),_pellet_diameter[_qp],_molefrac,&_mixed_gas[_qp]);
	success = calculate_properties(&_mixed_gas[_qp]);
	
	_gas_density[_qp] = (_total_pressure[_qp] * _gas_molecular_wieght[_qp]) / (8.3144621 * _temperature_old[_qp]);
	
	double kin_vis = _gas_viscosity[_qp] / _gas_density[_qp] * 3600.0;
	double ReNum = fabs(_velocity[_qp]) * _inner_dia[_qp] / kin_vis;
	double pre_coef = fabs(_velocity[_qp]) * _pellet_diameter[_qp] / _porosity[_qp];
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		_molecular_diffusion[_qp][i] = _mixed_gas[_qp].species_dat[i].molecular_diffusion * 3600.0;
		double ScNum = kin_vis / _molecular_diffusion[_qp][i] / _porosity[_qp];
		
		if (fabs(ReNum) < sqrt(DBL_EPSILON))
			_dispersion[_qp][i] = 20.0 * _pellet_diameter[_qp] * _molecular_diffusion[_qp][i] / _porosity[_qp] / _inner_dia[_qp];
		else
			_dispersion[_qp][i] = (1.0-exp(-4.02e-8*std::pow(_flow_rate,1.8)))*pre_coef * ((20.0/ReNum/ScNum) + 0.5) + (1.0-exp(-4.02e-8*std::pow(_flow_rate,1.8)))*fabs(_velocity[_qp])*_inner_dia[_qp];
	}
	
	double gas_cond = 4.21756e-5 + (7.19974e-7 * _temperature_old[_qp]); //J/s/cm/K
	double PrNum = _gas_heat_capacity[_qp] * _gas_viscosity[_qp] / gas_cond;
	double PeNum = ((0.73 * _porosity[_qp]) / (ReNum * PrNum)) + (0.5/(1.0+((9.7 * _porosity[_qp]) / (ReNum * PrNum))));
	_conductivity[_qp] = 40.0* PeNum * _gas_heat_capacity[_qp] * _gas_density[_qp] * fabs(_velocity[_qp]) * _pellet_diameter[_qp];
	_heat_retardation[_qp] = ((_gas_heat_capacity[_qp]*_gas_density[_qp]*_porosity[_qp]) + (_pellet_heat_capacity[_qp]*_pellet_density[_qp]*(1.0-_porosity[_qp])));
	
	// Calculate Film Mass Transfer and Pore Diffusion
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		_film_transfer[_qp][i] = FilmMTCoeff(_mixed_gas[_qp].species_dat[i].molecular_diffusion,_pellet_diameter[_qp],_mixed_gas[_qp].Reynolds,_mixed_gas[_qp].species_dat[i].Schmidt) * 3600.0;
		_pore_diffusion[_qp][i] = Dp(_mixed_gas[_qp].species_dat[i].molecular_diffusion,_binder_porosity[_qp]);
		_pore_diffusion[_qp][i] = avgDp(_pore_diffusion[_qp][i],Dk(_pore_size[_qp],_temperature_old[_qp],_mixed_gas[_qp].species_dat[i].molecular_weight)) * 3600.0;
				
	}
	
}
