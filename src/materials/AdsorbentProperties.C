/*!
 *  \file AdsorbentProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with the adsorbent
 *	\warning THIS KERNEL IS INCOMPLETE! ONLY USED FOR DATA STORAGE FOR PELLET DENSITY AND HEAT CAPACITY!
 *  \author Austin Ladshaw
 *	\date 11/20/2015
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

#include "AdsorbentProperties.h"
registerMooseObject("dgospreyApp", AdsorbentProperties);

template<>
// input parameters are the parameters that are constant and not calculated from other parameters
InputParameters validParams<AdsorbentProperties>()
{
  InputParameters params = validParams<Material>();
  
	params.addParam<Real>("binder_fraction",1.0,"Binder fraction of the adsorbent pellet");
	params.addParam<Real>("binder_porosity","Porosity of the binder material in the adsorbent pellet");
	params.addParam<Real>("crystal_radius",1.0,"Radius of the adsorbent crystals in the binder matrix (um)");
	params.addParam<Real>("pellet_diameter","Diameter of the adsorbent pellet (cm)");
	params.addParam<Real>("macropore_radius","Nominal pore size of the macropores in the binder material (cm)");
	params.addParam<Real>("pellet_density","Density of the adsorbent pellet (g/cm^3)");
	params.addParam<Real>("pellet_heat_capacity","Pellet heat capacity (J/g/K)");
	
	params.addParam<std::vector<Real> >("ref_diffusion","Reference Surface Diffusivity (um^2/hr)");
	params.addParam<std::vector<Real> >("activation_energy","Activation Energy of Surface Diffusion (J/mol)");
	params.addParam<std::vector<Real> >("ref_temperature","Reference Temperature for Surface Diffusion (K)");
	params.addParam<std::vector<Real> >("affinity","Affinity coefficient for Surface Diffusion (-)");
	
	params.addCoupledVar("temperature","Coupled variable for temperature");
	params.addCoupledVar("coupled_gases", "Gas concentrations variables being coupled");
  
  return params;
}

AdsorbentProperties::AdsorbentProperties(const InputParameters & parameters)
:Material(parameters),

_binder_fraction(getParam<Real>("binder_fraction")),
_eps_binder(getParam<Real>("binder_porosity")),
_crystal_rad(getParam<Real>("crystal_radius")),
_pellet_dia(getParam<Real>("pellet_diameter")),
_macropore_radius(getParam<Real>("macropore_radius")),
_rhos(getParam<Real>("pellet_density")),
_hs(getParam<Real>("pellet_heat_capacity")),

_ref_diff(getParam<std::vector<Real> >("ref_diffusion")),
_act_energy(getParam<std::vector<Real> >("activation_energy")),
_ref_temp(getParam<std::vector<Real> >("ref_temperature")),
_affinity(getParam<std::vector<Real> >("affinity")),

_ref_diffusion(declareProperty<std::vector<Real> >("ref_diffusion")),
_activation_energy(declareProperty<std::vector<Real> >("activation_energy")),
_ref_temperature(declareProperty<std::vector<Real> >("ref_temperature")),
_affinity_coeff(declareProperty<std::vector<Real> >("affinity_coeff")),

_pellet_density(declareProperty<Real>("pellet_density")),
_pellet_heat_capacity(declareProperty<Real>("pellet_heat_capacity")),
_pellet_diameter(declareProperty<Real>("pellet_diameter")),
_crystal_radius(declareProperty<Real>("crystal_radius")),
_binder_porosity(declareProperty<Real>("binder_porosity")),
_binder_ratio(declareProperty<Real>("binder_ratio")),
_pore_size(declareProperty<Real>("pore_size")),
_surface_diffusion(declareProperty<std::vector<Real> >("surface_diffusion")),
_surf_diff_old(getMaterialPropertyOld<std::vector<Real> >("surface_diffusion")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data")),
_temperature(coupledValue("temperature")),
_temperature_old(coupledValueOld("temperature"))

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
	/*
	 Note: When using _gas_conc[i], it must be appropriately dereferenced as follows...
	 (*_gas_conc[i])[_qp] = concentration of species i at node _qp
	 */
}

void AdsorbentProperties::initQpStatefulProperties()
{
	_surface_diffusion[_qp].resize(_gas_conc.size());
	_ref_diffusion[_qp].resize(_gas_conc.size());
	_activation_energy[_qp].resize(_gas_conc.size());
	_ref_temperature[_qp].resize(_gas_conc.size());
	_affinity_coeff[_qp].resize(_gas_conc.size());
	
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		_ref_diffusion[_qp][i] = _ref_diff[i];
		_activation_energy[_qp][i] = _act_energy[i];
		_ref_temperature[_qp][i] = _ref_temp[i];
		_affinity_coeff[_qp][i] = _affinity[i];
	}
	
	_pellet_density[_qp] = _rhos;
	_pellet_heat_capacity[_qp] = _hs;
	_pellet_diameter[_qp] = _pellet_dia;
	_crystal_radius[_qp] = _crystal_rad;
	_binder_porosity[_qp] = _eps_binder;
	_binder_ratio[_qp] = _binder_fraction;
	_pore_size[_qp] = _macropore_radius;
}

void AdsorbentProperties::computeQpProperties()
{
	
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		Real p = _magpie_dat[_qp].gpast_dat[i].y * _magpie_dat[_qp].sys_dat.PT;
		
		if (p < 0.0) p = 0.0;
		
		if (_affinity.size() != _gas_conc.size() || _ref_temp.size() != _gas_conc.size() || _act_energy.size() != _gas_conc.size() || _ref_diff.size() != _gas_conc.size() )
		{
			_surface_diffusion[_qp][i] = 0.0;
		}
		else
		{
			_surface_diffusion[_qp][i] = D_o(_ref_diff[i],_act_energy[i],_temperature_old[_qp]);
			_surface_diffusion[_qp][i] = D_inf(_surface_diffusion[_qp][i],_ref_temp[i],_affinity[i],p,_temperature_old[_qp]);
		}
		
		_ref_diffusion[_qp][i] = _ref_diff[i];
		_activation_energy[_qp][i] = _act_energy[i];
		_ref_temperature[_qp][i] = _ref_temp[i];
		_affinity_coeff[_qp][i] = _affinity[i];
	}
	
	//For constant bed properties...
	_pellet_density[_qp] = _rhos;
	_pellet_heat_capacity[_qp] = _hs;
	_pellet_diameter[_qp] = _pellet_dia;
	_crystal_radius[_qp] = _crystal_rad;
	_binder_porosity[_qp] = _eps_binder;
	_binder_ratio[_qp] = _binder_fraction;
	_pore_size[_qp] = _macropore_radius;
	//Note: some of these may vary with space or temperature or concentration, but for now we assume constant
}
