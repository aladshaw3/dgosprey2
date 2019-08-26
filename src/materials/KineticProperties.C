/*!
 *  \file KineticProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with SCOPSOWL simulations
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

#include "KineticProperties.h"
registerMooseObject("dgospreyApp", KineticProperties);

template<>
// input parameters are the parameters that are constant and not calculated from other parameters
InputParameters validParams<KineticProperties>()
{
	InputParameters params = validParams<Material>();
	
	params.addParam<bool>("dirichlet_bc",false,"Type of boundary condition to use in SCOPSOWL simulations");
	params.addParam<bool>("heterogeneous","Type of adsorbent configuration to use in SCOPSOWL simulations");
	params.addParam<bool>("surface_diffusion","Determines whether or not surface diffusion is included in SCOPSOWL");
	params.addParam<bool>("macro_spheres",true,"True = spherical pellets; False = cylindrical pellets");
	params.addParam<bool>("micro_spheres",true,"True = spherical crystals; False = cylindrical crystals");
	params.addParam<Real>("macro_length",1.0,"Length of the cylindrical pellets (cm)");
	params.addParam<Real>("micro_length",1.0,"Length of the cylindrical crystals (um)");
	params.addCoupledVar("coupled_adsorption", "Adsorption concentration variables being coupled");
	
	return params;
}

KineticProperties::KineticProperties(const InputParameters & parameters)
:Material(parameters),

DirichletBC(getParam<bool>("dirichlet_bc")),
Heterogeneous(getParam<bool>("heterogeneous")),
SurfaceDiffusion(getParam<bool>("surface_diffusion")),
MacroSpheres(getParam<bool>("macro_spheres")),
MicroSpheres(getParam<bool>("micro_spheres")),
MacroLength(getParam<Real>("macro_length")),
MicroLength(getParam<Real>("micro_length")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data")),
_coupled_gas(getMaterialProperty< MIXED_GAS >("mixed_gas_data")),
_velocity(getMaterialProperty<Real>("velocity")),
_pellet_diameter(getMaterialProperty<Real>("pellet_diameter")),
_pellet_density(getMaterialProperty<Real>("pellet_density")),
_binder_porosity(getMaterialProperty<Real>("binder_porosity")),
_pore_size(getMaterialProperty<Real>("pore_size")),
_crystal_radius(getMaterialProperty<Real>("crystal_radius")),
_binder_ratio(getMaterialProperty<Real>("binder_ratio")),
_ref_diffusion(getMaterialProperty<std::vector<Real> >("ref_diffusion")),
_activation_energy(getMaterialProperty<std::vector<Real> >("activation_energy")),
_ref_temperature(getMaterialProperty<std::vector<Real> >("ref_temperature")),
_affinity_coeff(getMaterialProperty<std::vector<Real> >("affinity_coeff")),
_owl_dat(declareProperty< SCOPSOWL_DATA >("owl_data")),
_owl_dat_old(getMaterialPropertyOld< SCOPSOWL_DATA >("owl_data")),
_gas_dat(declareProperty< MIXED_GAS >("gas_data")),
_gas_dat_old(getMaterialPropertyOld< MIXED_GAS >("gas_data"))

{
	unsigned int n = coupledComponents("coupled_adsorption");
	_solid_conc.resize(n);
	
	for (unsigned int i = 0; i<_solid_conc.size(); ++i)
	{
		_solid_conc[i] = &coupledValue("coupled_adsorption",i);
	}
}

void KineticProperties::initQpStatefulProperties()
{
	_owl_dat[_qp].Print2File = false;
	_owl_dat[_qp].NonLinear = true;
	_owl_dat[_qp].level = 1;
	_owl_dat[_qp].Print2Console = false;
	_gas_dat[_qp].CheckMolefractions = true;
	_owl_dat[_qp].total_steps = 0;
	
	_owl_dat[_qp].magpie_dat = _magpie_dat[_qp];
	
	_owl_dat[_qp].param_dat.resize(_magpie_dat[_qp].sys_dat.N);
	_owl_dat[_qp].finch_dat.resize(_magpie_dat[_qp].sys_dat.N);
	_owl_dat[_qp].y.resize(_magpie_dat[_qp].sys_dat.N);
	
	for (int i=0; i<_magpie_dat[_qp].sys_dat.N; i++)
	{
		if (_magpie_dat[_qp].gsta_dat[i].qmax <= 0.0)
			_owl_dat[_qp].param_dat[i].Adsorbable = false;
		else
			_owl_dat[_qp].param_dat[i].Adsorbable = true;
		
	}
	
	int success = 0;
	success = initialize_data(_owl_dat[_qp].magpie_dat.sys_dat.N, &_gas_dat[_qp]);
	if (success != 0) {mError(simulation_fail); return;}
	
	_owl_dat[_qp].DirichletBC = DirichletBC;
	_owl_dat[_qp].Heterogeneous = Heterogeneous;
	_owl_dat[_qp].SurfDiff = SurfaceDiffusion;
	
	if (MacroSpheres == true)
	{
		_owl_dat[_qp].coord_macro = 2;
		_owl_dat[_qp].char_macro = 1.0;
	}
	else
	{
		_owl_dat[_qp].coord_macro = 1;
		_owl_dat[_qp].char_macro = MacroLength;
	}
	
	if (MicroSpheres == true)
	{
		_owl_dat[_qp].coord_micro = 2;
		_owl_dat[_qp].char_micro = 1.0;
	}
	else
	{
		_owl_dat[_qp].coord_micro = 1;
		_owl_dat[_qp].char_micro = MicroLength;
	}
	
	if (_owl_dat[_qp].Heterogeneous == true && _owl_dat[_qp].SurfDiff == true)
	{
		_owl_dat[_qp].eval_surfDiff = (*default_Dc);
	}
	else
	{
		_owl_dat[_qp].eval_surfDiff = (*default_surf_diffusion);
		_owl_dat[_qp].coord_micro = 2;
		_owl_dat[_qp].char_micro = 1.0;
	}
	
	_owl_dat[_qp].pellet_radius = _pellet_diameter[_qp]/2.0;
	_owl_dat[_qp].pellet_density = _pellet_density[_qp];
	_owl_dat[_qp].binder_porosity = _binder_porosity[_qp];
	_owl_dat[_qp].binder_poresize = _pore_size[_qp];
	
	_gas_dat[_qp] = _coupled_gas[_qp];
	
	if (Heterogeneous == true)
	{
		_owl_dat[_qp].crystal_radius = _crystal_radius[_qp];
		_owl_dat[_qp].binder_fraction = _binder_ratio[_qp];
	}
	else
	{
		_owl_dat[_qp].binder_fraction = 1.0;
		_owl_dat[_qp].crystal_radius = 1.0;
	}
	
	success = setup_SCOPSOWL_DATA(NULL, default_adsorption, default_retardation, default_pore_diffusion, default_filmMassTransfer, _owl_dat[_qp].eval_surfDiff, (void *)&_owl_dat[_qp], &_gas_dat[_qp], &_owl_dat[_qp]);
	if (success != 0) {mError(simulation_fail); return;}
	
}

void KineticProperties::computeQpProperties()
{
	
	_owl_dat[_qp].total_pressure = _magpie_dat[_qp].sys_dat.PT;
	_owl_dat[_qp].gas_temperature = _magpie_dat[_qp].sys_dat.T;
	_owl_dat[_qp].gas_velocity = _velocity[_qp] / 3600.0;
	
	_owl_dat[_qp].pellet_radius = _pellet_diameter[_qp]/2.0;
	_owl_dat[_qp].pellet_density = _pellet_density[_qp];
	_owl_dat[_qp].binder_porosity = _binder_porosity[_qp];
	_owl_dat[_qp].binder_poresize = _pore_size[_qp];
	
	_gas_dat[_qp] = _coupled_gas[_qp];
	_owl_dat[_qp].magpie_dat = _magpie_dat[_qp];
	
	if (Heterogeneous == true)
	{
		_owl_dat[_qp].crystal_radius = _crystal_radius[_qp];
		_owl_dat[_qp].binder_fraction = _binder_ratio[_qp];
	}
	else
	{
		_owl_dat[_qp].binder_fraction = 1.0;
		_owl_dat[_qp].crystal_radius = 1.0;
	}
	
	for (int i=0; i<_owl_dat[_qp].magpie_dat.sys_dat.N; i++)
	{
		_owl_dat[_qp].param_dat[i].ref_diffusion = _ref_diffusion[_qp][i];
		_owl_dat[_qp].param_dat[i].activation_energy = _activation_energy[_qp][i];
		_owl_dat[_qp].param_dat[i].ref_temperature = _ref_temperature[_qp][i];
		_owl_dat[_qp].param_dat[i].affinity = _affinity_coeff[_qp][i];
		_owl_dat[_qp].y[i] = _magpie_dat[_qp].gpast_dat[i].y;
		
		if (_dt == 0.0)
			_owl_dat[_qp].param_dat[i].qIntegralAvg_old = (*_solid_conc[i])[_qp];
	}
	
}


