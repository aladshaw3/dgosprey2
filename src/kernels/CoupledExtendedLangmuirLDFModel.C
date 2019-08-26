/*!
 *  \file CoupledExtendedLangmuirLDFModel.h
 *	\brief Standard kernel for coupling non-linear variables via the extended langmuir model with LDF kinetics
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via the extended langmuir model, and applies linear driving force kinetics for the rate
 *			of adsorption.
 *
 *			This kernel extends the CoupledExtendedLangmuirModel kernel by calculating the model parameters from
 *			information in the input file for enthalpies and entropies. In addition, it calculates the Linear Driving Force
 *			parameter by estimating the overall LDF rate coefficient from material properties.
 *			The langmuir parameter (K) (described below) are to be estimated from the site enthalpies (dH)
 *			and entropies (dS) using the van't Hoff expression (shown below). Thus, this model
 *			is inherently a function of temperature and will require a different form of coupling
 *			with the temperature parameter.
 *
 *			In addition, the linear driving force parameter (k) is estimated using the Resistance-in-Series
 *			model, which couples film mass transfer (kf), pore diffusion (Dp), and surface diffusion (Dc)
 *			into a single lumped rate parameter (k). Those parameters all come from material properties
 *			files in the DGOSPREY framework.
 *
 *			Resistance-in-series: k = 1/(rhop*q*rp/(3*kf*C)) + (1/((rhop*q*rp*rp/(15*ep*Dp)) + (rc*rc/(15*Dc)))
 *			where rhop is the particle density, q is the adsorption, rp is the particle radius, kf is the film
 *			mass transfer parameter, C is the concentration in the gas phase, ep is the particle porosity,
 *			Dp is the pore diffusion parameter, rc is the adsorbent crystal radius, and Dc is the surface
 *			diffusion parameter.
 *
 *			van't Hoff: ln(K_i) = -dH_i/(R*T) + dS_i/R
 *			where R is the gas law constant and T is the column temperature.
 *
 *			Extended Langmuir isotherm: q_i = q_max_i * SUM(K_i*c_i)/(1+SUM(K_j*c_j))
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 10/10/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
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

#include "CoupledExtendedLangmuirLDFModel.h"
registerMooseObject("dgospreyApp", CoupledExtendedLangmuirLDFModel);

template<>
InputParameters validParams<CoupledExtendedLangmuirLDFModel>()
{
	InputParameters params = validParams<CoupledExtendedLangmuirModel>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	params.addParam<Real>("alpha",15.0,"Scaling parameter for maximum LDF parameter");
	params.addParam<Real>("beta",15.0,"Scaling parameter for minimum LDF parameter");
	return params;
}

CoupledExtendedLangmuirLDFModel::CoupledExtendedLangmuirLDFModel(const InputParameters & parameters)
: CoupledExtendedLangmuirModel(parameters),
_index(getParam<unsigned int>("index")),
_alpha(getParam<Real>("alpha")),
_beta(getParam<Real>("beta")),
_pellet_density(getMaterialProperty<Real>("pellet_density")),
_pellet_diameter(getMaterialProperty<Real>("pellet_diameter")),
_crystal_radius(getMaterialProperty<Real>("crystal_radius")),
_binder_porosity(getMaterialProperty<Real>("binder_porosity")),
_binder_fraction(getMaterialProperty<Real>("binder_ratio")),
_film_transfer(getMaterialProperty<std::vector<Real> >("film_transfer")),
_pore_diff(getMaterialProperty<std::vector<Real> >("pore_diffusion")),
_surf_diff(getMaterialProperty<std::vector<Real> >("surface_diffusion"))
{
	if (_alpha <= 0.0)
		_alpha = 15.0;
	if (_beta <= 0.0)
		_beta = 15.0;
	
	if (_alpha >= 30.0)
		_alpha = 30.0;
	if (_beta >= 30.0)
		_beta = 30.0;
	
	if (_beta/_alpha > 2.0)
		_beta = 2.0*_alpha;
}

void CoupledExtendedLangmuirLDFModel::computeLDFcoeff()
{
	CoupledExtendedLangmuirModel::computeAllLangmuirCoeffs();
	
	double part_coef = 1.0;
	double Henry = (_maxcap * _langmuircoef[_lang_index] );
	part_coef = _pellet_density[_qp]*Henry;
	
	double filmres, poreres, surfres;
	filmres = part_coef * _pellet_diameter[_qp] / (6.0*_film_transfer[_qp][_index]);
	if (_binder_fraction[_qp] <= 0.0)
		poreres = part_coef * _pellet_diameter[_qp] * _pellet_diameter[_qp] / (60.0*_pore_diff[_qp][_index]*_binder_porosity[_qp]);
	else
		poreres = part_coef * _pellet_diameter[_qp] * _pellet_diameter[_qp] / (60.0*_pore_diff[_qp][_index]*_binder_porosity[_qp]*_binder_fraction[_qp]);
	double surfdiff = _surf_diff[_qp][_index];
	if (surfdiff <= 0.0)
		surfres = 0.0;
	else
	{
		if (surfdiff < 0.001)
			surfdiff = 0.001;
		surfres = _crystal_radius[_qp] * _crystal_radius[_qp] / (15.0 * surfdiff);
	}
	
	_ldf_coeff = ((1.0/filmres) + (1.0/(poreres + surfres)));
}

void CoupledExtendedLangmuirLDFModel::computeScalingFactor()
{
	_scaling_factor = (_alpha*(1.0-(_u[_qp]/_maxcap)) + _beta*(_u[_qp]/_maxcap));
}

Real CoupledExtendedLangmuirLDFModel::computePartCoeffTempDerivative()
{
	return _maxcap*_pellet_density[_qp]*_langmuircoef[_lang_index]*_enthalpies[_lang_index]/Rstd/_coupled_temp[_qp]/_coupled_temp[_qp];
}

Real CoupledExtendedLangmuirLDFModel::computeLDFoffdiag()
{
	double part_coef = 1.0;
	double Henry = (_maxcap * _langmuircoef[_lang_index] );
	part_coef = _pellet_density[_qp]*Henry;
	
	double filmres, poreres, surfres;
	filmres = part_coef * _pellet_diameter[_qp] / (6.0*_film_transfer[_qp][_index]);
	if (_binder_fraction[_qp] <= 0.0)
		poreres = part_coef * _pellet_diameter[_qp] * _pellet_diameter[_qp] / (60.0*_pore_diff[_qp][_index]*_binder_porosity[_qp]);
	else
		poreres = part_coef * _pellet_diameter[_qp] * _pellet_diameter[_qp] / (60.0*_pore_diff[_qp][_index]*_binder_porosity[_qp]*_binder_fraction[_qp]);
	double surfdiff = _surf_diff[_qp][_index];
	if (surfdiff <= 0.0)
		surfres = 0.0;
	else
	{
		if (surfdiff < 0.001)
			surfdiff = 0.001;
		surfres = _crystal_radius[_qp] * _crystal_radius[_qp] / (15.0 * surfdiff);
	}
	
	if (_binder_fraction[_qp] <= 0.0)
	{
		return -(computePartCoeffTempDerivative()*6.0*_film_transfer[_qp][_index]/(part_coef * part_coef * _pellet_diameter[_qp])) - (computePartCoeffTempDerivative()*(_pellet_diameter[_qp] * _pellet_diameter[_qp] / (60.0*_pore_diff[_qp][_index]*_binder_porosity[_qp]))/((poreres + surfres)*(poreres + surfres)));
	}
	else
	{
		return -(computePartCoeffTempDerivative()*6.0*_film_transfer[_qp][_index]/(part_coef * part_coef * _pellet_diameter[_qp])) - (computePartCoeffTempDerivative()*(_pellet_diameter[_qp] * _pellet_diameter[_qp] / (60.0*_pore_diff[_qp][_index]*_binder_porosity[_qp]*_binder_fraction[_qp]))/((poreres + surfres)*(poreres + surfres)));
	}
}

Real CoupledExtendedLangmuirLDFModel::computeLDFjacobian()
{
	return _ldf_coeff*_phi[_j][_qp]*( (_beta/_maxcap) - (_alpha/_maxcap) );
}

Real CoupledExtendedLangmuirLDFModel::computeQpResidual()
{
	computeLDFcoeff();
	computeScalingFactor();
	return _ldf_coeff*_scaling_factor*CoupledExtendedLangmuirFunction::computeQpResidual() + _test[_i][_qp] * _u_dot[_qp];
}

Real CoupledExtendedLangmuirLDFModel::computeQpJacobian()
{
	computeLDFcoeff();
	computeScalingFactor();
	return _test[_i][_qp]*_ldf_coeff*_scaling_factor*_phi[_j][_qp] + _test[_i][_qp]*_u[_qp]*computeLDFjacobian() + _test[_i][_qp] * _phi[_j][_qp] * _du_dot_du[_qp];
}

Real CoupledExtendedLangmuirLDFModel::computeQpOffDiagJacobian(unsigned int jvar)
{
	computeLDFcoeff();
	computeScalingFactor();
	
	//Off-diagonals for non-main concentrations
	for (unsigned int i = 0; i<_coupled.size(); ++i)
	{
		if (jvar == _coupled_vars[i] && jvar != _coupled_var_i)
		{
			return _test[_i][_qp]*_ldf_coeff*_scaling_factor*CoupledExtendedLangmuirFunction::computeExtLangmuirOffJacobi(i);
		}
	}
	// Off-diagonal element for main coupled gas
	if (jvar == _coupled_var_i)
	{
		return -_test[_i][_qp]*_ldf_coeff*_scaling_factor*CoupledExtendedLangmuirFunction::computeExtLangmuirConcJacobi();
	}
	
	// Off-diagonal element for coupled temperature
	if (jvar == _coupled_var_temp)
	{
		return -_test[_i][_qp]*_ldf_coeff*_scaling_factor*CoupledExtendedLangmuirModel::computeExtLangmuirTempJacobi() + _test[_i][_qp]*_scaling_factor*computeLDFoffdiag()*(_u[_qp] - CoupledExtendedLangmuirFunction::computeExtLangmuirEquilibrium());
	}
	
	
	return 0.0;
}




