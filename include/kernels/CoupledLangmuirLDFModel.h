/*!
 *  \file CoupledLangmuirLDFModel.h
 *	\brief Standard kernel for coupling non-linear variables via the langmuir model with LDF kinetics
 *	\details This file creates a standard MOOSE kernel for the coupling of non-linear variables
 *			together via the langmuir model, and applies linear driving force kinetics for the rate
 *			of adsorption.
 *
 *			This kernel extends the CoupledLangmuirModel kernel by calculating the model parameters from
 *			information in the input file for enthalpy and entropy. In addition, it calculates the Linear Driving Force
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
 *			van't Hoff: ln(K) = -dH/(R*T) + dS/R
 *			where R is the gas law constant and T is the column temperature.
 *
 *			Langmuir isotherm: q = q_max * (K*c)/(1+(K*c))
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 10/09/2017
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

#pragma once

#include "CoupledLangmuirModel.h"

/// CoupledLangmuirLDFModel class object forward declarationss
class CoupledLangmuirLDFModel;

template<>
InputParameters validParams<CoupledLangmuirLDFModel>();

/// CoupledLangmuirLDFModel class object inherits from CoupledLangmuirModel object
/** This class object inherits from the CoupledLangmuirModel
	in the DGOSPREY framework. All public and protected members of this class are required
	function overrides. The kernel interfaces the two non-linear variables to couple the Langmuir
	isotherm model with concentration and temperature. Parameters of the Langmuir model are determined
	through the input parameters of enthalpy and entropy. Linear Driving Force parameters are
	determined through GasFlowProperties, AdsorbentProperties, and ThermodynamicProperties.*/
class CoupledLangmuirLDFModel : public CoupledLangmuirModel
{
public:
	/// Required constructor for objects in MOOSE
	CoupledLangmuirLDFModel(const InputParameters & parameters);
	
protected:
	/// Function to compute the linear driving force coefficient from material properties
	void computeLDFcoeff();
	
	/// Function to compute the scaling factor for the linear driving force coefficient
	void computeScalingFactor();
	
	/// Function to compute the temperature derivative of the partition coefficient
	Real computePartCoeffTempDerivative();
	
	/// Function to compute the temperature derivative of the LDF coefficient
	Real computeLDFoffdiag();
	
	/// Function to compute the jacobian for the linear driving force coefficient
	Real computeLDFjacobian();
	
	/// Required residual function for standard kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	
	/// Required Jacobian function for standard kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
	 computed is the associated diagonal element in the overall Jacobian matrix for the
	 system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
	
	/// Not Required, but aids in the preconditioning step
	/** This function returns the off diagonal Jacobian contribution for this object. By
	 returning a non-zero value we will hopefully improve the convergence rate for the
	 cross coupling of the variables. */
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);
	
	unsigned int _index;											///< Index of the gaseous species to calculate adsorbtion for
	Real _ldf_coeff;												///< Parameter place holder for LDF parameter (1/hr)
	Real _scaling_factor;											///< Parameter used to scale the LDF parameter
	Real _alpha, _beta;												///< Parameters for the scaling factor
	const MaterialProperty< Real > & _pellet_density;				///< Material Property for pellet density (kg/L)
	const MaterialProperty< Real > & _pellet_diameter;				///< Material Property for pellet diameter (cm)
	const MaterialProperty< Real > & _crystal_radius;				///< Material Property for crystal radius (um)
	const MaterialProperty< Real > & _binder_porosity;				///< Material Property for binder porosity (-)
	const MaterialProperty< Real > & _binder_fraction;				///< Material Property for binder fraction (-)
	const MaterialProperty<std::vector<Real> > & _film_transfer;	///< Material Property for film transfer coef (cm/hr)
	const MaterialProperty<std::vector<Real> > & _pore_diff;		///< Material Property for pore diffusion coef (cm^2/hr)
	const MaterialProperty<std::vector<Real> > & _surf_diff;		///< Material Property for pore diffusion coef (um^2/hr)
	
private:
	
};
