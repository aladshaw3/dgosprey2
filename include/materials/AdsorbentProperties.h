/*!
 *  \file AdsorbentProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with the adsorbent
 *	\details This file creates a material property object for various properties of a given adsorbent. These
 *			properties are then used in other material property files and/or kernels to calculate information
 *			such as linear velocity, mechanical dispersion, or any adsorption kinetic parameters. 
 *
 *	\note Currently, we do not couple with adsorption kinetics, so this file is only used in conjunction with
 *		the linear velocity and mechanical dispersion properties.
 *
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

#pragma once

#include "Material.h"
#include "flock.h"
#include "DataStruct_StoreLoad.h"

/// AdsorbentProperties class object forward declaration
class AdsorbentProperties;

template<>
InputParameters validParams<AdsorbentProperties>();

/// AdsorbentProperties class object inherits from Material object
/** This class object inherits from the Material object in the MOOSE framework.
	All public and protected members of this class are required function overrides. The object
	will set up the structural information about adsorbent in the system that will be used when
	determining flow properties, as well as some kinetic properties for adsorption dynamics. */
class AdsorbentProperties : public Material
{
public:
	/// Required constructor for objects in MOOSE
	AdsorbentProperties(const InputParameters & parameters);
	
protected:
	/// Required function override for Material objects in MOOSE
	/** This function computes the material properties when they are needed by other MOOSE objects.*/
	virtual void computeQpProperties();
	
	/// Required function override for Stateful Material objects in MOOSE
	/** This function is needed because we have to properly initialize our custom objects without
		having to reinitialize at each compute step. It takes more memory this way, but also prevents
		segfault errors and helps the kernel run faster after initialization. */
	virtual void initQpStatefulProperties();
	
private:
	Real _binder_fraction;		///< Binder fraction in the biporous adsorbent pellet (0 means no binder material)
	Real _eps_binder;			///< Macro-porosity of the binder material in the adsorbent pellet
	Real _crystal_rad;			///< Nominal radius of the adsorbent crystals suspended in the binder (um)
	Real _pellet_dia;			///< Nominal diameter of the adsorbent pellets in the system (cm)
	Real _macropore_radius;		///< Nominal size of the macro-pores in the pellet (cm)
	Real _rhos;					///< Density of the adsorbent pellet (g/cm^3)
	Real _hs;					///< Heat capacity of the adsorbent pellet (J/g/K)
	
	std::vector<Real> _ref_diff;	///< Reference Surface Diffusivity (um^2/hr)
	std::vector<Real> _act_energy;	///< Activation Energy of Surface Diffusion (J/mol)
	std::vector<Real> _ref_temp;	///< Reference Temperature for Surface Diffusion (K)
	std::vector<Real> _affinity;	///< Affinity coefficient for Surface Diffusion (-)
	
	MaterialProperty<std::vector<Real> > & _ref_diffusion;		///< MaterialProperty for reference diffusivity (um^2/hr)
	MaterialProperty<std::vector<Real> > & _activation_energy;	///< MaterialProperty for activation energy of surface diffusion (J/mol)
	MaterialProperty<std::vector<Real> > & _ref_temperature;	///< MaterialProperty for reference temperature (K)
	MaterialProperty<std::vector<Real> > & _affinity_coeff;		///< MaterialProperty for affinity coefficient of surface diffusion (-)
	
	MaterialProperty<Real> & _pellet_density;			///< MaterialProperty for the pellet density (g/cm^3)
	MaterialProperty<Real> & _pellet_heat_capacity;		///< MaterialProperty for the pellet heat capacity (J/g/K)
	MaterialProperty<Real> & _pellet_diameter;			///< MaterialProperty for pellet diameter (cm)
	MaterialProperty<Real> & _crystal_radius;			///< MaterialProperty for the crystal radius (um)
	MaterialProperty<Real> & _binder_porosity;			///< MaterialProperty for the binder porosity
	MaterialProperty<Real> & _binder_ratio;				///< MaterialProperty for the ratio of binder to pellet volumes
	MaterialProperty<Real> & _pore_size;				///< MaterialProperty for the macropore radius (cm)
	
	MaterialProperty<std::vector<Real> > & _surface_diffusion;		///< MaterialProperty for the surface diffusion (um^2/hr)
	const MaterialProperty<std::vector<Real> > & _surf_diff_old;	///< MaterialProperty for stateful surface diffusion
	
	const MaterialProperty< MAGPIE_DATA > & _magpie_dat;			///< Pointer to MAGPIE_DATA material property
	
	const VariableValue & _temperature;						///< Reference to the coupled column temperature
	const VariableValue & _temperature_old;					///< Reference to the coupled column temperature at previous time
	std::vector<unsigned int> _index;						///< List of indices for the coupled gases
	std::vector<const VariableValue *> _gas_conc;			///< Pointer list for the coupled gases
	std::vector<const VariableValue *> _gas_conc_old;		///< Pointer list for the coupled gases at previous time
	
};
