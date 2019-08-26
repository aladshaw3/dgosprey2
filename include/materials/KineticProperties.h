/*!
 *  \file KineticProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with SCOPSOWL simulations
 *	\details This file creates a material property object for the SCOPSOWL data structure and associated constants.
 *			That information is used in conjunction with the SCOPSOWL simulation functions (see scopsowl.h) in order
 *			to approximate the adsorption capacities and rates of each gas species in a given system
 *			for a given adsorbent.
 *
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

#pragma once

#include "Material.h"
#include "flock.h"
#include "DataStruct_StoreLoad.h"

/// KineticProperties class object forward declaration
class KineticProperties;

template<>
InputParameters validParams<KineticProperties>();

/// KineticProperties class object inherits from Material object
/** This class object inherits from the Material object in the MOOSE framework.
	All public and protected members of this class are required function overrides. The object
	will set up the SCOPSOWL_DATA structure (see scopsowl.h) based on user provied input from the
	input file and requires the MAGPIE_DATA structure from the MagpieAdsorbateProperties object.
 That information will be used to estimate the adsorption and rate of adsorption of each species in
	the system based on temperature, pressure, and concentrations of each species.
 */
class KineticProperties : public Material
{
public:
	/// Required constructor for objects in MOOSE
	KineticProperties(const InputParameters & parameters);
	
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
	
	bool DirichletBC;									///< Input parameter to determine type of BC for SCOPSOWL
	bool Heterogeneous;									///< Input parameter to determine type of adsorbent configuration
	bool SurfaceDiffusion;								///< Input parameter to determine whether or not to include surface diffusion
	bool MacroSpheres;									///< Input parameter to determine if macro pellets are spherical or cylindrical
	bool MicroSpheres;									///< Input parameter to determine if micro cystrals are spherical or cylindrical
	Real MacroLength;									///< Length of cylindrical pellets (cm)
	Real MicroLength;									///< Length of cylindrical crystals (um)
	
	const MaterialProperty< MAGPIE_DATA > & _magpie_dat;	///< Coupled material property for MAGPIE_DATA from MagpieProperties
	const MaterialProperty< MIXED_GAS > & _coupled_gas;		///< Coupled material property for MIXED_GAS from FlowProperties
	const MaterialProperty<Real> & _velocity;				///< Linear velocity in the column (cm/hr)
	const MaterialProperty<Real> & _pellet_diameter;		///< Diameter of the adsobent pellet (cm)
	const MaterialProperty<Real> & _pellet_density;			///< Pellet material density (g/cm^3)
	const MaterialProperty<Real> & _binder_porosity;		///< Porosity of the pellet binder material
	const MaterialProperty<Real> & _pore_size;				///< Nominal pore size of the binder material (cm)
	const MaterialProperty<Real> & _crystal_radius;			///< Nominal size of crystals in adsorbent pellet (um)
	const MaterialProperty<Real> & _binder_ratio;			///< Ratio of binder to crystal in the adsorbent
	
	const MaterialProperty<std::vector<Real> > & _ref_diffusion;		///< MaterialProperty for reference diffusivity (um^2/hr)
	const MaterialProperty<std::vector<Real> > & _activation_energy;	///< MaterialProperty for activation energy (J/mol)
	const MaterialProperty<std::vector<Real> > & _ref_temperature;		///< MaterialProperty for reference temperature (K)
	const MaterialProperty<std::vector<Real> > & _affinity_coeff;		///< MaterialProperty for affinity coefficient (-)
	
	/// MaterialProperty object to hold the SCOPSOWL_DATA structure and all relavent information
	/** This is the object that needs to interface with the MAGPIE functions in order to solve for
		variable information such as adsorption capacities, mixed gas adsorption equilibria, and
		heats of adsorption. */
	MaterialProperty< SCOPSOWL_DATA > & _owl_dat;
	
	/// Old MaterialProperty object to hold the SCOPSOWL_DATA structure and all relavent information
	/** This object is required to be created in order to use the stateful properties, which is how
		we intialize our custom objects in MOOSE correctly. */
	const MaterialProperty< SCOPSOWL_DATA > & _owl_dat_old;
	
	/// MaterialProperty object to hold the MIXED_GAS structure and all relavent information
	/** This is the object that needs to interface with the MAGPIE functions in order to solve for
		variable information such as adsorption capacities, mixed gas adsorption equilibria, and
		heats of adsorption. */
	MaterialProperty< MIXED_GAS > & _gas_dat;
	
	/// Old MaterialProperty object to hold the MIXED_GAS structure and all relavent information
	/** This object is required to be created in order to use the stateful properties, which is how
		we intialize our custom objects in MOOSE correctly. */
	const MaterialProperty< MIXED_GAS > & _gas_dat_old;
	
	std::vector<const VariableValue *> _solid_conc;		///< Pointer list to the coupled adsorption concentrations
	
};
