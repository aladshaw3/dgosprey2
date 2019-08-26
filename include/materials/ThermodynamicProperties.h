/*!
 *  \file ThermodynamicProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with equilibria simulations
 *	\details This file creates a material property object for the MAGPIE data structure and associated constants.
 *			That information is used in conjunction with the MAGPIE simulation functions (see magpie.h) in order
 *			to approximate the adsorption capacities and adsorbed amounts of each gas species in a given system
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

/// ThermodynamicProperties class object forward declaration
class ThermodynamicProperties;

template<>
InputParameters validParams<ThermodynamicProperties>();

/// ThermodynamicProperties class object inherits from Material object
/** This class object inherits from the Material object in the MOOSE framework.
	All public and protected members of this class are required function overrides. The object
	will set up the MAGPIE_DATA structure (see magpie.h) based on user provied input from the
	input file. That information will be used to estimate the adsorption of each species in
	the system based on temperature, pressure, and concentrations of each species.
 
	\note The GSTA isotherm model for each species allows upto 6 energetically distinct adsorption
	sites. If those sites are not used by a particular species, then those energies should be left
	as zeros in the input files and the number of relavent sites for each species needs to be recorded
	in the input file. Each species is allowed to have a different number of adsorption sites in a
	particular adsorbent. */
class ThermodynamicProperties : public Material
{
public:
	/// Required constructor for objects in MOOSE
	ThermodynamicProperties(const InputParameters & parameters);
	
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
	
	std::vector<unsigned int> _index;					///< Indices for the gas species in the system
	const VariableValue & _temperature;					///< Reference to the coupled column temperature
	const VariableValue & _temperature_old;				///< Reference to the coupled column temperature at previous time step
	const VariableValue & _total_pressure;				///< Reference to the coupled column pressure
	std::vector<const VariableValue *> _gas_conc;		///< Pointer list to the coupled gases
	std::vector<const VariableValue *> _gas_conc_old;	///< Pointer list to the coupled gases at previous time step
	
	std::vector<int> _num_sites;				///< List of the number of sites each gas species' isotherm contains
	std::vector<Real> _max_capacity;			///< List of the maximum adsorption capacities of each gas species
	std::vector<Real> _molar_volume;			///< List of the van der Waal's molar volumes of each species
	
	std::vector<Real> _enthalpy_1;				///< List of the site 1 enthalpies for each gas species
	std::vector<Real> _enthalpy_2;				///< List of the site 2 enthalpies for each gas species
	std::vector<Real> _enthalpy_3;				///< List of the site 3 enthalpies for each gas species
	std::vector<Real> _enthalpy_4;				///< List of the site 4 enthalpies for each gas species
	std::vector<Real> _enthalpy_5;				///< List of the site 5 enthalpies for each gas species
	std::vector<Real> _enthalpy_6;				///< List of the site 6 enthalpies for each gas species
	
	std::vector<Real> _entropy_1;				///< List of the site 1 entropies for each gas species
	std::vector<Real> _entropy_2;				///< List of the site 2 entropies for each gas species
	std::vector<Real> _entropy_3;				///< List of the site 3 entropies for each gas species
	std::vector<Real> _entropy_4;				///< List of the site 4 entropies for each gas species
	std::vector<Real> _entropy_5;				///< List of the site 5 entropies for each gas species
	std::vector<Real> _entropy_6;				///< List of the site 6 entropies for each gas species
	
	/// MaterialProperty object to hold the MAGPIE_DATA structure and all relavent information
	/** This is the object that needs to interface with the MAGPIE functions in order to solve for
		variable information such as adsorption capacities, mixed gas adsorption equilibria, and
		heats of adsorption. */
	MaterialProperty< MAGPIE_DATA > & _magpie_dat;
	
	/// Old MaterialProperty object to hold the MAGPIE_DATA structure and all relavent information
	/** This object is required to be created in order to use the stateful properties, which is how
		we intialize our custom objects in MOOSE correctly. */
	const MaterialProperty< MAGPIE_DATA > & _magpie_dat_old;
	
	MaterialProperty<std::vector<Real> > & _ads_heat;	///< MaterialProperty for each species' heat of adsorption coefficient
};
