/*!
 *  \file BedProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with the fixed-bed
 *	\details This file creates a material property object for various properties of the fixed bed. Those properties
 *			are used in conjunction with other kernels and materials to establish information such as heat transfer
 *			coefficients, conductivities, and size parameters.
 *
 *
 *	\warning THIS KERNEL IS INCOMPLETE! ONLY USED FOR DATA STORAGE FOR VARIOUS BED PARAMETERS!
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

/// BedProperties class object forward declaration
class BedProperties;

template<>
InputParameters validParams<BedProperties>();

/// BedProperties class object inherits from Material object
/** This class object inherits from the Material object in the MOOSE framework.
	All public and protected members of this class are required function overrides. The object
	will set up the parameters of the fixed-bed column. Those parameters include: length, diameter,
	thermal conductivity, heat transfer coefficients, bulk porosity, etc. */
class BedProperties : public Material
{
public:
	/// Required constructor for objects in MOOSE
	BedProperties(const InputParameters & parameters);

protected:
	/// Required function override for Material objects in MOOSE
	/** This function computes the material properties when they are needed by other MOOSE objects.*/
	virtual void computeQpProperties();

private:
	Real _length;			///< Bed length (cm)
	Real _din;				///< Column inner diameter (cm)
	Real _dout;				///< Column outer diameter (cm)
	Real _eb;				///< Bulk porosity of the bed
	Real _rhow;				///< Density of the column wall (g/cm^3)
	Real _hw;				///< Heat capacity of the column wall (J/g/K)
	Real _Uw;				///< Bed-Wall heat transfer coefficient (J/hr/cm^2/K)
	Real _Ua;				///< External-Wall heat transfer coefficient (J/hr/cm^2/K)
	
	MaterialProperty<Real> & _inner_dia;						///< MaterialProperty for column inner diameter
	MaterialProperty<Real> & _outer_dia;						///< MaterialProperty for column outer diameter
	MaterialProperty<Real> & _bed_length;							///< MaterialProperty for column length (cm)
	MaterialProperty<Real> & _porosity;							///< MaterialProperty for bulk porosity of the bed
	MaterialProperty<Real> & _wall_density;						///< MaterialProperty for column wall density
	MaterialProperty<Real> & _wall_heat_capacity;				///< MaterialProperty for column wall heat capacity
	MaterialProperty<Real> & _bed_wall_transfer_coeff;			///< MaterialProperty for bed-wall heat transfer coefficient
	MaterialProperty<Real> & _wall_exterior_transfer_coeff;		///< MaterialProperty for exterior-wall heat transfer coefficient

};
