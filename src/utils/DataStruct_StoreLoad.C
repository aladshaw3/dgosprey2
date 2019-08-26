/*!
 *  \file DataStruct_StoreLoad.h
 *	\brief Utilities kernel to perform serial loading and storing of custom data structures and classes
 *  \author Austin Ladshaw
 *	\date 01/17/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
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

#include "DataStruct_StoreLoad.h"

/// dataStore function for MAGPIE_DATA
template<>
void
dataStore(std::ostream & stream, MAGPIE_DATA & data, void * context)
{
	
}

/// dataLoad function for MAGPIE_DATA
template<>
void
dataLoad(std::istream & stream, MAGPIE_DATA & data, void * context)
{
	
}

/// dataStore function for SCOPSOWL_DATA
template<>
void
dataStore(std::ostream & stream, SCOPSOWL_DATA & data, void * context)
{
	
}

/// dataLoad function for SCOPSOWL_DATA
template<>
void
dataLoad(std::istream & stream, SCOPSOWL_DATA & data, void * context)
{
	
}

/// dataStore function for MIXED_GAS
template<>
void
dataStore(std::ostream & stream, MIXED_GAS & data, void * context)
{
	
}

/// dataLoad function for MIXED_GAS
template<>
void
dataLoad(std::istream & stream, MIXED_GAS & data, void * context)
{
	
}
