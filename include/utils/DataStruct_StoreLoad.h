/*!
 *  \file DataStruct_StoreLoad.h
 *	\brief Utilities kernel to perform serial loading and storing of custom data structures and classes
 *	\details This file is responsible for defining the functions that are used by MOOSE to Store and Load
 *				custom user information used in the MOOSE system. Storing and Loading of data is used in 
 *				MOOSE as part of the restart and multiapps functionality. Any custom user objects that are
 *				defined or used in MOOSE are now required to have specialized dataStore and dataLoad sub-
 *				routines, regardless of whether or not these features are actually used by the code.
 *
 *	\warning PLEASE NOTE: These functions are currently left blank just to appease the compiler gods. You
 *				cannot use restarts or multiapps until these functions are properly filled in.
 *
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

#pragma once

#include "GeneralUserObject.h"
#include "flock.h"

/// dataStore function for MAGPIE_DATA
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataStore functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
void
dataStore(std::ostream & stream, MAGPIE_DATA & data, void * context);

/// dataLoad function for MAGPIE_DATA
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataLoad functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
void
dataLoad(std::istream & stream, MAGPIE_DATA & data, void * context);

/// dataStore function for SCOPSOWL_DATA
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataStore functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
void
dataStore(std::ostream & stream, SCOPSOWL_DATA & data, void * context);

/// dataLoad function for SCOPSOWL_DATA
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataLoad functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
void
dataLoad(std::istream & stream, SCOPSOWL_DATA & data, void * context);

/// dataStore function for MIXED_GAS
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataStore functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
void
dataStore(std::ostream & stream, MIXED_GAS & data, void * context);

/// dataLoad function for MIXED_GAS
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataLoad functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
void
dataLoad(std::istream & stream, MIXED_GAS & data, void * context);

