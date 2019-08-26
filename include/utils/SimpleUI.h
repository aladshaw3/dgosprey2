/*!
 *  \file SimpleUI.h
 *	\brief Utilities kernel to provide a simpler user interface through yaml input files
 *	\details This file is responsible for providing users with a very simple, optional
 *				user interface that utilizes the yaml file structure. It is coupled with
 *				other data base files, stored in the project directory, that holds parameter
 *				information for common systems that the user may wish to model. The intent
 *				is to give an easier way for the uneducated user to be capable of utilizing
 *				this software without any advanced understanding or manual. 
 *
 *				The simple user interface will read input from a yaml (.yml) file and then
 *				read in the necessary data base information to create a standard MOOSE
 *				input file (.i) that the DGOSPREY program will run. As such, this interface
 *				will have significantly less control over the options that the lower level
 *				interface will have.
 *
 *
 *  \author Austin Ladshaw
 *	\date 01/19/2017
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

#pragma once

#include "yaml_wrapper.h"
#include <cstring>

/// Function to return true if the file extension of the given data is .yml
bool isYamlFile(char argv[]);

/// Execute the simple user interface to read yaml files and create DGOSPREY input files
int exec_SimpleUI(char *file);

/// Class structure for using the simple user interface
/** The following class object is used for reading in the simple yaml input files and/or
	data base files in order to determine the type of DGOSPREY simulation to run. After
	reading in the appropriate information, this class must create and save a MOOSE
	input file (.i) that is to be read in by DGOSPREY and executed. 
 
	\warning THIS OBJECT IS STILL UNDER ACTIVE DEVELOPMENT*/
class SimpleUI
{
public:
	SimpleUI();								///< Default Constructor
	~SimpleUI();							///< Default Destructor
	
	int readInputFile(const char *file);	///< Function to read the given yaml input file
	int writeOutputFile(const char *file);	///< Function to write an MOOSE input file based on the yaml input file
	void writeOutBlock(FILE *file,std::string name);	///< Function to write the specific block of the input file
	void createMooseBlank();				///< Function to create a blank MOOSE input file based on necessary arguments
	void createExample();					///< Function to create an example MOOSE input file (Kr_Xe_253K_OWL.i)
	
	YamlWrapper& getYamlInput();			///< Returns reference to the YamlWrapper Object for the input file
	YamlWrapper& getMooseInput();			///< Returns reference to the YamlWrapper that is used to build the MOOSE file
	
	void DisplayInput();					///< Function to display the contents of the input file
	void DisplayOutput();					///< Function to display the contents of the MOOSE input file we are making
	
private:
	yaml_cpp_class yaml_input;				///< Yaml object that will read and store the main input file
	YamlWrapper moose_input;				///< YamlWrapper object that will be used to create the moose input file
	
};
