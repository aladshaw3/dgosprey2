/*!
 *  \file SimpleUI.h
 *	\brief Utilities kernel to provide a simpler user interface through yaml input files
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

#include "SimpleUI.h"

/// Function to return true if the file extension of the given data is .yml
bool isYamlFile(char argv[])
{
	std::string arg = argv;
	if((arg.substr(arg.find_last_of(".") + 1) == "yml") || (arg.substr(arg.find_last_of(".") + 1) == "yaml"))
		return true;
	else
		return false;
}

/// Execute the simple user interface to read yaml files and create DGOSPREY input files
int exec_SimpleUI(char *file)
{
	//Declarations
	int success = 0;
	std::string arg = file;
	SimpleUI sui;
	
	// Create a blank digital MOOSE input file with default and required information
	sui.createMooseBlank();
	
	//Read the input file
	success = sui.readInputFile(file);
	if (success != 0)
	{
		mError(read_error);
		return -1;
	}
	
	//Check input file for success flags
	
	//Perform preprocessing for finalization of MOOSE YamlWrapper object
	
	//Create the MOOSE YamlWrapper Object
	
	//create example file (replace with real file)
	sui.createExample();
	
	//Create output file for MOOSE
	std::strcpy (file, ((arg.erase(arg.find_last_of(".") + 1)).append("i")).c_str());
	sui.writeOutputFile(file);
	
	return 1;
}

//Default constructor
SimpleUI::SimpleUI()
{
}

//Default destructor
SimpleUI::~SimpleUI()
{
}

//Function to read a given input file
int SimpleUI::readInputFile(const char *file)
{
	return this->yaml_input.executeYamlRead(file);
}

// Function to write MOOSE input file
int SimpleUI::writeOutputFile(const char *file)
{
	//Declarations
	int success = 0;
	FILE *Output;
	
	//Open the output file for editing
	Output = fopen(file,"w+");
	
	/** Will not iterate through documents because we want a specific order */
	this->writeOutBlock(Output,"GlobalParams");
	this->writeOutBlock(Output,"Problem");
	this->writeOutBlock(Output,"Mesh");
	this->writeOutBlock(Output,"Variables");
	this->writeOutBlock(Output,"AuxVariables");
	this->writeOutBlock(Output,"ICs");
	this->writeOutBlock(Output,"Kernels");
	this->writeOutBlock(Output,"DGKernels");
	this->writeOutBlock(Output,"AuxKernels");
	this->writeOutBlock(Output,"BCs");
	this->writeOutBlock(Output,"Materials");
	this->writeOutBlock(Output,"Postprocessors");
	this->writeOutBlock(Output,"Executioner");
	this->writeOutBlock(Output,"Preconditioning");
	this->writeOutBlock(Output,"Outputs");
	
	//Close the file and return
	fclose(Output);
	return success;
}

//Function to write an specific block to the file
void SimpleUI::writeOutBlock(FILE *file,std::string name)
{
	fprintf(file,"[%s]\n",this->moose_input.getDocument(name).getName().c_str());
	for (auto &x: this->moose_input.getDocument(name).getDataMap())
		fprintf(file,"\t%s = %s\n",x.first.c_str(),x.second.getValue().c_str());
	for (auto &x: this->moose_input.getDocument(name).getHeadMap())
	{
		fprintf(file,"\t[./%s]\n",x.second.getName().c_str());
		for (auto &y: x.second.getDataMap())
			fprintf(file,"\t\t%s = %s\n",y.first.c_str(),y.second.getValue().c_str());
		fprintf(file,"\t[../]\n");
	}
	fprintf(file,"[] #END %s\n\n",name.c_str());
}

//Function to create a blank MOOSE input file
void SimpleUI::createMooseBlank()
{
	this->moose_input.addDocKey("GlobalParams");
	
	this->moose_input.addDocKey("Problem");
	this->moose_input.getDocument("Problem").addPair("coord_type","RZ");
	
	this->moose_input.addDocKey("Mesh");
	this->moose_input.getDocument("Mesh").addPair("type","GeneratedMesh");
	this->moose_input.getDocument("Mesh").addPair("dim","2");
	this->moose_input.getDocument("Mesh").addPair("nx","10");
	this->moose_input.getDocument("Mesh").addPair("ny","40");
	this->moose_input.getDocument("Mesh").addPair("xmin","0.0");
	this->moose_input.getDocument("Mesh").addPair("ymin","0.0");
	
	this->moose_input.addDocKey("Variables");
	this->moose_input.getDocument("Variables").addHeadKey("wall_temp");
	this->moose_input.getDocument("Variables").getHeader("wall_temp").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("wall_temp").addPair("family","MONOMIAL");
	this->moose_input.getDocument("Variables").addHeadKey("column_temp");
	this->moose_input.getDocument("Variables").getHeader("column_temp").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("column_temp").addPair("family","MONOMIAL");
	
	this->moose_input.addDocKey("AuxVariables");
	this->moose_input.getDocument("AuxVariables").addHeadKey("total_pressure");
	this->moose_input.getDocument("AuxVariables").getHeader("total_pressure").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("total_pressure").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("total_pressure").addPair("initial_condition","101.35");
	this->moose_input.getDocument("AuxVariables").addHeadKey("ambient_temp");
	this->moose_input.getDocument("AuxVariables").getHeader("ambient_temp").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("ambient_temp").addPair("family","MONOMIAL");
	
	this->moose_input.addDocKey("ICs");
	
	this->moose_input.addDocKey("Kernels");
	this->moose_input.getDocument("Kernels").addHeadKey("wallAccum");
	this->moose_input.getDocument("Kernels").getHeader("wallAccum").addPair("type","WallHeatAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("wallAccum").addPair("variable","wall_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("wall_bed_trans");
	this->moose_input.getDocument("Kernels").getHeader("wall_bed_trans").addPair("type","BedWallHeatTransfer");
	this->moose_input.getDocument("Kernels").getHeader("wall_bed_trans").addPair("variable","wall_temp");
	this->moose_input.getDocument("Kernels").getHeader("wall_bed_trans").addPair("coupled","column_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("wall_amb_trans");
	this->moose_input.getDocument("Kernels").getHeader("wall_amb_trans").addPair("type","WallAmbientHeatTransfer");
	this->moose_input.getDocument("Kernels").getHeader("wall_amb_trans").addPair("variable","wall_temp");
	this->moose_input.getDocument("Kernels").getHeader("wall_amb_trans").addPair("coupled","ambient_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("columnAccum");
	this->moose_input.getDocument("Kernels").getHeader("columnAccum").addPair("type","BedHeatAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("columnAccum").addPair("variable","column_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("columnConduction");
	this->moose_input.getDocument("Kernels").getHeader("columnConduction").addPair("type","GColumnHeatDispersion");
	this->moose_input.getDocument("Kernels").getHeader("columnConduction").addPair("variable","column_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("columnAdvection");
	this->moose_input.getDocument("Kernels").getHeader("columnAdvection").addPair("type","GColumnHeatAdvection");
	this->moose_input.getDocument("Kernels").getHeader("columnAdvection").addPair("variable","column_temp");
	
	this->moose_input.addDocKey("DGKernels");
	this->moose_input.getDocument("DGKernels").addHeadKey("DGcolumnConduction");
	this->moose_input.getDocument("DGKernels").getHeader("DGcolumnConduction").addPair("type","DGColumnHeatDispersion");
	this->moose_input.getDocument("DGKernels").getHeader("DGcolumnConduction").addPair("variable","column_temp");
	this->moose_input.getDocument("DGKernels").addHeadKey("DGcolumnAdvection");
	this->moose_input.getDocument("DGKernels").getHeader("DGcolumnAdvection").addPair("type","DGColumnHeatAdvection");
	this->moose_input.getDocument("DGKernels").getHeader("DGcolumnAdvection").addPair("variable","column_temp");
	
	this->moose_input.addDocKey("AuxKernels");
	this->moose_input.getDocument("AuxKernels").addHeadKey("column_pressure");
	this->moose_input.getDocument("AuxKernels").getHeader("column_pressure").addPair("type","TotalColumnPressure");
	this->moose_input.getDocument("AuxKernels").getHeader("column_pressure").addPair("variable","total_pressure");
	this->moose_input.getDocument("AuxKernels").getHeader("column_pressure").addPair("temperature","column_temp");
	
	this->moose_input.addDocKey("BCs");
	this->moose_input.getDocument("BCs").addHeadKey("Heat_Gas_Flux");
	this->moose_input.getDocument("BCs").getHeader("Heat_Gas_Flux").addPair("variable","column_temp");
	this->moose_input.getDocument("BCs").getHeader("Heat_Gas_Flux").addPair("boundary","'top bottom'");
	this->moose_input.getDocument("BCs").addHeadKey("Heat_Wall_Flux");
	this->moose_input.getDocument("BCs").getHeader("Heat_Wall_Flux").addPair("variable","column_temp");
	this->moose_input.getDocument("BCs").getHeader("Heat_Wall_Flux").addPair("boundary","'right left'");
	this->moose_input.getDocument("BCs").getHeader("Heat_Wall_Flux").addPair("wall_temp","wall_temp");
	
	
	
	this->moose_input.addDocKey("Materials");
	this->moose_input.getDocument("Materials").addHeadKey("BedMaterials");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("type","BedProperties");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("temperature","column_temp");
	
	this->moose_input.getDocument("Materials").addHeadKey("FlowMaterials");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("type","FlowProperties");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("temperature","column_temp");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("total_pressure","total_pressure");
	
	this->moose_input.getDocument("Materials").addHeadKey("AdsorbentMaterials");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("type","AdsorbentProperties");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("temperature","column_temp");
	
	this->moose_input.getDocument("Materials").addHeadKey("AdsorbateMaterials");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("type","MagpieAdsorbateProperties");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("temperature","column_temp");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("total_pressure","total_pressure");
	
	this->moose_input.addDocKey("Postprocessors");
	this->moose_input.getDocument("Postprocessors").addHeadKey("wall_temp");
	this->moose_input.getDocument("Postprocessors").getHeader("wall_temp").addPair("type","SideAverageValue");
	this->moose_input.getDocument("Postprocessors").getHeader("wall_temp").addPair("boundary","'right'");
	this->moose_input.getDocument("Postprocessors").getHeader("wall_temp").addPair("variable","wall_temp");
	this->moose_input.getDocument("Postprocessors").getHeader("wall_temp").addPair("execute_on","'initial timestep_end'");
	
	this->moose_input.addDocKey("Executioner");
	this->moose_input.getDocument("Executioner").addPair("type","Transient");
	this->moose_input.getDocument("Executioner").addPair("scheme","implicit-euler");
	this->moose_input.getDocument("Executioner").addPair("nl_rel_tol","1e-6");
	this->moose_input.getDocument("Executioner").addPair("nl_abs_tol","1e-6");
	this->moose_input.getDocument("Executioner").addPair("nl_rel_step_tol","1e-10");
	this->moose_input.getDocument("Executioner").addPair("nl_abs_step_tol","1e-10");
	this->moose_input.getDocument("Executioner").addPair("l_tol","1e-6");
	this->moose_input.getDocument("Executioner").addPair("l_max_its","100");
	this->moose_input.getDocument("Executioner").addPair("nl_max_its","10");
	this->moose_input.getDocument("Executioner").addPair("solve_type","newton");
	this->moose_input.getDocument("Executioner").addPair("line_search","none");
	this->moose_input.getDocument("Executioner").addPair("start_time","0.0");
	this->moose_input.getDocument("Executioner").addPair("petsc_options_iname","'-pc_type -pc_hypre_type -ksp_gmres_restart'");
	this->moose_input.getDocument("Executioner").addPair("petsc_options_value","'hypre boomeramg 100'");
	
	this->moose_input.getDocument("Executioner").addHeadKey("TimeStepper");
	
	this->moose_input.addDocKey("Preconditioning");
	
	this->moose_input.addDocKey("Outputs");
	this->moose_input.getDocument("Outputs").addPair("exodus","true");
	this->moose_input.getDocument("Outputs").addPair("csv","true");
	this->moose_input.getDocument("Outputs").addPair("print_linear_residuals","true");
}

//Function to create a MOOSE example input file
void SimpleUI::createExample()
{
	this->moose_input.getDocument("GlobalParams").addPair("initial_dt","0.01");
	this->moose_input.getDocument("GlobalParams").addPair("length","22.86");
	
	this->moose_input.getDocument("Mesh").addPair("xmax","0.8636");
	this->moose_input.getDocument("Mesh").addPair("ymax","22.86");
	
	this->moose_input.getDocument("Variables").addHeadKey("Kr");
	this->moose_input.getDocument("Variables").getHeader("Kr").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("Kr").addPair("family","MONOMIAL");
	
	this->moose_input.getDocument("Variables").addHeadKey("Xe");
	this->moose_input.getDocument("Variables").getHeader("Xe").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("Xe").addPair("family","MONOMIAL");
	
	this->moose_input.getDocument("Variables").addHeadKey("He");
	this->moose_input.getDocument("Variables").getHeader("He").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("He").addPair("family","MONOMIAL");
	
	this->moose_input.getDocument("Variables").getHeader("wall_temp").addPair("initial_condition","253.15");
	this->moose_input.getDocument("Variables").getHeader("column_temp").addPair("initial_condition","253.15");
	
	this->moose_input.getDocument("AuxVariables").getHeader("ambient_temp").addPair("initial_condition","253.15");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("He_Adsorbed");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Adsorbed").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Adsorbed").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Adsorbed").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Kr_Adsorbed");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Adsorbed").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Adsorbed").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Adsorbed").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Xe_Adsorbed");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Adsorbed").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Adsorbed").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Adsorbed").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("He_Perturb");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Perturb").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Perturb").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Perturb").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Kr_Perturb");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Perturb").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Perturb").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Perturb").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Xe_Perturb");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Perturb").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Perturb").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Perturb").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("He_AdsorbedHeat");
	this->moose_input.getDocument("AuxVariables").getHeader("He_AdsorbedHeat").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("He_AdsorbedHeat").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("He_AdsorbedHeat").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Kr_AdsorbedHeat");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_AdsorbedHeat").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_AdsorbedHeat").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_AdsorbedHeat").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Xe_AdsorbedHeat");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_AdsorbedHeat").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_AdsorbedHeat").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_AdsorbedHeat").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("ICs").addHeadKey("Kr_IC");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("type","ConcentrationIC");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("variable","Kr");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("initial_mole_frac","0.0");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("initial_press","101.35");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("initial_temp","253.15");
	
	this->moose_input.getDocument("ICs").addHeadKey("Xe_IC");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("type","ConcentrationIC");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("variable","Xe");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("initial_mole_frac","0.0");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("initial_press","101.35");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("initial_temp","253.15");
	
	this->moose_input.getDocument("ICs").addHeadKey("He_IC");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("type","ConcentrationIC");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("variable","He");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("initial_mole_frac","1.0");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("initial_press","101.35");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("initial_temp","253.15");
	
	this->moose_input.getDocument("Kernels").addHeadKey("Kr_accum");
	this->moose_input.getDocument("Kernels").getHeader("Kr_accum").addPair("type","BedMassAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("Kr_accum").addPair("variable","Kr");
	this->moose_input.getDocument("Kernels").getHeader("Kr_accum").addPair("index","0");
	this->moose_input.getDocument("Kernels").addHeadKey("Kr_masstrans");
	this->moose_input.getDocument("Kernels").getHeader("Kr_masstrans").addPair("type","AdsorptionMassTransfer");
	this->moose_input.getDocument("Kernels").getHeader("Kr_masstrans").addPair("variable","Kr");
	this->moose_input.getDocument("Kernels").getHeader("Kr_masstrans").addPair("solid_conc","Kr_Adsorbed");
	this->moose_input.getDocument("Kernels").addHeadKey("Kr_diff");
	this->moose_input.getDocument("Kernels").getHeader("Kr_diff").addPair("type","GColumnMassDispersion");
	this->moose_input.getDocument("Kernels").getHeader("Kr_diff").addPair("variable","Kr");
	this->moose_input.getDocument("Kernels").getHeader("Kr_diff").addPair("index","0");
	this->moose_input.getDocument("Kernels").addHeadKey("Kr_adv");
	this->moose_input.getDocument("Kernels").getHeader("Kr_adv").addPair("type","GColumnMassAdvection");
	this->moose_input.getDocument("Kernels").getHeader("Kr_adv").addPair("variable","Kr");
	
	this->moose_input.getDocument("Kernels").addHeadKey("Xe_accum");
	this->moose_input.getDocument("Kernels").getHeader("Xe_accum").addPair("type","BedMassAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("Xe_accum").addPair("variable","Xe");
	this->moose_input.getDocument("Kernels").getHeader("Xe_accum").addPair("index","1");
	this->moose_input.getDocument("Kernels").addHeadKey("Xe_masstrans");
	this->moose_input.getDocument("Kernels").getHeader("Xe_masstrans").addPair("type","AdsorptionMassTransfer");
	this->moose_input.getDocument("Kernels").getHeader("Xe_masstrans").addPair("variable","Xe");
	this->moose_input.getDocument("Kernels").getHeader("Xe_masstrans").addPair("solid_conc","Xe_Adsorbed");
	this->moose_input.getDocument("Kernels").addHeadKey("Xe_diff");
	this->moose_input.getDocument("Kernels").getHeader("Xe_diff").addPair("type","GColumnMassDispersion");
	this->moose_input.getDocument("Kernels").getHeader("Xe_diff").addPair("variable","Xe");
	this->moose_input.getDocument("Kernels").getHeader("Xe_diff").addPair("index","1");
	this->moose_input.getDocument("Kernels").addHeadKey("Xe_adv");
	this->moose_input.getDocument("Kernels").getHeader("Xe_adv").addPair("type","GColumnMassAdvection");
	this->moose_input.getDocument("Kernels").getHeader("Xe_adv").addPair("variable","Xe");
	
	this->moose_input.getDocument("Kernels").addHeadKey("He_accum");
	this->moose_input.getDocument("Kernels").getHeader("He_accum").addPair("type","BedMassAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("He_accum").addPair("variable","He");
	this->moose_input.getDocument("Kernels").getHeader("He_accum").addPair("index","2");
	this->moose_input.getDocument("Kernels").addHeadKey("He_masstrans");
	this->moose_input.getDocument("Kernels").getHeader("He_masstrans").addPair("type","AdsorptionMassTransfer");
	this->moose_input.getDocument("Kernels").getHeader("He_masstrans").addPair("variable","He");
	this->moose_input.getDocument("Kernels").getHeader("He_masstrans").addPair("solid_conc","He_Adsorbed");
	this->moose_input.getDocument("Kernels").addHeadKey("He_diff");
	this->moose_input.getDocument("Kernels").getHeader("He_diff").addPair("type","GColumnMassDispersion");
	this->moose_input.getDocument("Kernels").getHeader("He_diff").addPair("variable","He");
	this->moose_input.getDocument("Kernels").getHeader("He_diff").addPair("index","2");
	this->moose_input.getDocument("Kernels").addHeadKey("He_adv");
	this->moose_input.getDocument("Kernels").getHeader("He_adv").addPair("type","GColumnMassAdvection");
	this->moose_input.getDocument("Kernels").getHeader("He_adv").addPair("variable","He");
	
	this->moose_input.getDocument("Kernels").addHeadKey("column_AdsHeat");
	this->moose_input.getDocument("Kernels").getHeader("column_AdsHeat").addPair("type","AdsorptionHeatAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("column_AdsHeat").addPair("variable","column_temp");
	this->moose_input.getDocument("Kernels").getHeader("column_AdsHeat").addPair("solid_heats","'Kr_AdsorbedHeat Xe_AdsorbedHeat He_AdsorbedHeat'");
	
	this->moose_input.getDocument("DGKernels").addHeadKey("Kr_DGdiff");
	this->moose_input.getDocument("DGKernels").getHeader("Kr_DGdiff").addPair("type","DGColumnMassDispersion");
	this->moose_input.getDocument("DGKernels").getHeader("Kr_DGdiff").addPair("variable","Kr");
	this->moose_input.getDocument("DGKernels").getHeader("Kr_DGdiff").addPair("index","0");
	this->moose_input.getDocument("DGKernels").addHeadKey("Kr_DGadv");
	this->moose_input.getDocument("DGKernels").getHeader("Kr_DGadv").addPair("type","DGColumnMassAdvection");
	this->moose_input.getDocument("DGKernels").getHeader("Kr_DGadv").addPair("variable","Kr");
	
	this->moose_input.getDocument("DGKernels").addHeadKey("Xe_DGdiff");
	this->moose_input.getDocument("DGKernels").getHeader("Xe_DGdiff").addPair("type","DGColumnMassDispersion");
	this->moose_input.getDocument("DGKernels").getHeader("Xe_DGdiff").addPair("variable","Xe");
	this->moose_input.getDocument("DGKernels").getHeader("Xe_DGdiff").addPair("index","1");
	this->moose_input.getDocument("DGKernels").addHeadKey("Xe_DGadv");
	this->moose_input.getDocument("DGKernels").getHeader("Xe_DGadv").addPair("type","DGColumnMassAdvection");
	this->moose_input.getDocument("DGKernels").getHeader("Xe_DGadv").addPair("variable","Xe");
	
	this->moose_input.getDocument("DGKernels").addHeadKey("He_DGdiff");
	this->moose_input.getDocument("DGKernels").getHeader("He_DGdiff").addPair("type","DGColumnMassDispersion");
	this->moose_input.getDocument("DGKernels").getHeader("He_DGdiff").addPair("variable","He");
	this->moose_input.getDocument("DGKernels").getHeader("He_DGdiff").addPair("index","2");
	this->moose_input.getDocument("DGKernels").addHeadKey("He_DGadv");
	this->moose_input.getDocument("DGKernels").getHeader("He_DGadv").addPair("type","DGColumnMassAdvection");
	this->moose_input.getDocument("DGKernels").getHeader("He_DGadv").addPair("variable","He");
	
	this->moose_input.getDocument("AuxKernels").getHeader("column_pressure").addPair("coupled_gases","'Kr Xe He'");
	
	this->moose_input.getDocument("AuxKernels").addHeadKey("Kr_ads");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_ads").addPair("type","Scopsowl_Adsorption");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_ads").addPair("variable","Kr_Adsorbed");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_ads").addPair("index","0");
	this->moose_input.getDocument("AuxKernels").addHeadKey("Xe_ads");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_ads").addPair("type","Scopsowl_Adsorption");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_ads").addPair("variable","Xe_Adsorbed");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_ads").addPair("index","1");
	this->moose_input.getDocument("AuxKernels").addHeadKey("He_ads");
	this->moose_input.getDocument("AuxKernels").getHeader("He_ads").addPair("type","Scopsowl_Adsorption");
	this->moose_input.getDocument("AuxKernels").getHeader("He_ads").addPair("variable","He_Adsorbed");
	this->moose_input.getDocument("AuxKernels").getHeader("He_ads").addPair("index","2");
	
	this->moose_input.getDocument("AuxKernels").addHeadKey("Kr_pert");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_pert").addPair("type","Scopsowl_Perturbation");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_pert").addPair("variable","Kr_Perturb");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_pert").addPair("index","0");
	this->moose_input.getDocument("AuxKernels").addHeadKey("Xe_pert");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_pert").addPair("type","Scopsowl_Perturbation");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_pert").addPair("variable","Xe_Perturb");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_pert").addPair("index","1");
	this->moose_input.getDocument("AuxKernels").addHeadKey("He_pert");
	this->moose_input.getDocument("AuxKernels").getHeader("He_pert").addPair("type","Scopsowl_Perturbation");
	this->moose_input.getDocument("AuxKernels").getHeader("He_pert").addPair("variable","He_Perturb");
	this->moose_input.getDocument("AuxKernels").getHeader("He_pert").addPair("index","2");
	
	this->moose_input.getDocument("AuxKernels").addHeadKey("Kr_ads_heat");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_ads_heat").addPair("type","MAGPIE_AdsorptionHeat");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_ads_heat").addPair("variable","Kr_AdsorbedHeat");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_ads_heat").addPair("solid_conc","Kr_Adsorbed");
	this->moose_input.getDocument("AuxKernels").getHeader("Kr_ads_heat").addPair("index","0");
	this->moose_input.getDocument("AuxKernels").addHeadKey("Xe_ads_heat");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_ads_heat").addPair("type","MAGPIE_AdsorptionHeat");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_ads_heat").addPair("variable","Xe_AdsorbedHeat");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_ads_heat").addPair("solid_conc","Xe_Adsorbed");
	this->moose_input.getDocument("AuxKernels").getHeader("Xe_ads_heat").addPair("index","1");
	this->moose_input.getDocument("AuxKernels").addHeadKey("He_ads_heat");
	this->moose_input.getDocument("AuxKernels").getHeader("He_ads_heat").addPair("type","MAGPIE_AdsorptionHeat");
	this->moose_input.getDocument("AuxKernels").getHeader("He_ads_heat").addPair("variable","He_AdsorbedHeat");
	this->moose_input.getDocument("AuxKernels").getHeader("He_ads_heat").addPair("solid_conc","He_Adsorbed");
	this->moose_input.getDocument("AuxKernels").getHeader("He_ads_heat").addPair("index","2");

	this->moose_input.getDocument("BCs").getHeader("Heat_Gas_Flux").addPair("type","DGHeatFluxLimitedBC");
	this->moose_input.getDocument("BCs").getHeader("Heat_Gas_Flux").addPair("input_temperature","253.15");
	this->moose_input.getDocument("BCs").getHeader("Heat_Wall_Flux").addPair("type","DGColumnWallHeatFluxLimitedBC");
	
	this->moose_input.getDocument("BCs").addHeadKey("Kr_Flux");
	this->moose_input.getDocument("BCs").getHeader("Kr_Flux").addPair("type","DGMassFluxLimitedBC");
	this->moose_input.getDocument("BCs").getHeader("Kr_Flux").addPair("variable","Kr");
	this->moose_input.getDocument("BCs").getHeader("Kr_Flux").addPair("boundary","'top bottom'");
	this->moose_input.getDocument("BCs").getHeader("Kr_Flux").addPair("input_temperature","253.15");
	this->moose_input.getDocument("BCs").getHeader("Kr_Flux").addPair("input_pressure","101.35");
	this->moose_input.getDocument("BCs").getHeader("Kr_Flux").addPair("input_molefraction","0.000131792");
	this->moose_input.getDocument("BCs").getHeader("Kr_Flux").addPair("index","0");
	
	this->moose_input.getDocument("BCs").addHeadKey("Xe_Flux");
	this->moose_input.getDocument("BCs").getHeader("Xe_Flux").addPair("type","DGMassFluxLimitedBC");
	this->moose_input.getDocument("BCs").getHeader("Xe_Flux").addPair("variable","Xe");
	this->moose_input.getDocument("BCs").getHeader("Xe_Flux").addPair("boundary","'top bottom'");
	this->moose_input.getDocument("BCs").getHeader("Xe_Flux").addPair("input_temperature","253.15");
	this->moose_input.getDocument("BCs").getHeader("Xe_Flux").addPair("input_pressure","101.35");
	this->moose_input.getDocument("BCs").getHeader("Xe_Flux").addPair("input_molefraction","0.000863107");
	this->moose_input.getDocument("BCs").getHeader("Xe_Flux").addPair("index","1");
	
	this->moose_input.getDocument("BCs").addHeadKey("He_Flux");
	this->moose_input.getDocument("BCs").getHeader("He_Flux").addPair("type","DGMassFluxLimitedBC");
	this->moose_input.getDocument("BCs").getHeader("He_Flux").addPair("variable","He");
	this->moose_input.getDocument("BCs").getHeader("He_Flux").addPair("boundary","'top bottom'");
	this->moose_input.getDocument("BCs").getHeader("He_Flux").addPair("input_temperature","253.15");
	this->moose_input.getDocument("BCs").getHeader("He_Flux").addPair("input_pressure","101.35");
	this->moose_input.getDocument("BCs").getHeader("He_Flux").addPair("input_molefraction","0.999005101");
	this->moose_input.getDocument("BCs").getHeader("He_Flux").addPair("index","2");
	
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("inner_diameter","1.7272");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("outer_diameter","1.905");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("bulk_porosity","0.885");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("axial_conductivity","6.292E-05");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("wall_density","7.7");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("wall_heat_capacity","0.5");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("wall_heat_trans_coef","9.0");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("extern_heat_trans_coef","90.0");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("coupled_gases","'Kr Xe He'");
	
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("molecular_weight","'83.8 131.29 4.0026'");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("comp_heat_capacity","'0.25 0.16 5.1916'");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("comp_ref_viscosity","'0.00023219 0.00021216 0.0001885'");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("comp_ref_temp","'273.15 273.15 273.15'");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("comp_Sutherland_const","'266.505 232.746 80.0'");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("flow_rate","2994.06");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("coupled_gases","'Kr Xe He'");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("coupled_adsorption","'Kr_Adsorbed Xe_Adsorbed He_Adsorbed'");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("coupled_perturbation","'Kr_Perturb Xe_Perturb He_Perturb'");
	
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("binder_porosity","0.384");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("pellet_diameter","0.16");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("macropore_radius","1.5e-4");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("pellet_density","3.06");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("pellet_heat_capacity","1.045");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("ref_diffusion","'0 0 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("activation_energy","'0 0 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("ref_temperature","'0 0 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("affinity","'0 0 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("coupled_gases","'Kr Xe He'");
	
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("coupled_gases","'Kr Xe He'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("number_sites","'2 3 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("maximum_capacity","'1.716 1.479 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("molar_volume","'20.785 25.412 0'");
	
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("enthalpy_site_1","'-44696.86 -18455.18 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("enthalpy_site_2","'-65465.52 -35511.74 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("enthalpy_site_3","'0 -53315.13 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("enthalpy_site_4","'0 0 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("enthalpy_site_5","'0 0 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("enthalpy_site_6","'0 0 0'");
	
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("entropy_site_1","'-170.45 -23.25 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("entropy_site_2","'-248.55 -62.45 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("entropy_site_3","'0 -100.10 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("entropy_site_4","'0 0 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("entropy_site_5","'0 0 0'");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("entropy_site_6","'0 0 0'");
	
	this->moose_input.getDocument("Materials").addHeadKey("KineticMaterials");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("type","ScopsowlProperties");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("dirichlet_bc","false");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("heterogeneous","false");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("surface_diffusion","true");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("macro_spheres","false");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("macro_length","0.4");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("coupled_gases","'Kr Xe He'");
	this->moose_input.getDocument("Materials").getHeader("KineticMaterials").addPair("coupled_adsorption","'Kr_Adsorbed Xe_Adsorbed He_Adsorbed'");
	
	this->moose_input.getDocument("Postprocessors").addHeadKey("Kr_exit");
	this->moose_input.getDocument("Postprocessors").getHeader("Kr_exit").addPair("type","SideAverageValue");
	this->moose_input.getDocument("Postprocessors").getHeader("Kr_exit").addPair("boundary","'top'");
	this->moose_input.getDocument("Postprocessors").getHeader("Kr_exit").addPair("variable","Kr");
	this->moose_input.getDocument("Postprocessors").getHeader("Kr_exit").addPair("execute_on","'initial timestep_end'");
	
	this->moose_input.getDocument("Postprocessors").addHeadKey("Xe_exit");
	this->moose_input.getDocument("Postprocessors").getHeader("Xe_exit").addPair("type","SideAverageValue");
	this->moose_input.getDocument("Postprocessors").getHeader("Xe_exit").addPair("boundary","'top'");
	this->moose_input.getDocument("Postprocessors").getHeader("Xe_exit").addPair("variable","Xe");
	this->moose_input.getDocument("Postprocessors").getHeader("Xe_exit").addPair("execute_on","'initial timestep_end'");
	
	this->moose_input.getDocument("Postprocessors").addHeadKey("He_exit");
	this->moose_input.getDocument("Postprocessors").getHeader("He_exit").addPair("type","SideAverageValue");
	this->moose_input.getDocument("Postprocessors").getHeader("He_exit").addPair("boundary","'top'");
	this->moose_input.getDocument("Postprocessors").getHeader("He_exit").addPair("variable","He");
	this->moose_input.getDocument("Postprocessors").getHeader("He_exit").addPair("execute_on","'initial timestep_end'");
	
	this->moose_input.getDocument("Postprocessors").addHeadKey("temp_exit");
	this->moose_input.getDocument("Postprocessors").getHeader("temp_exit").addPair("type","SideAverageValue");
	this->moose_input.getDocument("Postprocessors").getHeader("temp_exit").addPair("boundary","'top'");
	this->moose_input.getDocument("Postprocessors").getHeader("temp_exit").addPair("variable","column_temp");
	this->moose_input.getDocument("Postprocessors").getHeader("temp_exit").addPair("execute_on","'initial timestep_end'");
	
	this->moose_input.getDocument("Postprocessors").addHeadKey("pressure_exit");
	this->moose_input.getDocument("Postprocessors").getHeader("pressure_exit").addPair("type","SideAverageValue");
	this->moose_input.getDocument("Postprocessors").getHeader("pressure_exit").addPair("boundary","'top'");
	this->moose_input.getDocument("Postprocessors").getHeader("pressure_exit").addPair("variable","total_pressure");
	this->moose_input.getDocument("Postprocessors").getHeader("pressure_exit").addPair("execute_on","'initial timestep_end'");
	
	this->moose_input.getDocument("Postprocessors").addHeadKey("Kr_adsorbed");
	this->moose_input.getDocument("Postprocessors").getHeader("Kr_adsorbed").addPair("type","ElementAverageValue");
	this->moose_input.getDocument("Postprocessors").getHeader("Kr_adsorbed").addPair("variable","Kr_Adsorbed");
	this->moose_input.getDocument("Postprocessors").getHeader("Kr_adsorbed").addPair("execute_on","'initial timestep_end'");
	
	this->moose_input.getDocument("Postprocessors").addHeadKey("Xe_adsorbed");
	this->moose_input.getDocument("Postprocessors").getHeader("Xe_adsorbed").addPair("type","ElementAverageValue");
	this->moose_input.getDocument("Postprocessors").getHeader("Xe_adsorbed").addPair("variable","Xe_Adsorbed");
	this->moose_input.getDocument("Postprocessors").getHeader("Xe_adsorbed").addPair("execute_on","'initial timestep_end'");
	
	this->moose_input.getDocument("Executioner").addPair("end_time","50.0");
	this->moose_input.getDocument("Executioner").addPair("dtmax","1.0");
	
	this->moose_input.getDocument("Executioner").getHeader("TimeStepper").addPair("type","SolutionTimeAdaptiveDT");
	this->moose_input.getDocument("Executioner").getHeader("TimeStepper").addPair("dt","0.01");
}

//Return reference to the YamlWrapper object for the input file
YamlWrapper& SimpleUI::getYamlInput()
{
	return this->yaml_input.getYamlWrapper();
}

//Return reference to the YamlWrapper object to create the MOOSE input file
YamlWrapper& SimpleUI::getMooseInput()
{
	return this->moose_input;
}

//Function to display input file contents
void SimpleUI::DisplayInput()
{
	this->yaml_input.DisplayContents();
}

//Function to display MOOSE input file contents that we are making
void SimpleUI::DisplayOutput()
{
	this->moose_input.DisplayContents();
}
