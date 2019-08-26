/*!
 *  \file skua.h skua.cpp
 *	\brief Surface Kinetics for Uptake by Adsorption
 *	\details This file contains structures and functions associated with solving the surface diffusion
 *			partial differential equations for adsorption kinetics in spherical and/or cylindrical
 *			adsorbents. For this system, it is assumed that the pore size is so small that all molecules
 *			are confined to movement exclusively on the surface area of the adsorbent. The total amount
 *			of adsorption for each species is drive by the MAGPIE model for non-ideal mixed gas adsorption.
 *			Spatial and temporal varience in adsorption is caused by a combination of different kinetics
 *			between adsorbing species and different adsorption affinities for the surface.
 *
 *			The function for surface diffusion involves four parameters, although not all of these parameters
 *			are required to be used. Surface diffusion theoretically varies with temperature according to the
 *			Arrhenius rate expression, but we also add in an empirical correction term to account for variations
 *			in diffusivity with the partial pressure of the species in the gas phase. \n
 *
 *			D_surf = D_ref * exp(-E / (R*T) ) * pow(p , (T_ref/T) - B ) \n
 *
 *			D_ref is the Reference Diffusivity (um^2/hr), E is the activation energy for adsorption
 *			(J/mol), R is the gas law constant (J/K/mol), T is the system temperature (K), p is the
 *			partial pressure of the adsorbing species (kPa), T_ref is the Reference Temperature (K), and
 *			B is the Affinity constant. \n
 *
 *  \author Austin Ladshaw
 *	\date 01/26/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#pragma once

#include "finch.h"				//FINCH handles the physics solver and discretization of the problem
#include "magpie.h"				//MAGPIE handles the adsorption equilibria equations
#include "egret.h"				//EGRET handles the parameter estimation for gas phase properties

#ifndef SKUA_HPP_
#define SKUA_HPP_

#ifndef D_inf
#define D_inf(Dref,Tref,B,p,T) ( Dref * pow(p+sqrt(DBL_EPSILON),(Tref/T)-B) ) ///< Empirical correction of diffusivity (um^2/hr)
#endif

#ifndef D_o
#define D_o(Diff,E,T) ( Diff * exp(-E/(Rstd*T)) ) ///< Arrhenius Rate Expression for Diffusivity (um^2/hr)
#endif

#ifndef D_c
#define D_c(Diff,phi) ( Diff * (1.0/((1.0+1.1E-6)-phi) ) ) ///< Approximate Darken Diffusivity Equation (um^2/hr)
#endif

/// Data structure for species' parameters in SKUA
/** C-style object holding data and parameters associated with the gas/solid species in the overall SKUA system.
	These parameters are used in to modify surface diffusivity with temperature, establish film mass transfer
	coefficients, formulate the initial conditions, and store solution results for heat of adsorption and adsorbed
	mole fractions. One of these objects will be created for each species in the gas system.*/
typedef struct
{
	double activation_energy;				//Activation energy for surface diffusion (J/mol)
	double ref_diffusion;					//Reference state diffusivity (um^2/hr)
	double ref_temperature;					//Reference temperature for empirical adjustments (K)
	double affinity;						//Affinity parameter used in empirical adjustments (-)
	double ref_pressure;					//Reference pressure for empirical adjustments (kPa)
	double film_transfer;					//Film mass transfer coeff (um/hr)
	
	double xIC;								//Inside initial mole fractions of each component (-)
	double y_eff;							//Effective interior gas mole fraction based on adsorption (-)
	
	double Qstn;							//Old heat of adsorption (J/mol)
	double Qstnp1;							//New heat of adsorption (J/mol)
	double xn;								//Old adsorbed mole fraction (-)
	double xnp1;							//New adsorbed mole fraction (-)
	
	bool Adsorbable;						//Boolean to identify with components are adsorbable
	
	std::string speciesName;				//String to hold the name of each species
}SKUA_PARAM;

/// Data structure for all simulation information in SKUA
/** C-style object holding all data, functions, and other objects needed to successfully run a SKUA simulation.
	This object holds system information, such as boundary condition type, adsorbent size, and total adsorption,
	and also contains structure for EGRET (egret.h), FINCH (finch.h), and MAGPIE (magpie.h) calculations. Function
	pointers for evaluation of the surface diffusivity and film mass transfer coefficients can be overriden by the
	user to change the behavior of the SKUA simulation. However, defaults are also provided for these functions.*/
typedef struct
{
	unsigned long int total_steps;			///< Running total of all calculation steps
	int coord;								///< Used to determine the coordinates of the problem
	double sim_time;						///< Stopping time for the simulation (hrs)
	double t_old;							///< Old time of the simulations (hrs)
	double t;								///< Current time of the simulations (hrs)
	double t_counter = 0.0;					///< Counts for print times for output (hrs)
	double t_print;							///< Prints out every t_print time (hrs)
	double qTn;								///< Old total amounts adsorbed (mol/kg)
	double qTnp1;							///< New total amounts adsorbed (mol/kg)
	bool Print2File = true;					///< True = results to .txt; False = no printing
	bool Print2Console = true;				///< True = results to console; False = no printing
	
	double gas_velocity;					///< Superficial Gas Velocity arount pellet (cm/s)
	double pellet_radius;					///< Nominal radius of the pellet/crystal (um)
	double char_measure;					///< Length or Area if in Cylindrical or Cartesian coordinates (um or um^2)
	bool DirichletBC = true;				///< True = Dirichlet BC; False = Neumann BC
	bool NonLinear = true;					///< True = Non-linear solver; False = Linear solver
	std::vector<double> y;					///< Outside mole fractions of each component (-)
	
	FILE *OutputFile;						///< Output file pointer to the output file
	double (*eval_diff) (int i, int l, const void *user_data);	///< Function pointer for evaluating surface diffusivity
	double (*eval_kf) (int i, const void *user_data);			///< Function pointer for evaluating film mass transfer
	const void *user_data;					///< Data structure for user's information needed in parameter functions
	MAGPIE_DATA magpie_dat;					///< Data structure for adsorption equilibria (see magpie.h)
	MIXED_GAS *gas_dat;						///< Pointer to the MIXED_GAS data structure (see egret.h)
	std::vector<FINCH_DATA> finch_dat;		///< Data structure for adsorption kinetics (see finch.h)
	std::vector<SKUA_PARAM> param_dat;		///< Data structure for SKUA specific parameters
}SKUA_DATA;

/// Function to print out the species' headers to output file
void print2file_species_header(FILE *Output, SKUA_DATA *skua_dat, int i);

/// Function to print out time and space headers to output file
void print2file_SKUA_time_header(FILE *Output, SKUA_DATA *skua_dat, int i);

/// Function calls the other header functions to establish output file structure
void print2file_SKUA_header(SKUA_DATA *skua_dat);

/// Function to print out the old time step simulation results to the output file
void print2file_SKUA_results_old(SKUA_DATA *skua_dat);

/// Function to print out the new time step simulation results to the output file
void print2file_SKUA_results_new(SKUA_DATA *skua_dat);

//--------- Default Parameter Functions -------------

/// Default function for surface diffusivity
/** This is the default function provided by SKUA for the calculation of the surface diffusivity
	parameter. The diffusivity is calculated based on the Arrhenius rate expression, then corrected
	for using the empirical correction term with the outside partial pressure of the gas species.
 
	\param i index of the gas/adsorbed phase species that this function acts on
	\param l index of the node in the spatial discretization that this function acts on
	\param data pointer to the SKUA_DATA structure*/
double default_Dc(int i, int l, const void *data);

/// Default function for film mass transfer coefficent
/** This is the default function provided by SKUA for the calculation of the film mass transfer
	parameter. By default, we are usually going to couple the SKUA model with a pore diffusion
	model (see scopsowl.h). Therefore, the film mass transfer coefficient would be zero, because
	we would only consider a Dirichlet boundary condition for this sub-problem.
 
	\param i index of the gas/adsorbed phase species that this function acts on
	\param data pointer to the SKUA_DATA structure*/
double default_kf(int i, const void *data);
//---------------------------------------------------

//-------- Optional Parameter Functions -----------------

/// Constant surface diffusivity function
/** This function allows the user to specify just a single constant value for surface diffusivity. The
	value of diffusivity applied at all nodes will be the ref_diffusion parameter in SKUA_PARAM.
 
	\param i index of the gas/adsorbed phase species that this function acts on
	\param l index of the node in the spatial discretization that this function acts on
	\param data pointer to the SKUA_DATA structure*/
double const_Dc(int i, int l, const void *data);

/// Simple Darken model for surface diffusivity
/** This function uses an approximation to Darken's model for surface diffusion. The approximation is
	exact if the isotherm for adsorption takes the form of the Langmuir model, but is only approximate
	if the isotherm is heterogeneous. Forming the approximation in this manner is significantly cheaper
	than forming the true Darken model expression for the GSTA isotherm.
 
	\param i index of the gas/adsorbed phase species that this function acts on
	\param l index of the node in the spatial discretization that this function acts on
	\param data pointer to the SKUA_DATA structure*/
double simple_darken_Dc(int i, int l, const void *data);

/// Theoretical Darken model for surface diffusivity
/** This function uses the full theoretical expression of the Darken's diffusion model to calculate
	the surface diffusivity. This calculation involves formulating the reference state pressures for
	the adsorbed amount at every node, then calculating derivatives of the adsorption isotherm for
	each species. It is more accurate than the simple Darken model function, but costs significantly
	more computational time.
 
	\param i index of the gas/adsorbed phase species that this function acts on
	\param l index of the node in the spatial discretization that this function acts on
	\param data pointer to the SKUA_DATA structure*/
double theoretical_darken_Dc(int i, int l, const void *data);

/// Empirical function for film mass transfer coefficent
/** This function provides an empirical estimate of the mass transfer coefficient using the gas
	velocity, molecular diffusivities, and dimensionless numbers (see egret.h). It is used as
	the default film mass transfer function IF the boundary condition is specified to be a
	Neumann type boundary by the user.
 
	\param i index of the gas/adsorbed phase species that this function acts on
	\param data pointer to the SKUA_DATA structure*/
double empirical_kf(int i, const void *data);

/// Constant function for film mass transfer coefficent
/** This function allows the user to specify a constant value for the film mass transfer coefficient.
	The value of the film mass transfer coefficient will be the value of film_transfer given in the
	SKUA_PARAM data structure.
 
	\param i index of the gas/adsorbed phase species that this function acts on
	\param data pointer to the SKUA_DATA structure*/
double const_kf(int i, const void *data);
//---------------------------------------------------

/// Function to check mole fractions in gas and solid phases for errors
/** This function is called after reading input and before calling the primary solution routines. It will
	force and error and quit the program if their are inconsistencies in the mole fractions it was given.
	All mole fractions must sum to 1, otherwise there is missing information.*/
int molefractionCheck(SKUA_DATA *skua_dat);

/// Function to setup the function pointers and vector objects in memory to setup the SKUA simulation
/** This function is called to setup the SKUA problem in memory and set function pointers to either
	defaults or user specified functions. It must be called prior to calling any other SKUA function and
	will report an error if the object was not setup properly.
 
	\param file pointer to the output file for SKUA simulations
	\param eval_Dc pointer to the function to evaluate the surface diffusivity
	\param eval_Kf pointer to the function to evaluate the film mass transfer coefficient
	\param user_data pointer to a user defined data structure used in the calculation the the parameters
	\param gas_data pointer to the MIXED_GAS data structure for egret.h calculations
	\param skua_dat pointer to the SKUA_DATA data structure */
int setup_SKUA_DATA(FILE *file, double (*eval_Dc) (int i, int l, const void *user_data),
					double (*eval_Kf) (int i, const void *user_data), const void *user_data,
					MIXED_GAS *gas_data, SKUA_DATA *skua_dat);

/// Function to execute preprocesses, solvers, and postprocesses for a SKUA simulation
/** This function calls the preprocess, solver, and postprocess functions to complete a single
	time step in a SKUA simulation. User's will want to call this function whenever a time step
	simulation result is needed. This is used primarily when coupling with other models (see
	scopsowl.h). */
int SKUA_Executioner(SKUA_DATA *skua_dat);

/// Function to establish the initial conditions of adsorption in the adsorbent
/** This function needs to be called before doing any simulation or execution of a time step, but
	only once per simulation. It sets the value of adsorption for each adsorbable species to the
	specified initial values given via qT and xIC in SKUA_DATA. */
int set_SKUA_ICs(SKUA_DATA *skua_dat);

/// Function to establish the time step for the current simulation
/** This function is called to set a time step value for a particular simulation step. By default,
	the time step is set to (1/4)x space step size. If you need to change the step size, you must
	do so manually. */
int set_SKUA_timestep(SKUA_DATA *skua_dat);

/// Function to perform the necessary preprocess operations before a solve
/** This function performs preprocess operations prior to calling the solver routine. Those preprocesses
	include establishing boundary conditions and performing a MAGPIE simulation for the adsorption on the
	surface (see magpie.h). */
int SKUA_preprocesses(SKUA_DATA *skua_dat);

/// Function to call the diffusivity function during the solve
/** This is the function passed into FINCH to be called during the FINCH solver (see finch.h). It will
	call the diffusion functions set by the user in the setup function above. This is not overridable. */
int set_SKUA_params(const void *user_data);

/// Function to perform the necessary postprocess operations after a solve
/** This function performs postprocess operations after a solve was completed successfully. Those
	operations include estimating average total adsorption, average adsorbed mole fractions, and
	heat of adsorption for each species. Results are then printed to the output file. */
int SKUA_postprocesses(SKUA_DATA *skua_dat);

/// Function to reset the stateful information in SKUA after a simulation
/** This function sets all the old state data to the newly formed state data. It needs to be called
	after a successful execution of the simulation step and before calling for the next time step
	to be solved. Do not call out of turn, otherwise information will be lost. */
int SKUA_reset(SKUA_DATA *skua_dat);

/// Function to iteratively call all execution steps to evolve a simulation through time
/** This function is used in conjunction with the scenario call from the UI to numerically solve
	the adsorption kinetics problem in time. It will call the initial conditions function once,
	then iteratively call the reset, time step, and executioner functions for SKUA to push the
	simulation forward in time. This function will be called from the SKUA_SCENARIOS function.*/
int SKUA(SKUA_DATA *skua_dat);

//-------- Running Specific Tests ------------
/// \cond

int SKUA_CYCLE_TEST01(SKUA_DATA *skua_dat);

int SKUA_CYCLE_TEST02(SKUA_DATA *skua_dat);

int SKUA_LOW_TEST03(SKUA_DATA *skua_dat);

int SKUA_MID_TEST04(SKUA_DATA *skua_dat);

/// \endcond
//--------------------------------------------

/// Function callable from the UI to perform a SKUA simulation based on user supplied input files
/** This is the primary function to be called when running a stand-alone SKUA simulation. Parameters
	and system information for the simulation are given in a series of input files that come in as
	character arrays. These inputs are all required to call this function.
 
	\param scene Sceneario Input File
	\param sorbent Adsorbent Input File
	\param comp Component Input File
	\param sorbate Adsorbate Input File*/

/** \note Each input file has a particular format that must be strictly adhered to in order for the
	simulation to be carried out correctly. The format for each input file, and an example, is provided
	below...
 
	Scenario Input Format
	---------------------
	System Temperature (K) [tab] Total Pressure (kPa) [tab] Gas Velocity (cm/s) \n
	Simulation Time (hrs) [tab] Print Out Time (hrs) \n
	BC Type (0 = Neumann, 1 = Dirichlet) \n
	Number of Gas Species \n
	Initial Total Adsorption (mol/kg) \n
	Name of ith Species [tab] Adsorbable? (0 = false, 1 = true) [tab] Gas Phase Molefraction [tab] Initial Sorbed Molefraction \n
	(repeat above for all species) \n
 
	Example Scenario Input
	----------------------
	353.15	101.35	0.36 \n
	4.0		0.05 \n
	0 \n
	5 \n
	0.0 \n
	N2	0	0.7634	0.0 \n
	O2	0	0.2081	0.0 \n
	Ar	0	0.009	0.0 \n
	CO2	0	0.0004	0.0 \n
	H2O	1	0.0191	0.0 \n
 
	Above example is for a 5-component mixture of N2, O2, Ar, CO2, and H2O, but we are only considering the H2O as adsorbable. \n
	
	Adsorbent Input File
	--------------------
	Domain Coord. (2 = spherical, 1 = cylindrical) { [tab] Char. Length (um) (i.e., cylinder length) } \n
	(NOTE: Char. Length is only needed if problem is not spherical) \n
	Pellet Radius (um)  \n
 
	Example Adsorbent Input
	-----------------------
	1	6.0	\n
	2.0	\n
 
	Above example is for a cylindrical adsorbent with a length of 5 um and radius of 2 um.  \n
 
	Component Input File
	--------------------
	Molar Weight of ith species (g/mol) [tab] Specific Heat of ith species (J/g/K) \n
	Sutherland Viscosity (g/cm/s) [tab] Sutherland Temperature (K) [tab] Sutherland Constant (K) of ith species \n
	(repeat above for all species in same order they appeared in the Scenario Input File) \n
 
	Example Component Input
	-----------------------
	28.016	1.04	\n
	0.0001781	300.55	111.0	\n
	32.0	0.919	\n
	0.0002018	292.25	127.0	\n
	39.948	0.522	\n
	0.0002125	273.11	144.4	\n
	44.009	0.846	\n
	0.000148	293.15	240.0	\n
	18.0	1.97	\n
	0.0001043	298.16	784.72	\n
 
	Above example is a continuation of the Scenario Input example wherein each grouping represents parameters that
	are associated with N2, O2, Ar, CO2, and H2O, respectively. The order is VERY important! \n
 
	Adsorbate Input File
	--------------------
	Type of Surface Diffusion Function (0 = constant, 1 = simple Darken, 2 = theoretical Darken) \n
	Reference Diffusivity (um^2/hr) [tab] Activation Energy (J/mol) of ith adsorbable species\n
	Reference Temperature (K) [tab] Affinity Constant (-) of ith adsorbable species \n
	van der Waals Volume (cm^3/mol) of ith species \n
	GSTA adsorption capacity (mol/kg) of ith species \n
	Number of GSTA parameters of ith species \n
	Enthalpy (J/mol) of nth site       [tab]       Entropy of nth site (J/K/mol)       of ith species \n
	(repeat enthalpy and entropy for all n sites in species i) \n
	(repeat above for all species i) \n
 
	Example Adsorbate Input
	-----------------------
	0	\n
	0.8814	0.0	\n
	267.999	0.0	\n
	13.91	\n
	11.67	\n
	4	\n
	-46597.5	-53.6994	\n
	-125024		-221.073	\n
	-193619		-356.728	\n
	-272228		-567.459	\n
	1.28 	540.1	\n
	374.99	0.01	\n
	3.01	\n
	1.27	\n
	2	\n
	-46597.5	-53.6994	\n
	-125024		-221.073	\n
	
	Above example would be for a simulation involving two adsorbable species using a constant surface diffusion
	function. Each adsorbable species has it's own set of kinetic and equilibrium parameters that must be given
	in the same order as the species appeared in the Scenario Input. Note: we do not need to supply this
	information for non-adsorbable species. \n
 */
int SKUA_SCENARIOS(const char *scene, const char *sorbent, const char *comp, const char *sorbate);

/// Function to perform a test of the SKUA functions and routines
/** This function is callable from the UI and will perform a test simulation of the SKUA system of
	equations. Results from that test are output into a sub-directory called output and named
	SKUA_Test_Output.txt. */
int SKUA_TESTS();

#endif
