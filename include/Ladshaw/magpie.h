/*!
 *  \file magpie.h magpie.cpp
 *	\brief Multicomponent Adsorption Generalized Procedure for Isothermal Equilibria
 *	\details This file contains all functions and routines associated with predicting isothermal
 *			adsorption equilibria from only single component isotherm information. The basis of
 *			the model is the Adsorbed Solution Theory developed by Myers and Prausnitz (1965).
 *			Added to that base model is a procedure by which we can predict the non-idealities
 *			present at the surface phase by solving a closed system of equations involving
 *			the activity model.
 *
 *			For more details on this procedure, check out our publication in AIChE where we
 *			give a fully feature explaination of our Generalized Predictive Adsorbed Solution
 *			Theory (GPAST).
 *
 *			Reference: Ladshaw, A., Yiacoumi, S., and Tsouris, C., "A generalized procedure for
 *				the prediction of multicomponent adsorption equilibria", AIChE J., vol. 61, No. 8,
 *				p. 2600-2610, 2015.
 *
 *			MAGPIE represents a special case of the more general GPAST procedure, wherin the isotherm
 *			for each species is respresent by the GSTA isotherm (see gsta_opt.h) and the activity
 *			model for non-ideality at the adsorbent surface is a Modified Spreading Pressure Dependent
 *			(MSPD) model. See the above paper reference for more details.
 *
 *  \author Austin Ladshaw
 *	\date 12/17/2013
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#pragma once

#ifndef MAGPIE_HPP_
#define MAGPIE_HPP_

#include "lmcurve.h"			  //Main include to use the lmfit solver library
#include <stdio.h>				  //Line to allow for printf functions
#include <math.h>                 //Line added to allow usage of the pow (e, x) function
#include <iostream>				  //Line to allow for read/write to the console using cpp functions
#include <fstream>				  //Line to allow for read/write to and from .txt files
#include <stdlib.h>				  //Line need to convert strings to doubles
#include <vector>				  //Line needed to use dynamic arrays called vectors
#include <time.h>				  //Line needed to display program runtime
#include <float.h>				  //Line to allow use of machine constants
#include <string>    			  //Line to allow use of c++ strings
#include "error.h"				  //Line to allow use of custom error reporting
#include "lark.h"				  //Line to allow use of custom numerical methods

#ifndef DBL_EPSILON
#define DBL_EPSILON    2.2204460492503131e-016 		///< Machine precision value used for approximating gradients
#endif

#ifndef	Z
#define Z 10.0	///< Surface coordination number used in the MSPD activity model
#endif

#ifndef	A
#define A 3.13E+09				///< Corresponding van der Waals standard area for our coordination number (cm^2/mol)
#endif

#ifndef	V
#define V 18.92					///< Corresponding van der Waals standard volume for our coordination number (cm^3/mol)
#endif

#ifndef	Po
#define Po 100.0				///< Standard State Pressure - Units: kPa
#endif

#ifndef	R
#define R 8.3144621				///< Gas Constant - Units: J/(K*mol) = kB * Na
#endif

#ifndef	Na
#define Na 6.0221413E+23		///< Avagadro's Number - Units: molecules/mol
#endif

#ifndef	kB
#define kB 1.3806488E-23		///< Boltzmann's Constant - Units: J/K
#endif

/// This macro replaces all instances of shapeFactor(#) with the following single line calculation
#ifndef shapeFactor
#define shapeFactor(v_i) ( ( (Z - 2) * v_i ) / ( Z * V ) ) + ( 2 / Z )
#endif

/// This macro calculates the natural log of the dimensionless isotherm parameter
#ifndef lnKo
#define lnKo(H,S,T)	-( H / ( R * T ) ) + ( S / R )
#endif

/// This macro calculates the Henry's Coefficient for the ith component
#ifndef He
#define He(qm,K1,m) ( qm * K1 ) / ( m * Po )
#endif

/// GSTA Data Structure
/** C-style object holding all parameter information associated with the Generalized
	Statistical Thermodynamic Adsorption (GSTA) isotherm model. Each species in the gas phase
	will have one of these objects. */
typedef struct
{
	double qmax;						///< Theoretical maximum capacity of adsorbate-adsorbent pair (mol/kg)
	int m;								///< Number of parameters in the GSTA isotherm
	std::vector<double> dHo;			///< Enthalpies for each site (J/mol)
	std::vector<double> dSo;			///< Entropies for each site (J/(K*mol))
}GSTA_DATA;

/// MSPD Data Structure
/** C-Style object holding all parameter information associated with the Modified
	Spreading Pressure Dependent (SPD) activity model. Each species in the gas phase will have one
	of these objects. */
typedef struct
{
	double s;							///< Area shape factor
	double v;							///< van der Waals Volume (cm^3/mol)
	double eMax;						///< Maximum lateral interaction energy (J/mol)
	std::vector<double> eta;			///< Binary interaction parameter matrix (i,j)
	double gama;						///< Activity coefficient calculated from mSPD
}mSPD_DATA;

/// GPAST Data Structure
/** C-style object holding all parameter information associated with the Generalized
	Predictive Adsorbed Solution Theory (GPAST) system of equations. Each species in the
	gas phase will have one of these objects. */
typedef struct
{
	double x;							///< Adsorbed mole fraction
	double y;							///< Gas phase mole fraction
	double He;							///< Henry's Coefficient (mol/kg/kPa)
	double q;							///< Amount adsorbed for each component (mol/kg)
	std::vector<double> gama_inf;		///< Infinite dilution activities
	double qo;							///< Pure component capacities (mol/kg)
	double PIo;							///< Pure component spreading pressures (mol/kg)
	std::vector<double> po;				///< Pure component reference state pressures (kPa)
	double poi;							///< Reference state pressures solved for using Recover eval GPAST
	bool present;						///< If true, then the component is present; if false, then the component is not present
}GPAST_DATA;

/// System Data Structure
/** C-style object holding all the data associated with the overall system to be modeled. */
typedef struct
{
	double T;							///< System Temperature (K)
	double PT;							///< Total Pressure (kPa)
	double qT;							///< Total Amount adsorbed (mol/kg)
	double PI;							///< Total Lumped Spreading Pressure (mol/kg)
	double pi;							///< Actual Spreading pressure (J/m^2)
	double As;							///< Specific surface area of adsorbent (m^2/kg)
	int N;								///< Total Number of Components
	int I,J,K;							///< Special indices used to keep track of sub-systems
	unsigned long int total_eval;		///< Counter to keep track of total number of non-linear steps
	double avg_norm;					///< Used to store all norms from evaluations then average at end of run
	double max_norm;					///< Used to store the maximum e.norm calculated from non-linear iterations
	int Sys;							///< Number of sub-systems to solve
	int Par;							///< Number of binary parameters to solve for
	bool Recover;						///< If Recover == false, standard GPAST using y's as knowns
	bool Carrier;						///< If there is an inert carrier gas, Carrier == true
	bool Ideal;							///< If the behavior of the system is determined to be ideal, then Ideal == true
	bool Output;						///< Boolean to suppress output if desired (true = display, false = no display
}SYSTEM_DATA;

/// MAGPIE Data Structure
/** C-style object holding all information necessary to run a MAGPIE simulation. This is the data
	structure that will be used in other sub-routines when a mixed gas adsorption simulation needs
	to be run.*/
typedef struct
{
	std::vector<GSTA_DATA> gsta_dat;
	std::vector<mSPD_DATA> mspd_dat;
	std::vector<GPAST_DATA> gpast_dat;
	SYSTEM_DATA sys_dat;
}MAGPIE_DATA;

/// Function computes the result of the GSTA isotherm for the ith species
/** This function just computes the result of the GSTA isotherm model for the ith
	species given the partial pressure po.
 
	\param po partial pressure in kPa at which to evaluate the GSTA model
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double qo(double po, const void *data, int i);

/// Function computes the derivative of the GSTA model with respect to partial pressure
/** This function just computes the result of the derivative of GSTA isotherm model
	for the ith species at the given the partial pressure p.
 
	\param p partial pressure in kPa at which to evaluate the GSTA model
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double dq_dp(double p, const void *data, int i);

/// Function computes the ratio between the adsorbed amount and partial pressure for the GSTA isotherm
/** This function just computes the ratio between the adsorbed amount q (mol/kg) and the
	partial pressure p (kPa) at the given partial pressure. If p == 0, then this function
	returns the Henry's Law constant for the isotherm of the ith species.
 
	\param p partial pressure in kPa at which to evaluate the GSTA model
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double q_p(double p, const void *data, int i);

/// Function computes the spreading pressure integral of the ith species
/** This function uses an analytical solution to the spreading pressure integral with
	the GSTA isotherm to evaluate and return the value computed by that integral equation.
 
	\param po partial pressure in kPa at which to evaluate the lumped spreading pressure
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double PI(double po, const void *data, int i);

/// Function computes the heat of adsorption based on the ith species GSTA parameters
/** This function computes the isosteric heat of adsorption (J/mol) for the GSTA parameters
	of the ith species.
 
	\param po partial pressure in kPa at which to evaluate the heat of adsorption
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double Qst(double po, const void *data, int i);

/// Function to approximate the maximum lateral energy term for the ith species
/** The function attempts to approximate the maximum lateral energy term for the ith
	species. This is not a true maximum, but a cheaper estimate. Value being computed
	is used to shift the geometric mean and formulate the average cross-lateral
	energy term between species i and j. */
double eMax(const void *data, int i);

/// Function to evaluate the MSPD activity coefficient for the ith species
/** This function will return the natural log of the ith species activity coefficient
	using the Modified Spreading Pressure Dependent (MSPD) activity model. The par
	argument holds the variable values being solved for by GPAST and their contents
	will change depending on whether we are doing a forward or reverse evaluation.
	This function should not be called by the user and will only be called when needed
	in the GPAST routine.
	
	\param par list of parameters representing variables to be solved for in GPAST
	\param data void pointer for the MAGPIE_DATA data structure
	\param i ith species that we want to calculate the activity coefficient for
	\param PI lumped spreading pressure term used in gradient estimations */
double lnact_mSPD(const double *par, const void *data, int i, volatile double PI);

/// Function to approximate the derivative of the MSPD activity model with spreading pressure
/** This function returns a 2nd order, finite different approximation of the derivative of the
	MSPD activity model with the spreading pressure. The par argument will either hold the current
	iterates estimate of spreading pressure or should be passed as null. User does not need to call
	this function. GPAST will call automatically when needed.
 
	\param par list of parameters representing variables to be solved for in GPAST
	\param data void pointer for the MAGPIE_DATA data structure
	\param i ith species for which we will approximate the activty model gradient*/
double grad_mSPD(const double *par, const void *data, int i);

/// Function to calculate the total adsorbed amount (mol/kg) for the mixed surface phase
/** This function will uses the obtained system parameters from par and estimate the total
	amount of gases adsorbed to the surface in mol/kg. The user does not need to call this
	function, since this result will be stored in the SYSTEM_DATA structure.
 
	\param par list of parameters representing variables to be solved for in GPAST
	\param data void pointer for the MAGPIE_DATA data structure*/
double qT(const double *par, const void *data);

/// Function to provide an initial guess to the unknown parameters being solved for in GPAST
/** This function intends to provide an initial guess for the unknown values being solved for
	in the GPAST system. Depending on what type of solve is requested, this algorithm will
	provide a guess for the adsorbed or gas phase composition.
 
	\param par list of parameters representing variables to be solved for in GPAST
	\param data void pointer for the MAGPIE_DATA data structure*/
void initialGuess_mSPD(double *par, const void *data);

/// Function used with lmfit to evaluate the reference state pressure of a species based on spreading pressure
/** This function is used inside of the MSPD activity model to calculate the reference state pressure of a
	particular species at a given spreading pressure for the system. User does not need to call
	this function. GPAST will call automatically when needed.
 
	\param par list of parameters representing variables to be solved for in GPAST
	\param m_dat number of functions/variables in the GPAST system of equations
	\param data void pointer for the MAGPIE_DATA data structure
	\param fvec list of residuals formed by the functions in GPAST
	\param info integer flag variable used in the lmfit routine*/
void eval_po_PI(const double *par, int m_dat, const void *data, double *fvec, int *info);

/// Function used with lmfit to evaluate the reference state pressure of a species based on that species isotherm
/** This function is used to evaluate the partial pressure or reference state pressure for a particular species
	given single-component adsorbed amount. User does not need to call this function. GPAST will call automatically
	when needed.
 
	\param par list of parameters representing variables to be solved for in GPAST
	\param m_dat number of functions/variables in the GPAST system of equations
	\param data void pointer for the MAGPIE_DATA data structure
	\param fvec list of residuals formed by the functions in GPAST
	\param info integer flag variable used in the lmfit routine*/
void eval_po_qo(const double *par, int m_dat, const void *data, double *fvec, int *info);

/// Function used with lmfit to evaluate the reference state pressure of a species based on a sub-system
/** This function is used to approximate reference state pressures based on the spreading pressure of a
	sub-system in GPAST. The sub-system will be one of the unique binary systems that exist in the overall
	mixed gas system. User does not need to call this function. GPAST will call automatically when needed.
 
	\param par list of parameters representing variables to be solved for in GPAST
	\param m_dat number of functions/variables in the GPAST system of equations
	\param data void pointer for the MAGPIE_DATA data structure
	\param fvec list of residuals formed by the functions in GPAST
	\param info integer flag variable used in the lmfit routine*/
void eval_po(const double *par, int m_dat, const void *data, double *fvec, int *info);

/// Function used with lmfit to evaluate the binary interaction parameters for each unique species pair
/** This function is used to estimate the binary interaction parameters for all species pairs in a given
	sub-system. Those parameters are then stored for later used when evaluating the activity coefficients
	for the overall mixture. User does not need to call this function. GPAST will call automatically when needed.
 
	\param par list of parameters representing variables to be solved for in GPAST
	\param m_dat number of functions/variables in the GPAST system of equations
	\param data void pointer for the MAGPIE_DATA data structure
	\param fvec list of residuals formed by the functions in GPAST
	\param info integer flag variable used in the lmfit routine*/
void eval_eta(const double *par, int m_dat, const void *data, double *fvec, int *info);

/// Function used with lmfit to solve the GPAST system of equations
/** This function is used after having calculated and stored all necessary information to solve a closed form
	GPAST system of equations. User does not need to call this function. GPAST will call automatically when needed.
 
	\param par list of parameters representing variables to be solved for in GPAST
	\param m_dat number of functions/variables in the GPAST system of equations
	\param data void pointer for the MAGPIE_DATA data structure
	\param fvec list of residuals formed by the functions in GPAST
	\param info integer flag variable used in the lmfit routine*/
void eval_GPAST(const double *par, int m_dat, const void *data, double *fvec, int *info);

/// Function to call all sub-routines to solve a MAGPIE/GPAST problem at a given temperature and pressure
/** This is the function that a typical user will want to incorporate into their own codes when evaluating
	adsorption of a gas mixture. Prior to calling this function, all required structures and information in
	the MAGPIE_DATA structure must have been properly initialized. After this function has completed it's
	operations, it will return an integer used to denote a success or failure of the routine. Integers
	0, 1, 2, and 3 all denote success. Anything else is considered a failure.
 
	To setup the MAGPIE_DATA structure correctly, you must reserve space for all vector objects based on the
	number of gas species in the mixture. In general, you only need to reserve space for the adsorbing species.
	However, you can also reserve space for non-adsorbing species, but you MUST give a gas/adsorbed mole fraction
	of the non-adsorbing species 0.0 so that the routine knows to ignore them (very important)!
 
	After setting up the memory for the vector objects, you can intialize information specific to the simulation
	you want to request. The number of species (N), total pressure (PT) and gas temperature (T) must always be given.
	You can neglect the non-idealities of the surface phase by setting the Ideal bool to true. This will result in
	faster calculations, because MAGPIE will just revert down to the Ideal Adsorbed Solution Theory (IAST).
 
	The Recover bool will denote whether we are doing a forward or reverse GPAST evaluation. Forward evaluation is for
	solving for the composition of the adsorbed phase given the composition of the gas phase (Recover = false).
	Reverse evaluation is for solve for the composition of the gas phase given the composition of the adsorbed phase
	(Recover = true).
 
	For a reverse evaluation (Recover = true) you will also need to stipulate whether or not there is a carrier gas
	(Carrier = true or false). A carrier gas is considered any non-adsorbing species that may be present in the gas
	phase and contributing to the total pressure in the system.
 
 
	The parameters that must be initialized for all species include all GSTA_DATA parameters and the van der Waals
	volume parameter (v) in the mSPD_DATA structure. For non-adsorbing species, you can ignore these parameters, but
	need to set the sites (m) from GSTA_DATA to 1. GPAST cannot run any evaluations without these parameters being
	set properly AND set in the same order for all species (i.e., make sure that gpast_dat[i].qmax corresponds to
	mspd_dat[i].v and so on).
 
	Lastly, you need to give either the gas phase or adsorbed phase mole fractions, depending on whether you are
	going to run a forward or reverse evaluation, respectively. For a forward evaluation, provide the gas mole
	fractions (y) in GPAST_DATA for each species (non-adsorbing species should have this value set to 0.0). For a
	reverse evaluation, provide the adsorbed mole fractions (x) in GPAST_DATA for each species, as well as the total
	adsorbed amount (qT) in SYSTEM_DATA. Again, non-adsorbing species should have their respective phase mole fractions
	set to 0.0 to exclude them from the simulation. Additionally, if there are non-adsorbing species present, then the
	Carrier bool in SYSTEM_DATA must be set to true.
 
	\param data void pointer for the MAGPIE_DATA data structure holding all necessary information*/
int MAGPIE(const void *data);

/// Function to perform a series of MAGPIE simulations based on given input files
/** This function is callable from the UI and is used to perform a series of isothermal equilibria evaluations
	using the MAGPIE routines. There are two input files that must be provided: (i) inputFileName - containing
	parameter information for the species and (ii) sceneFileName - containing information for each MAGPIE
	simulation. Each of these files have a specific structure (see below). NOTE: this may change in future versions.
 
	inputFileName Text File Structure:
	---------------------------------
	Integer for Number of Adsorbing Species \n
	van der Waals Volume (cm^3/mol) of ith species \n
	GSTA adsorption capacity (mol/kg) of ith species \n
	Number of GSTA parameters of ith species \n
	Enthalpy (J/mol) of nth site       [tab]       Entropy of nth site (J/K/mol)       of ith species \n
 (repeat above for all n sites in species i) \n
	(repeat above for all species i) \n
 
	Example Input File:
	-------------------
	5 \n
	17.1 \n
	5.8797 \n
	1 \n
	-20351.9	-81.8369 \n
	16.2 \n
	5.14934 \n
	1 \n
	-16662.7	-74.4766 \n
	19.7 \n
	9.27339 \n
	4 \n
	-46597.5	-53.6994	\n
	-125024	-221.073	\n
	-193619	-356.728	\n
	-272228	-567.459 \n
	13.25 \n
	4.59144 \n
	1 \n
	-13418.5	-84.888 \n
	18.0 \n
	10.0348 \n
	1 \n
	-20640.4	-72.6119 \n
 
	(The above input file gives the parameter information for 5 adsorbing species) \n
 
	sceneFileName Text File Structure:
	---------------------------------
	Integer Flag to mark Forward (0) or { Reverse (1) evaluations } \n
	Number of Simulations to Run \n
	Total Pressure (kPa) [tab] Temperature (K) { [tab] Total Adsorption (mol/kg) [tab] Carrier Gas Flag (0=false, 1=true) } \n
	Gas/Adsorbed Mole Fractions for each species in the order given in prior file (tab separated) \n
	(repeat above for all simulations desired) \n
	NOTE: only provide the Total Adsorption and Carrier Flag if doing Reverse evaluations! \n
 
	Example Scenario File 1:
	------------------------
	0 \n
	4 \n
	0.65	303.15 \n
	0.364	0.318	0.318 \n
	3.25	303.15 \n
	0.371	0.32	0.309 \n
	6.85	303.15 \n
	0.388	0.299	0.313 \n
	13.42	303.15 \n
	0.349	0.326	0.325 \n
 
	(The above scenario file is for 4 forward evaluations/simulations for a 3-adsorbing species system) \n
 
	Example Scenario File 2:
	------------------------
	1 \n
	4 \n
	0.65	303.15  5.4   0\n
	0.364	0.318	0.318 \n
	3.25	303.15  7.7   0\n
	0.371	0.32	0.309 \n
	6.85	303.15  9.8   0\n
	0.388	0.299	0.313 \n
	13.42	303.15  10.4  0\n
	0.349	0.326	0.325 \n
 
	(The above scenario file is for 4 reverse evaluations/simulations for a 3-adsorbing species system and no carrier gas) */
int MAGPIE_SCENARIOS(const char *inputFileName, const char *sceneFileName);

#endif
