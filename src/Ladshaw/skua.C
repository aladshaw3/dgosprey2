/*!
 *  \file skua.h skua.cpp
 *	\brief Surface Kinetics for Uptake by Adsorption
 *
 *  \author Austin Ladshaw
 *	\date 01/26/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "skua.h"

//Function to print out the species header
void print2file_species_header(FILE *Output, SKUA_DATA *skua_dat, int i)
{
	const char *name = skua_dat->param_dat[i].speciesName.c_str();
	fprintf(Output, "%s\t", name);
	if (skua_dat->finch_dat[i].Dirichlet == true)
		fprintf(Output,"-\t");
	for (int l=0; l<skua_dat->finch_dat[i].LN+4; l++)
	{
		fprintf(Output,"-\t");
	}
}

//Print time header for skua output
void print2file_SKUA_time_header(FILE *Output, SKUA_DATA *skua_dat, int i)
{
	fprintf(Output,"Time\t");
	fprintf(Output,"q[0] (mol/kg)\t");
	for (int l=0; l<skua_dat->finch_dat[i].LN; l++)
	{
		if (skua_dat->finch_dat[i].Dirichlet == false && l == skua_dat->finch_dat[i].LN-1)
			break;
		fprintf(Output,"q[%i] (mol/kg)\t",l+1);
	}
	fprintf(Output,"q Total (mol/kg)\tq Average (mol/kg)\t");
	fprintf(Output,"x[%i] (-)\tQst[%i] (J/mol)\t",i,i);
}

//Print out the skua header to file
void print2file_SKUA_header(SKUA_DATA *skua_dat)
{
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_species_header(skua_dat->OutputFile, skua_dat, i);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_dim_header(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_SKUA_time_header(skua_dat->OutputFile, skua_dat, i);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
}

//Print old results to file
void print2file_SKUA_results_old(SKUA_DATA *skua_dat)
{
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_result_old(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			fprintf(skua_dat->OutputFile, "%.6g\t%.6g\t",skua_dat->param_dat[i].xn,skua_dat->param_dat[i].Qstn);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
}

//Print new results to file
void print2file_SKUA_results_new(SKUA_DATA *skua_dat)
{
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_result_new(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			fprintf(skua_dat->OutputFile, "%.6g\t%.6g\t",skua_dat->param_dat[i].xnp1,skua_dat->param_dat[i].Qstnp1);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
}

//Default function for computing Dc
double default_Dc(int i, int l, const void *data)
{
	double Dc = 0.0;
	SKUA_DATA *dat = (SKUA_DATA *) data;
	
	if (i < 0)
	{
		mError(indexing_error);
		return -1;
	}
	Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
	Dc = D_inf(Dc, dat->param_dat[i].ref_temperature, dat->param_dat[i].affinity, dat->param_dat[i].ref_pressure, dat->magpie_dat.sys_dat.T);
	
	return Dc;
}

//Default function for computing kf (Film Mass Transfer: um/hr)
double default_kf(int i, const void *data)
{
	return 0.0;
}

//Function to set Dc to a constant based on ref_diffusion
double const_Dc(int i, int l, const void *data)
{
	SKUA_DATA *dat = (SKUA_DATA *) data;
	return dat->param_dat[i].ref_diffusion;
}

//Simple Darken Dc
double simple_darken_Dc(int i, int l, const void *data)
{
	double Dc = 0.0;
	SKUA_DATA *dat = (SKUA_DATA *) data;
	
	if (i < 0)
	{
		mError(indexing_error);
		return -1;
	}
	
	//Boundary Dc
	if (l < 0)
	{
		Dc = D_c(D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T),dat->finch_dat[i].uo/(dat->magpie_dat.gsta_dat[i].qmax+1.0E-6));
	}
	//Interior Dc
	else
	{
		Dc = D_c(D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T),dat->finch_dat[i].unp1(l,0)/(dat->magpie_dat.gsta_dat[i].qmax+1.0E-6));
	}
	Dc = D_inf(Dc, dat->param_dat[i].ref_temperature, dat->param_dat[i].affinity, dat->param_dat[i].ref_pressure, dat->magpie_dat.sys_dat.T);
	
	return Dc;
}

//Theoretical Darken Dc (Much more computationally intense)
double theoretical_darken_Dc(int i, int l, const void *data)
{
	double Dc = 0.0;
	SKUA_DATA *dat = (SKUA_DATA *) data;
	
	if (i < 0)
	{
		mError(indexing_error);
		return -1;
	}
	
	if (dat->param_dat[i].Adsorbable == true)
	{
		//Boundary Dc
		if (l < 0)
		{
			Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
			Dc = Dc * q_p(dat->y[i]*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i) / dq_dp(dat->y[i]*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i);
		}
		//Interior Dc
		else
		{
			
			dat->magpie_dat.sys_dat.Recover = true;
			dat->magpie_dat.sys_dat.Carrier = true;
			dat->magpie_dat.gpast_dat[i].x = 1.0;
			dat->magpie_dat.sys_dat.qT = dat->finch_dat[i].unp1(l,0);
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
			{
				if (j!=i)
					dat->magpie_dat.gpast_dat[j].x = 0.0;
			}
			
			//Call MAGPIE to determine equilibrium adsorption
			double success = 0;
			if (dat->magpie_dat.sys_dat.qT > 0.0)
			{
				success = MAGPIE((void *)&dat->magpie_dat);
				if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
				else success = 0;
			}
			else
			{
				dat->magpie_dat.gpast_dat[i].y = 0.0;
			}
			dat->magpie_dat.sys_dat.Recover = false;
			dat->magpie_dat.sys_dat.Carrier = false;
			dat->total_steps = dat->total_steps + dat->magpie_dat.sys_dat.total_eval;
			
			
			Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
			Dc = Dc * q_p(dat->magpie_dat.gpast_dat[i].y*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i) / dq_dp(dat->magpie_dat.gpast_dat[i].y*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i);
		}
	}
	else
	{
		Dc = 0.0;
	}
	Dc = D_inf(Dc, dat->param_dat[i].ref_temperature, dat->param_dat[i].affinity, dat->param_dat[i].ref_pressure, dat->magpie_dat.sys_dat.T);
	
	return Dc;
}

//Estimation of Film Mass Transfer based on gas phase properties
double empirical_kf(int i, const void *data)
{
	int success = 0;
	SKUA_DATA *dat = (SKUA_DATA *) data;
	if (dat->gas_dat != NULL)
	{
		success = set_variables(dat->magpie_dat.sys_dat.PT, dat->magpie_dat.sys_dat.T, dat->gas_velocity, dat->pellet_radius*2.0E-4, dat->y, dat->gas_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		//Calculate Properties
		success = calculate_properties(dat->gas_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		return FilmMTCoeff(dat->gas_dat->species_dat[i].molecular_diffusion, dat->gas_dat->char_length, dat->gas_dat->Reynolds, dat->gas_dat->species_dat[i].Schmidt) * 3.6E+7;
	}
	else
	{
		std::cout << "\nWARNING!!! Calling a Function without the correct pointer initialized!\n" << std::endl;
		return default_kf(i, data);
	}
}

//Function to give the constant film mass transfer rate {NOT USED}
double const_kf(int i, const void *data)
{
	SKUA_DATA *dat = (SKUA_DATA *) data;
	return dat->param_dat[i].film_transfer;
}

//Check the input mole fractions from files in SKUA for errors
int molefractionCheck(SKUA_DATA *skua_dat)
{
	int success = 0;
	double sumx = 0.0, sumy = 0.0;
	
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		sumx = sumx + skua_dat->param_dat[i].xIC;
		sumy = sumy + skua_dat->y[i];
	}
	if (sumy > (1.0 + 1.0E-6) || sumx > (1.0 + 1.0E-6))
	{
		mError(invalid_molefraction);
		success = -1;
	}
	
	return success;
}

//Function to setup remaining memory and structures after initializing parameter information
int setup_SKUA_DATA(FILE *file, double (*eval_Dc) (int i, int l, const void *user_data),
					double (*eval_Kf) (int i, const void *user_data), const void *user_data,
					MIXED_GAS *gas_data, SKUA_DATA *skua_dat)
{
	int success = 0;
	
	if (file == NULL)
		skua_dat->Print2File = false;
	else
	{
		skua_dat->Print2File = true;
		skua_dat->OutputFile = file;
	}
	if ((*eval_Dc) == NULL)
		skua_dat->eval_diff = (*default_Dc);
	else
		skua_dat->eval_diff = (*eval_Dc);
	if ((*eval_Kf) == NULL && skua_dat->DirichletBC == true)
		skua_dat->eval_kf = (*default_kf);
	else if ((*eval_Kf) == NULL && skua_dat->DirichletBC == false)
		skua_dat->eval_kf = (*empirical_kf);
	else
		skua_dat->eval_kf = (*eval_Kf);
	skua_dat->user_data = user_data;
	if (gas_data == NULL)
		skua_dat->DirichletBC = true;
	else
		skua_dat->gas_dat = gas_data;
	
	//Loop to setup the FINCH data
	skua_dat->t = 0.0;
	skua_dat->t_old	= 0.0;
	skua_dat->t_counter = 0.0;
	skua_dat->magpie_dat.sys_dat.Output = false;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->finch_dat[i].t = skua_dat->t;
		skua_dat->finch_dat[i].t_old = skua_dat->t_old;
		skua_dat->finch_dat[i].L = skua_dat->pellet_radius;			//um
		skua_dat->finch_dat[i].T = skua_dat->sim_time;				//hrs
		skua_dat->finch_dat[i].Dirichlet = skua_dat->DirichletBC;	//-
		skua_dat->finch_dat[i].LN = 10;						//-
		skua_dat->finch_dat[i].SteadyState = false;			//-
		skua_dat->finch_dat[i].CheckMass = true;				//-
		skua_dat->finch_dat[i].Iterative = skua_dat->NonLinear;		//-
		skua_dat->finch_dat[i].NormTrack = skua_dat->Print2Console;	//-
		skua_dat->finch_dat[i].nl_method = FINCH_Picard;					//-
		skua_dat->finch_dat[i].d = skua_dat->coord;							//-
		skua_dat->finch_dat[i].s = skua_dat->char_measure;						//-
		skua_dat->finch_dat[i].t = 0.0;						//hrs
		skua_dat->finch_dat[i].t_old = 0.0;					//hrs
		skua_dat->finch_dat[i].Ro = 1.0;						//-
		skua_dat->finch_dat[i].RIC = 1.0;						//-
		skua_dat->finch_dat[i].vo = 0.0;						//um/hr
		skua_dat->finch_dat[i].vIC = 0.0;						//um/hr
		skua_dat->finch_dat[i].kIC = 0.0;						//1/hr
		skua_dat->finch_dat[i].ko = 0.0;						//1/hr
		skua_dat->param_dat[i].film_transfer = 0.0;
		
		success = setup_FINCH_DATA(NULL,NULL,NULL,NULL,NULL,set_SKUA_params,NULL,NULL,NULL,NULL,NULL,NULL,&skua_dat->finch_dat[i],user_data);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	
	return success;
}
//Execution order for the SKUA routine (Note: the time step must be set prior to calling this routine)
int SKUA_Executioner(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//Perform Preprocess Actions
	success = SKUA_preprocesses(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Loop for all components
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			//Solve the system
			success = (*skua_dat->finch_dat[i].solve) ((void *)&skua_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Check for negative concentrations
			success = check_Mass(&skua_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Form Totals
			success = uTotal(&skua_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Form Averages
			success = uAverage(&skua_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		
	}
	
	//Perform Postprocess Actions
	success = SKUA_postprocesses(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	return success;
}

//Set the initial conditions
int set_SKUA_ICs(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//Set MAGPIE parameters to perform reverse evaluations
	skua_dat->magpie_dat.sys_dat.Recover = true;
	skua_dat->magpie_dat.sys_dat.Carrier = false;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == false)
		{
			skua_dat->magpie_dat.sys_dat.Carrier = true;
			break;
		}
	}
	
	//Set the initial conditions
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->param_dat[i].ref_pressure = skua_dat->y[i]*skua_dat->magpie_dat.sys_dat.PT;
		if (skua_dat->param_dat[i].Adsorbable == true)
			skua_dat->magpie_dat.gpast_dat[i].x = skua_dat->param_dat[i].xIC;
		else
			skua_dat->magpie_dat.gpast_dat[i].x = 0.0;
		skua_dat->param_dat[i].xn = skua_dat->magpie_dat.gpast_dat[i].x;
		skua_dat->param_dat[i].xnp1 = skua_dat->param_dat[i].xn;
		skua_dat->finch_dat[i].uAvg = skua_dat->magpie_dat.sys_dat.qT * skua_dat->magpie_dat.gpast_dat[i].x;
	}
	
	//Call MAGPIE to determine equilibrium adsorption
	success = MAGPIE((void *)&skua_dat->magpie_dat);
	if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
	else success = 0;
	skua_dat->total_steps = skua_dat->total_steps + skua_dat->magpie_dat.sys_dat.total_eval;
	
	//Use MAGPIE result to form initial condition in adsorbent
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->finch_dat[i].Sn.ConstantICFill(0.0);
		skua_dat->finch_dat[i].Snp1.ConstantICFill(0.0);
		skua_dat->finch_dat[i].vn.ConstantICFill(0.0);
		skua_dat->finch_dat[i].vnp1.ConstantICFill(0.0);
		skua_dat->finch_dat[i].kn.ConstantICFill(0.0);
		skua_dat->finch_dat[i].knp1.ConstantICFill(0.0);
		skua_dat->finch_dat[i].Rn.ConstantICFill(1.0);
		skua_dat->finch_dat[i].Rnp1.ConstantICFill(1.0);
		skua_dat->finch_dat[i].un.ConstantICFill(skua_dat->magpie_dat.sys_dat.qT*skua_dat->param_dat[i].xIC);
		skua_dat->finch_dat[i].unm1 = skua_dat->finch_dat[i].un;
		skua_dat->finch_dat[i].unp1 = skua_dat->finch_dat[i].un;
		skua_dat->param_dat[i].y_eff = skua_dat->magpie_dat.gpast_dat[i].y;
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			skua_dat->param_dat[i].Qstn = Qst((skua_dat->magpie_dat.gpast_dat[i].y*skua_dat->magpie_dat.sys_dat.PT),(void *)&skua_dat->magpie_dat,i);
		}
		else
		{
			skua_dat->param_dat[i].Qstn = 0.0;
		}
		skua_dat->param_dat[i].Qstnp1 = skua_dat->param_dat[i].Qstn;
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			for (int l=0; l<skua_dat->finch_dat[i].LN; l++)
			{
				skua_dat->finch_dat[i].Dn(l,0) = (*skua_dat->eval_diff) (i,l,skua_dat->user_data);
				skua_dat->finch_dat[i].Dnp1(l,0) = skua_dat->finch_dat[i].Dn(l,0);
			}
		}
		else
		{
			skua_dat->finch_dat[i].Dn.ConstantICFill(0.0);
			skua_dat->finch_dat[i].Dnp1.ConstantICFill(0.0);
		}
		skua_dat->finch_dat[i].Do = skua_dat->finch_dat[i].Dn(skua_dat->finch_dat[i].LN-1,0);
		skua_dat->finch_dat[i].uo = skua_dat->finch_dat[i].un(skua_dat->finch_dat[i].LN-1,0);
		skua_dat->finch_dat[i].kfn = (*skua_dat->eval_kf) (i,skua_dat->user_data);
		if (skua_dat->finch_dat[i].kfn < 0.0)
		{
			mError(simulation_fail);
			return -1;
		}
		
		//Form Totals
		success = uTotal(&skua_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Form Averages
		success = uAverage(&skua_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Update averages and totals
		skua_dat->finch_dat[i].uT_old = skua_dat->finch_dat[i].uT;
		skua_dat->finch_dat[i].uAvg_old = skua_dat->finch_dat[i].uAvg;
	}
	skua_dat->qTnp1 = skua_dat->qTn;
	
	//Reset the MAGPIE parameters before continuing
	skua_dat->magpie_dat.sys_dat.Recover = false;
	skua_dat->magpie_dat.sys_dat.Carrier = false;
	
	return success;
}

//Setup the time step
int set_SKUA_timestep(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->finch_dat[i].dt = skua_dat->finch_dat[i].dz / 4.0;
		skua_dat->finch_dat[i].t = skua_dat->finch_dat[i].dt + skua_dat->finch_dat[i].t_old;
	}
	skua_dat->t_old = skua_dat->finch_dat[0].t_old;
	skua_dat->t = skua_dat->finch_dat[0].t;
	
	return success;
}

//Perform Preprocess Actions
int SKUA_preprocesses(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//First, establish the boundary conditions by solving a MAGPIE problem
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		double p_out = skua_dat->y[i]*skua_dat->magpie_dat.sys_dat.PT;
		if (p_out >= skua_dat->param_dat[i].ref_pressure) skua_dat->param_dat[i].ref_pressure = p_out;
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			skua_dat->magpie_dat.gpast_dat[i].y = skua_dat->y[i];
		}
		else
		{
			skua_dat->magpie_dat.gpast_dat[i].y = 0.0;
		}
	}
	
	//Call MAGPIE to determine equilibrium adsorption
	success = MAGPIE((void *)&skua_dat->magpie_dat);
	if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
	else success = 0;
	skua_dat->total_steps = skua_dat->total_steps + skua_dat->magpie_dat.sys_dat.total_eval;
	
	//Loop to set bounary values and parameters
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		//Bounary value is adsorbed amount
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			skua_dat->finch_dat[i].uo = skua_dat->magpie_dat.gpast_dat[i].q;
			skua_dat->finch_dat[i].kfnp1 = (*skua_dat->eval_kf) (i,skua_dat->user_data);
			if (skua_dat->finch_dat[i].kfnp1 < 0.0)
			{
				mError(simulation_fail);
				return -1;
			}
		}
		else
		{
			skua_dat->finch_dat[i].uo = 0.0;
			skua_dat->finch_dat[i].kfnp1 = 0.0;
		}
	}
	
	return success;
}

//Set the parameters for the FINCH problem
int set_SKUA_params(const void *user_data)
{
	int success = 0;
	
	//In this file, we force an Implicit Only FINCH solution
	SKUA_DATA *dat = (SKUA_DATA *) user_data;
	for (int i=0; i<dat->magpie_dat.sys_dat.N; i++)
	{
		dat->finch_dat[i].CN = false;
		dat->finch_dat[i].beta = 1.0;
		if (dat->param_dat[i].Adsorbable == true)
		{
			dat->finch_dat[i].Do = (*dat->eval_diff) (i,-1,dat->user_data);
			//Loop over the nodes
			for (int l=0; l<dat->finch_dat[i].LN; l++)
			{
				dat->finch_dat[i].Dnp1(l,0) = (*dat->eval_diff) (i,l,dat->user_data);
			}
		}
		else
		{
			dat->finch_dat[i].Do = 0.0;
			dat->finch_dat[i].Dnp1.ConstantICFill(0.0);
		}
	}
	
	return success;
}

//Perform Postprocess Actions
int SKUA_postprocesses(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//Update some SKUA information
	skua_dat->qTnp1 = 0.0;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->qTnp1 = skua_dat->qTnp1 + skua_dat->finch_dat[i].uAvg;
		skua_dat->total_steps = skua_dat->total_steps + skua_dat->finch_dat[i].total_iter;
		skua_dat->finch_dat[i].total_iter = 0;
	}
	if (skua_dat->qTnp1 >=0.0)
		skua_dat->magpie_dat.sys_dat.qT = skua_dat->qTnp1;
	else
		skua_dat->magpie_dat.sys_dat.qT = 0.0;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->qTnp1 > 0.0)
			skua_dat->param_dat[i].xnp1 = skua_dat->finch_dat[i].uAvg / skua_dat->qTnp1;
		else
			skua_dat->param_dat[i].xnp1 = 0.0;
		skua_dat->magpie_dat.gpast_dat[i].x = skua_dat->param_dat[i].xnp1;
	}
	
	//Setup for another reverse MAGPIE evaluation
	skua_dat->magpie_dat.sys_dat.Recover = true;
	skua_dat->magpie_dat.sys_dat.Carrier = false;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == false)
		{
			skua_dat->magpie_dat.sys_dat.Carrier = true;
			break;
		}
	}
	
	//Call MAGPIE to determine equilibrium adsorption
	success = MAGPIE((void *)&skua_dat->magpie_dat);
	if (success < 0 || success > 5)
	{
		mError(simulation_fail);
		std::cout << success << std::endl;
		for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
		{
			std::cout << "x[" << i << "] = " << skua_dat->magpie_dat.gpast_dat[i].x << "\t";
			std::cout << "q[" << i << "] = " << skua_dat->magpie_dat.gpast_dat[i].q << "\t";
			std::cout << "y[" << i << "] = " << skua_dat->magpie_dat.gpast_dat[i].y << "\n";
		}
		return -1;
	}
	else success = 0;
	skua_dat->total_steps = skua_dat->total_steps + skua_dat->magpie_dat.sys_dat.total_eval;
	
	//Calculate the new heat of adsorption
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			skua_dat->param_dat[i].Qstnp1 = Qst((skua_dat->magpie_dat.gpast_dat[i].y*skua_dat->magpie_dat.sys_dat.PT),(void *)&skua_dat->magpie_dat,i);
		}
		else
		{
			skua_dat->param_dat[i].Qstnp1 = 0.0;
		}
		skua_dat->param_dat[i].y_eff = skua_dat->magpie_dat.gpast_dat[i].y;
	}
	
	//Print results
	if (skua_dat->Print2File == true)
	{
		skua_dat->t_counter = skua_dat->t_counter + skua_dat->finch_dat[0].dt;
		if (skua_dat->t_counter > skua_dat->t_print)
		{
			print2file_SKUA_results_new(skua_dat);
			skua_dat->t_counter = 0.0;
		}
	}
	
	//Reset the MAGPIE parameters before continuing
	skua_dat->magpie_dat.sys_dat.Recover = false;
	skua_dat->magpie_dat.sys_dat.Carrier = false;
	
	return success;
}

//Reset the values for next time step
int SKUA_reset(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		success = (*skua_dat->finch_dat[i].resettime) ((void *)&skua_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
		skua_dat->param_dat[i].Qstn = skua_dat->param_dat[i].Qstnp1;
		skua_dat->param_dat[i].xn = skua_dat->param_dat[i].xnp1;
	}
	skua_dat->qTn = skua_dat->qTnp1;
	if (skua_dat->t_old == 0.0)
		skua_dat->t_old = 0.0;
	else
		skua_dat->t_old = skua_dat->t;
	
	return success;
}

//Calling the SKUA routine
int SKUA(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//Print to file
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_header(skua_dat);
	}
	
	//Set Initial Conditions
	success = set_SKUA_ICs(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Print ICs
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_results_old(skua_dat);
	}
	
	//Solve a series of time steps implicitly
	do
	{
		if (skua_dat->finch_dat[0].Update == true)
		{
			success = SKUA_reset(skua_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_SKUA_timestep(skua_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		std::cout << "Evaluating time: " << skua_dat->t << " hrs..." << std::endl;
		success = SKUA_Executioner(skua_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
				skua_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); skua_dat->finch_dat[0].Update = false; return -1;}
		skua_dat->total_steps++;
	} while (skua_dat->t < skua_dat->sim_time);
	
	return success;
}

//Calling the standard SKUA CYCLE Test
int SKUA_CYCLE_TEST01(SKUA_DATA *skua_dat)
{
	int success = 0;
	skua_dat->magpie_dat.sys_dat.T = 313.15; 			//K
	skua_dat->magpie_dat.sys_dat.PT = 101.35;			//kPa
	skua_dat->sim_time = 525.0;							//hrs
	skua_dat->t_print = skua_dat->sim_time / 10000.0;	//hrs
	skua_dat->y[0] = 0.787289;							//-
	skua_dat->y[1] = 0.209386;							//-
	skua_dat->y[2] = 0.003325;							//-
	if (skua_dat->gas_dat != NULL)
		skua_dat->gas_dat->CheckMolefractions = false;
	
	//Print to file
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_header(skua_dat);
	}
	
	//Set Initial Conditions
	success = set_SKUA_ICs(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Print ICs
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_results_old(skua_dat);
	}
	
	do
	{
		if (skua_dat->finch_dat[0].Update == true)
		{
			success = SKUA_reset(skua_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_SKUA_timestep(skua_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Cycling Info Here
		if (skua_dat->t > 18.3)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 0.000434;					//-
			
			double m = (0.000128 - 0.000434) / (18.8 - 18.3);
			skua_dat->y[2] = (m * (skua_dat->t - 18.3)) + 0.000434;
		}
		if (skua_dat->t > 18.8)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 0.000128;					//-
			
			double m = (3.6E-5 - 0.000128) / (22.1 - 18.8);
			skua_dat->y[2] = (m * (skua_dat->t - 18.8)) + 0.000128;
		}
		if (skua_dat->t > 22.1)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 3.6E-5;					//-
			
			double m = (2.9E-5 - 3.6E-5) / (25.3 - 22.1);
			skua_dat->y[2] = (m * (skua_dat->t - 22.1)) + 3.6E-5;
		}
		if (skua_dat->t > 25.3)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 2.9E-5;					//-
			
			double m = (2.0E-5 - 2.9E-5) / (42.1 - 25.3);
			skua_dat->y[2] = (m * (skua_dat->t - 25.3)) + 2.9E-5;
		}
		if (skua_dat->t > 42.1)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 2.0E-5;					//-
			
			double m = (8.2E-6 - 2.0E-5) / (97.1 - 42.1);
			skua_dat->y[2] = (m * (skua_dat->t - 42.1)) + 2.0E-5;
		}
		if (skua_dat->t > 97.1)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 8.2E-6;					//-
			
			double m = (6.0E-6 - 8.2E-6) / (114.5 - 97.1);
			skua_dat->y[2] = (m * (skua_dat->t - 97.1)) + 8.2E-6;
		}
		if (skua_dat->t > 114.5)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 6.0E-6;					//-
		}
		if (skua_dat->t > 354.0)
		{
			skua_dat->y[0] = 0.787289;				//-
			skua_dat->y[1] = 0.209386;				//-
			skua_dat->y[2] = 0.003325;				//-
		}
		if (skua_dat->t > 378.35)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 6.0E-6;					//-
		}
		if (skua_dat->t > 498.2)
		{
			skua_dat->y[0] = 0.787289;				//-
			skua_dat->y[1] = 0.209386;				//-
			skua_dat->y[2] = 0.003325;				//-
		}
		if (skua_dat->t > 522.3)
		{
			skua_dat->y[0] = 0.789994;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 6.0E-6;					//-
		}
		
		std::cout << "Evaluating time: " << skua_dat->t << " hrs..." << std::endl;
		success = SKUA_Executioner(skua_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
				skua_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); skua_dat->finch_dat[0].Update = false; return -1;}
		skua_dat->total_steps++;
	} while (skua_dat->t < skua_dat->sim_time);
	
	return success;
}

//Calling the SKUA Cycle Test 02
int SKUA_CYCLE_TEST02(SKUA_DATA *skua_dat)
{
	int success = 0;
	skua_dat->magpie_dat.sys_dat.T = 353.15; 			//K
	skua_dat->magpie_dat.sys_dat.PT = 101.35;			//kPa
	skua_dat->sim_time = 6.0;							//hrs
	skua_dat->t_print = skua_dat->sim_time / 1000.0;	//hrs
	skua_dat->y[0] = 0.7748;					//-
	skua_dat->y[1] = 0.2061;					//-
	skua_dat->y[2] = 0.0191;					//-
	if (skua_dat->gas_dat != NULL)
		skua_dat->gas_dat->CheckMolefractions = false;
	
	//Print to file
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_header(skua_dat);
	}
	
	//Set Initial Conditions
	success = set_SKUA_ICs(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Print ICs
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_results_old(skua_dat);
	}
	
	//Solve a series of time steps implicitly
	do
	{
		if (skua_dat->finch_dat[0].Update == true)
		{
			success = SKUA_reset(skua_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_SKUA_timestep(skua_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Cycling Info Here
		if (skua_dat->t > 2.0)
		{
			/*
			 skua_dat->param_dat[0].y = 0.79;					//-
			 skua_dat->param_dat[1].y = 0.21;					//-
			 skua_dat->param_dat[2].y = 0.0;					//-
			 */
			skua_dat->y[0] = 0.789995;					//-
			skua_dat->y[1] = 0.21;					//-
			skua_dat->y[2] = 5.0E-6;					//-
		}
		
		std::cout << "Evaluating time: " << skua_dat->t << " hrs..." << std::endl;
		success = SKUA_Executioner(skua_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
				skua_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); skua_dat->finch_dat[0].Update = false; return -1;}
		skua_dat->total_steps++;
	} while (skua_dat->t < skua_dat->sim_time);
	
	return success;
}

//Calling the routine to test low pressure, low temperature SKUA
int SKUA_LOW_TEST03(SKUA_DATA *skua_dat)
{
	int success = 0;
	skua_dat->magpie_dat.sys_dat.T = 298.15; 			//K
	skua_dat->magpie_dat.sys_dat.PT = 101.35;			//kPa
	skua_dat->sim_time = 500.0;							//hrs
	skua_dat->t_print = skua_dat->sim_time / 10000.0;	//hrs
	skua_dat->y[0] = 0.79;					//-
	skua_dat->y[1] = 0.21;					//-
	skua_dat->y[2] = 3.07E-6;					//-
	if (skua_dat->gas_dat != NULL)
		skua_dat->gas_dat->CheckMolefractions = false;
	
	//Print to file
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_header(skua_dat);
	}
	
	//Set Initial Conditions
	success = set_SKUA_ICs(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Print ICs
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_results_old(skua_dat);
	}
	
	//Solve a series of time steps implicitly
	do
	{
		if (skua_dat->finch_dat[0].Update == true)
		{
			success = SKUA_reset(skua_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_SKUA_timestep(skua_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		std::cout << "Evaluating time: " << skua_dat->t << " hrs..." << std::endl;
		success = SKUA_Executioner(skua_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
				skua_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); skua_dat->finch_dat[0].Update = false; return -1;}
		skua_dat->total_steps++;
	} while (skua_dat->t < skua_dat->sim_time);
	
	return success;
}

//Mid pressure and temp test
int SKUA_MID_TEST04(SKUA_DATA *skua_dat)
{
	int success = 0;
	skua_dat->magpie_dat.sys_dat.T = 333.15; 			//K
	skua_dat->magpie_dat.sys_dat.PT = 101.35;			//kPa
	skua_dat->sim_time = 23.0;							//hrs
	skua_dat->t_print = skua_dat->sim_time / 1000.0;	//hrs
	skua_dat->y[0] = 0.79;					//-
	skua_dat->y[1] = 0.21;					//-
	skua_dat->y[2] = 1.34E-3;					//-
	if (skua_dat->gas_dat != NULL)
		skua_dat->gas_dat->CheckMolefractions = false;
	
	//Print to file
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_header(skua_dat);
	}
	
	//Set Initial Conditions
	success = set_SKUA_ICs(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Print ICs
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_results_old(skua_dat);
	}
	
	//Solve a series of time steps implicitly
	do
	{
		if (skua_dat->finch_dat[0].Update == true)
		{
			success = SKUA_reset(skua_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_SKUA_timestep(skua_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		std::cout << "Evaluating time: " << skua_dat->t << " hrs..." << std::endl;
		success = SKUA_Executioner(skua_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
				skua_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); skua_dat->finch_dat[0].Update = false; return -1;}
		skua_dat->total_steps++;
	} while (skua_dat->t < skua_dat->sim_time);
	
	return success;
}

//Running a SKUA scenario (UNFINISHED)
int SKUA_SCENARIOS(const char *scene, const char *sorbent, const char *comp, const char *sorbate)
{
	int success = 0;
	std::string sceneName, sorbentName, compName, sorbateName;
	
	//Check the input files
	if (scene == NULL || sorbent == NULL || comp == NULL || sorbate == NULL)
	{
		std::cout << "Enter name of Scenario File: ";
		std::cin >> sceneName;
		std::cout << "Enter name of Adsorbent Property File: ";
		std:: cin >> sorbentName;
		std::cout << "Enter name of Component Property File: ";
		std::cin >> compName;
		std::cout << "Enter name of Adsorbate Property File: ";
		std::cin >> sorbateName;
		std::cout << "\n";
		
		scene = sceneName.c_str();
		sorbent = sorbentName.c_str();
		comp = compName.c_str();
		sorbate = sorbateName.c_str();
	}
	
	std::ifstream sceneFile ( scene );
	std::ifstream sorbentFile ( sorbent );
	std::ifstream compFile ( comp );
	std::ifstream sorbateFile ( sorbate );
	
	if (sceneFile.good()==false || sorbentFile.good()==false || compFile.good()==false || sorbateFile.good()==false)
	{
		mError(file_dne);
		return -1;
	}
	
	//Declarations
	SKUA_DATA dat;
	MIXED_GAS mixture;
	int i_read;
	double d_read;
	std::string s_read;
	double time;
	FILE *Output;
	
	//Initializations
	time = clock();
	Output = fopen("output/SKUA_Scenario_Output.txt","w+");
	if (Output == nullptr)
	{
		system("mkdir output");
		Output = fopen("output/SKUA_Scenario_Output.txt","w+");
	}
	dat.total_steps = 0;
	
	//	(1) - Read the Scenario File
	sceneFile >> d_read; dat.magpie_dat.sys_dat.T = d_read;
	sceneFile >> d_read; dat.magpie_dat.sys_dat.PT = d_read;
	sceneFile >> d_read; dat.gas_velocity = d_read;
	sceneFile >> d_read; dat.sim_time = d_read;
	sceneFile >> d_read; dat.t_print = d_read;
	sceneFile >> i_read;
	if (i_read == 0) dat.DirichletBC = false;
	else if (i_read == 1) dat.DirichletBC = true;
	else {mError(invalid_boolean); return -1;}
	sceneFile >> i_read; dat.magpie_dat.sys_dat.N = i_read;
	dat.magpie_dat.gpast_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.magpie_dat.gsta_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.magpie_dat.mspd_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.param_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.finch_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.y.resize(dat.magpie_dat.sys_dat.N);
	sceneFile >> d_read; dat.magpie_dat.sys_dat.qT = d_read;
	dat.qTn = dat.magpie_dat.sys_dat.qT;
	for (int i=0; i<dat.magpie_dat.sys_dat.N; i++)
	{
		sceneFile >> s_read; dat.param_dat[i].speciesName = s_read;
		sceneFile >> i_read;
		if (i_read == 0) dat.param_dat[i].Adsorbable = false;
		else if (i_read == 1) dat.param_dat[i].Adsorbable = true;
		else {mError(invalid_boolean); return -1;}
		sceneFile >> d_read; dat.y[i] = d_read;
		sceneFile >> d_read; dat.param_dat[i].xIC = d_read;
	}
	
	//Initialize gas mixture data
	success = initialize_data(dat.magpie_dat.sys_dat.N, &mixture);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	std::cout << "Scenario file to simulate surface kinetic adsorption at " << dat.magpie_dat.sys_dat.T << " K and " << dat.magpie_dat.sys_dat.PT << " kPa for " << dat.sim_time << " hours with the following gas species...\n";
	for (int i=0; i<dat.magpie_dat.sys_dat.N; i++)
	{
		//Additional MAGPIE initializations
		dat.magpie_dat.mspd_dat[i].eta.resize(dat.magpie_dat.sys_dat.N);
		dat.magpie_dat.gpast_dat[i].gama_inf.resize(dat.magpie_dat.sys_dat.N);
		dat.magpie_dat.gpast_dat[i].po.resize(dat.magpie_dat.sys_dat.N);
		
		std::cout << "\tSpecies[" << i << "]: " << dat.param_dat[i].speciesName << " at molefraction of " << dat.y[i] << "\t";
		if (dat.param_dat[i].Adsorbable == true)
			std::cout << "Adsorbable species \n";
		else
			std::cout << "Carrier gas component \n";
	}
	std::cout << "\n";
	success = molefractionCheck(&dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	sceneFile.close();
	
	//	(2) - Read the Adsorbent File
	sorbentFile >> i_read; dat.coord = i_read;
	if (dat.coord > 2 || dat.coord < 0)
	{
		mError(invalid_boolean);
		return -1;
	}
	if (dat.coord == 0 || dat.coord == 1)
	{
		sorbentFile >> d_read; dat.char_measure = d_read;
	}
	else
		dat.char_measure = 1.0;
	sorbentFile >> d_read; dat.pellet_radius = d_read;
	sorbentFile.close();
	
	//	(3) - Read the Component File
	for (int i=0; i<dat.magpie_dat.sys_dat.N; i++)
	{
		compFile >> d_read; mixture.species_dat[i].molecular_weight = d_read;
		compFile >> d_read; mixture.species_dat[i].specific_heat = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Viscosity = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Temp = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Const = d_read;
	}
	compFile.close();
	
	//	(4) - Read the Adsorbate File
	sorbateFile >> i_read; if (i_read < 0) {mError(invalid_boolean); return -1;}
	if (i_read == 0) dat.eval_diff = (*default_Dc);
	else if (i_read == 1) dat.eval_diff = (*simple_darken_Dc);
	else if (i_read == 2) dat.eval_diff = (*theoretical_darken_Dc);
	else {mError(invalid_boolean); return -1;}
	for (int i=0; i<dat.magpie_dat.sys_dat.N; i++)
	{
		//Only read in data for adsorbable components in order of appearence
		if (dat.param_dat[i].Adsorbable == true)
		{
			sorbateFile >> d_read; dat.param_dat[i].ref_diffusion = d_read;			//um^2/hr
			sorbateFile >> d_read; dat.param_dat[i].activation_energy = d_read;		//J/mol
			sorbateFile >> d_read; dat.param_dat[i].ref_temperature = d_read;		//K
			sorbateFile >> d_read; dat.param_dat[i].affinity = d_read;				//-
			sorbateFile >> d_read; dat.magpie_dat.mspd_dat[i].v = d_read;				//cm^3/mol
			sorbateFile >> d_read; dat.magpie_dat.gsta_dat[i].qmax = d_read;			//mol/kg
			sorbateFile >> i_read; dat.magpie_dat.gsta_dat[i].m = i_read;				//-
			dat.magpie_dat.gsta_dat[i].dHo.resize(dat.magpie_dat.gsta_dat[i].m);
			dat.magpie_dat.gsta_dat[i].dSo.resize(dat.magpie_dat.gsta_dat[i].m);
			for (int n=0; n<dat.magpie_dat.gsta_dat[i].m; n++)
			{
				sorbateFile >> d_read; dat.magpie_dat.gsta_dat[i].dHo[n] = d_read;	//J/mol
				sorbateFile >> d_read; dat.magpie_dat.gsta_dat[i].dSo[n] = d_read;	//J/K/mol
			}
		}
		//Otherwise, set values to zeros
		else
		{
			dat.param_dat[i].ref_diffusion = 0.0;			//um^2/hr
			dat.param_dat[i].activation_energy = 0.0;		//J/mol
			dat.param_dat[i].ref_temperature = 0.0;			//K
			dat.param_dat[i].affinity = 0.0;				//-
			dat.magpie_dat.mspd_dat[i].v = 0.0;				//cm^3/mol
			dat.magpie_dat.gsta_dat[i].qmax = 0.0;			//mol/kg
			dat.magpie_dat.gsta_dat[i].m = 1;				//-
			dat.magpie_dat.gsta_dat[i].dHo.resize(dat.magpie_dat.gsta_dat[i].m);
			dat.magpie_dat.gsta_dat[i].dSo.resize(dat.magpie_dat.gsta_dat[i].m);
			for (int n=0; n<dat.magpie_dat.gsta_dat[i].m; n++)
			{
				dat.magpie_dat.gsta_dat[i].dHo[n] = 0.0;	//J/mol
				dat.magpie_dat.gsta_dat[i].dSo[n] = 0.0;	//J/K/mol
			}
		}
	}
	sorbateFile.close();
	
	//Setup the SKUA Problem with the correct pointer arguments
	success = setup_SKUA_DATA(Output, dat.eval_diff, NULL, (void *)&dat, &mixture, &dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Call the standard simulation routine
	success = SKUA(&dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Close files
	fclose(Output);
	
	//END execution
	time = clock() - time;
	std::cout << "Simulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds.\n";
	std::cout << "Total Evaluations: " << dat.total_steps << "\n";
	std::cout << "Evaluations/sec: " << dat.total_steps/(time / CLOCKS_PER_SEC) << "\n";
	return success;
}

//Running the SKUA tests
int SKUA_TESTS()
{
	int success = 0;
	
	//Declarations
	SKUA_DATA dat;
	MIXED_GAS mixture;
	double time;
	FILE *TestOutput;
	
	//Initializations
	time = clock();
	TestOutput = fopen("output/SKUA_Test_Output.txt","w+");
	if (TestOutput == nullptr)
	{
		system("mkdir output");
		TestOutput = fopen("output/SKUA_Test_Output.txt","w+");
	}
	dat.Print2File = true;
	dat.NonLinear = true;						//-
	dat.Print2Console = true;					//-
	dat.magpie_dat.sys_dat.Output = false;
	dat.total_steps = 0;
	
	//	(1) - Scenario to Test (with initializations)
	dat.magpie_dat.sys_dat.T = 313.15; 			//K
	dat.magpie_dat.sys_dat.PT = 101.35;			//kPa
	dat.gas_velocity = 0.36;					//cm/s
	dat.sim_time = 10.0;						//hrs
	dat.DirichletBC = true;						//-
	dat.t = 0.0;								//hrs
	dat.t_old = 0.0;							//hrs
	dat.t_print = dat.sim_time / 1000.0;			//hrs
	dat.magpie_dat.sys_dat.N = 3;				//-
	dat.magpie_dat.sys_dat.qT = 0.0;			//mol/kg
	dat.qTn = dat.magpie_dat.sys_dat.qT;		//mol/kg
	dat.magpie_dat.gpast_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.magpie_dat.gsta_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.magpie_dat.mspd_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.param_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.finch_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.y.resize(dat.magpie_dat.sys_dat.N);
	
	dat.param_dat[0].speciesName = "N2";
	dat.param_dat[1].speciesName = "O2";
	dat.param_dat[2].speciesName = "H2O";
	dat.y[0] = 0.787289;				//-
	dat.y[1] = 0.209386;				//-
	dat.y[2] = 0.003325;				//-
	dat.param_dat[0].xIC = 0.0;					//-
	dat.param_dat[1].xIC = 0.0;					//-
	dat.param_dat[2].xIC = 0.0;					//-
	
	success = molefractionCheck(&dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//	(2) - Sorbent to test (with more initializations)
	dat.coord = 2;
	dat.char_measure = 1.0;
	dat.pellet_radius = 1.5;							//um
	
	//	(3) - Components to test
	dat.param_dat[0].Adsorbable = false;		//-
	dat.param_dat[1].Adsorbable = false;		//-
	dat.param_dat[2].Adsorbable = true;			//-
	
	//Initialize gas mixture data
	success = initialize_data(dat.magpie_dat.sys_dat.N, &mixture);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Set the Constants
	mixture.species_dat[0].molecular_weight = 28.016;
	mixture.species_dat[1].molecular_weight = 32.0;
	mixture.species_dat[2].molecular_weight = 18.0;
	
	mixture.species_dat[0].specific_heat = 1.04;
	mixture.species_dat[1].specific_heat = 0.919;
	mixture.species_dat[2].specific_heat = 1.97;
	
	mixture.species_dat[0].Sutherland_Viscosity = 0.0001781;
	mixture.species_dat[1].Sutherland_Viscosity = 0.0002018;
	mixture.species_dat[2].Sutherland_Viscosity = 0.0001043;
	
	mixture.species_dat[0].Sutherland_Temp = 300.55;
	mixture.species_dat[1].Sutherland_Temp = 292.25;
	mixture.species_dat[2].Sutherland_Temp = 298.16;
	
	mixture.species_dat[0].Sutherland_Const = 111.0;
	mixture.species_dat[1].Sutherland_Const = 127.0;
	mixture.species_dat[2].Sutherland_Const = 784.72;
	
	
	//	(4) - Adsorbate to test
	for (int i=0; i<dat.magpie_dat.sys_dat.N; i++)
	{
		dat.magpie_dat.mspd_dat[i].eta.resize(dat.magpie_dat.sys_dat.N);
		dat.magpie_dat.gpast_dat[i].gama_inf.resize(dat.magpie_dat.sys_dat.N);
		dat.magpie_dat.gpast_dat[i].po.resize(dat.magpie_dat.sys_dat.N);
		
		if (i==0)
		{
			dat.param_dat[i].ref_diffusion = 0.0;			//um^2/hr
			dat.param_dat[i].activation_energy = 0.0;		//J/mol
			dat.param_dat[i].ref_temperature = 0.0;			//K
			dat.param_dat[i].affinity = 0.0;				//-
			dat.magpie_dat.mspd_dat[i].v = 0.0;				//cm^3/mol
			dat.magpie_dat.gsta_dat[i].qmax = 0.0;			//mol/kg
			dat.magpie_dat.gsta_dat[i].m = 1;				//-
			dat.magpie_dat.gsta_dat[i].dHo.resize(dat.magpie_dat.gsta_dat[i].m);
			dat.magpie_dat.gsta_dat[i].dSo.resize(dat.magpie_dat.gsta_dat[i].m);
			for (int n=0; n<dat.magpie_dat.gsta_dat[i].m; n++)
			{
				dat.magpie_dat.gsta_dat[i].dHo[n] = 0.0;	//J/mol
				dat.magpie_dat.gsta_dat[i].dSo[n] = 0.0;	//J/K/mol
			}
		}
		else if (i==1)
		{
			dat.param_dat[i].ref_diffusion = 0.0;			//um^2/hr
			dat.param_dat[i].activation_energy = 0.0;		//J/mol
			dat.param_dat[i].ref_temperature = 0.0;			//K
			dat.param_dat[i].affinity = 0.0;				//-
			dat.magpie_dat.mspd_dat[i].v = 0.0;				//cm^3/mol
			dat.magpie_dat.gsta_dat[i].qmax = 0.0;			//mol/kg
			dat.magpie_dat.gsta_dat[i].m = 1;				//-
			dat.magpie_dat.gsta_dat[i].dHo.resize(dat.magpie_dat.gsta_dat[i].m);
			dat.magpie_dat.gsta_dat[i].dSo.resize(dat.magpie_dat.gsta_dat[i].m);
			for (int n=0; n<dat.magpie_dat.gsta_dat[i].m; n++)
			{
				dat.magpie_dat.gsta_dat[i].dHo[n] = 0.0;	//J/mol
				dat.magpie_dat.gsta_dat[i].dSo[n] = 0.0;	//J/K/mol
			}
		}
		else if (i==2)
		{
			//dat.param_dat[i].ref_diffusion = 1.416;			//um^2/hr		(Old Value)
			
			//--------- Diffusion Values for CYCLE TEST 01---------------------------------
			//dat.param_dat[i].ref_diffusion = 0.7;			//um^2/hr			(Const Dc)
			//dat.param_dat[i].ref_diffusion = 0.2;			//um^2/hr			(Simple Darken)
			//dat.param_dat[i].ref_diffusion = 0.3;			//um^2/hr			(Average Darken)
			//dat.param_dat[i].ref_diffusion = 0.05;			//um^2/hr		(Theoretical Darken)
			//dat.param_dat[i].ref_diffusion = 0.25;			//um^2/hr			(Average Theoretical Darken)
			
			
			//--------- Diffusion Values for CYCLE TEST 02---------------------------------
			dat.param_dat[i].ref_diffusion = 4.4;			//um^2/hr			(Const Dc)
			//dat.param_dat[i].ref_diffusion = 1.8;			//um^2/hr			(Simple Darken)
			//dat.param_dat[i].ref_diffusion = 2.5;			//um^2/hr			(Average Darken)
			//dat.param_dat[i].ref_diffusion = 0.75;			//um^2/hr		(Theoretical Darken)
			//dat.param_dat[i].ref_diffusion = 2.0;			//um^2/hr			(Average Theoretical Darken)
			
			//--------- Diffusion Values for LOW TEST 03---------------------------------
			//dat.param_dat[i].ref_diffusion = 0.00981;			//um^2/hr			(Const Dc)
			//dat.param_dat[i].ref_diffusion = 0.00981;			//um^2/hr			(Simple Darken)
			//dat.param_dat[i].ref_diffusion = 0.00981;			//um^2/hr			(Average Darken)
			//dat.param_dat[i].ref_diffusion = 0.00981;			//um^2/hr		(Theoretical Darken)
			//dat.param_dat[i].ref_diffusion = 0.00981;			//um^2/hr			(Average Theoretical Darken)
			
			//--------- Diffusion Values for MID TEST 04---------------------------------
			//dat.param_dat[i].ref_diffusion = 0.4194;			//um^2/hr			(Const Dc)
			//dat.param_dat[i].ref_diffusion = 0.2194;			//um^2/hr			(Simple Darken)
			//dat.param_dat[i].ref_diffusion = 0.3194;			//um^2/hr			(Average Darken)
			//dat.param_dat[i].ref_diffusion = 0.1194;			//um^2/hr		(Theoretical Darken)
			//dat.param_dat[i].ref_diffusion = 0.2694;			//um^2/hr			(Average Theoretical Darken)
			
			
			dat.param_dat[i].activation_energy = 5657.18;	//J/mol
			dat.param_dat[i].ref_temperature = 0.0;			//K
			dat.param_dat[i].affinity = 0.0;				//-
			dat.magpie_dat.mspd_dat[i].v = 13.91;			//cm^3/mol
			dat.magpie_dat.gsta_dat[i].qmax = 11.67;		//mol/kg
			dat.magpie_dat.gsta_dat[i].m = 4;				//-
			dat.magpie_dat.gsta_dat[i].dHo.resize(dat.magpie_dat.gsta_dat[i].m);
			dat.magpie_dat.gsta_dat[i].dSo.resize(dat.magpie_dat.gsta_dat[i].m);
			
			dat.magpie_dat.gsta_dat[i].dHo[0] = -46597.5;	//J/mol
			dat.magpie_dat.gsta_dat[i].dSo[0] = -53.6994;	//J/K/mol
			dat.magpie_dat.gsta_dat[i].dHo[1] = -125024;	//J/mol
			dat.magpie_dat.gsta_dat[i].dSo[1] = -221.073;	//J/K/mol
			dat.magpie_dat.gsta_dat[i].dHo[2] = -193619;	//J/mol
			dat.magpie_dat.gsta_dat[i].dSo[2] = -356.728;	//J/K/mol
			dat.magpie_dat.gsta_dat[i].dHo[3] = -272228;	//J/mol
			dat.magpie_dat.gsta_dat[i].dSo[3] = -567.459;	//J/K/mol
		}
		else {mError(indexing_error); std::cout << "Test case invalid!" << std::endl; return -1;}
		
	}
	
	//Setup the SKUA Problem with the correct pointer arguments
	success = setup_SKUA_DATA(TestOutput, default_Dc, default_kf, (void *)&dat, &mixture, &dat);
	//success = setup_SKUA_DATA(TestOutput, simple_darken_Dc, default_kf, (void *)&dat, &mixture, &dat);
	//success = setup_SKUA_DATA(TestOutput, average_darken_Dc, default_kf, (void *)&dat, &mixture, &dat);
	//success = setup_SKUA_DATA(TestOutput, theoretical_darken_Dc, default_kf, (void *)&dat, &mixture, &dat);
	//success = setup_SKUA_DATA(TestOutput, average_theoretical_darken_Dc, default_kf, (void *)&dat, &mixture, &dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Call the simulation routine
	success = SKUA(&dat);
	//success = SKUA_CYCLE_TEST01(&dat);
	//success = SKUA_CYCLE_TEST02(&dat);
	//success = SKUA_LOW_TEST03(&dat);
	//success = SKUA_MID_TEST04(&dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//END execution
	fclose(TestOutput);
	time = clock() - time;
	std::cout << "Simulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Evaluations: " << dat.total_steps << "\n";
	std::cout << "Evaluations/sec: " << dat.total_steps/(time / CLOCKS_PER_SEC) << "\n";
	return success;
}
