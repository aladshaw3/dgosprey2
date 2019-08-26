/*!
 *  \file scopsowl.cpp scopsowl.h
 *	\brief Simultaneously Coupled Objects for Pore and Surface diffusion Operations With Linear systems
 *
 *  \author Austin Ladshaw
 *	\date 01/29/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "scopsowl.h"

//Creating the header file for SCOPSOWL Output
void print2file_species_header(FILE *Output, SCOPSOWL_DATA *owl_dat, int i)
{
	const char *name = owl_dat->param_dat[i].speciesName.c_str();
	fprintf(Output, "%s\t", name);
	if (owl_dat->finch_dat[i].Dirichlet == true)
		fprintf(Output,"-\t");
	for (int l=0; l<owl_dat->finch_dat[i].LN+2; l++)
	{
		fprintf(Output,"-\t");
	}
}

//Print out the time header for output
void print2file_SCOPSOWL_time_header(FILE *Output, SCOPSOWL_DATA *owl_dat, int i)
{
	fprintf(Output,"Time\t");
	fprintf(Output,"q[0] (mol/kg)\t");
	for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
	{
		if (owl_dat->finch_dat[i].Dirichlet == false && l == owl_dat->finch_dat[i].LN-1)
			break;
		fprintf(Output,"q[%i] (mol/kg)\t",l+1);
	}
	fprintf(Output,"q Average (mol/kg)\tQst Average (J/mol)\t");
}

//Print out the header file for output
void print2file_SCOPSOWL_header(SCOPSOWL_DATA *owl_dat)
{
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		if (owl_dat->param_dat[i].Adsorbable == true)
		{
			print2file_species_header(owl_dat->OutputFile, owl_dat, i);
			print2file_tab(owl_dat->OutputFile, &owl_dat->finch_dat[i]);
		}
	}
	print2file_newline(owl_dat->OutputFile, NULL);
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		if (owl_dat->param_dat[i].Adsorbable == true)
		{
			print2file_dim_header(owl_dat->OutputFile, &owl_dat->finch_dat[i]);
			print2file_tab(owl_dat->OutputFile, &owl_dat->finch_dat[i]);
			print2file_tab(owl_dat->OutputFile, &owl_dat->finch_dat[i]);
			print2file_tab(owl_dat->OutputFile, &owl_dat->finch_dat[i]);
		}
	}
	print2file_newline(owl_dat->OutputFile, NULL);
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		if (owl_dat->param_dat[i].Adsorbable == true)
		{
			print2file_SCOPSOWL_time_header(owl_dat->OutputFile, owl_dat, i);
			print2file_tab(owl_dat->OutputFile, &owl_dat->finch_dat[i]);
		}
	}
	print2file_newline(owl_dat->OutputFile, NULL);
}

//Print out old results
void print2file_SCOPSOWL_result_old(SCOPSOWL_DATA *owl_dat)
{
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		if (owl_dat->param_dat[i].Adsorbable == true)
		{
			fprintf(owl_dat->OutputFile,"%.6g\t",owl_dat->t_old);
			for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
			{
				fprintf(owl_dat->OutputFile,"%.6g\t",owl_dat->param_dat[i].qAvg_old(l,0));
			}
			if (owl_dat->finch_dat[i].Dirichlet == true)
				fprintf(owl_dat->OutputFile,"%.6g\t",owl_dat->param_dat[i].qo);
			fprintf(owl_dat->OutputFile,"%.6g\t%.6g\t",owl_dat->param_dat[i].qIntegralAvg_old,owl_dat->param_dat[i].QstAvg_old);
			print2file_tab(owl_dat->OutputFile, &owl_dat->finch_dat[i]);
		}
	}
	print2file_newline(owl_dat->OutputFile, NULL);
}

//Print out new results
void print2file_SCOPSOWL_result_new(SCOPSOWL_DATA *owl_dat)
{
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		if (owl_dat->param_dat[i].Adsorbable == true)
		{
			fprintf(owl_dat->OutputFile,"%.6g\t",owl_dat->t);
			for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
			{
				fprintf(owl_dat->OutputFile,"%.6g\t",owl_dat->param_dat[i].qAvg(l,0));
			}
			if (owl_dat->finch_dat[i].Dirichlet == true)
				fprintf(owl_dat->OutputFile,"%.6g\t",owl_dat->param_dat[i].qo);
			fprintf(owl_dat->OutputFile,"%.6g\t%.6g\t",owl_dat->param_dat[i].qIntegralAvg,owl_dat->param_dat[i].QstAvg);
			print2file_tab(owl_dat->OutputFile, &owl_dat->finch_dat[i]);
		}
	}
	print2file_newline(owl_dat->OutputFile, NULL);
}

//Default Adsorption Strength function is based on input args
double default_adsorption(int i, int l, const void *user_data)
{
	SCOPSOWL_DATA *dat = (SCOPSOWL_DATA *) user_data;
	double dq_dc = 0.0;
	double eps = sqrt(DBL_EPSILON);
	double qEPS = 0.0;
	double uo = 0.0;
	int success = 0;
	
	//Check if any work needs to be done
	if (dat->param_dat[i].Adsorbable == false)
	{
		if (l < 0)
			dat->param_dat[i].dq_dco = 0.0;
		else
			dat->param_dat[i].dq_dc.edit(l, 0, 0.0);
		return 0.0;
	}
	
	//Boundary value of dq/dc
	if (l < 0)
	{
		//Based on MAGPIE simulations ONLY!
		uo = Cstd(dat->magpie_dat.sys_dat.PT*dat->y[i], dat->magpie_dat.sys_dat.T);
		
		//Perturb the ith concentration
		for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
		{
			if (j == i)
				dat->magpie_dat.gpast_dat[j].y = Pstd(uo+eps,dat->magpie_dat.sys_dat.T)/dat->magpie_dat.sys_dat.PT;
			else if (dat->param_dat[j].Adsorbable == true)
				dat->magpie_dat.gpast_dat[j].y = Pstd(uo,dat->magpie_dat.sys_dat.T)/dat->magpie_dat.sys_dat.PT;
			else
				dat->magpie_dat.gpast_dat[j].y = 0.0;
		}
		
		//Call MAGPIE to evaluate perturbed adsorption
		success = MAGPIE((void *)&dat->magpie_dat);
		if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
		else success = 0;
		dat->total_steps = dat->total_steps + dat->magpie_dat.sys_dat.total_eval;
		
		//Store perturbed result in temp location
		qEPS = dat->magpie_dat.sys_dat.qT * dat->magpie_dat.gpast_dat[i].x;
		
		//Restore original info
		for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
		{
			if (dat->param_dat[j].Adsorbable == true)
				dat->magpie_dat.gpast_dat[j].y = dat->y[j];
			else
				dat->magpie_dat.gpast_dat[j].y = 0.0;
		}
		
		//Call MAGPIE to evaluate adsorption
		success = MAGPIE((void *)&dat->magpie_dat);
		if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
		else success = 0;
		dat->total_steps = dat->total_steps + dat->magpie_dat.sys_dat.total_eval;
		
		//Store the solution
		if (dat->param_dat[i].Adsorbable == true)
		{
			dat->param_dat[i].qo = dat->magpie_dat.sys_dat.qT * dat->magpie_dat.gpast_dat[i].x;
			dat->param_dat[i].Qsto = Qst(Pstd(uo,dat->magpie_dat.sys_dat.T), (void *)&dat->magpie_dat, i);
		}
		else
		{
			dat->param_dat[i].qo = 0.0;
			dat->param_dat[i].Qsto = 0.0;
		}
		
		dq_dc = ( qEPS - dat->param_dat[i].qo ) / eps;
		dat->param_dat[i].dq_dco = dq_dc;
		
	}
	//Nodal values of dq/dc
	else
	{
		if (dat->SurfDiff == true && dat->Heterogeneous == true)
		{
			//Adjust interior pressure with perturbation
			double tempPT = 0.0;
			for (int j=0; j<dat->skua_dat[l].magpie_dat.sys_dat.N; j++)
			{
				if (i == j)
					tempPT = tempPT + Pstd(dat->finch_dat[j].unp1(l,0)+eps, dat->skua_dat[l].magpie_dat.sys_dat.T);
				else
					tempPT = tempPT + Pstd(dat->finch_dat[j].unp1(l,0), dat->skua_dat[l].magpie_dat.sys_dat.T);
			}
			dat->skua_dat[l].magpie_dat.sys_dat.PT = tempPT;
			
			//Perform SKUA simulation
			dat->skua_dat[l].finch_dat[i].dt = dat->finch_dat[i].dt;
			dat->skua_dat[l].finch_dat[i].dt_old = dat->finch_dat[i].dt_old;
			dat->skua_dat[l].finch_dat[i].t = dat->finch_dat[i].t;
			dat->skua_dat[l].finch_dat[i].t_old = dat->finch_dat[i].t_old;
			dat->skua_dat[l].t = dat->skua_dat[l].finch_dat[i].t;
			dat->skua_dat[l].t_old = dat->skua_dat[l].finch_dat[i].t_old;
			
			//Set up SKUA data
			dat->skua_dat[l].magpie_dat.sys_dat.T = dat->magpie_dat.sys_dat.T;
			for (int j=0; j<dat->skua_dat[l].magpie_dat.sys_dat.N; j++)
			{
				if (i == j)
					dat->skua_dat[l].y[j] = Pstd(dat->finch_dat[j].unp1(l,0)+eps,dat->skua_dat[l].magpie_dat.sys_dat.T)/dat->skua_dat[l].magpie_dat.sys_dat.PT;
				else
					dat->skua_dat[l].y[j] = Pstd(dat->finch_dat[j].unp1(l,0),dat->skua_dat[l].magpie_dat.sys_dat.T)/dat->skua_dat[l].magpie_dat.sys_dat.PT;
			}
			
			//Call SKUA executioner to evaluate adsorption
			dat->skua_dat[l].total_steps = 0;
			success = SKUA_Executioner(&dat->skua_dat[l]);
			if (success != 0) {mError(simulation_fail); return -1;}
			dat->total_steps = dat->total_steps + dat->skua_dat[l].total_steps;
			
			//Store the results temporarily
			qEPS = dat->skua_dat[l].finch_dat[i].uAvg;
			
			//Adjust again without perturbation
			tempPT = 0.0;
			for (int j=0; j<dat->skua_dat[l].magpie_dat.sys_dat.N; j++)
				tempPT = tempPT + Pstd(dat->finch_dat[j].unp1(l,0), dat->skua_dat[l].magpie_dat.sys_dat.T);
			dat->skua_dat[l].magpie_dat.sys_dat.PT = tempPT;
			
			//Undo perturbation
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
			{
				dat->skua_dat[l].y[j] = Pstd(dat->finch_dat[j].unp1(l,0),dat->magpie_dat.sys_dat.T)/dat->skua_dat[l].magpie_dat.sys_dat.PT;
			}
			
			//Call SKUA again to evaluate perturbation
			dat->skua_dat[l].total_steps = 0;
			success = SKUA_Executioner(&dat->skua_dat[l]);
			if (success != 0) {mError(simulation_fail); return -1;}
			dat->total_steps = dat->total_steps + dat->skua_dat[l].total_steps;
			
			//Store perturbed result in temp location
			dat->param_dat[i].qAvg.edit(l, 0, dat->skua_dat[l].finch_dat[i].uAvg);
			dat->param_dat[i].Qst.edit(l, 0, dat->skua_dat[l].param_dat[i].Qstnp1);
			dq_dc = ( qEPS - dat->param_dat[i].qAvg(l,0) ) / eps;
			dat->param_dat[i].dq_dc.edit(l, 0, dq_dc);
			
		}
		else
		{
			//Perturb the ith value and form the total pressure
			double tempPT = 0.0;
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
			{
				if (i == j)
					tempPT = tempPT + Pstd(dat->finch_dat[j].unp1(l,0)+eps, dat->magpie_dat.sys_dat.T);
				else
					tempPT = tempPT + Pstd(dat->finch_dat[j].unp1(l,0), dat->magpie_dat.sys_dat.T);
			}
			dat->magpie_dat.sys_dat.PT = tempPT;
			
			//Perturb the ith concentration
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
			{
				if (j == i)
					dat->magpie_dat.gpast_dat[j].y = Pstd(dat->finch_dat[j].unp1(l,0)+eps,dat->magpie_dat.sys_dat.T)/dat->magpie_dat.sys_dat.PT;
				else if (dat->param_dat[j].Adsorbable == true)
					dat->magpie_dat.gpast_dat[j].y = Pstd(dat->finch_dat[j].unp1(l,0),dat->magpie_dat.sys_dat.T)/dat->magpie_dat.sys_dat.PT;
				else
					dat->magpie_dat.gpast_dat[j].y = 0.0;
			}
			
			//Call MAGPIE to evaluate perturbed adsorption
			success = MAGPIE((void *)&dat->magpie_dat);
			if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
			else success = 0;
			dat->total_steps = dat->total_steps + dat->magpie_dat.sys_dat.total_eval;
			
			//Store perturbed result in temp location
			qEPS = dat->magpie_dat.sys_dat.qT * dat->magpie_dat.gpast_dat[i].x;
			
			
			//Adjust interior pressure back to true values
			tempPT = 0.0;
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
				tempPT = tempPT + Pstd(dat->finch_dat[j].unp1(l,0), dat->magpie_dat.sys_dat.T);
			dat->magpie_dat.sys_dat.PT = tempPT;
			
			//Perform MAGPIE simulation
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
			{
				if (dat->param_dat[j].Adsorbable == true)
					dat->magpie_dat.gpast_dat[j].y = Pstd(dat->finch_dat[j].unp1(l,0),dat->magpie_dat.sys_dat.T)/dat->magpie_dat.sys_dat.PT;
				else
					dat->magpie_dat.gpast_dat[j].y = 0.0;
			}
			
			//Call MAGPIE to evaluate adsorption
			success = MAGPIE((void *)&dat->magpie_dat);
			if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
			else success = 0;
			dat->total_steps = dat->total_steps + dat->magpie_dat.sys_dat.total_eval;
			
			//Store the solutions
			if (dat->param_dat[i].Adsorbable == true)
			{
				dat->param_dat[i].qAvg.edit(l, 0, dat->magpie_dat.sys_dat.qT * dat->magpie_dat.gpast_dat[i].x);
				dat->param_dat[i].Qst.edit(l,0,Qst(Pstd(dat->finch_dat[i].unp1(l,0),dat->magpie_dat.sys_dat.T), (void *)&dat->magpie_dat, i));
			}
			else
			{
				dat->param_dat[i].qAvg.edit(l, 0, 0.0);
				dat->param_dat[i].Qst.edit(l,0,0.0);
			}
			
			
			dq_dc = ( qEPS - dat->param_dat[i].qAvg(l,0) ) / eps;
			dat->param_dat[i].dq_dc.edit(l, 0, dq_dc);
			
		}
	}
	
	return dq_dc;
}

//Default Retardation
double default_retardation(int i, int l, const void *user_data)
{
	double Ret = 1.0;
	SCOPSOWL_DATA *dat = (SCOPSOWL_DATA *) user_data;
	double ads_coeff = (dat->pellet_density*(*dat->eval_ads)(i,l,user_data));
	if (l < 0 && dat->DirichletBC == false)
		Ret = (dat->binder_porosity*dat->binder_fraction);
	else
	{
		if (ads_coeff > 0.0)
			Ret = (dat->binder_porosity*dat->binder_fraction) + ads_coeff;
		else
			Ret = (dat->binder_porosity*dat->binder_fraction);
	}
	if (Ret < 0.0)
	{
		mError(unstable_matrix);
		std::cout << "\nNode: " << l << "\tSpecies: " << i << "\tName: " << dat->param_dat[i].speciesName << std::endl;
		std::cout << "NOTE: This error typically indicates a bad gradient approximation for adsorption! Try using a different diffusion function or check the isotherm for negative slopes...\n";
		Ret = -1.0;
	}
	return Ret;
}

//Default Pore Diffusion
double default_pore_diffusion(int i, int l, const void *user_data)
{
	double Dp = 0.0;
	int success = 0;
	SCOPSOWL_DATA *dat = (SCOPSOWL_DATA *) user_data;
	//Check for existance of MIXED GAS DATA
	if (l < 0)
	{
		if (dat->gas_dat != NULL)
		{
			dat->gas_dat->CheckMolefractions = false;
			success = set_variables(dat->magpie_dat.sys_dat.PT, dat->magpie_dat.sys_dat.T, dat->gas_velocity, dat->pellet_radius*2.0, dat->y, dat->gas_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
			//Calculate Properties
			success = calculate_properties(dat->gas_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
			Dp = avgDp(Dp(dat->gas_dat->species_dat[i].molecular_diffusion, dat->binder_porosity), Dk(dat->binder_poresize, dat->magpie_dat.sys_dat.T, dat->gas_dat->species_dat[i].molecular_weight))*3600.0/dat->binder_porosity/dat->binder_porosity;
		}
		else
		{
			Dp = Dp(dat->gas_velocity*dat->pellet_radius, dat->binder_porosity)*3600.0/dat->binder_porosity/dat->binder_porosity;
		}
	}
	else
	{
		if (dat->gas_dat != NULL && dat->SurfDiff == true && dat->Heterogeneous == true)
		{
			double tempPT = 0.0;
			for (int j=0; j<dat->skua_dat[l].magpie_dat.sys_dat.N; j++)
				tempPT = tempPT + Pstd(dat->finch_dat[j].unp1(l,0), dat->skua_dat[l].magpie_dat.sys_dat.T);
			dat->skua_dat[l].magpie_dat.sys_dat.PT = tempPT;
			for (int j=0; j<dat->skua_dat[l].magpie_dat.sys_dat.N; j++)
			{
				dat->skua_dat[l].y[j] = Pstd(dat->finch_dat[j].unp1(l,0), dat->skua_dat[l].magpie_dat.sys_dat.T) / dat->skua_dat[l].magpie_dat.sys_dat.PT;
			}
			dat->gas_dat->CheckMolefractions = false;
			success = set_variables(dat->skua_dat[l].magpie_dat.sys_dat.PT, dat->skua_dat[l].magpie_dat.sys_dat.T, dat->gas_velocity, dat->pellet_radius*2.0, dat->skua_dat[l].y, dat->gas_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
			//Calculate Properties
			success = calculate_properties(dat->gas_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
			Dp = avgDp(Dp(dat->gas_dat->species_dat[i].molecular_diffusion, dat->binder_porosity), Dk(dat->binder_poresize, dat->skua_dat[l].magpie_dat.sys_dat.T, dat->gas_dat->species_dat[i].molecular_weight))*3600.0/dat->binder_porosity/dat->binder_porosity;
		}
		else if (dat->gas_dat != NULL)
		{
			double tempPT = 0.0;
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
				tempPT = tempPT + Pstd(dat->finch_dat[j].unp1(l,0), dat->magpie_dat.sys_dat.T);
			dat->magpie_dat.sys_dat.PT = tempPT;
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
			{
				dat->tempy[j] = Pstd(dat->finch_dat[j].unp1(l,0), dat->magpie_dat.sys_dat.T) / dat->magpie_dat.sys_dat.PT;
			}
			dat->gas_dat->CheckMolefractions = false;
			success = set_variables(dat->magpie_dat.sys_dat.PT, dat->magpie_dat.sys_dat.T, dat->gas_velocity, dat->pellet_radius*2.0, dat->tempy, dat->gas_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
			//Calculate Properties
			success = calculate_properties(dat->gas_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
			Dp = avgDp(Dp(dat->gas_dat->species_dat[i].molecular_diffusion, dat->binder_porosity), Dk(dat->binder_poresize, dat->magpie_dat.sys_dat.T, dat->gas_dat->species_dat[i].molecular_weight))*3600.0/dat->binder_porosity/dat->binder_porosity;
		}
		else
		{
			Dp = Dp(dat->gas_velocity*dat->pellet_radius, dat->binder_porosity)*3600.0/dat->binder_porosity/dat->binder_porosity;
		}
	}
	
	return dat->binder_porosity*dat->binder_fraction*Dp;
}

//Function to evaluate the default surface diffusion
double default_surf_diffusion(int i, int l, const void *user_data)
{
	double Dc = 0.0;
	SCOPSOWL_DATA *dat = (SCOPSOWL_DATA *) user_data;
	dat->magpie_dat.sys_dat.T = dat->gas_temperature;
	Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
	Dc = D_inf(Dc, dat->param_dat[i].ref_temperature, dat->param_dat[i].affinity, dat->param_dat[i].ref_pressure, dat->magpie_dat.sys_dat.T);
	return Dc;
}

//Function to return 0 as surface diffusion
double zero_surf_diffusion(int i, int l, const void *user_data)
{
	return 0.0;
}

//Function to evaluate an effective diffusivity based on pore and surface diffusion
double default_effective_diffusion(int i, int l, const void *user_data)
{
	double De = 0.0;
	double Dc = 0.0, Dp = 0.0;
	SCOPSOWL_DATA *dat = (SCOPSOWL_DATA *) user_data;
	
	Dp = default_pore_diffusion(i, l, user_data);
	Dc = (*dat->eval_surfDiff) (i,l,user_data);
	if (l < 0)
		De = (Dp*dat->binder_porosity*dat->binder_porosity) + (dat->pellet_density*Dc*dat->param_dat[i].dq_dco*1.0E-8);
	else
		De = (Dp*dat->binder_porosity*dat->binder_porosity) + (dat->pellet_density*Dc*dat->param_dat[i].dq_dc(l,0)*1.0E-8);
	
	return De;
}

//Function to set pore diffusion to a constant
double const_pore_diffusion(int i, int l, const void *user_data)
{
	SCOPSOWL_DATA *dat = (SCOPSOWL_DATA *) user_data;
	return dat->param_dat[i].pore_diffusion;
}

//Default Film Mass Transfer
double default_filmMassTransfer(int i, const void *user_data)
{
	int success = 0;
	SCOPSOWL_DATA *dat = (SCOPSOWL_DATA *) user_data;
	if (dat->gas_dat != NULL)
	{
		success = set_variables(dat->magpie_dat.sys_dat.PT, dat->magpie_dat.sys_dat.T, dat->gas_velocity, dat->pellet_radius, dat->y, dat->gas_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		//Calculate Properties
		success = calculate_properties(dat->gas_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		return FilmMTCoeff(dat->gas_dat->species_dat[i].molecular_diffusion, dat->gas_dat->char_length, dat->gas_dat->Reynolds, dat->gas_dat->species_dat[i].Schmidt) * 3600.0;
	}
	else
	{
		std::cout << "\nWARNING!!! Calling a Function without the correct pointer initialized!\n" << std::endl;
		return default_kf(i, user_data); //This will return 0.0 and result in errors
	}
}

//Function to set film mass transfer to a constant
double const_filmMassTransfer(int i, const void *user_data)
{
	SCOPSOWL_DATA *dat = (SCOPSOWL_DATA *) user_data;
	return dat->param_dat[i].film_transfer;
}

//Setup of the SCOPSOWL Object
int setup_SCOPSOWL_DATA(FILE *file,
						double (*eval_sorption) (int i, int l, const void *user_data),
						double (*eval_retardation) (int i, int l, const void *user_data),
						double (*eval_pore_diff) (int i, int l, const void *user_data),
						double (*eval_filmMT) (int i, const void *user_data),
						double (*eval_surface_diff) (int i, int l, const void *user_data),
						const void *user_data,MIXED_GAS *gas_data,SCOPSOWL_DATA *owl_data)
{
	int success = 0;
	
	//Check input args
	if (file == NULL)
		owl_data->Print2File = false;
	else
	{
		owl_data->Print2File = true;
		owl_data->OutputFile = file;
	}
	if ((*eval_sorption) == NULL)
	{
		owl_data->eval_ads = (*default_adsorption);
	}
	else
	{
		owl_data->eval_ads = (*eval_sorption);
	}
	if ((*eval_retardation) == NULL)
		owl_data->eval_retard = (*default_retardation);
	else
		owl_data->eval_retard = (*eval_retardation);
	if ((*eval_pore_diff) == NULL)
		owl_data->eval_diff = (*default_pore_diffusion);
	else
		owl_data->eval_diff = (*eval_pore_diff);
	if ((*eval_filmMT) == NULL)
		owl_data->eval_kf = (*default_filmMassTransfer);
	else
		owl_data->eval_kf = (*eval_filmMT);
	if (gas_data == NULL)
	{
		owl_data->DirichletBC = true;
	}
	else
	{
		owl_data->gas_dat = gas_data;
	}
	owl_data->user_data = user_data;
	
	//Setup the finch problems for SCOPSOWL
	owl_data->t = 0.0; owl_data->t_old = 0.0;
	owl_data->t_counter = 0.0;
	owl_data->total_steps = 0;
	owl_data->tempy.resize(owl_data->magpie_dat.sys_dat.N);
	owl_data->magpie_dat.sys_dat.T = owl_data->gas_temperature;
	owl_data->magpie_dat.sys_dat.PT = owl_data->total_pressure;
	if (owl_data->NonLinear == false)
		owl_data->level = 1;
	if (owl_data->level <= 0)
		owl_data->level = 2;
	if (owl_data->level > 5)
		owl_data->level = 5;
	for (int i=0; i<owl_data->magpie_dat.sys_dat.N; i++)
	{
		owl_data->finch_dat[i].L = owl_data->pellet_radius;
		owl_data->finch_dat[i].T = owl_data->sim_time;
		owl_data->finch_dat[i].Dirichlet = owl_data->DirichletBC;
		owl_data->finch_dat[i].LN = 10;
		owl_data->finch_dat[i].SteadyState = false;
		owl_data->finch_dat[i].CheckMass = true;
		if (owl_data->param_dat[i].Adsorbable == true)
			owl_data->finch_dat[i].Iterative = owl_data->NonLinear;
		else
			owl_data->finch_dat[i].Iterative = false;
		owl_data->finch_dat[i].NormTrack = owl_data->Print2Console;
		owl_data->finch_dat[i].nl_method = FINCH_Picard;
		owl_data->finch_dat[i].d = owl_data->coord_macro;
		owl_data->finch_dat[i].s = owl_data->char_macro;
		owl_data->finch_dat[i].t = 0.0;
		owl_data->finch_dat[i].t_old = 0.0;
		owl_data->finch_dat[i].vo = 0.0;						//cm/hr
		owl_data->finch_dat[i].vIC = 0.0;						//cm/hr
		owl_data->finch_dat[i].kIC = 0.0;						//1/hr
		owl_data->finch_dat[i].ko = 0.0;						//1/hr
		
		success=setup_FINCH_DATA(NULL,NULL,NULL,NULL,NULL,set_SCOPSOWL_params,NULL,NULL,NULL,NULL,NULL,NULL,&owl_data->finch_dat[i],user_data);
		if (success != 0) {mError(simulation_fail); return -1;}
		owl_data->finch_dat[i].max_iter = 300;
		owl_data->finch_dat[i].nl_method = FINCH_Picard;
		owl_data->finch_dat[i].pjfnk_dat.LineSearch = false;
		owl_data->finch_dat[i].pjfnk_dat.Bounce = false;
		owl_data->finch_dat[i].pjfnk_dat.linear_solver = GMRESRP;
		owl_data->finch_dat[i].pjfnk_dat.nl_maxit = 300;
		
		owl_data->param_dat[i].qAvg.set_size(owl_data->finch_dat[i].LN, 1);
		owl_data->param_dat[i].qAvg_old.set_size(owl_data->finch_dat[i].LN, 1);
		owl_data->param_dat[i].Qst.set_size(owl_data->finch_dat[i].LN, 1);
		owl_data->param_dat[i].Qst_old.set_size(owl_data->finch_dat[i].LN, 1);
		owl_data->param_dat[i].dq_dc.set_size(owl_data->finch_dat[i].LN, 1);
	}
	
	//Setup the SKUA problem
	if (owl_data->SurfDiff == true && owl_data->Heterogeneous == true)
	{
		owl_data->skua_dat.resize(owl_data->finch_dat[0].LN);
		for (int l=0; l<owl_data->finch_dat[0].LN; l++)
		{
			owl_data->skua_dat[l].param_dat.resize(owl_data->magpie_dat.sys_dat.N);
			owl_data->skua_dat[l].finch_dat.resize(owl_data->magpie_dat.sys_dat.N);
			owl_data->skua_dat[l].y.resize(owl_data->magpie_dat.sys_dat.N);
			owl_data->skua_dat[l].magpie_dat = owl_data->magpie_dat;
			owl_data->skua_dat[l].pellet_radius = owl_data->crystal_radius;
			owl_data->skua_dat[l].sim_time = owl_data->sim_time;
			owl_data->skua_dat[l].DirichletBC = true;
			owl_data->skua_dat[l].NonLinear = owl_data->NonLinear;
			owl_data->skua_dat[l].Print2Console = false;
			owl_data->skua_dat[l].Print2File = false;
			owl_data->skua_dat[l].coord = owl_data->coord_micro;
			owl_data->skua_dat[l].char_measure = owl_data->char_micro;
			
			success = setup_SKUA_DATA(NULL,eval_surface_diff, NULL, (void *)&owl_data->skua_dat[l], gas_data, &owl_data->skua_dat[l]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			for (int i=0; i<owl_data->magpie_dat.sys_dat.N; i++)
			{
				owl_data->skua_dat[l].param_dat[i].Adsorbable = owl_data->param_dat[i].Adsorbable;
				owl_data->skua_dat[l].param_dat[i].speciesName = owl_data->param_dat[i].speciesName;
				owl_data->skua_dat[l].param_dat[i].ref_diffusion = owl_data->param_dat[i].ref_diffusion;
				owl_data->skua_dat[l].param_dat[i].activation_energy = owl_data->param_dat[i].activation_energy;
				owl_data->skua_dat[l].param_dat[i].ref_temperature = owl_data->param_dat[i].ref_temperature;
				owl_data->skua_dat[l].param_dat[i].affinity = owl_data->param_dat[i].affinity;
			}
		}
	}
	else if (owl_data->SurfDiff == true && owl_data->Heterogeneous == false)
	{
		owl_data->eval_diff = (*default_effective_diffusion);
		owl_data->eval_surfDiff = (*eval_surface_diff);
	}
	else if (owl_data->SurfDiff == false && owl_data->Heterogeneous == false)
	{
		owl_data->eval_diff = (*default_effective_diffusion);
		owl_data->eval_surfDiff = (*zero_surf_diffusion);
	}
	else
	{
		//No Action
	}
	
	return success;
}

//Specify the order of execution of each function
int SCOPSOWL_Executioner(SCOPSOWL_DATA *owl_dat)
{
	int success = 0;
	
	//Perform Preprocess Actions
	success = SCOPSOWL_preprocesses(owl_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Loop for each coupling level
	for (int n=0; n<owl_dat->level; n++)
	{
		//Loop for all components
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			
			//Solve the system
			success = (*owl_dat->finch_dat[i].solve) ((void *)&owl_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Check for negative concentrations
			success = check_Mass(&owl_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Form Totals
			success = uTotal(&owl_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Form Averages
			success = uAverage(&owl_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
		}
		
		if (n < owl_dat->level-1 && owl_dat->Print2File == true)
		{
			owl_dat->Print2File = false;
			success = SCOPSOWL_postprocesses(owl_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
			owl_dat->Print2File = true;
		}
		else
		{
			//Perform Postprocess Actions
			success = SCOPSOWL_postprocesses(owl_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
	}
	
	return success;
}

//Set up the ICs for the problem
int set_SCOPSOWL_ICs(SCOPSOWL_DATA *owl_dat)
{
	int success = 0;
	
	//Print file header
	if (owl_dat->Print2File == true)
		print2file_SCOPSOWL_header(owl_dat);
	
	owl_dat->magpie_dat.sys_dat.PT = owl_dat->total_pressure;
	owl_dat->magpie_dat.sys_dat.T = owl_dat->gas_temperature;
	if (owl_dat->SurfDiff == true && owl_dat->Heterogeneous == true)
	{
		//Call SKUA ICs to establish SCOPSOWL ICs
		for (int l=0; l<owl_dat->finch_dat[0].LN; l++)
		{
			owl_dat->skua_dat[l].total_steps = 0;
			owl_dat->skua_dat[l].magpie_dat.sys_dat.T = owl_dat->gas_temperature;
			owl_dat->skua_dat[l].magpie_dat.sys_dat.PT = owl_dat->total_pressure;
			owl_dat->skua_dat[l].magpie_dat.sys_dat.qT = owl_dat->magpie_dat.sys_dat.qT;
			for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
			{
				owl_dat->skua_dat[l].param_dat[i].ref_diffusion = D_inf(owl_dat->param_dat[i].ref_diffusion, owl_dat->param_dat[i].ref_temperature, owl_dat->param_dat[i].affinity, owl_dat->y[i]*owl_dat->total_pressure, owl_dat->gas_temperature);
				owl_dat->skua_dat[l].param_dat[i].activation_energy = owl_dat->param_dat[i].activation_energy;
				owl_dat->skua_dat[l].param_dat[i].ref_temperature = 0.0;
				owl_dat->skua_dat[l].param_dat[i].affinity = 0.0;
				
				owl_dat->skua_dat[l].param_dat[i].xIC = owl_dat->param_dat[i].xIC;
			}
			success = set_SKUA_ICs(&owl_dat->skua_dat[l]);
			if (success != 0) {mError(simulation_fail); return -1;}
			owl_dat->total_steps = owl_dat->total_steps + owl_dat->skua_dat[l].total_steps;
			
			//Loop for species
			for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
			{
				owl_dat->param_dat[i].qAvg.edit(l, 0, owl_dat->skua_dat[l].finch_dat[i].uAvg);
				owl_dat->param_dat[i].Qst.edit(l, 0, owl_dat->skua_dat[l].param_dat[i].Qstn);
				owl_dat->param_dat[i].qAvg_old.edit(l, 0, owl_dat->param_dat[i].qAvg(l,0));
				owl_dat->param_dat[i].Qst_old.edit(l, 0, owl_dat->param_dat[i].Qst(l,0));
				owl_dat->magpie_dat.gpast_dat[i].y = owl_dat->skua_dat[l].magpie_dat.gpast_dat[i].y;
			}
		}
		
		//Outer species loop
		double yad_sum = 0.0;
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			//Now form averages and totals
			owl_dat->param_dat[i].Qsto = owl_dat->param_dat[i].Qst(0,0);
			owl_dat->param_dat[i].qo = owl_dat->param_dat[i].qAvg(owl_dat->finch_dat[i].LN-1,0);
			owl_dat->param_dat[i].qIntegralAvg = owl_dat->param_dat[i].qAvg.sphericalAvg(owl_dat->finch_dat[i].L,owl_dat->finch_dat[i].dz,owl_dat->param_dat[i].qo,owl_dat->DirichletBC);
			owl_dat->param_dat[i].QstAvg = owl_dat->param_dat[i].Qst.sphericalAvg(owl_dat->finch_dat[i].L, owl_dat->finch_dat[i].dz,owl_dat->param_dat[i].Qsto,owl_dat->DirichletBC);
			owl_dat->param_dat[i].QstAvg_old = owl_dat->param_dat[i].QstAvg;
			owl_dat->param_dat[i].qIntegralAvg_old = owl_dat->param_dat[i].qIntegralAvg;
			
			//SCOPSOWL ICs
			owl_dat->finch_dat[i].Sn.ConstantICFill(0.0);
			owl_dat->finch_dat[i].Snp1.ConstantICFill(0.0);
			owl_dat->finch_dat[i].vn.ConstantICFill(0.0);
			owl_dat->finch_dat[i].vnp1.ConstantICFill(0.0);
			owl_dat->finch_dat[i].kn.ConstantICFill(0.0);
			owl_dat->finch_dat[i].knp1.ConstantICFill(0.0);
			owl_dat->finch_dat[i].ko = 0.0;
			owl_dat->finch_dat[i].vo = 0.0;
			
			owl_dat->finch_dat[i].kfn = (*owl_dat->eval_kf) (i,owl_dat);
			owl_dat->finch_dat[i].kfnp1 = owl_dat->finch_dat[i].kfn;
			if (owl_dat->finch_dat[i].kfn < 0.0)
			{
				mError(simulation_fail);
				return -1;
			}
			
			if (owl_dat->param_dat[i].Adsorbable == true)
			{
				owl_dat->finch_dat[i].un.ConstantICFill(Cstd(owl_dat->magpie_dat.sys_dat.PT*owl_dat->magpie_dat.gpast_dat[i].y,owl_dat->magpie_dat.sys_dat.T));
			}
			
			if (owl_dat->param_dat[i].Adsorbable == true)
				yad_sum = yad_sum + owl_dat->magpie_dat.gpast_dat[i].y;
		}
		
		//Need to re-adjust the initial conditions for non-adsorbing species
		double rat_frac = 0.0;
		int first = 0;
		bool changed = false;
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			if (owl_dat->param_dat[i].Adsorbable == false)
			{
				if (changed == false)
				{
					first = i; changed = true;
					rat_frac = rat_frac + (owl_dat->y[i] / owl_dat->y[first]);
				}
				else
				{
					rat_frac = rat_frac + (owl_dat->y[i] / owl_dat->y[first]);
				}
			}
		}
		double yFirst = (1.0 - yad_sum) / rat_frac;
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			double yi = 0.0;
			if (owl_dat->param_dat[i].Adsorbable == false)
			{
				yi = yFirst * (owl_dat->y[i] / owl_dat->y[first]);
				owl_dat->finch_dat[i].un.ConstantICFill(Cstd(owl_dat->magpie_dat.sys_dat.PT*yi,owl_dat->magpie_dat.sys_dat.T));
			}
			owl_dat->finch_dat[i].uo = owl_dat->finch_dat[i].un(0,0);
			owl_dat->finch_dat[i].unp1 = owl_dat->finch_dat[i].un;
			owl_dat->finch_dat[i].unm1 = owl_dat->finch_dat[i].un;
			
			//Form Totals
			success = uTotal(&owl_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Form Averages
			success = uAverage(&owl_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		
		//Print out ICs
		if (owl_dat->Print2File == true)
			print2file_SCOPSOWL_result_old(owl_dat);
		
		//Call parameter functions
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
			{
				double p_out = owl_dat->y[i]*owl_dat->total_pressure;
				double p_in = Pstd(owl_dat->finch_dat[i].un(l,0), owl_dat->gas_temperature);
				if (p_out >= p_in)	owl_dat->param_dat[i].ref_pressure = p_out;
				else owl_dat->param_dat[i].ref_pressure = p_in;
				owl_dat->skua_dat[l].param_dat[i].ref_diffusion = D_inf(owl_dat->param_dat[i].ref_diffusion, owl_dat->param_dat[i].ref_temperature, owl_dat->param_dat[i].affinity, owl_dat->param_dat[i].ref_pressure, owl_dat->gas_temperature);
				owl_dat->skua_dat[l].param_dat[i].activation_energy = owl_dat->param_dat[i].activation_energy;
				owl_dat->skua_dat[l].param_dat[i].ref_temperature = 0.0;
				owl_dat->skua_dat[l].param_dat[i].affinity = 0.0;
				
				//Now Recall skua parameter functions to apply corrections
				owl_dat->skua_dat[l].finch_dat[i].Do = (*owl_dat->skua_dat[l].eval_diff) (i, -1, (void *)&owl_dat->skua_dat[l]);
				owl_dat->skua_dat[l].finch_dat[i].Dn.ConstantICFill(owl_dat->skua_dat[l].finch_dat[i].Do);
				owl_dat->skua_dat[l].finch_dat[i].Dnp1.ConstantICFill(owl_dat->skua_dat[l].finch_dat[i].Do);
				
				owl_dat->finch_dat[i].Rn.edit(l, 0, (*owl_dat->eval_retard) (i,l,owl_dat));
				owl_dat->finch_dat[i].Dn.edit(l, 0, (*owl_dat->eval_diff) (i,l,owl_dat));
			}
			owl_dat->finch_dat[i].Ro = (*owl_dat->eval_retard) (i,-1,owl_dat);
			owl_dat->finch_dat[i].Do = (*owl_dat->eval_diff) (i,-1,owl_dat);
			owl_dat->finch_dat[i].Rnp1 = owl_dat->finch_dat[i].Rn;
			owl_dat->finch_dat[i].Dnp1 = owl_dat->finch_dat[i].Dn;
			
			//Update averages and totals
			owl_dat->finch_dat[i].uT_old = owl_dat->finch_dat[i].uT;
			owl_dat->finch_dat[i].uAvg_old = owl_dat->finch_dat[i].uAvg;
		}
		
	}
	else
	{
		//Call MAGPIE to establish SCOPSOWL ICs
		owl_dat->magpie_dat.sys_dat.Recover = true;
		owl_dat->magpie_dat.sys_dat.Carrier = false;
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			if (owl_dat->param_dat[i].Adsorbable == false)
			{
				owl_dat->magpie_dat.sys_dat.Carrier = true;
				break;
			}
		}
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
			owl_dat->magpie_dat.gpast_dat[i].x = owl_dat->param_dat[i].xIC;
		
		//Temperature, pressure, total adsorption, and x's are set: Call MAGPIE
		success = MAGPIE((void *)&owl_dat->magpie_dat);
		if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
		else success = 0;
		owl_dat->total_steps = owl_dat->total_steps + owl_dat->magpie_dat.sys_dat.total_eval;
		
		//Reset the MAGPIE parameters before continuing
		owl_dat->magpie_dat.sys_dat.Recover = false;
		owl_dat->magpie_dat.sys_dat.Carrier = false;
		
		//The GPAST structure of MAGPIE now holds the y's which establish the SCOPSOWL ICs
		double yad_sum = 0.0;
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			owl_dat->param_dat[i].qAvg.ConstantICFill(owl_dat->magpie_dat.sys_dat.qT * owl_dat->magpie_dat.gpast_dat[i].x);
			owl_dat->param_dat[i].qAvg_old = owl_dat->param_dat[i].qAvg;
			owl_dat->param_dat[i].qo = owl_dat->magpie_dat.sys_dat.qT * owl_dat->magpie_dat.gpast_dat[i].x;
			
			owl_dat->finch_dat[i].Sn.ConstantICFill(0.0);
			owl_dat->finch_dat[i].Snp1.ConstantICFill(0.0);
			owl_dat->finch_dat[i].vn.ConstantICFill(0.0);
			owl_dat->finch_dat[i].vnp1.ConstantICFill(0.0);
			owl_dat->finch_dat[i].kn.ConstantICFill(0.0);
			owl_dat->finch_dat[i].knp1.ConstantICFill(0.0);
			owl_dat->finch_dat[i].ko = 0.0;
			owl_dat->finch_dat[i].vo = 0.0;
			
			owl_dat->finch_dat[i].kfn = (*owl_dat->eval_kf) (i,owl_dat);
			owl_dat->finch_dat[i].kfnp1 = owl_dat->finch_dat[i].kfn;
			if (owl_dat->finch_dat[i].kfn < 0.0)
			{
				mError(simulation_fail);
				return -1;
			}
			
			if (owl_dat->param_dat[i].Adsorbable == true)
			{
				owl_dat->finch_dat[i].un.ConstantICFill(Cstd(owl_dat->magpie_dat.sys_dat.PT*owl_dat->magpie_dat.gpast_dat[i].y,owl_dat->magpie_dat.sys_dat.T));
			}
			
			if (owl_dat->param_dat[i].Adsorbable == true)
				yad_sum = yad_sum + owl_dat->magpie_dat.gpast_dat[i].y;
		}
		
		//Need to re-adjust the initial conditions for non-adsorbing species
		double rat_frac = 0.0;
		int first = 0;
		bool changed = false;
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			if (owl_dat->param_dat[i].Adsorbable == false)
			{
				if (changed == false)
				{
					first = i; changed = true;
					rat_frac = rat_frac + (owl_dat->y[i] / owl_dat->y[first]);
				}
				else
				{
					rat_frac = rat_frac + (owl_dat->y[i] / owl_dat->y[first]);
				}
			}
		}
		double yFirst = (1.0 - yad_sum) / rat_frac;
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			double yi = 0.0;
			if (owl_dat->param_dat[i].Adsorbable == false)
			{
				yi = yFirst * (owl_dat->y[i] / owl_dat->y[first]);
				owl_dat->finch_dat[i].un.ConstantICFill(Cstd(owl_dat->magpie_dat.sys_dat.PT*yi,owl_dat->magpie_dat.sys_dat.T));
			}
			owl_dat->finch_dat[i].uo = owl_dat->finch_dat[i].un(0,0);
			owl_dat->finch_dat[i].unp1 = owl_dat->finch_dat[i].un;
			owl_dat->finch_dat[i].unm1 = owl_dat->finch_dat[i].un;
			
			//Form Totals
			success = uTotal(&owl_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Form Averages
			success = uAverage(&owl_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		
		//Use established info to call the retardation and diffusion functions as well as heats of adsorption
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			if (owl_dat->param_dat[i].Adsorbable == true)
			{
				owl_dat->param_dat[i].Qst.ConstantICFill(Qst(Pstd(owl_dat->finch_dat[i].un(0,0),owl_dat->magpie_dat.sys_dat.T), (void *)&owl_dat->magpie_dat, i));
			}
			else
			{
				owl_dat->param_dat[i].Qst.ConstantICFill(0.0);
			}
			owl_dat->param_dat[i].Qst_old = owl_dat->param_dat[i].Qst;
			owl_dat->param_dat[i].Qsto = owl_dat->param_dat[i].Qst(0,0);
			
			//Now form averages and totals
			owl_dat->param_dat[i].qIntegralAvg = owl_dat->param_dat[i].qAvg.sphericalAvg(owl_dat->finch_dat[i].L,owl_dat->finch_dat[i].dz,owl_dat->param_dat[i].qo,owl_dat->DirichletBC);
			owl_dat->param_dat[i].QstAvg = owl_dat->param_dat[i].Qst.sphericalAvg(owl_dat->finch_dat[i].L, owl_dat->finch_dat[i].dz,owl_dat->param_dat[i].Qsto,owl_dat->DirichletBC);
			
			//Update adsorption info
			owl_dat->param_dat[i].QstAvg_old = owl_dat->param_dat[i].QstAvg;
			owl_dat->param_dat[i].qIntegralAvg_old = owl_dat->param_dat[i].qIntegralAvg;
		}
		
		//Print out ICs
		if (owl_dat->Print2File == true)
			print2file_SCOPSOWL_result_old(owl_dat);
		
		//Call the parameter functions
		for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
		{
			//Set the reference pressure
			for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
			{
				double p_out = owl_dat->y[i]*owl_dat->total_pressure;
				double p_in = Pstd(owl_dat->finch_dat[i].un(l,0), owl_dat->gas_temperature);
				if (p_out >= p_in)	owl_dat->param_dat[i].ref_pressure = p_out;
				else owl_dat->param_dat[i].ref_pressure = p_in;
			}
			for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
			{
				owl_dat->finch_dat[i].Rn.edit(l, 0, (*owl_dat->eval_retard) (i,l,owl_dat));
				owl_dat->finch_dat[i].Dn.edit(l, 0, (*owl_dat->eval_diff) (i,l,owl_dat));
			}
			owl_dat->finch_dat[i].Ro = (*owl_dat->eval_retard) (i,-1,owl_dat);
			owl_dat->finch_dat[i].Do = (*owl_dat->eval_diff) (i,-1,owl_dat);
			owl_dat->finch_dat[i].Rnp1 = owl_dat->finch_dat[i].Rn;
			owl_dat->finch_dat[i].Dnp1 = owl_dat->finch_dat[i].Dn;
			
			//Update averages and totals
			owl_dat->finch_dat[i].uT_old = owl_dat->finch_dat[i].uT;
			owl_dat->finch_dat[i].uAvg_old = owl_dat->finch_dat[i].uAvg;
		}
		
	}
	
	return success;
}

//Determine the current time step size
int set_SCOPSOWL_timestep(SCOPSOWL_DATA *owl_dat)
{
	int success = 0;
	
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		owl_dat->finch_dat[i].dt = owl_dat->finch_dat[i].dz / 4.0;
		owl_dat->finch_dat[i].t = owl_dat->finch_dat[i].dt + owl_dat->finch_dat[i].t_old;
		
		if (owl_dat->SurfDiff == true && owl_dat->Heterogeneous == true)
		{
			for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
			{
				owl_dat->skua_dat[l].finch_dat[i].dt = owl_dat->finch_dat[i].dt;
				owl_dat->skua_dat[l].finch_dat[i].t = owl_dat->finch_dat[i].t;
				owl_dat->skua_dat[l].t_old = owl_dat->finch_dat[i].t_old;
				owl_dat->skua_dat[l].t = owl_dat->finch_dat[i].t;
			}
		}
	}
	owl_dat->t_old = owl_dat->finch_dat[0].t_old;
	owl_dat->t = owl_dat->finch_dat[0].t;
	
	return success;
}

//Evaluate any routines prior to solving the system
int SCOPSOWL_preprocesses(SCOPSOWL_DATA *owl_dat)
{
	int success = 0;
	
	//Establish boundary conditions before simulation starts
	owl_dat->magpie_dat.sys_dat.T = owl_dat->gas_temperature;
	owl_dat->magpie_dat.sys_dat.PT = owl_dat->total_pressure;
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		if (owl_dat->SurfDiff == true && owl_dat->Heterogeneous == true)
		{
			for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
			{
				double p_out = owl_dat->y[i]*owl_dat->total_pressure;
				if (p_out >= owl_dat->param_dat[i].ref_pressure) owl_dat->param_dat[i].ref_pressure = p_out;
				owl_dat->skua_dat[l].param_dat[i].ref_diffusion = D_inf(owl_dat->param_dat[i].ref_diffusion, owl_dat->param_dat[i].ref_temperature, owl_dat->param_dat[i].affinity, owl_dat->param_dat[i].ref_pressure, owl_dat->gas_temperature);
				owl_dat->skua_dat[l].param_dat[i].activation_energy = owl_dat->param_dat[i].activation_energy;
				owl_dat->skua_dat[l].param_dat[i].ref_temperature = 0.0;
				owl_dat->skua_dat[l].param_dat[i].affinity = 0.0;
			}
		}
		owl_dat->finch_dat[i].uo = Cstd(owl_dat->magpie_dat.sys_dat.PT*owl_dat->y[i], owl_dat->magpie_dat.sys_dat.T);
		owl_dat->finch_dat[i].Ro = (*owl_dat->eval_retard) (i,-1,owl_dat);
		if (owl_dat->finch_dat[i].Ro < 0.0)
		{
			mError(simulation_fail);
			return -1;
		}
		owl_dat->finch_dat[i].Do = (*owl_dat->eval_diff) (i,-1,owl_dat);
		owl_dat->finch_dat[i].kfnp1 = (*owl_dat->eval_kf) (i,owl_dat);
	}
	
	return success;
}

//Evaluate all parameters for SCOPSOWL
int set_SCOPSOWL_params(const void *user_data)
{
	int success = 0;
	SCOPSOWL_DATA *owl_dat = (SCOPSOWL_DATA *) user_data;
	
	//Faster methods use explicit coupling for multiscale problems
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		owl_dat->finch_dat[i].CN = false;
		owl_dat->finch_dat[i].beta = 1.0;
		//Loop over all nodes
		for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
		{
			owl_dat->finch_dat[i].Dnp1.edit(l, 0, (*owl_dat->eval_diff) (i,l,owl_dat));
		}
	}
	
	return success;
}

//Evaluate any routines after successfully finding a solution
int SCOPSOWL_postprocesses(SCOPSOWL_DATA *owl_dat)
{
	int success = 0;
	
	//Faster methods use explicit coupling for multiscale problems
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		//Loop over all nodes
		for (int l=0; l<owl_dat->finch_dat[i].LN; l++)
		{
			owl_dat->finch_dat[i].Rnp1.edit(l, 0, (*owl_dat->eval_retard) (i,l,owl_dat));
			if (owl_dat->finch_dat[i].Rnp1(l,0) < 0.0)
			{
				mError(simulation_fail);
				return -1;
			}
		}
	}
	
	//Store current results in correct locations then print to a file
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		//Now form averages and totals
		owl_dat->param_dat[i].qIntegralAvg = owl_dat->param_dat[i].qAvg.sphericalAvg(owl_dat->finch_dat[i].L,owl_dat->finch_dat[i].dz,owl_dat->param_dat[i].qo,owl_dat->DirichletBC);
		owl_dat->param_dat[i].QstAvg = owl_dat->param_dat[i].Qst.sphericalAvg(owl_dat->finch_dat[i].L, owl_dat->finch_dat[i].dz,owl_dat->param_dat[i].Qsto,owl_dat->DirichletBC);
	}
	
	//Print results
	if (owl_dat->Print2File == true)
	{
		owl_dat->t_counter = owl_dat->t_counter + owl_dat->finch_dat[0].dt;
		if (owl_dat->t_counter >= owl_dat->t_print)
		{
			print2file_SCOPSOWL_result_new(owl_dat);
			owl_dat->t_counter = 0.0;
		}
	}
	
	return success;
}

//Reset the current time level to prepare for next step
int SCOPSOWL_reset(SCOPSOWL_DATA *owl_dat)
{
	int success = 0;
	
	//Reset all time dependent information
	for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
	{
		owl_dat->param_dat[i].qAvg_old = owl_dat->param_dat[i].qAvg;
		owl_dat->param_dat[i].Qst_old = owl_dat->param_dat[i].Qst;
		owl_dat->param_dat[i].qIntegralAvg_old = owl_dat->param_dat[i].qIntegralAvg;
		owl_dat->param_dat[i].QstAvg_old = owl_dat->param_dat[i].QstAvg;
		success = (*owl_dat->finch_dat[i].resettime) ((void *)&owl_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	
	//Update SKUA if available
	if (owl_dat->SurfDiff == true && owl_dat->Heterogeneous == true)
	{
		for (int l=0; l<owl_dat->finch_dat[0].LN; l++)
		{
			success = SKUA_reset(&owl_dat->skua_dat[l]);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
	}
	
	//Reset time counter
	if (owl_dat->t_old == 0.0)
		owl_dat->t_old = 0.0;
	else
		owl_dat->t_old = owl_dat->t;
	
	return success;
}

//Standard simulation execution for self running
int SCOPSOWL(SCOPSOWL_DATA *owl_dat)
{
	int success = 0;
	
	//Set Initial Conditions
	success = set_SCOPSOWL_ICs(owl_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Loop till simulation complete
	do
	{
		if (owl_dat->finch_dat[0].Update == true)
		{
			success = SCOPSOWL_reset(owl_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_SCOPSOWL_timestep(owl_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		std::cout << "Evaluating time: " << owl_dat->t << " hrs..." << std::endl;
		success = SCOPSOWL_Executioner(owl_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
				owl_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); owl_dat->finch_dat[0].Update = false; return -1;}
		owl_dat->total_steps++;
	} while (owl_dat->t < owl_dat->sim_time);
	
	return success;
}

//Run the large scale cycle test
int LARGE_CYCLE_TEST(SCOPSOWL_DATA *owl_dat)
{
	int success = 0;
	
	//Set specific values for this test (H2O on MS3A)
	/*
	 owl_dat->total_pressure = 101.35;//kPa
	 owl_dat->gas_temperature = 313.15;//K
	 owl_dat->gas_velocity = 0.36;//cm/s
	 owl_dat->sim_time = 530.0;//hr
	 owl_dat->t_print = 0.5;
	 owl_dat->y[0] = 0.779736;
	 owl_dat->y[1] = 0.20736;
	 owl_dat->y[2] = 0.00934;
	 owl_dat->y[3] = 0.000314;
	 owl_dat->y[4] = 0.0030325;
	 */
	
	// H2O and I2 on Ag0Z
	owl_dat->total_pressure = 101.35;//kPa
	owl_dat->gas_temperature = 423.15;//K
	owl_dat->gas_velocity = 1.833;//cm/s
	owl_dat->sim_time = 200.0;//hr
	owl_dat->t_print = 0.5;
	owl_dat->y[0] = 0.779876175;				//-
	owl_dat->y[1] = 0.209476;				//-
	owl_dat->y[2] = 0.00934;				//-
	//owl_dat->y[3] = 7.6425e-5; //I2
	owl_dat->y[3] = 0.0; //I2
	owl_dat->y[4] = 1.2314e-3; //H2O
	
	//Set Initial Conditions
	success = set_SCOPSOWL_ICs(owl_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Loop till simulation complete
	do
	{
		if (owl_dat->finch_dat[0].Update == true)
		{
			success = SCOPSOWL_reset(owl_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_SCOPSOWL_timestep(owl_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		std::cout << "Evaluating time: " << owl_dat->t << " hrs..." << std::endl;
		
		//I2 and H2O on Ag0Z
		if (owl_dat->t > 25.0)
		{
			owl_dat->y[3] = 7.6425e-5;
		}
		if (owl_dat->t > 125.0)
		{
			owl_dat->y[4] = 0.0;
		}
		
		//Cycling Info Here (H2O on MS3A)
		/*
		 if (owl_dat->t > 18.3)
		 {
			owl_dat->y[0] = 0.789994;					//-
			owl_dat->y[1] = 0.21;					//-
			owl_dat->y[4] = 0.000434;					//-
			
			double m = (0.000128 - 0.000434) / (18.8 - 18.3);
			owl_dat->y[4] = (m * (owl_dat->t - 18.3)) + 0.000434;
		 }
		 if (owl_dat->t > 18.8)
		 {
			owl_dat->y[0] = 0.789994;					//-
			owl_dat->y[1] = 0.21;					//-
			owl_dat->y[4] = 0.000128;					//-
			
			double m = (3.6E-5 - 0.000128) / (22.1 - 18.8);
			owl_dat->y[4] = (m * (owl_dat->t - 18.8)) + 0.000128;
		 }
		 if (owl_dat->t > 22.1)
		 {
			owl_dat->y[0] = 0.789994;					//-
			owl_dat->y[1] = 0.21;					//-
			owl_dat->y[4] = 3.6E-5;					//-
			
			double m = (2.9E-5 - 3.6E-5) / (25.3 - 22.1);
			owl_dat->y[4] = (m * (owl_dat->t - 22.1)) + 3.6E-5;
		 }
		 if (owl_dat->t > 25.3)
		 {
			owl_dat->y[0] = 0.789994;					//-
			owl_dat->y[1] = 0.21;					//-
			owl_dat->y[4] = 2.9E-5;					//-
			
			double m = (2.0E-5 - 2.9E-5) / (42.1 - 25.3);
			owl_dat->y[4] = (m * (owl_dat->t - 25.3)) + 2.9E-5;
		 }
		 if (owl_dat->t > 42.1)
		 {
			owl_dat->y[0] = 0.789994;					//-
			owl_dat->y[1] = 0.21;					//-
			owl_dat->y[4] = 2.0E-5;					//-
			
			double m = (1.0E-5 - 2.0E-5) / (97.1 - 42.1);
			owl_dat->y[4] = (m * (owl_dat->t - 42.1)) + 2.0E-5;
		 }
		 if (owl_dat->t > 97.1)
		 {
			owl_dat->y[4] = 1.0E-5;
			
			double m = (7.2E-6 - 1.0E-5) / (150.0 - 97.1);
			owl_dat->y[4] = (m * (owl_dat->t - 97.1)) + 1.0E-5;
		 }
		 if (owl_dat->t > 150.0)
		 {
			owl_dat->y[4] = 7.2E-6;
		 }
		 if (owl_dat->t > 353.4)
		 {
			owl_dat->y[4] = 7.2E-6;
			
			double m = (7.5299E-4 - 7.2E-6) / (353.8 - 353.4);
			owl_dat->y[4] = (m * (owl_dat->t - 353.4)) + 7.2E-6;
		 }
		 if (owl_dat->t > 353.8)
		 {
			owl_dat->y[4] = 7.5299E-4;
			
			double m = (1.85358E-3 - 7.5299E-4) / (353.9 - 353.8);
			owl_dat->y[4] = (m * (owl_dat->t - 353.8)) + 7.5299E-4;
		 }
		 if (owl_dat->t > 353.9)
		 {
			owl_dat->y[4] = 1.85358E-3;
			
			double m = (3.03265E-3 - 1.85358E-3) / (355.3 - 353.9);
			owl_dat->y[4] = (m * (owl_dat->t - 353.9)) + 1.85358E-3;
		 }
		 if (owl_dat->t > 355.3)
		 {
			owl_dat->y[4] = 3.03265E-3;
			
			double m = (2.98526E-3 - 3.03265E-3) / (359.1 - 355.3);
			owl_dat->y[4] = (m * (owl_dat->t - 355.3)) + 3.03265E-3;
		 }
		 if (owl_dat->t > 359.1)
		 {
			owl_dat->y[4] = 2.98526E-3;
			
			double m = (2.89248E-3 - 2.98526E-3) / (378.2 - 359.1);
			owl_dat->y[4] = (m * (owl_dat->t - 359.1)) + 2.98526E-3;
		 }
		 if (owl_dat->t > 378.2)
		 {
			owl_dat->y[4] = 2.89248E-3;
			
			double m = (3.09567E-4 - 2.89248E-3) / (378.4 - 378.2);
			owl_dat->y[4] = (m * (owl_dat->t - 378.2)) + 2.89248E-3;
		 }
		 if (owl_dat->t > 378.4)
		 {
			owl_dat->y[4] = 3.09567E-4;
			
			double m = (1.77397E-4 - 3.09567E-4) / (378.7 - 378.4);
			owl_dat->y[4] = (m * (owl_dat->t - 378.4)) + 3.09567E-4;
		 }
		 if (owl_dat->t > 378.7)
		 {
			owl_dat->y[4] = 1.77397E-4;
			
			double m = (5.13752E-5 - 1.77397E-4) / (380.2 - 378.7);
			owl_dat->y[4] = (m * (owl_dat->t - 378.7)) + 1.77397E-4;
		 }
		 if (owl_dat->t > 380.2)
		 {
			owl_dat->y[4] = 5.13752E-5;
			
			double m = (3.24215E-5 - 5.13752E-5) / (382.0 - 380.2);
			owl_dat->y[4] = (m * (owl_dat->t - 380.2)) + 5.13752E-5;
		 }
		 if (owl_dat->t > 382.0)
		 {
			owl_dat->y[4] = 3.24215E-5;
			
			double m = (2.57389E-5 - 3.24215E-5) / (384.9 - 382.0);
			owl_dat->y[4] = (m * (owl_dat->t - 382.0)) + 3.24215E-5;
		 }
		 if (owl_dat->t > 384.9)
		 {
			owl_dat->y[4] = 2.57389E-5;
			
			double m = (1.20E-5 - 2.57389E-5) / (450.0 - 384.9);
			owl_dat->y[4] = (m * (owl_dat->t - 384.9)) + 2.57389E-5;
		 }
		 if (owl_dat->t > 450.0)
		 {
			owl_dat->y[4] = 1.20E-5;
			
			double m = (9.8E-6 - 1.20E-5) / (478.2 - 450.0);
			owl_dat->y[4] = (m * (owl_dat->t - 450.0)) + 1.20E-5;
		 }
		 if (owl_dat->t > 478.2)
		 {
			owl_dat->y[4] = 9.8E-6;
		 }
		 if (owl_dat->t > 498.1)
		 {
			owl_dat->y[4] = 9.8E-6;
			
			double m = (9.934E-4 - 9.8E-6) / (498.2 - 498.1);
			owl_dat->y[4] = (m * (owl_dat->t - 498.1)) + 9.8E-6;
		 }
		 if (owl_dat->t > 498.2)
		 {
			owl_dat->y[4] = 9.934E-4;
			
			double m = (2.331E-3 - 9.934E-4) / (498.3 - 498.2);
			owl_dat->y[4] = (m * (owl_dat->t - 498.2)) + 9.934E-4;
		 }
		 if (owl_dat->t > 498.3)
		 {
			owl_dat->y[4] = 2.331E-3;
			
			double m = (2.8925E-3 - 2.331E-3) / (498.7 - 498.3);
			owl_dat->y[4] = (m * (owl_dat->t - 498.3)) + 2.331E-3;
		 }
		 if (owl_dat->t > 498.7)
		 {
			owl_dat->y[4] = 2.8925E-3;
			
			double m = (2.9385E-3 - 2.8925E-3) / (499.7 - 498.7);
			owl_dat->y[4] = (m * (owl_dat->t - 498.7)) + 2.8925E-3;
		 }
		 if (owl_dat->t > 499.7)
		 {
			owl_dat->y[4] = 2.9385E-3;
			
			double m = (2.8697E-3 - 2.9385E-3) / (503.0 - 499.7);
			owl_dat->y[4] = (m * (owl_dat->t - 499.7)) + 2.9385E-3;
		 }
		 if (owl_dat->t > 503.0)
		 {
			owl_dat->y[4] = 2.8697E-3;
			
			double m = (2.9385E-3 - 2.8697E-3) / (522.1 - 503.0);
			owl_dat->y[4] = (m * (owl_dat->t - 503.0)) + 2.8697E-3;
		 }
		 if (owl_dat->t > 522.1)
		 {
			owl_dat->y[4] = 2.9385E-3;
			
			double m = (3.551E-4 - 2.9385E-3) / (522.3 - 522.1);
			owl_dat->y[4] = (m * (owl_dat->t - 522.1)) + 2.9385E-3;
		 }
		 if (owl_dat->t > 522.3)
		 {
			owl_dat->y[4] = 3.551E-4;
			
			double m = (2.1572E-4 - 3.551E-4) / (522.5 - 522.3);
			owl_dat->y[4] = (m * (owl_dat->t - 522.3)) + 3.551E-4;
		 }
		 if (owl_dat->t > 522.5)
		 {
			owl_dat->y[4] = 2.1572E-4;
			
			double m = (3.321E-5 - 2.1572E-4) / (525.8 - 522.5);
			owl_dat->y[4] = (m * (owl_dat->t - 522.5)) + 2.1572E-4;
		 }
		 if (owl_dat->t > 525.8)
		 {
			owl_dat->y[4] = 3.321E-5;
			
			double m = (2.943E-5 - 3.321E-5) / (527.0 - 525.8);
			owl_dat->y[4] = (m * (owl_dat->t - 525.8)) + 3.321E-5;
		 }
		 */
		
		success = SCOPSOWL_Executioner(owl_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<owl_dat->magpie_dat.sys_dat.N; i++)
				owl_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); owl_dat->finch_dat[0].Update = false; return -1;}
		owl_dat->total_steps++;
	} while (owl_dat->t < owl_dat->sim_time);
	
	return success;
}

//Run a particular scenario based on input files
int SCOPSOWL_SCENARIOS(const char *scene, const char *sorbent, const char *comp, const char *sorbate)
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
	SCOPSOWL_DATA dat;
	MIXED_GAS mixture;
	double time;
	FILE *Output;
	int i_read;
	double d_read;
	std::string s_read;
	
	//Initializations
	time = clock();
	Output = fopen("output/SCOPSOWL_Output.txt","w+");
	if (Output == nullptr)
	{
		system("mkdir output");
		Output = fopen("output/SCOPSOWL_Output.txt","w+");
	}
	dat.Print2File = true;
	dat.NonLinear = true;						//-
	dat.level = 1;
	dat.Print2Console = true;					//-
	dat.magpie_dat.sys_dat.Output = false;
	mixture.CheckMolefractions = true;
	dat.total_steps = 0;
	
	//	(1) - Read the Scenario File
	sceneFile >> d_read; dat.magpie_dat.sys_dat.T = d_read;
	sceneFile >> d_read; dat.magpie_dat.sys_dat.PT = d_read;
	dat.total_pressure = dat.magpie_dat.sys_dat.PT;
	dat.gas_temperature = dat.magpie_dat.sys_dat.T;
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
	
	for (int i=0; i<dat.magpie_dat.sys_dat.N; i++)
	{
		sceneFile >> s_read; dat.param_dat[i].speciesName = s_read;
		sceneFile >> i_read;
		if (i_read == 0) dat.param_dat[i].Adsorbable = false;
		else if (i_read == 1) dat.param_dat[i].Adsorbable = true;
		else {mError(invalid_boolean); return -1;}
		sceneFile >> d_read; dat.y[i] = d_read;
		sceneFile >> d_read; dat.param_dat[i].xIC = d_read;
		dat.param_dat[i].qIntegralAvg_old = dat.param_dat[i].xIC * dat.magpie_dat.sys_dat.qT;
		
		//Additional magpie initializations
		dat.magpie_dat.mspd_dat[i].eta.resize(dat.magpie_dat.sys_dat.N);
		dat.magpie_dat.gpast_dat[i].gama_inf.resize(dat.magpie_dat.sys_dat.N);
		dat.magpie_dat.gpast_dat[i].po.resize(dat.magpie_dat.sys_dat.N);
	}
	sceneFile.close();
	
	//Initialize gas mixture data
	success = initialize_data(dat.magpie_dat.sys_dat.N, &mixture);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	// Read (2) sorbentFile
	sorbentFile >> i_read;
	if (i_read == 0) dat.Heterogeneous = false;
	else if (i_read == 1) dat.Heterogeneous = true;
	else {mError(invalid_boolean); return -1;}
	
	sorbentFile >> i_read;
	if (i_read == 0) dat.SurfDiff = false;
	else if (i_read == 1) dat.SurfDiff = true;
	else {mError(invalid_boolean); return -1;}
	
	sorbentFile >> i_read; dat.coord_macro = i_read;
	if (i_read > 2 || i_read < 0) {mError(invalid_boolean); return -1;}
	if (i_read == 0 || i_read == 1)
	{
		sorbentFile >> d_read; dat.char_macro = d_read;
	}
	else
		dat.char_macro = 1.0;
	sorbentFile	>> d_read; dat.pellet_radius = d_read;
	sorbentFile >> d_read; dat.pellet_density = d_read;
	sorbentFile >> d_read; dat.binder_porosity = d_read;
	sorbentFile >> d_read; dat.binder_poresize = d_read;
	if (dat.Heterogeneous == true)
	{
		sorbentFile >> i_read; dat.coord_micro = i_read;
		if (i_read > 2 || i_read < 0) {mError(invalid_boolean); return -1;}
		if (i_read == 0 || i_read == 1)
		{
			sorbentFile >> d_read; dat.char_micro = d_read;
		}
		else
			dat.char_micro = 1.0;
		sorbentFile >> d_read; dat.crystal_radius = d_read;
		sorbentFile >> d_read; dat.binder_fraction = d_read;
	}
	else
	{
		dat.binder_fraction = 1.0;
		dat.crystal_radius = 1.0;
		dat.coord_micro = 2;
		dat.char_micro = 1.0;
	}
	sorbentFile.close();
	
	// Read (3) compFile
	for (int i=0; i<dat.magpie_dat.sys_dat.N; i++)
	{
		compFile >> d_read; mixture.species_dat[i].molecular_weight = d_read;
		compFile >> d_read; mixture.species_dat[i].specific_heat = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Viscosity = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Temp = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Const = d_read;
	}
	compFile.close();
	
	// Read (4) sorbateFile
	if (dat.Heterogeneous == true)
	{
		sorbateFile >> i_read;
		if (i_read == 0) dat.eval_surfDiff = (*default_Dc);
		else if (i_read == 1) dat.eval_surfDiff = (*simple_darken_Dc);
		else if (i_read == 2) dat.eval_surfDiff = (*theoretical_darken_Dc);
		else {mError(invalid_boolean); return -1;}
	}
	else
	{
		dat.eval_surfDiff = (*default_surf_diffusion);
	}
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
			dat.param_dat[i].ref_diffusion = 0.0;				//um^2/hr
			dat.param_dat[i].activation_energy = 0.0;			//J/mol
			dat.param_dat[i].ref_temperature = 0.0;				//K
			dat.param_dat[i].affinity = 0.0;					//-
			dat.magpie_dat.mspd_dat[i].v = 0.0;				//cm^3/mol
			dat.magpie_dat.gsta_dat[i].qmax = 0.0;			//mol/kg
			dat.magpie_dat.gsta_dat[i].m = 1;					//-
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
	
	//Call the setup function
	success = setup_SCOPSOWL_DATA(Output, default_adsorption, default_retardation, default_pore_diffusion, default_filmMassTransfer, dat.eval_surfDiff, (void *)&dat, &mixture, &dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Call Routine
	success = SCOPSOWL(&dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//END execution
	fclose(Output);
	time = clock() - time;
	std::cout << "Simulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Evaluations: " << dat.total_steps << "\n";
	std::cout << "Evaluations/sec: " << dat.total_steps/(time / CLOCKS_PER_SEC) << "\n";
	
	return success;
}

//Running tests for SCOPSOWL
int SCOPSOWL_TESTS()
{
	int success = 0;
	
	//Declarations
	SCOPSOWL_DATA dat;
	MIXED_GAS mixture;
	double time;
	FILE *TestOutput;
	
	//Initializations
	time = clock();
	TestOutput = fopen("output/SCOPSOWL_Test_Output.txt","w+");
	if (TestOutput == nullptr)
	{
		system("mkdir output");
		TestOutput = fopen("output/SCOPSOWL_Test_Output.txt","w+");
	}
	dat.Print2File = true;
	dat.NonLinear = true;						//-
	dat.level = 1;
	dat.Print2Console = true;					//-
	dat.magpie_dat.sys_dat.Output = false;
	mixture.CheckMolefractions = false;
	dat.total_steps = 0;
	
	//	(1) - Scenario to Test (with initializations)
	dat.gas_temperature = 423.15;				//K
	dat.total_pressure = 101.35;				//kPa
	dat.magpie_dat.sys_dat.T = dat.gas_temperature; 			//K
	dat.magpie_dat.sys_dat.PT = dat.total_pressure;			//kPa
	//dat.gas_velocity = 0.36;
	dat.gas_velocity = 1.833;					//cm/s
	dat.sim_time = 1.0;						//hrs
	dat.t_print = 0.1;
	dat.DirichletBC = false;					//-
	dat.SurfDiff = false;						//-
	dat.Heterogeneous = false;
	dat.t = 0.0;								//hrs
	dat.t_old = 0.0;							//hrs
	dat.magpie_dat.sys_dat.N = 5;				//-
	dat.magpie_dat.sys_dat.qT = 0.0;			//mol/kg
	dat.magpie_dat.gpast_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.magpie_dat.gsta_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.magpie_dat.mspd_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.param_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.finch_dat.resize(dat.magpie_dat.sys_dat.N);
	dat.y.resize(dat.magpie_dat.sys_dat.N);
	
	dat.param_dat[0].speciesName = "N2";
	dat.param_dat[1].speciesName = "O2";
	dat.param_dat[2].speciesName = "Ar";
	dat.param_dat[3].speciesName = "I2";
	dat.param_dat[4].speciesName = "H2O";
	
	dat.param_dat[0].Adsorbable = false;		//-
	dat.param_dat[1].Adsorbable = false;		//-
	dat.param_dat[2].Adsorbable = false;		//-
	dat.param_dat[3].Adsorbable = true;		//-
	dat.param_dat[4].Adsorbable = true;			//-
	
	dat.y[0] = 0.779876175;				//-
	dat.y[1] = 0.209476;				//-
	dat.y[2] = 0.00934;				//-
	dat.y[3] = 7.6425e-5;
	dat.y[4] = 1.2314e-3;
	
	dat.param_dat[0].xIC = 0.0;
	dat.param_dat[1].xIC = 0.0;
	dat.param_dat[2].xIC = 0.0;
	dat.param_dat[3].xIC = 0.0;
	dat.param_dat[4].xIC = 0.0;
	
	//Initialize gas mixture data
	success = initialize_data(dat.magpie_dat.sys_dat.N, &mixture);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//NOTE: Should Check Molefractions here ------------------------------
	
	//	(2) - Sorbent to test (with more initializations)
	dat.coord_macro = 1;
	dat.char_macro = 0.4;									//cm
	dat.pellet_radius = 0.08;								//cm
	dat.pellet_density = 3.06;								//kg/L
	dat.binder_porosity = 0.384;							//-
	dat.binder_poresize = 1.5E-4;							//cm
	dat.coord_micro = 2;
	dat.char_micro = 1.0;
	dat.crystal_radius = 2.0;								//um
	dat.binder_fraction = 0.175;
	
	/*
	 dat.coord_macro = 1;
	 dat.char_macro = 0.335;
	 dat.pellet_radius = 0.08;								//cm
	 dat.pellet_density = 0.93;								//kg/L
	 dat.binder_porosity = 0.272;							//-
	 dat.binder_poresize = 3.5E-6;							//cm
	 dat.coord_micro = 2;
	 dat.char_micro = 1.0;
	 dat.crystal_radius = 8.0;								//um
	 dat.binder_fraction = 0.175;
	 */
	
	if (dat.Heterogeneous == false)
	{
		dat.binder_fraction = 1.0;
	}
	
	//	(3) - Components to test
	
	//Set the Constants in mixture data (Had to be initialized first)
	mixture.species_dat[0].molecular_weight = 28.016;
	mixture.species_dat[1].molecular_weight = 32.0;
	mixture.species_dat[2].molecular_weight = 39.948;
	mixture.species_dat[3].molecular_weight = 253.808;
	mixture.species_dat[4].molecular_weight = 18.0;
	
	mixture.species_dat[0].specific_heat = 1.04;
	mixture.species_dat[1].specific_heat = 0.919;
	mixture.species_dat[2].specific_heat = 0.522;
	mixture.species_dat[3].specific_heat = 0.214;
	mixture.species_dat[4].specific_heat = 1.97;
	
	mixture.species_dat[0].Sutherland_Viscosity = 0.0001781;
	mixture.species_dat[1].Sutherland_Viscosity = 0.0002018;
	mixture.species_dat[2].Sutherland_Viscosity = 0.0002125;
	mixture.species_dat[3].Sutherland_Viscosity = 0.00013283;
	mixture.species_dat[4].Sutherland_Viscosity = 0.0001043;
	
	mixture.species_dat[0].Sutherland_Temp = 300.55;
	mixture.species_dat[1].Sutherland_Temp = 292.25;
	mixture.species_dat[2].Sutherland_Temp = 273.11;
	mixture.species_dat[3].Sutherland_Temp = 295.496;
	mixture.species_dat[4].Sutherland_Temp = 298.16;
	
	mixture.species_dat[0].Sutherland_Const = 111.0;
	mixture.species_dat[1].Sutherland_Const = 127.0;
	mixture.species_dat[2].Sutherland_Const = 144.4;
	mixture.species_dat[3].Sutherland_Const = 573.474;
	mixture.species_dat[4].Sutherland_Const = 784.72;
	
	//	(4) - Adsorbate to test
	for (int i=0; i<dat.magpie_dat.sys_dat.N; i++)
	{
		dat.magpie_dat.mspd_dat[i].eta.resize(dat.magpie_dat.sys_dat.N);
		dat.magpie_dat.gpast_dat[i].gama_inf.resize(dat.magpie_dat.sys_dat.N);
		dat.magpie_dat.gpast_dat[i].po.resize(dat.magpie_dat.sys_dat.N);
		
		if (dat.param_dat[i].Adsorbable == false)
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
		else
		{
			if (i==3) //I2
			{
				//Const Dc Parameters
				dat.param_dat[i].ref_diffusion = 0.0;		//um^2/hr
				dat.param_dat[i].activation_energy = 0.0;	//J/mol
				dat.param_dat[i].ref_temperature = 0.0;			//K
				dat.param_dat[i].affinity = 0.0;				//-
				
				dat.magpie_dat.mspd_dat[i].v = 18.8;			//cm^3/mol
				dat.magpie_dat.gsta_dat[i].qmax = 1.062207984;		//mol/kg
				dat.magpie_dat.gsta_dat[i].m = 1;				//-
				dat.magpie_dat.gsta_dat[i].dHo.resize(dat.magpie_dat.gsta_dat[i].m);
				dat.magpie_dat.gsta_dat[i].dSo.resize(dat.magpie_dat.gsta_dat[i].m);
			 
				dat.magpie_dat.gsta_dat[i].dHo[0] = -8308.67985;	//J/mol
				dat.magpie_dat.gsta_dat[i].dSo[0] = 92.3491515;	//J/K/mol
				
			}
			
			if (i==4) //H2O
			{
				
				//Const Dc Parameters
				dat.param_dat[i].ref_diffusion = 0.0;		//um^2/hr
				dat.param_dat[i].activation_energy = 0.0;	//J/mol
				dat.param_dat[i].ref_temperature = 0.0;			//K
				dat.param_dat[i].affinity = 0.0;				//-
				
				dat.magpie_dat.mspd_dat[i].v = 13.91;			//cm^3/mol
				dat.magpie_dat.gsta_dat[i].qmax = 4.71;		//mol/kg
				dat.magpie_dat.gsta_dat[i].m = 3;				//-
				dat.magpie_dat.gsta_dat[i].dHo.resize(dat.magpie_dat.gsta_dat[i].m);
				dat.magpie_dat.gsta_dat[i].dSo.resize(dat.magpie_dat.gsta_dat[i].m);
			 
				dat.magpie_dat.gsta_dat[i].dHo[0] = -43260.13434;	//J/mol
				dat.magpie_dat.gsta_dat[i].dSo[0] = -38.88589704;	//J/K/mol
				
				dat.magpie_dat.gsta_dat[i].dHo[1] = -110183.2213;	//J/mol
				dat.magpie_dat.gsta_dat[i].dSo[1] = -179.30964;	//J/K/mol
				
				dat.magpie_dat.gsta_dat[i].dHo[2] = -125398.6827;	//J/mol
				dat.magpie_dat.gsta_dat[i].dSo[2] = -180.2574885;	//J/K/mol
			}
		}
		
	}
	
	//Setup the problem
	if (dat.Heterogeneous == true)
	{
		success = setup_SCOPSOWL_DATA(TestOutput, default_adsorption, default_retardation, default_pore_diffusion, default_filmMassTransfer, default_Dc, (void *)&dat, &mixture, &dat);
	}
	//Setups for Homo Pellet
	else if (dat.Heterogeneous == false && dat.SurfDiff == true)
	{
		success = setup_SCOPSOWL_DATA(TestOutput, default_adsorption, default_retardation, default_pore_diffusion, default_filmMassTransfer, default_surf_diffusion, (void *)&dat, &mixture, &dat);
	}
	else if (dat.Heterogeneous == false && dat.SurfDiff == false)
	{
		success = setup_SCOPSOWL_DATA(TestOutput, default_adsorption, default_retardation, default_pore_diffusion, default_filmMassTransfer, zero_surf_diffusion, (void *)&dat, &mixture, &dat);
	}
	else
	{
		//Nothing
	}
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Call Routine
	success = SCOPSOWL(&dat);
	//success = LARGE_CYCLE_TEST(&dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//END execution
	fclose(TestOutput);
	time = clock() - time;
	std::cout << "Simulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Evaluations: " << dat.total_steps << "\n";
	std::cout << "Evaluations/sec: " << dat.total_steps/(time / CLOCKS_PER_SEC) << "\n";
	
	return success;
}
