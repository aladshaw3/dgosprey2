/*!
 *  \file finch.h finch.cpp
 *	\brief Flux-limiting Implicit Non-oscillatory Conservative High-resolution scheme
 *  \author Austin Ladshaw
 *	\date 01/29/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "finch.h"

//Find the maximum of a vector of data
double max(std::vector<double> &values)
{
	double max = values[0];
	unsigned long int size = values.size();
	for (int i=0; i<size; i++)
	{
		if (values[i] >= max)
			max = values[i];
	}
	return max;
}

//Find the minium of a vector of data
double min(std::vector<double> &values)
{
	double min = values[0];
	unsigned long int size = values.size();
	for (int i=0; i<size; i++)
	{
		if (values[i] <= min)
			min = values[i];
	}
	return min;
}

//Evaluate the minmod function for a vector of values
double minmod(std::vector<double> &values)
{
	double minmod = 0.0;
	unsigned long int size = values.size();
	bool pos = true;
	bool change = false;
	for (int i=0; i<size; i++)
	{
		if (change == false)
		{
			if (values[i] > 0.0 && pos == true) {pos = true;}
			else if (values[i] < 0.0) {pos = false;}
			else change = true;
		}
		else
			break;
	}
	if (pos == true && change == false) minmod = min(values);
	else if (pos == false && change == false) minmod = max(values);
	else minmod = 0.0;
	
	return minmod;
}

//Function to return the Total value of the conserved quantity in the domain
int uTotal(FINCH_DATA *dat)
{
	int success = 0;
	double total = 0.0;
	
	for (int l=0; l<dat->LN; l++)
	{
		if (l == 0)
		{
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				total = total + (0.5 * (dat->unp1(l,0) +dat->uo) * (pow((double)(l+1),(double)dat->d+1) - pow((double)l,(double)dat->d+1)));
				total = total + (0.5 * (dat->unp1(l+1,0) +dat->unp1(l,0)) * (pow((double)(l+1),(double)dat->d+1) - pow((double)l,(double)dat->d+1)));
			}
			else
			{
				total = total + (0.5 * (dat->unp1(l+1,0) + dat->unp1(l,0)) * (pow((double)(l+1),(double)dat->d+1) - pow((double)l,(double)dat->d+1)));
			}
		}
		else if (l == dat->LN-1)
		{
			if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				total = total + (0.5 * (dat->uo + dat->unp1(l,0)) * (pow((double)(l+1),(double)dat->d+1) - pow((double)l,(double)dat->d+1)));
			}
			else {/*No Action*/}
		}
		else
		{
			total = total + (0.5 * (dat->unp1(l+1,0) + dat->unp1(l,0)) * (pow((double)(l+1),(double)dat->d+1) - pow((double)l,(double)dat->d+1)));
		}
	}
	dat->uT = pow(dat->dz,(double)dat->d+1)*total;
	
	if (dat->d == 0)
		dat->uT = dat->uT * dat->s;
	else if (dat->d == 1)
		dat->uT = M_PI * dat->uT * dat->s;
	else if (dat->d == 2)
		dat->uT = (4.0/3.0) * M_PI * dat->uT;
	else {mError(invalid_boolean); return -1;}
	
	return success;
}

//Function to return the Average value of the conserved quantity in the domain
int uAverage(FINCH_DATA *dat)
{
	int success = 0;
	
	if (dat->d == 0)
	{
		dat->uAvg = dat->uT / (dat->L * dat->s);
	}
	else if (dat->d == 1)
	{
		dat->uAvg = dat->uT / (M_PI * pow(dat->L,2.0) * dat->s);
	}
	else if (dat->d == 2)
	{
		dat->uAvg = dat->uT / ( (4.0/3.0) * M_PI * pow(dat->L,3.0) );
	}
	else {mError(invalid_boolean); return -1;}
	
	return success;
}

//Check solution vectors for negative mass
int check_Mass(FINCH_DATA *dat)
{
	int success = 0;
	
	//Small negative numbers will be treated as zero or given previously valid values
	if (dat->CheckMass == false)
		return 0;
	
	//Check for negative concentrations first
	for (int l=0; l<dat->LN; l++)
	{
		//Check for negative concentrations in u
		if (dat->u_star(l,0) < -1.0e-4 || dat->un(l,0) < -1.0e-4 || dat->unm1(l,0) < -1.0e-4 || dat->unp1(l,0) < -1.0e-4)
		{
			//Error Reporting is temporarilly disabled disabled
			/*
			 mError(negative_mass);
			 std::cout << "u_star at node " << l << ": " << dat->u_star(l,0) << std::endl;
			 std::cout << "unp1 at node " << l << ": " << dat->unp1(l,0) << std::endl;
			 std::cout << "un at node " << l << ": " << dat->u_star(l,0) << std::endl;
			 std::cout << "unm1 at node " << l << ": " << dat->u_star(l,0) << std::endl;
			 return -1;
			 */
			
			dat->unp1.edit(l,0,dat->un(l,0));
			dat->u_star.edit(l,0,dat->un(l,0));
			dat->unm1.edit(l,0,0.0);
			dat->un.edit(l,0,0.0);
			
		}
		if (dat->u_star(l,0) >= -1.0e-4 && dat->u_star(l,0) < 0.0)
		{
			dat->u_star.edit(l,0,dat->un(l,0));
		}
		if (dat->unp1(l,0) >= -1.0e-4 && dat->unp1(l,0) < 0.0)
		{
			dat->unp1.edit(l,0,dat->un(l,0));
		}
		if (dat->un(l,0) >= -1.0e-4 && dat->un(l,0) < 0.0)
		{
			dat->un.edit(l,0,0.0);
		}
		if (dat->unm1(l,0) >= -1.0e-4 && dat->unm1(l,0) < 0.0)
		{
			dat->unm1.edit(l,0,0.0);
		}
		else {/*No Action*/}
	}
	
	return success;
}

//Function to evaluate the RHS of the system (Default Supplied Solver called by default_solve)
int l_direct(FINCH_DATA *dat)
{
	int success = 0;
	
	//Set Up all parameters
	success = (*dat->setparams) (dat->param_data);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Set Up all Coefficients and slope limiters
	success = (*dat->discretize) (dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Set up the BCs for the system
	success = (*dat->setbcs) (dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	for (int l=0; l<dat->LN; l++)
	{
		if (l==0)
		{
			dat->res.edit(l, 0, dat->NE[l]*dat->un(l,0) + dat->OE[l]*dat->un(l+1,0) + dat->gE(l,0) - dat->gI(l,0) + dat->Fn(l,0) - dat->Fnp1(l,0) + dat->Sn(l,0) - dat->Snp1(l,0));
		}
		else if (l==dat->LN-1)
		{
			dat->res.edit(l, 0, dat->ME[l]*dat->un(l-1,0) + dat->NE[l]*dat->un(l,0) + dat->gE(l,0) - dat->gI(l,0) + dat->Fn(l,0) - dat->Fnp1(l,0) + dat->Sn(l,0) - dat->Snp1(l,0));
		}
		else
		{
			dat->res.edit(l, 0, dat->ME[l]*dat->un(l-1,0) + dat->NE[l]*dat->un(l,0) + dat->OE[l]*dat->un(l+1,0) + dat->gE(l,0) - dat->gI(l,0) + dat->Fn(l,0) - dat->Fnp1(l,0) + dat->Sn(l,0) - dat->Snp1(l,0));
		}
	}
	
	//Solve the system in one step
	success = (*dat->evalprecon) (dat->res,dat->unp1,dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	return success;
}

//Function to evaluate the Picard Step x = G(x) (Used by LARK for Method Testing)
int lark_picard_step(const Matrix<double> &x, Matrix<double> &G, const void *data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) data;
	
	//Call Residual Function
	success = (*dat->evalres) (x,dat->res,dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Form Search direction
	success = (*dat->evalprecon) (dat->res,dat->pres,dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Form guess
	for (int n=0; n<G.rows(); n++)
		G.edit(n, 0, G(n,0) + dat->pres(n,0));
	
	return success;
}

//Function to solve non-linear portion iteratively (Default Supplied Solver called by default_solve)
int nl_picard(FINCH_DATA *dat)
{
	int success = 0;
	double rel_res_base, res_norm, res_old, bestres;
	
	//Set up System for First Iterate
	dat->unp1 = dat->un;
	
	//Call Residual Function
	success = (*dat->evalres) (dat->unp1,dat->res,dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	res_norm = dat->res.norm();
	bestres = res_norm;
	dat->ubest = dat->unp1;
	res_old = res_norm;
	rel_res_base = res_norm;
	//Report Norms if applicable
	if (dat->NormTrack == true)
	{
		std::cout << "Begin FINCH Picard/Inexact-Newton iterations...\n" << std::endl;
		std::cout << "\tIterate:\tRes:\tRelRes:" << std::endl;
		std::cout << "\t" << 0 << "\t" << res_norm << "\t" << (res_norm/rel_res_base) << std::endl;
	}
	
	//Loop until convergence
	int k;
	for (k=0; k<dat->max_iter; k++)
	{
		dat->total_iter++;
		
		//Form Search direction
		success = (*dat->evalprecon) (dat->res,dat->pres,dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		for (int n=0; n<dat->unp1.rows(); n++)
			dat->unp1.edit(n, 0, dat->unp1(n,0) + dat->pres(n,0));
		
		//Call the Residual Function
		success = (*dat->evalres) (dat->unp1,dat->res,dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		res_norm = dat->res.norm();
		if (dat->NormTrack == true)
			std::cout << "\t" << k+1 << "\t" << res_norm << "\t" << (res_norm/rel_res_base) << std::endl;
		if (res_norm <= dat->tol_abs)
		{
			if (dat->NormTrack == true)
				std::cout << "\nSolution Converged in " << k+1 << " iteration(s) within Residual Tolerance!" << std::endl;
			return success;
		}
		else if ( (res_norm/rel_res_base) <= dat->tol_rel)
		{
			if (dat->NormTrack == true)
				std::cout << "\nSolution Converged in " << k+1 << " iteration(s) within Relative Residual Tolerance!" << std::endl;
			return success;
		}
		else
		{
			if (res_norm < bestres)
			{
				dat->ubest = dat->unp1;
				bestres = res_norm;
			}
		}
		
		res_old = res_norm;
	}
	
	if (k >= dat->max_iter)
	{
		if (dat->NormTrack == true)
		{
			std::cout << "\nReached Maximum Iterations without Convergence!" << std::endl;
			std::cout << "Best Reported Norm: " << bestres << std::endl;
		}
		dat->unp1 = dat->ubest;
	}
	if (success == -1)
	{
		if (dat->NormTrack == true)
		{
			std::cout << "Picard did not find exact solution within " << k << " iterations!" << std::endl;
			std::cout << "Best reported norm =\t" << bestres << std::endl;
		}
		dat->unp1 = dat->ubest;
		success = 0;
	}
	
	return success;
}

//Function to setup memory and default arguments based on input
int setup_FINCH_DATA( int (*user_callroutine) (const void *user_data),
					 int (*user_setic) (const void *user_data),
					 int (*user_timestep) (const void *user_data),
					 int (*user_preprocess) (const void *user_data),
					 int (*user_solve) (const void *user_data),
					 int (*user_setparams) (const void *user_data),
					 int (*user_discretize) (const void *user_data),
					 int (*user_bcs) (const void *user_data),
					 int (*user_res) (const Matrix<double>&x, Matrix<double>& res, const void *user_data),
					 int (*user_precon) (const Matrix<double>&b, Matrix<double>& p, const void *user_data),
					 int (*user_postprocess) (const void *user_data),
					 int (*user_reset) (const void *user_data),
					 FINCH_DATA *dat, const void *param_data)
{
	int success = 0;
	
	//Set up function pointers
	if ((*user_callroutine) == NULL)
		dat->callroutine = (*default_execution);
	else
		dat->callroutine = (*user_callroutine);
	if ((*user_setic) == NULL)
		dat->setic = (*default_ic);
	else
		dat->setic = (*user_setic);
	if ((*user_timestep) == NULL)
		dat->settime = (*default_timestep);
	else
		dat->settime = (*user_timestep);
	if ((*user_preprocess) == NULL)
		dat->setpreprocess = (*default_preprocess);
	else
		dat->setpreprocess = (*user_preprocess);
	if ((*user_solve) == NULL)
		dat->solve = (*default_solve);
	else
		dat->solve = (*user_solve);
	if ((*user_setparams) == NULL)
		dat->setparams = (*default_params);
	else
		dat->setparams = (*user_setparams);
	if ((*user_discretize) == NULL)
		dat->discretize = (*ospre_discretization);
	else
		dat->discretize = (*user_discretize);
	if ((*user_bcs) == NULL)
		dat->setbcs = (*default_bcs);
	else
		dat->setbcs = (*user_bcs);
	if ((*user_res) == NULL)
		dat->evalres = (*default_res);
	else
		dat->evalres = (*user_res);
	if ((*user_precon) == NULL)
		dat->evalprecon = (*default_precon);
	else
		dat->evalprecon = (*user_precon);
	if ((*user_postprocess) == NULL)
		dat->setpostprocess = (*default_postprocess);
	else
		dat->setpostprocess = (*user_postprocess);
	if ((*user_reset) == NULL)
		dat->resettime = (*default_reset);
	else
		dat->resettime = (*user_reset);
	dat->param_data = param_data;
	
	//Set default arguments
	if (dat->SteadyState == true)
		dat->CN = false;
	else
		dat->CN = true;
	if (dat->CN == true)
		dat->beta = 0.5;
	else
		dat->beta = 1.0;
	dat->dz = dat->L/dat->LN;
	if (dat->Dirichlet == false)
		dat->LN = dat->LN+1;
	if (dat->Iterative == true)
		dat->ExplicitFlux = false;
	
	//Set up method
	if (dat->SteadyState == true)
		dat->max_iter = std::min(30*dat->LN,1000);
	else
		dat->max_iter = std::min(3*dat->LN,1000);
	dat->picard_dat.maxit = dat->max_iter;
	dat->pjfnk_dat.nl_maxit = dat->max_iter;
	dat->total_iter = 0;
	dat->pjfnk_dat.NL_Output = dat->NormTrack;
	dat->picard_dat.Output = dat->NormTrack;
	
	//Set up Memory
	//Vectors
	dat->CL_I.resize(dat->LN);
	dat->CL_E.resize(dat->LN);
	dat->CC_I.resize(dat->LN);
	dat->CC_E.resize(dat->LN);
	dat->CR_I.resize(dat->LN);
	dat->CR_E.resize(dat->LN);
	
	dat->fL_I.resize(dat->LN);
	dat->fL_E.resize(dat->LN);
	dat->fC_I.resize(dat->LN);
	dat->fC_E.resize(dat->LN);
	dat->fR_I.resize(dat->LN);
	dat->fR_E.resize(dat->LN);
	
	dat->OI.resize(dat->LN);
	dat->OE.resize(dat->LN);
	dat->NI.resize(dat->LN);
	dat->NE.resize(dat->LN);
	dat->MI.resize(dat->LN);
	dat->ME.resize(dat->LN);
	
	dat->uz_l_E.resize(3);
	dat->uz_lm1_E.resize(3);
	dat->uz_lp1_E.resize(3);
	dat->uz_l_I.resize(3);
	dat->uz_lm1_I.resize(3);
	dat->uz_lp1_I.resize(3);
	
	//Matrices
	dat->unm1.set_size(dat->LN,1);
	dat->un.set_size(dat->LN,1);
	dat->unp1.set_size(dat->LN,1);
	dat->u_star.set_size(dat->LN,1);
	dat->ubest.set_size(dat->LN, 1);
	dat->vn.set_size(dat->LN,1);
	dat->vnp1.set_size(dat->LN,1);
	dat->Dn.set_size(dat->LN,1);
	dat->Dnp1.set_size(dat->LN,1);
	dat->kn.set_size(dat->LN, 1);
	dat->knp1.set_size(dat->LN, 1);
	dat->Sn.set_size(dat->LN, 1);
	dat->Snp1.set_size(dat->LN, 1);
	dat->Rn.set_size(dat->LN, 1);
	dat->Rnp1.set_size(dat->LN, 1);
	dat->Fn.set_size(dat->LN,1);
	dat->Fnp1.set_size(dat->LN,1);
	dat->gI.set_size(dat->LN,1);
	dat->gE.set_size(dat->LN,1);
	dat->res.set_size(dat->LN,1);
	dat->pres.set_size(dat->LN,1);
	
	return success;
}

//Function to print out the space header for the output file
void print2file_dim_header(FILE *Output, FINCH_DATA *dat)
{
	fprintf(Output,"z_dim\t%.6g",0.0);
	for (int l=0; l<dat->LN; l++)
	{
		if (dat->Dirichlet == false && l == dat->LN-1)
			break;
		fprintf(Output,"\t%.6g",(l+1)*dat->dz);
	}
}

//Function to print out the time header for the output file
void print2file_time_header(FILE *Output, FINCH_DATA *dat)
{
	fprintf(Output,"Time\t");
	fprintf(Output,"u[0]\t");
	for (int l=0; l<dat->LN; l++)
	{
		if (dat->Dirichlet == false && l == dat->LN-1)
			break;
		fprintf(Output,"u[%i]\t",l+1);
	}
	fprintf(Output,"uTotal\tuAverage");
}

//Function to print out the previous time step results to a file
void print2file_result_old(FILE *Output, FINCH_DATA *dat)
{
	
	if (dat->vo > 0.0)
	{
		fprintf(Output,"%.6g\t",dat->t_old);
		if (dat->Dirichlet == true)
			fprintf(Output,"%.6g\t",dat->uo);
		for (int l=0; l<dat->LN; l++)
		{
			fprintf(Output,"%.6g\t",dat->un(l,0));
		}
		fprintf(Output,"%.6g\t%.6g",dat->uT,dat->uAvg);
	}
	else
	{
		fprintf(Output,"%.6g\t",dat->t_old);
		for (int l=0; l<dat->LN; l++)
		{
			fprintf(Output,"%.6g\t",dat->un(l,0));
		}
		if (dat->Dirichlet == true)
			fprintf(Output,"%.6g\t",dat->uo);
		fprintf(Output,"%.6g\t%.6g",dat->uT,dat->uAvg);
	}
	
}

//Function to print out the new simulation results
void print2file_result_new(FILE *Output, FINCH_DATA *dat)
{
	
	if (dat->vo > 0.0)
	{
		fprintf(Output,"%.6g\t",dat->t);
		if (dat->Dirichlet == true)
			fprintf(Output,"%.6g\t",dat->uo);
		for (int l=0; l<dat->LN; l++)
		{
			fprintf(Output,"%.6g\t",dat->unp1(l,0));
		}
		fprintf(Output,"%.6g\t%.6g",dat->uT,dat->uAvg);
	}
	else
	{
		fprintf(Output,"%.6g\t",dat->t);
		for (int l=0; l<dat->LN; l++)
		{
			fprintf(Output,"%.6g\t",dat->unp1(l,0));
		}
		if (dat->Dirichlet == true)
			fprintf(Output,"%.6g\t",dat->uo);
		fprintf(Output,"%.6g\t%.6g",dat->uT,dat->uAvg);
	}
	
}

//Function to manually move down a line in the output file
void print2file_newline(FILE *Output, FINCH_DATA *dat)
{
	fprintf(Output, "\n");
}

//Function to manually move over one tab in output file
void print2file_tab(FILE *Output, FINCH_DATA *dat)
{
	fprintf(Output, "\t");
}

//Default execution function for FINCH
int default_execution(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	//Perform Preprocess Actions
	success = (*dat->setpreprocess) (dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Solve the system
	success = (*dat->solve) (dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Perform Postprocess Actions
	success = (*dat->setpostprocess) (dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Check for negative concentrations
	success = check_Mass(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Form Totals
	success = uTotal(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Form Averages
	success = uAverage(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	return success;
}

//Function for the default initial conditions
int default_ic(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	if (dat->SteadyState == true)
		dat->Rn.ConstantICFill(0.0);
	else
		dat->Rn.ConstantICFill(dat->RIC);
	
	dat->un.ConstantICFill(dat->uIC);
	dat->vn.ConstantICFill(dat->vIC);
	dat->Dn.ConstantICFill(dat->DIC);
	dat->kn.ConstantICFill(dat->kIC);
	dat->Sn.ConstantICFill(0.0);
	dat->unm1 = dat->un;
	dat->unp1 = dat->un;
	
	success = uTotal(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	success = uAverage(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	dat->uT_old = dat->uT;
	dat->uAvg_old = dat->uAvg;
	
	return success;
}

//Default time step
int default_timestep(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	double maxV = 0.0;
	
	maxV = fabs(dat->vo);
	for (int l=0; l<dat->LN; l++)
	{
		if (maxV < fabs(dat->vn(l,0)))
			maxV = fabs(dat->vn(l,0));
		if (maxV < fabs(dat->vnp1(l,0)))
			maxV = fabs(dat->vnp1(l,0));
	}
	if (fabs(maxV) <= 0.0 || dat->Iterative == true)
	{
		dat->dt = dat->dz / 2.0;
	}
	else
	{
		dat->dt = dat->dz / fabs(maxV) / 8.0;
	}
	
	return success;
}

//Default preprocess function
int default_preprocess(const void *user_data)
{
	std::cout << "Preprocess Actions..." << std::endl;
	return 0;
}

//Default function to solve the system either linearly or with picard iterations
int default_solve(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	if (dat->Iterative == false)
	{
		success = l_direct(dat);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	else
	{
		if (dat->nl_method == FINCH_Picard)
		{
			success = nl_picard(dat);
			if (success != 0)
			{
				mError(simulation_fail); return -1;
			}
		}
		else if (dat->nl_method == LARK_Picard)
		{
			success = picard(dat->evalres,lark_picard_step,dat->unp1,&dat->picard_dat,dat,dat);
			if (success !=0) {mError(simulation_fail); return -1;}
			dat->total_iter = dat->total_iter + dat->picard_dat.iter;
		}
		else if (dat->nl_method == LARK_PJFNK)
		{
			success = pjfnk(dat->evalres, dat->evalprecon, dat->unp1, &dat->pjfnk_dat, dat, dat);
			if (success !=0) {mError(simulation_fail); return -1;}
			dat->total_iter = dat->total_iter + dat->pjfnk_dat.nl_iter + dat->pjfnk_dat.l_iter;
		}
		else
		{
			std::cout << "\nInvalid Method Selection!\n" << std::endl;
			mError(simulation_fail);
			return -1;
		}
	}
	return success;
}

//Function for the default parameter values
int default_params(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	for (int l=0; l<dat->LN; l++)
	{
		dat->vnp1(l,0) = dat->vo;
		dat->Dnp1(l,0) = dat->Do;
		dat->knp1(l,0) = dat->ko;
		dat->Snp1(l,0) = 0.0;
		if (dat->SteadyState == false)
			dat->Rnp1(l,0) = dat->Ro;
		else
			dat->Rnp1(l,0) = 0.0;
	}
	
	return success;
}

//Default discretization of the FINCH
int minmod_discretization(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	double Dn_neg, Dn_pos, vn_neg, vn_pos;
	double Dnp1_neg, Dnp1_pos, vnp1_neg, vnp1_pos;
	double z_pos, z_neg, zl;
	double check = 0.0;
	double w = 0.0; //weight
	double un, u_star, uo;
	
	//Fill in the lambda parameters for Neumann BC
	if (dat->vo > 0.0 && dat->d == 0)
	{
		dat->lambda_I = (2.0 * dat->dz * dat->vnp1(0,0)) / dat->Dnp1(0,0);
		dat->lambda_E = (2.0 * dat->dz * dat->vn(0,0)) / dat->Dn(0,0);
	}
	else if (dat->vo <= 0.0 && dat->d == 0)
	{
		dat->lambda_I = (2.0 * dat->dz * dat->vnp1(dat->LN-1,0)) / dat->Dnp1(dat->LN-1,0);
		dat->lambda_E = (2.0 * dat->dz * dat->vn(dat->LN-1,0)) / dat->Dn(dat->LN-1,0);
	}
	else
	{
		dat->lambda_I = -(2.0 * dat->dz * dat->kfnp1) / dat->Dnp1(dat->LN-1,0);
		dat->lambda_E = -(2.0 * dat->dz * dat->kfn) / dat->Dn(dat->LN-1,0);
	}
	
	//Check the Neumann parameters for problems
	if (dat->Dirichlet == false)
	{
		if (isnan(dat->lambda_E) || isinf(dat->lambda_E) || isnan(dat->lambda_I) || isinf(dat->lambda_I))
		{
			std::cout << "\nBreakdown of the Neumann BCs!\n Restart with Dirichlet BCs!\n";
			mError(simulation_fail);
			return -1;
		}
	}
	
	//Set up the necessary information for the slope limiters
	if (dat->Iterative == false)
	{
		dat->u_star.columnProjection(dat->un, dat->unm1, dat->dt, dat->dt_old);
		uo = dat->uo;
		//Check for negative concentrations
		success = check_Mass(dat);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	else
	{
		dat->u_star = dat->unp1;
		uo = dat->uo;
	}
	
	//Determine the weight factor as some weighted average between 1 and 2
	if (dat->ExplicitFlux == true)
		w = w + 0.25;
	else
		w = w + 0.5;
	if (dat->Iterative == true)
		w = w + 0.25;
	else
		w = w + 0.5;
	if (dat->SteadyState == false)
		w = w + 0.25;
	else
		w = w + 0.5;
	if (dat->CN == false)
		w = w + 0.25;
	else
		w = w + 0.5;
	
	/*	NOTE:
	 If we are using the iterative approach, the implicit fluxes will be based on dat->unp1,
	 instead of u_star. Weight factor still depends on ExplicitFlux flag
	 
	 If we are using the Explicit only flux approach, then all fluxes are based on dat->un
	 and we use the smaller weight factor (1.0) for additional stability.
	 
	 As a consequence, if we are using the iterative approach AND the explicit flux approach,
	 then our preconditioning matrix M is exactly equal to A (in Ax=b), such that the iterative
	 method will converage after the first iteration.
	 
	 PROOF: (Ax = Mx + Fnp1), where M = GI, the coefficient matrix for the scheme. Subtract
	 the explicit side (b = GE*un + gE - gI + Fn) from the implicit side (Ax), such
	 that (b - Ax) = (b - Mx + Fnp1). If fluxes are determined explicitly, then
	 Fnp1 = 0, making (b - Ax) = (b - Mx). Therefore, A = M. If A = M, then our
	 preconditioned residual becomes r = A^-1 (b - Ax) = A^-1*b - x = 0.
	 */
	
	//Loop to fill in parameters for all nodes
	for (int l=0; l<dat->LN; l++)
	{
		//Boundary 0
		if (l==0)
		{
			//Positive Velocity
			if (dat->vo > 0.0)
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Do);
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Do);
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vo);
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vo);
				
			}
			//Negative Velocity
			else
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
			}
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(0)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		//Boundary 1
		else if (l==dat->LN-1)
		{
			//Positive Velocity
			if (dat->vo > 0.0)
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			}
			//Negative Velocity
			else
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Do);
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Do);
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vo);
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vo);
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			}
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(l-1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		//Interior Nodes
		else
		{
			Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
			Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
			
			Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
			Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
			
			vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
			vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
			
			vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
			vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(l-1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		zl = pow((double)(l)*dat->dz,(double)dat->d);
		
		//Fill in Coefficients--------------------------------------------------------------------
		dat->CL_I[l] = ((z_neg*Dnp1_neg)/pow(dat->dz,2.0)) + ( (z_neg*(vnp1_neg + fabs(vnp1_neg)))/(2.0*dat->dz) );
		dat->CL_E[l] = ((z_neg*Dn_neg)/pow(dat->dz,2.0)) + ( (z_neg*(vn_neg + fabs(vn_neg)))/(2.0*dat->dz) );
		
		dat->CC_I[l] = ((z_pos*Dnp1_pos)/pow(dat->dz,2.0)) + ((z_neg*Dnp1_neg)/pow(dat->dz,2.0)) + ( (z_pos*(vnp1_pos + fabs(vnp1_pos)))/(2.0*dat->dz) ) - ( (z_neg*(vnp1_neg - fabs(vnp1_neg)))/(2.0*dat->dz) ) + (zl*dat->knp1(l,0));
		dat->CC_E[l] = ((z_pos*Dn_pos)/pow(dat->dz,2.0)) + ((z_neg*Dn_neg)/pow(dat->dz,2.0)) + ( (z_pos*(vn_pos + fabs(vn_pos)))/(2.0*dat->dz) ) - ( (z_neg*(vn_neg - fabs(vn_neg)))/(2.0*dat->dz) ) + (zl*dat->kn(l,0));
		
		dat->CR_I[l] = ((z_pos*Dnp1_pos)/pow(dat->dz,2.0)) - ( (z_pos*(vnp1_pos - fabs(vnp1_pos)))/(2.0*dat->dz) );
		dat->CR_E[l] = ((z_pos*Dn_pos)/pow(dat->dz,2.0)) - ( (z_pos*(vn_pos - fabs(vn_pos)))/(2.0*dat->dz) );
		
		dat->fL_I[l] = ( (z_neg*(vnp1_neg + fabs(vnp1_neg)))/4.0 );
		dat->fL_E[l] = ( (z_neg*(vn_neg + fabs(vn_neg)))/4.0 );
		
		dat->fC_I[l] = ( (z_pos*(vnp1_pos + fabs(vnp1_pos)))/4.0 ) + ( (z_neg*(vnp1_neg - fabs(vnp1_neg)))/4.0 );
		dat->fC_E[l] = ( (z_pos*(vn_pos + fabs(vn_pos)))/4.0 ) + ( (z_neg*(vn_neg - fabs(vn_neg)))/4.0 );
		
		dat->fR_I[l] = ( (z_pos*(vnp1_pos - fabs(vnp1_pos)))/4.0 );
		dat->fR_E[l] = ( (z_pos*(vn_pos - fabs(vn_pos)))/4.0 );
		
		dat->Sn(l,0) = zl*dat->Sn(l,0);
		dat->Snp1(1,0) = zl*dat->Snp1(l,0);
		//----------------------------------------------------------------------------------------
		
		//Formulate Discretization----------------------------------------------------------------
		if (dat->SteadyState == true)
		{
			dat->dt = 1.0;
			dat->CN = false;
			dat->beta = 1.0;
		}
		
		//Implicit Side
		dat->MI[l] = -dat->beta*dat->dt*dat->CL_I[l];
		dat->NI[l] = (dat->Rnp1(l,0)*zl)+(dat->beta*dat->dt*dat->CC_I[l]);
		dat->OI[l] = -dat->beta*dat->dt*dat->CR_I[l];
		dat->Snp1(l,0) = dat->beta*dat->dt*dat->Snp1(l,0);
		
		//Explicit Side
		dat->ME[l] = (1.0-dat->beta)*dat->dt*dat->CL_E[l];
		dat->NE[l] = (dat->Rn(l,0)*zl)-((1.0-dat->beta)*dat->dt*dat->CC_E[l]);
		dat->OE[l] = (1.0-dat->beta)*dat->dt*dat->CR_E[l];
		dat->Sn(l,0) = -(1.0-dat->beta)*dat->dt*dat->Sn(l,0);
		
		if (l==0)
			check = dat->NE[l];
		else
		{
			if (dat->NE[l] < check)
				check = dat->NE[l];
		}
		//-------------------------------------------------------------------------------------
		
		//Evaluate Slope Limiter Functions-----------------------------------------------------
		//Node 0
		if (l == 0)
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - uo)/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - uo)/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(uo - uo)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - uo)/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - uo)/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - uo)/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(uo - uo)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - uo)/dat->dz;
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l+1,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l+1,0))/dat->dz;
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				un = (dat->lambda_E*uo) + dat->un(1,0) - (dat->lambda_E*dat->un(0,0));
				u_star = (dat->lambda_I*uo) + dat->u_star(1,0) - (dat->lambda_I*dat->u_star(0,0));;
				
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - u_star)/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - u_star)/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(u_star - u_star)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - u_star)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - u_star)/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - un)/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - un)/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(un - un)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - un)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - un)/dat->dz;
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l+1,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l+1,0))/dat->dz;
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node 1
		else if (l == 1)
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - uo)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - uo)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				un = (dat->lambda_E*uo) + dat->un(1,0) - (dat->lambda_E*dat->un(0,0));
				u_star = (dat->lambda_I*uo) + dat->u_star(1,0) - (dat->lambda_I*dat->u_star(0,0));
				
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - u_star)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - u_star)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - un)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - un)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node L-1
		else if (l == (dat->LN-2))
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (uo - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(uo - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (uo - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(uo - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l-1,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l-1,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				un = (dat->lambda_E*dat->un(dat->LN-1,0)) + dat->un(dat->LN-2,0) - (dat->lambda_E*uo);
				u_star = (dat->lambda_I*dat->u_star(dat->LN-1,0)) + dat->u_star(dat->LN-2,0) - (dat->lambda_I*uo);
				
				dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (u_star - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(u_star - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (un - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(un - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node L
		else if (l == (dat->LN-1))
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l-1,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l-1,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l-1,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l-1,0) - dat->un(l-1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = w*(uo - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (uo - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(uo - uo)/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (uo - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(uo - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(uo - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (uo - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(uo - uo)/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (uo - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(uo - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l-1,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(dat->u_star(l-1,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l-1,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(dat->un(l-1,0) - dat->un(l-1,0))/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				un = (dat->lambda_E*dat->un(dat->LN-1,0)) + dat->un(dat->LN-2,0) - (dat->lambda_E*uo);
				u_star = (dat->lambda_I*dat->u_star(dat->LN-1,0)) + dat->u_star(dat->LN-2,0) - (dat->lambda_I*uo);
				
				dat->uz_lp1_I[0] = w*(u_star - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (u_star - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = w*(u_star - u_star)/dat->dz;
				
				dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (u_star - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = w*(u_star - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = w*(un - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (un - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = w*(un - un)/dat->dz;
				
				dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (un - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = w*(un - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			}
			else {mError(simulation_fail); return -1;}
		}
		//Interior Nodes
		else
		{
			dat->uz_lp1_I[0] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
			dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
			dat->uz_lp1_I[2] = w*(dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
			
			dat->uz_l_I[0] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
			dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
			dat->uz_l_I[2] = w*(dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
			
			dat->uz_lm1_I[0] = w*(dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
			dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
			dat->uz_lm1_I[2] = w*(dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
			
			dat->uz_lp1_E[0] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
			dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
			dat->uz_lp1_E[2] = w*(dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
			
			dat->uz_l_E[0] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
			dat->uz_l_E[2] = w*(dat->un(l+1,0) - dat->un(l,0))/dat->dz;
			
			dat->uz_lm1_E[0] = w*(dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
			dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
			dat->uz_lm1_E[2] = w*(dat->un(l,0) - dat->un(l-1,0))/dat->dz;
		}
		
		//Fill in fluxes
		if (dat->ExplicitFlux == true)
		{
			dat->Fnp1.edit(l,0,0.0);
			dat->Fn.edit(l,0,dat->dt*( (dat->fL_E[l]*minmod(dat->uz_lm1_E)) - (dat->fC_E[l]*minmod(dat->uz_l_E)) + (dat->fR_E[l]*minmod(dat->uz_lp1_E)) ) );
		}
		else
		{
			dat->Fnp1.edit(l,0,-dat->beta*dat->dt*( (dat->fL_I[l]*minmod(dat->uz_lm1_I)) - (dat->fC_I[l]*minmod(dat->uz_l_I)) + (dat->fR_I[l]*minmod(dat->uz_lp1_I)) ) );
			
			dat->Fn.edit(l,0,(1-dat->beta)*dat->dt*( (dat->fL_E[l]*minmod(dat->uz_lm1_E)) - (dat->fC_E[l]*minmod(dat->uz_l_E)) + (dat->fR_E[l]*minmod(dat->uz_lp1_E)) ) );
		}
		//-------------------------------------------------------------------------------------
		
	}//END LOOP
	
	//Check Stability of Explicit Side
	if (dat->CN == true)
	{
		if (check <= 0.0 && dat->d == 0)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		if (dat->vo > 0.0 && (dat->NE[0] - (dat->lambda_E*dat->ME[0])) <= 0.0 && dat->Dirichlet == false)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		else if (dat->vo <= 0.0 && (dat->NE[dat->LN-1] + (dat->lambda_E*dat->OE[dat->LN-1])) <= 0.0 && dat->Dirichlet == false)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		else
		{
			//No action
		}
	}
	
	return success;
	
}

//Function to perform discretization using a smooth van Albada Slope Limiter function
int vanAlbada_discretization(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	double Dn_neg, Dn_pos, vn_neg, vn_pos;
	double Dnp1_neg, Dnp1_pos, vnp1_neg, vnp1_pos;
	double z_pos, z_neg, zl;
	double check = 0.0;
	double un, u_star, uo;
	double rl_I, rlm1_I, rlp1_I, rl_E, rlm1_E, rlp1_E;
	
	//Fill in the lambda parameters for Neumann BC
	if (dat->vo > 0.0 && dat->d == 0)
	{
		dat->lambda_I = (2.0 * dat->dz * dat->vnp1(0,0)) / dat->Dnp1(0,0);
		dat->lambda_E = (2.0 * dat->dz * dat->vn(0,0)) / dat->Dn(0,0);
	}
	else if (dat->vo <= 0.0 && dat->d == 0)
	{
		dat->lambda_I = (2.0 * dat->dz * dat->vnp1(dat->LN-1,0)) / dat->Dnp1(dat->LN-1,0);
		dat->lambda_E = (2.0 * dat->dz * dat->vn(dat->LN-1,0)) / dat->Dn(dat->LN-1,0);
	}
	else
	{
		dat->lambda_I = -(2.0 * dat->dz * dat->kfnp1) / dat->Dnp1(dat->LN-1,0);
		dat->lambda_E = -(2.0 * dat->dz * dat->kfn) / dat->Dn(dat->LN-1,0);
	}
	
	//Check the Neumann parameters for problems
	if (dat->Dirichlet == false)
	{
		if (isnan(dat->lambda_E) || isinf(dat->lambda_E) || isnan(dat->lambda_I) || isinf(dat->lambda_I))
		{
			std::cout << "\nBreakdown of the Neumann BCs!\n Restart with Dirichlet BCs!\n";
			mError(simulation_fail);
			return -1;
		}
	}
	
	//Set up the necessary information for the slope limiters
	if (dat->Iterative == false)
	{
		dat->u_star.columnProjection(dat->un, dat->unm1, dat->dt, dat->dt_old);
		uo = dat->uo;
		//Check for negative concentrations
		success = check_Mass(dat);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	else
	{
		dat->u_star = dat->unp1;
		uo = dat->uo;
	}
	
	/*	NOTE:
	 If we are using the iterative approach, the implicit fluxes will be based on dat->unp1,
	 instead of u_star. Weight factor still depends on ExplicitFlux flag
	 
	 If we are using the Explicit only flux approach, then all fluxes are based on dat->un
	 and we use the smaller weight factor (1.0) for additional stability.
	 
	 As a consequence, if we are using the iterative approach AND the explicit flux approach,
	 then our preconditioning matrix M is exactly equal to A (in Ax=b), such that the iterative
	 method will converage after the first iteration.
	 
	 PROOF: (Ax = Mx + Fnp1), where M = GI, the coefficient matrix for the scheme. Subtract
	 the explicit side (b = GE*un + gE - gI + Fn) from the implicit side (Ax), such
	 that (b - Ax) = (b - Mx + Fnp1). If fluxes are determined explicitly, then
	 Fnp1 = 0, making (b - Ax) = (b - Mx). Therefore, A = M. If A = M, then our
	 preconditioned residual becomes r = A^-1 (b - Ax) = A^-1*b - x = 0.
	 */
	
	//Loop to fill in parameters for all nodes
	for (int l=0; l<dat->LN; l++)
	{
		//Boundary 0
		if (l==0)
		{
			//Positive Velocity
			if (dat->vo > 0.0)
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Do);
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Do);
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vo);
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vo);
				
			}
			//Negative Velocity
			else
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
			}
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(0)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		//Boundary 1
		else if (l==dat->LN-1)
		{
			//Positive Velocity
			if (dat->vo > 0.0)
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			}
			//Negative Velocity
			else
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Do);
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Do);
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vo);
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vo);
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			}
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(l-1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		//Interior Nodes
		else
		{
			Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
			Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
			
			Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
			Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
			
			vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
			vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
			
			vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
			vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(l-1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		zl = pow((double)(l)*dat->dz,(double)dat->d);
		
		//Fill in Coefficients--------------------------------------------------------------------
		dat->CL_I[l] = ((z_neg*Dnp1_neg)/pow(dat->dz,2.0)) + ( (z_neg*(vnp1_neg + fabs(vnp1_neg)))/(2.0*dat->dz) );
		dat->CL_E[l] = ((z_neg*Dn_neg)/pow(dat->dz,2.0)) + ( (z_neg*(vn_neg + fabs(vn_neg)))/(2.0*dat->dz) );
		
		dat->CC_I[l] = ((z_pos*Dnp1_pos)/pow(dat->dz,2.0)) + ((z_neg*Dnp1_neg)/pow(dat->dz,2.0)) + ( (z_pos*(vnp1_pos + fabs(vnp1_pos)))/(2.0*dat->dz) ) - ( (z_neg*(vnp1_neg - fabs(vnp1_neg)))/(2.0*dat->dz) ) + (zl*dat->knp1(l,0));
		dat->CC_E[l] = ((z_pos*Dn_pos)/pow(dat->dz,2.0)) + ((z_neg*Dn_neg)/pow(dat->dz,2.0)) + ( (z_pos*(vn_pos + fabs(vn_pos)))/(2.0*dat->dz) ) - ( (z_neg*(vn_neg - fabs(vn_neg)))/(2.0*dat->dz) ) + (zl*dat->kn(l,0));
		
		dat->CR_I[l] = ((z_pos*Dnp1_pos)/pow(dat->dz,2.0)) - ( (z_pos*(vnp1_pos - fabs(vnp1_pos)))/(2.0*dat->dz) );
		dat->CR_E[l] = ((z_pos*Dn_pos)/pow(dat->dz,2.0)) - ( (z_pos*(vn_pos - fabs(vn_pos)))/(2.0*dat->dz) );
		
		dat->fL_I[l] = ( (z_neg*(vnp1_neg + fabs(vnp1_neg)))/4.0 );
		dat->fL_E[l] = ( (z_neg*(vn_neg + fabs(vn_neg)))/4.0 );
		
		dat->fC_I[l] = ( (z_pos*(vnp1_pos + fabs(vnp1_pos)))/4.0 ) + ( (z_neg*(vnp1_neg - fabs(vnp1_neg)))/4.0 );
		dat->fC_E[l] = ( (z_pos*(vn_pos + fabs(vn_pos)))/4.0 ) + ( (z_neg*(vn_neg - fabs(vn_neg)))/4.0 );
		
		dat->fR_I[l] = ( (z_pos*(vnp1_pos - fabs(vnp1_pos)))/4.0 );
		dat->fR_E[l] = ( (z_pos*(vn_pos - fabs(vn_pos)))/4.0 );
		
		dat->Sn(l,0) = zl*dat->Sn(l,0);
		dat->Snp1(1,0) = zl*dat->Snp1(l,0);
		//----------------------------------------------------------------------------------------
		
		//Formulate Discretization----------------------------------------------------------------
		if (dat->SteadyState == true)
		{
			dat->dt = 1.0;
			dat->CN = false;
			dat->beta = 1.0;
		}
		
		//Implicit Side
		dat->MI[l] = -dat->beta*dat->dt*dat->CL_I[l];
		dat->NI[l] = (dat->Rnp1(l,0)*zl)+(dat->beta*dat->dt*dat->CC_I[l]);
		dat->OI[l] = -dat->beta*dat->dt*dat->CR_I[l];
		dat->Snp1(l,0) = dat->beta*dat->dt*dat->Snp1(l,0);
		
		//Explicit Side
		dat->ME[l] = (1.0-dat->beta)*dat->dt*dat->CL_E[l];
		dat->NE[l] = (dat->Rn(l,0)*zl)-((1.0-dat->beta)*dat->dt*dat->CC_E[l]);
		dat->OE[l] = (1.0-dat->beta)*dat->dt*dat->CR_E[l];
		dat->Sn(l,0) = -(1.0-dat->beta)*dat->dt*dat->Sn(l,0);
		
		if (l==0)
			check = dat->NE[l];
		else
		{
			if (dat->NE[l] < check)
				check = dat->NE[l];
		}
		//-------------------------------------------------------------------------------------
		
		//Evaluate Slope Limiter Functions-----------------------------------------------------
		//Node 0
		if (l == 0)
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - uo)/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - uo)/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (uo - uo)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - uo)/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - uo)/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - uo)/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (uo - uo)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - uo)/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - uo) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - uo) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (uo - uo) / (dat->un(l,0) - uo);
				rlm1_I = (uo - uo) / (dat->u_star(l,0) - uo);
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l+1,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l+1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l+1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l+1,0) - dat->un(l+1,0)) / (dat->un(l,0) - dat->un(l+1,0));
				rlm1_I = (dat->u_star(l+1,0) - dat->u_star(l+1,0)) / (dat->u_star(l,0) - dat->u_star(l+1,0));
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				un = (dat->lambda_E*uo) + dat->un(1,0) - (dat->lambda_E*dat->un(0,0));
				u_star = (dat->lambda_I*uo) + dat->u_star(1,0) - (dat->lambda_I*dat->u_star(0,0));;
				
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - u_star)/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - u_star)/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (u_star - u_star)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - u_star)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - u_star)/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - un)/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - un)/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (un - un)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - un)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - un)/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - un) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - u_star) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (un - un) / (dat->un(l,0) - un);
				rlm1_I = (u_star - u_star) / (dat->u_star(l,0) - u_star);
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l+1,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l+1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l+1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l+1,0) - dat->un(l+1,0)) / (dat->un(l,0) - dat->un(l+1,0));
				rlm1_I = (dat->u_star(l+1,0) - dat->u_star(l+1,0)) / (dat->u_star(l,0) - dat->u_star(l+1,0));
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node 1
		else if (l == 1)
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - uo)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - uo)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - uo) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - uo) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				un = (dat->lambda_E*uo) + dat->un(1,0) - (dat->lambda_E*dat->un(0,0));
				u_star = (dat->lambda_I*uo) + dat->u_star(1,0) - (dat->lambda_I*dat->u_star(0,0));
				
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - u_star)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - u_star)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - un)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - un)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - un) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - u_star) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node L-1
		else if (l == (dat->LN-2))
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (uo - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (uo - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (uo - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (uo - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (uo - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (uo - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l-1,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l-1,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				un = (dat->lambda_E*dat->un(dat->LN-1,0)) + dat->un(dat->LN-2,0) - (dat->lambda_E*uo);
				u_star = (dat->lambda_I*dat->u_star(dat->LN-1,0)) + dat->u_star(dat->LN-2,0) - (dat->lambda_I*uo);
				
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (u_star - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (u_star - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (un - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (un - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (un - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (u_star - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node L
		else if (l == (dat->LN-1))
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l-1,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l-1,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l-1,0) - dat->un(l-1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l-1,0) - dat->un(l,0)) / (dat->un(l-1,0) - dat->un(l-1,0));
				rlp1_I = (dat->u_star(l-1,0) - dat->u_star(l,0)) / (dat->u_star(l-1,0) - dat->u_star(l-1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l-1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l-1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (uo - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (uo - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (uo - uo)/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (uo - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (uo - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (uo - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (uo - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (uo - uo)/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (uo - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (uo - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (uo - dat->un(l,0)) / (uo - uo);
				rlp1_I = (uo - dat->u_star(l,0)) / (uo - uo);
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (uo - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (uo - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l-1,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l-1,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l-1,0) - dat->un(l-1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l-1,0) - dat->un(l,0)) / (dat->un(l-1,0) - dat->un(l-1,0));
				rlp1_I = (dat->u_star(l-1,0) - dat->u_star(l,0)) / (dat->u_star(l-1,0) - dat->u_star(l-1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l-1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l-1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				un = (dat->lambda_E*dat->un(dat->LN-1,0)) + dat->un(dat->LN-2,0) - (dat->lambda_E*uo);
				u_star = (dat->lambda_I*dat->u_star(dat->LN-1,0)) + dat->u_star(dat->LN-2,0) - (dat->lambda_I*uo);
				
				dat->uz_lp1_I[0] = (u_star - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (u_star - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (u_star - u_star)/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (u_star - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (u_star - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (un - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (un - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (un - un)/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (un - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (un - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (un - dat->un(l,0)) / (un - un);
				rlp1_I = (u_star - dat->u_star(l,0)) / (u_star - u_star);
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (un - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (u_star - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else {mError(simulation_fail); return -1;}
		}
		//Interior Nodes
		else
		{
			dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
			dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
			dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
			
			dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
			dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
			dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
			
			dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
			dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
			dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
			
			dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
			dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
			dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
			
			dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
			dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
			
			dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
			dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
			dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			
			rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
			rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
			
			rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
			rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
			
			rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
			rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
		}
		
		//Evaluate the van Albada Slope Limiter
		if (isinf(rlp1_E))
		{
			rlp1_E = 1.0;
		}
		else if (isnan(rlp1_E))
		{
			rlp1_E = 0.0;
		}
		else
		{
			rlp1_E = ( (rlp1_E*rlp1_E) + rlp1_E ) / ( (rlp1_E*rlp1_E) + 1.0 );
		}
		
		if (isinf(rlp1_I))
		{
			rlp1_I = 1.0;
		}
		else if (isnan(rlp1_I))
		{
			rlp1_I = 0.0;
		}
		else
		{
			rlp1_I = ( (rlp1_I*rlp1_I) + rlp1_I ) / ( (rlp1_I*rlp1_I) + 1.0 );
		}
		
		if (isinf(rl_E))
		{
			rl_E = 1.0;
		}
		else if (isnan(rl_E))
		{
			rl_E = 0.0;
		}
		else
		{
			rl_E = ( (rl_E*rl_E) + rl_E ) / ( (rl_E*rl_E) + 1.0 );
		}
		
		if (isinf(rl_I))
		{
			rl_I = 1.0;
		}
		else if (isnan(rl_I))
		{
			rl_I = 0.0;
		}
		else
		{
			rl_I = ( (rl_I*rl_I) + rl_I ) / ( (rl_I*rl_I) + 1.0 );
		}
		
		if (isinf(rlm1_E))
		{
			rlm1_E = 1.0;
		}
		else if (isnan(rlm1_E))
		{
			rlm1_E = 0.0;
		}
		else
		{
			rlm1_E = ( (rlm1_E*rlm1_E) + rlm1_E ) / ( (rlm1_E*rlm1_E) + 1.0 );
		}
		
		if (isinf(rlm1_I))
		{
			rlm1_I = 1.0;
		}
		else if (isnan(rlm1_I))
		{
			rlm1_I = 0.0;
		}
		else
		{
			rlm1_I = ( (rlm1_I*rlm1_I) + rlm1_I ) / ( (rlm1_I*rlm1_I) + 1.0 );
		}
		
		//Fill in fluxes
		if (dat->ExplicitFlux == true)
		{
			dat->Fnp1.edit(l,0,0.0);
			
			if (dat->vn(l,0) > 0.0)
			{
				dat->Fn.edit(l,0,dat->dt*( (dat->fL_E[l]*((1.0-rlm1_E)*dat->uz_lm1_E[0]+rlm1_E*dat->uz_lm1_E[1]) ) - (dat->fC_E[l]*((1.0-rl_E)*dat->uz_l_E[0]+rl_E*dat->uz_l_E[1]) ) + (dat->fR_E[l]*((1.0-rlp1_E)*dat->uz_lp1_E[0]+rlp1_E*dat->uz_lp1_E[1]) ) ) );
			}
			else
			{
				dat->Fn.edit(l,0,dat->dt*( (dat->fL_E[l]*((1.0-rlm1_E)*dat->uz_lm1_E[2]+rlm1_E*dat->uz_lm1_E[1]) ) - (dat->fC_E[l]*((1.0-rl_E)*dat->uz_l_E[2]+rl_E*dat->uz_l_E[1]) ) + (dat->fR_E[l]*((1.0-rlp1_E)*dat->uz_lp1_E[2]+rlp1_E*dat->uz_lp1_E[1]) ) ) );
			}
		}
		else
		{
			if (dat->vnp1(l,0) > 0.0)
			{
				
				dat->Fnp1.edit(l,0,-dat->beta*dat->dt*( (dat->fL_I[l]*((1.0-rlm1_I)*dat->uz_lm1_I[0]+rlm1_I*dat->uz_lm1_I[1]) ) - (dat->fC_I[l]*((1.0-rl_I)*dat->uz_l_I[0]+rl_I*dat->uz_l_I[1]) ) + (dat->fR_I[l]*((1.0-rlp1_I)*dat->uz_lp1_I[0]+rlp1_I*dat->uz_lp1_I[1]) ) ) );
			}
			else
			{
				dat->Fnp1.edit(l,0,-dat->beta*dat->dt*( (dat->fL_I[l]*((1.0-rlm1_I)*dat->uz_lm1_I[2]+rlm1_I*dat->uz_lm1_I[1]) ) - (dat->fC_I[l]*((1.0-rl_I)*dat->uz_l_I[2]+rl_I*dat->uz_l_I[1]) ) + (dat->fR_I[l]*((1.0-rlp1_I)*dat->uz_lp1_I[2]+rlp1_I*dat->uz_lp1_I[1]) ) ) );
			}
			if (dat->vn(l,0) > 0.0)
			{
				dat->Fn.edit(l,0,(1-dat->beta)*dat->dt*( (dat->fL_E[l]*((1.0-rlm1_E)*dat->uz_lm1_E[0]+rlm1_E*dat->uz_lm1_E[1]) ) - (dat->fC_E[l]*((1.0-rl_E)*dat->uz_l_E[0]+rl_E*dat->uz_l_E[1]) ) + (dat->fR_E[l]*((1.0-rlp1_E)*dat->uz_lp1_E[0]+rlp1_E*dat->uz_lp1_E[1]) ) ) );
			}
			else
			{
				dat->Fn.edit(l,0,(1-dat->beta)*dat->dt*( (dat->fL_E[l]*((1.0-rlm1_E)*dat->uz_lm1_E[2]+rlm1_E*dat->uz_lm1_E[1]) ) - (dat->fC_E[l]*((1.0-rl_E)*dat->uz_l_E[2]+rl_E*dat->uz_l_E[1]) ) + (dat->fR_E[l]*((1.0-rlp1_E)*dat->uz_lp1_E[2]+rlp1_E*dat->uz_lp1_E[1]) ) ) );
			}
		}
		//-------------------------------------------------------------------------------------
		
	}//END LOOP
	
	//Check Stability of Explicit Side
	if (dat->CN == true)
	{
		if (check <= 0.0 && dat->d == 0)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		if (dat->vo > 0.0 && (dat->NE[0] - (dat->lambda_E*dat->ME[0])) <= 0.0 && dat->Dirichlet == false)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		else if (dat->vo <= 0.0 && (dat->NE[dat->LN-1] + (dat->lambda_E*dat->OE[dat->LN-1])) <= 0.0 && dat->Dirichlet == false)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		else
		{
			//No action
		}
	}
	
	return success;
}

//Function to perform discretization using the Ospre Slope Limiter
int ospre_discretization(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	double Dn_neg, Dn_pos, vn_neg, vn_pos;
	double Dnp1_neg, Dnp1_pos, vnp1_neg, vnp1_pos;
	double z_pos, z_neg, zl;
	double check = 0.0;
	double un, u_star, uo;
	double rl_I, rlm1_I, rlp1_I, rl_E, rlm1_E, rlp1_E;
	
	//Fill in the lambda parameters for Neumann BC
	if (dat->vo > 0.0 && dat->d == 0)
	{
		dat->lambda_I = (2.0 * dat->dz * dat->vnp1(0,0)) / dat->Dnp1(0,0);
		dat->lambda_E = (2.0 * dat->dz * dat->vn(0,0)) / dat->Dn(0,0);
	}
	else if (dat->vo <= 0.0 && dat->d == 0)
	{
		dat->lambda_I = (2.0 * dat->dz * dat->vnp1(dat->LN-1,0)) / dat->Dnp1(dat->LN-1,0);
		dat->lambda_E = (2.0 * dat->dz * dat->vn(dat->LN-1,0)) / dat->Dn(dat->LN-1,0);
	}
	else
	{
		dat->lambda_I = -(2.0 * dat->dz * dat->kfnp1) / dat->Dnp1(dat->LN-1,0);
		dat->lambda_E = -(2.0 * dat->dz * dat->kfn) / dat->Dn(dat->LN-1,0);
	}
	
	//Check the Neumann parameters for problems
	if (dat->Dirichlet == false)
	{
		if (isnan(dat->lambda_E) || isinf(dat->lambda_E) || isnan(dat->lambda_I) || isinf(dat->lambda_I))
		{
			std::cout << "\nBreakdown of the Neumann BCs!\n Restart with Dirichlet BCs!\n";
			mError(simulation_fail);
			return -1;
		}
	}
	
	//Set up the necessary information for the slope limiters
	if (dat->Iterative == false)
	{
		dat->u_star.columnProjection(dat->un, dat->unm1, dat->dt, dat->dt_old);
		uo = dat->uo;
		//Check for negative concentrations
		success = check_Mass(dat);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	else
	{
		dat->u_star = dat->unp1;
		uo = dat->uo;
	}
	
	/*	NOTE:
	 If we are using the iterative approach, the implicit fluxes will be based on dat->unp1,
	 instead of u_star. Weight factor still depends on ExplicitFlux flag
	 
	 If we are using the Explicit only flux approach, then all fluxes are based on dat->un
	 and we use the smaller weight factor (1.0) for additional stability.
	 
	 As a consequence, if we are using the iterative approach AND the explicit flux approach,
	 then our preconditioning matrix M is exactly equal to A (in Ax=b), such that the iterative
	 method will converage after the first iteration.
	 
	 PROOF: (Ax = Mx + Fnp1), where M = GI, the coefficient matrix for the scheme. Subtract
	 the explicit side (b = GE*un + gE - gI + Fn) from the implicit side (Ax), such
	 that (b - Ax) = (b - Mx + Fnp1). If fluxes are determined explicitly, then
	 Fnp1 = 0, making (b - Ax) = (b - Mx). Therefore, A = M. If A = M, then our
	 preconditioned residual becomes r = A^-1 (b - Ax) = A^-1*b - x = 0.
	 */
	
	//Loop to fill in parameters for all nodes
	for (int l=0; l<dat->LN; l++)
	{
		//Boundary 0
		if (l==0)
		{
			//Positive Velocity
			if (dat->vo > 0.0)
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Do);
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Do);
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vo);
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vo);
				
			}
			//Negative Velocity
			else
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
				
			}
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(0)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		//Boundary 1
		else if (l==dat->LN-1)
		{
			//Positive Velocity
			if (dat->vo > 0.0)
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			}
			//Negative Velocity
			else
			{
				Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Do);
				Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Do);
				
				Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
				Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
				
				vn_pos = 0.5 * (dat->vn(l,0) + dat->vo);
				vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vo);
				
				vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
				vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			}
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(l-1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		//Interior Nodes
		else
		{
			Dn_pos = 0.5 * (dat->Dn(l,0) + dat->Dn(l+1,0));
			Dnp1_pos = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l+1,0));
			
			Dn_neg = 0.5 * (dat->Dn(l,0) + dat->Dn(l-1,0));
			Dnp1_neg = 0.5 * (dat->Dnp1(l,0) + dat->Dnp1(l-1,0));
			
			vn_pos = 0.5 * (dat->vn(l,0) + dat->vn(l+1,0));
			vnp1_pos = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l+1,0));
			
			vn_neg = 0.5 * (dat->vn(l,0) + dat->vn(l-1,0));
			vnp1_neg = 0.5 * (dat->vnp1(l,0) + dat->vnp1(l-1,0));
			
			z_pos = 0.5*(pow((double)(l+1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
			z_neg = 0.5*(pow((double)(l-1)*dat->dz,(double)dat->d)+pow((double)(l)*dat->dz,(double)dat->d));
		}
		zl = pow((double)(l)*dat->dz,(double)dat->d);
		
		//Fill in Coefficients--------------------------------------------------------------------
		dat->CL_I[l] = ((z_neg*Dnp1_neg)/pow(dat->dz,2.0)) + ( (z_neg*(vnp1_neg + fabs(vnp1_neg)))/(2.0*dat->dz) );
		dat->CL_E[l] = ((z_neg*Dn_neg)/pow(dat->dz,2.0)) + ( (z_neg*(vn_neg + fabs(vn_neg)))/(2.0*dat->dz) );
		
		dat->CC_I[l] = ((z_pos*Dnp1_pos)/pow(dat->dz,2.0)) + ((z_neg*Dnp1_neg)/pow(dat->dz,2.0)) + ( (z_pos*(vnp1_pos + fabs(vnp1_pos)))/(2.0*dat->dz) ) - ( (z_neg*(vnp1_neg - fabs(vnp1_neg)))/(2.0*dat->dz) ) + (zl*dat->knp1(l,0));
		dat->CC_E[l] = ((z_pos*Dn_pos)/pow(dat->dz,2.0)) + ((z_neg*Dn_neg)/pow(dat->dz,2.0)) + ( (z_pos*(vn_pos + fabs(vn_pos)))/(2.0*dat->dz) ) - ( (z_neg*(vn_neg - fabs(vn_neg)))/(2.0*dat->dz) ) + (zl*dat->kn(l,0));
		
		dat->CR_I[l] = ((z_pos*Dnp1_pos)/pow(dat->dz,2.0)) - ( (z_pos*(vnp1_pos - fabs(vnp1_pos)))/(2.0*dat->dz) );
		dat->CR_E[l] = ((z_pos*Dn_pos)/pow(dat->dz,2.0)) - ( (z_pos*(vn_pos - fabs(vn_pos)))/(2.0*dat->dz) );
		
		dat->fL_I[l] = ( (z_neg*(vnp1_neg + fabs(vnp1_neg)))/4.0 );
		dat->fL_E[l] = ( (z_neg*(vn_neg + fabs(vn_neg)))/4.0 );
		
		dat->fC_I[l] = ( (z_pos*(vnp1_pos + fabs(vnp1_pos)))/4.0 ) + ( (z_neg*(vnp1_neg - fabs(vnp1_neg)))/4.0 );
		dat->fC_E[l] = ( (z_pos*(vn_pos + fabs(vn_pos)))/4.0 ) + ( (z_neg*(vn_neg - fabs(vn_neg)))/4.0 );
		
		dat->fR_I[l] = ( (z_pos*(vnp1_pos - fabs(vnp1_pos)))/4.0 );
		dat->fR_E[l] = ( (z_pos*(vn_pos - fabs(vn_pos)))/4.0 );
		
		dat->Sn(l,0) = zl*dat->Sn(l,0);
		dat->Snp1(1,0) = zl*dat->Snp1(l,0);
		//----------------------------------------------------------------------------------------
		
		//Formulate Discretization----------------------------------------------------------------
		if (dat->SteadyState == true)
		{
			dat->dt = 1.0;
			dat->CN = false;
			dat->beta = 1.0;
		}
		
		//Implicit Side
		dat->MI[l] = -dat->beta*dat->dt*dat->CL_I[l];
		dat->NI[l] = (dat->Rnp1(l,0)*zl)+(dat->beta*dat->dt*dat->CC_I[l]);
		dat->OI[l] = -dat->beta*dat->dt*dat->CR_I[l];
		dat->Snp1(l,0) = dat->beta*dat->dt*dat->Snp1(l,0);
		
		//Explicit Side
		dat->ME[l] = (1.0-dat->beta)*dat->dt*dat->CL_E[l];
		dat->NE[l] = (dat->Rn(l,0)*zl)-((1.0-dat->beta)*dat->dt*dat->CC_E[l]);
		dat->OE[l] = (1.0-dat->beta)*dat->dt*dat->CR_E[l];
		dat->Sn(l,0) = -(1.0-dat->beta)*dat->dt*dat->Sn(l,0);
		
		if (l==0)
			check = dat->NE[l];
		else
		{
			if (dat->NE[l] < check)
				check = dat->NE[l];
		}
		//-------------------------------------------------------------------------------------
		
		//Evaluate Slope Limiter Functions-----------------------------------------------------
		//Node 0
		if (l == 0)
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - uo)/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - uo)/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (uo - uo)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - uo)/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - uo)/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - uo)/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (uo - uo)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - uo)/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - uo) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - uo) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (uo - uo) / (dat->un(l,0) - uo);
				rlm1_I = (uo - uo) / (dat->u_star(l,0) - uo);
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l+1,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l+1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l+1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l+1,0) - dat->un(l+1,0)) / (dat->un(l,0) - dat->un(l+1,0));
				rlm1_I = (dat->u_star(l+1,0) - dat->u_star(l+1,0)) / (dat->u_star(l,0) - dat->u_star(l+1,0));
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				un = (dat->lambda_E*uo) + dat->un(1,0) - (dat->lambda_E*dat->un(0,0));
				u_star = (dat->lambda_I*uo) + dat->u_star(1,0) - (dat->lambda_I*dat->u_star(0,0));;
				
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - u_star)/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - u_star)/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (u_star - u_star)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - u_star)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - u_star)/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - un)/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - un)/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (un - un)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - un)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - un)/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - un) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - u_star) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (un - un) / (dat->un(l,0) - un);
				rlm1_I = (u_star - u_star) / (dat->u_star(l,0) - u_star);
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l+1,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l+1,0) - dat->un(l+1,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l+1,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l+1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l+1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l+1,0) - dat->un(l+1,0)) / (dat->un(l,0) - dat->un(l+1,0));
				rlm1_I = (dat->u_star(l+1,0) - dat->u_star(l+1,0)) / (dat->u_star(l,0) - dat->u_star(l+1,0));
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node 1
		else if (l == 1)
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - uo)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - uo)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - uo)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - uo) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - uo) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				un = (dat->lambda_E*uo) + dat->un(1,0) - (dat->lambda_E*dat->un(0,0));
				u_star = (dat->lambda_I*uo) + dat->u_star(1,0) - (dat->lambda_I*dat->u_star(0,0));
				
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - u_star)/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - u_star)/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - un)/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - un)/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - un) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - u_star) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node L-1
		else if (l == (dat->LN-2))
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (uo - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (uo - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (uo - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (uo - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (uo - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (uo - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l-1,0) - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l-1,0) - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l,0) - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l,0) - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				un = (dat->lambda_E*dat->un(dat->LN-1,0)) + dat->un(dat->LN-2,0) - (dat->lambda_E*uo);
				u_star = (dat->lambda_I*dat->u_star(dat->LN-1,0)) + dat->u_star(dat->LN-2,0) - (dat->lambda_I*uo);
				
				dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (u_star - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (u_star - dat->u_star(l+1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (un - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (un - dat->un(l+1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (un - dat->un(l+1,0));
				rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (u_star - dat->u_star(l+1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else {mError(simulation_fail); return -1;}
		}
		//Node L
		else if (l == (dat->LN-1))
		{
			
			if (dat->Dirichlet == true && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l-1,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l-1,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l-1,0) - dat->un(l-1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l-1,0) - dat->un(l,0)) / (dat->un(l-1,0) - dat->un(l-1,0));
				rlp1_I = (dat->u_star(l-1,0) - dat->u_star(l,0)) / (dat->u_star(l-1,0) - dat->u_star(l-1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l-1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l-1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == true && dat->vo <= 0.0)
			{
				dat->uz_lp1_I[0] = (uo - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (uo - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (uo - uo)/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (uo - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (uo - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (uo - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (uo - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (uo - uo)/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (uo - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (uo - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (uo - dat->un(l,0)) / (uo - uo);
				rlp1_I = (uo - dat->u_star(l,0)) / (uo - uo);
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (uo - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (uo - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo > 0.0)
			{
				dat->uz_lp1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (dat->u_star(l-1,0) - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (dat->u_star(l-1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (dat->u_star(l-1,0) - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (dat->un(l-1,0) - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (dat->un(l-1,0) - dat->un(l-1,0))/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (dat->un(l-1,0) - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (dat->un(l-1,0) - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (dat->un(l-1,0) - dat->un(l,0)) / (dat->un(l-1,0) - dat->un(l-1,0));
				rlp1_I = (dat->u_star(l-1,0) - dat->u_star(l,0)) / (dat->u_star(l-1,0) - dat->u_star(l-1,0));
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l-1,0) - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l-1,0) - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else if (dat->Dirichlet == false && dat->vo <= 0.0)
			{
				un = (dat->lambda_E*dat->un(dat->LN-1,0)) + dat->un(dat->LN-2,0) - (dat->lambda_E*uo);
				u_star = (dat->lambda_I*dat->u_star(dat->LN-1,0)) + dat->u_star(dat->LN-2,0) - (dat->lambda_I*uo);
				
				dat->uz_lp1_I[0] = (u_star - dat->u_star(l,0))/dat->dz;
				dat->uz_lp1_I[1] = (u_star - dat->u_star(l,0))/(2.0*dat->dz);
				dat->uz_lp1_I[2] = (u_star - u_star)/dat->dz;
				
				dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				dat->uz_l_I[1] = (u_star - dat->u_star(l-1,0))/(2.0*dat->dz);
				dat->uz_l_I[2] = (u_star - dat->u_star(l,0))/dat->dz;
				
				dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
				dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
				
				dat->uz_lp1_E[0] = (un - dat->un(l,0))/dat->dz;
				dat->uz_lp1_E[1] = (un - dat->un(l,0))/(2.0*dat->dz);
				dat->uz_lp1_E[2] = (un - un)/dat->dz;
				
				dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				dat->uz_l_E[1] = (un - dat->un(l-1,0))/(2.0*dat->dz);
				dat->uz_l_E[2] = (un - dat->un(l,0))/dat->dz;
				
				dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
				dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
				dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
				
				rlp1_E = (un - dat->un(l,0)) / (un - un);
				rlp1_I = (u_star - dat->u_star(l,0)) / (u_star - u_star);
				
				rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (un - dat->un(l,0));
				rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (u_star - dat->u_star(l,0));
				
				rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
				rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
			}
			else {mError(simulation_fail); return -1;}
		}
		//Interior Nodes
		else
		{
			dat->uz_lp1_I[0] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
			dat->uz_lp1_I[1] = (dat->u_star(l+2,0) - dat->u_star(l,0))/(2.0*dat->dz);
			dat->uz_lp1_I[2] = (dat->u_star(l+2,0) - dat->u_star(l+1,0))/dat->dz;
			
			dat->uz_l_I[0] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
			dat->uz_l_I[1] = (dat->u_star(l+1,0) - dat->u_star(l-1,0))/(2.0*dat->dz);
			dat->uz_l_I[2] = (dat->u_star(l+1,0) - dat->u_star(l,0))/dat->dz;
			
			dat->uz_lm1_I[0] = (dat->u_star(l-1,0) - dat->u_star(l-2,0))/dat->dz;
			dat->uz_lm1_I[1] = (dat->u_star(l,0) - dat->u_star(l-2,0))/(2.0*dat->dz);
			dat->uz_lm1_I[2] = (dat->u_star(l,0) - dat->u_star(l-1,0))/dat->dz;
			
			dat->uz_lp1_E[0] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
			dat->uz_lp1_E[1] = (dat->un(l+2,0) - dat->un(l,0))/(2.0*dat->dz);
			dat->uz_lp1_E[2] = (dat->un(l+2,0) - dat->un(l+1,0))/dat->dz;
			
			dat->uz_l_E[0] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			dat->uz_l_E[1] = (dat->un(l+1,0) - dat->un(l-1,0))/(2.0*dat->dz);
			dat->uz_l_E[2] = (dat->un(l+1,0) - dat->un(l,0))/dat->dz;
			
			dat->uz_lm1_E[0] = (dat->un(l-1,0) - dat->un(l-2,0))/dat->dz;
			dat->uz_lm1_E[1] = (dat->un(l,0) - dat->un(l-2,0))/(2.0*dat->dz);
			dat->uz_lm1_E[2] = (dat->un(l,0) - dat->un(l-1,0))/dat->dz;
			
			rlp1_E = (dat->un(l+1,0) - dat->un(l,0)) / (dat->un(l+2,0) - dat->un(l+1,0));
			rlp1_I = (dat->u_star(l+1,0) - dat->u_star(l,0)) / (dat->u_star(l+2,0) - dat->u_star(l+1,0));
			
			rl_E = (dat->un(l,0) - dat->un(l-1,0)) / (dat->un(l+1,0) - dat->un(l,0));
			rl_I = (dat->u_star(l,0) - dat->u_star(l-1,0)) / (dat->u_star(l+1,0) - dat->u_star(l,0));
			
			rlm1_E = (dat->un(l-1,0) - dat->un(l-2,0)) / (dat->un(l,0) - dat->un(l-1,0));
			rlm1_I = (dat->u_star(l-1,0) - dat->u_star(l-2,0)) / (dat->u_star(l,0) - dat->u_star(l-1,0));
		}
		
		//Evaluate the Ospre Slope Limiter
		if (isinf(rlp1_E))
		{
			rlp1_E = 1.5;
		}
		else if (isnan(rlp1_E))
		{
			rlp1_E = 0.0;
		}
		else
		{
			rlp1_E = ( 1.5 * ((rlp1_E*rlp1_E) + rlp1_E )) / ( (rlp1_E*rlp1_E) + rlp1_E+ 1.0 );
		}
		
		if (isinf(rlp1_I))
		{
			rlp1_I = 1.5;
		}
		else if (isnan(rlp1_I))
		{
			rlp1_I = 0.0;
		}
		else
		{
			rlp1_I = ( 1.5 * ((rlp1_I*rlp1_I) + rlp1_I )) / ( (rlp1_I*rlp1_I) + rlp1_I+ 1.0 );
		}
		
		if (isinf(rl_E))
		{
			rl_E = 1.5;
		}
		else if (isnan(rl_E))
		{
			rl_E = 0.0;
		}
		else
		{
			rl_E = ( 1.5 * ((rl_E*rl_E) + rl_E )) / ( (rl_E*rl_E) + rl_E+ 1.0 );
		}
		
		if (isinf(rl_I))
		{
			rl_I = 1.5;
		}
		else if (isnan(rl_I))
		{
			rl_I = 0.0;
		}
		else
		{
			rl_I = ( 1.5 * ((rl_I*rl_I) + rl_I )) / ( (rl_I*rl_I) + rl_I+ 1.0 );
		}
		
		if (isinf(rlm1_E))
		{
			rlm1_E = 1.5;
		}
		else if (isnan(rlm1_E))
		{
			rlm1_E = 0.0;
		}
		else
		{
			rlm1_E = ( 1.5 * ((rlm1_E*rlm1_E) + rlm1_E )) / ( (rlm1_E*rlm1_E) + rlm1_E+ 1.0 );
		}
		
		if (isinf(rlm1_I))
		{
			rlm1_I = 1.5;
		}
		else if (isnan(rlm1_I))
		{
			rlm1_I = 0.0;
		}
		else
		{
			rlm1_I = ( 1.5 * ((rlm1_I*rlm1_I) + rlm1_I )) / ( (rlm1_I*rlm1_I) + rlm1_I+ 1.0 );
		}
		
		//Fill in fluxes
		if (dat->ExplicitFlux == true)
		{
			dat->Fnp1.edit(l,0,0.0);
			
			if (dat->vn(l,0) > 0.0)
			{
				dat->Fn.edit(l,0,dat->dt*( (dat->fL_E[l]*((1.0-rlm1_E)*dat->uz_lm1_E[0]+rlm1_E*dat->uz_lm1_E[1]) ) - (dat->fC_E[l]*((1.0-rl_E)*dat->uz_l_E[0]+rl_E*dat->uz_l_E[1]) ) + (dat->fR_E[l]*((1.0-rlp1_E)*dat->uz_lp1_E[0]+rlp1_E*dat->uz_lp1_E[1]) ) ) );
			}
			else
			{
				dat->Fn.edit(l,0,dat->dt*( (dat->fL_E[l]*((1.0-rlm1_E)*dat->uz_lm1_E[2]+rlm1_E*dat->uz_lm1_E[1]) ) - (dat->fC_E[l]*((1.0-rl_E)*dat->uz_l_E[2]+rl_E*dat->uz_l_E[1]) ) + (dat->fR_E[l]*((1.0-rlp1_E)*dat->uz_lp1_E[2]+rlp1_E*dat->uz_lp1_E[1]) ) ) );
			}
		}
		else
		{
			if (dat->vnp1(l,0) > 0.0)
			{
				
				dat->Fnp1.edit(l,0,-dat->beta*dat->dt*( (dat->fL_I[l]*((1.0-rlm1_I)*dat->uz_lm1_I[0]+rlm1_I*dat->uz_lm1_I[1]) ) - (dat->fC_I[l]*((1.0-rl_I)*dat->uz_l_I[0]+rl_I*dat->uz_l_I[1]) ) + (dat->fR_I[l]*((1.0-rlp1_I)*dat->uz_lp1_I[0]+rlp1_I*dat->uz_lp1_I[1]) ) ) );
			}
			else
			{
				dat->Fnp1.edit(l,0,-dat->beta*dat->dt*( (dat->fL_I[l]*((1.0-rlm1_I)*dat->uz_lm1_I[2]+rlm1_I*dat->uz_lm1_I[1]) ) - (dat->fC_I[l]*((1.0-rl_I)*dat->uz_l_I[2]+rl_I*dat->uz_l_I[1]) ) + (dat->fR_I[l]*((1.0-rlp1_I)*dat->uz_lp1_I[2]+rlp1_I*dat->uz_lp1_I[1]) ) ) );
			}
			if (dat->vn(l,0) > 0.0)
			{
				dat->Fn.edit(l,0,(1-dat->beta)*dat->dt*( (dat->fL_E[l]*((1.0-rlm1_E)*dat->uz_lm1_E[0]+rlm1_E*dat->uz_lm1_E[1]) ) - (dat->fC_E[l]*((1.0-rl_E)*dat->uz_l_E[0]+rl_E*dat->uz_l_E[1]) ) + (dat->fR_E[l]*((1.0-rlp1_E)*dat->uz_lp1_E[0]+rlp1_E*dat->uz_lp1_E[1]) ) ) );
			}
			else
			{
				dat->Fn.edit(l,0,(1-dat->beta)*dat->dt*( (dat->fL_E[l]*((1.0-rlm1_E)*dat->uz_lm1_E[2]+rlm1_E*dat->uz_lm1_E[1]) ) - (dat->fC_E[l]*((1.0-rl_E)*dat->uz_l_E[2]+rl_E*dat->uz_l_E[1]) ) + (dat->fR_E[l]*((1.0-rlp1_E)*dat->uz_lp1_E[2]+rlp1_E*dat->uz_lp1_E[1]) ) ) );
			}
		}
		//-------------------------------------------------------------------------------------
		
	}//END LOOP
	
	//Check Stability of Explicit Side
	if (dat->CN == true)
	{
		if (check <= 0.0 && dat->d == 0)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		if (dat->vo > 0.0 && (dat->NE[0] - (dat->lambda_E*dat->ME[0])) <= 0.0 && dat->Dirichlet == false)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		else if (dat->vo <= 0.0 && (dat->NE[dat->LN-1] + (dat->lambda_E*dat->OE[dat->LN-1])) <= 0.0 && dat->Dirichlet == false)
		{
			if (dat->NormTrack == true)
				std::cout << "Instability in Explicit Side... Converting to full Implicit..." << std::endl;
			dat->beta = 1.0;
			dat->CN = false;
			success = (*dat->discretize) (dat);
		}
		else
		{
			//No action
		}
	}
	
	return success;
}

//Function to apply BCs (either Dirichlet or Neumann based on boolean)
int default_bcs(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	//Positive Velocity
	if (dat->vo > 0.0 && dat->d == 0)
	{
		//Fill in Input Dirichlet BC
		if (dat->Dirichlet == true)
		{
			dat->gI.dirichletBCFill(0, dat->MI[0], dat->uo);
			dat->gE.dirichletBCFill(0, dat->ME[0], dat->uo);
		}
		//Fill in Input Neumann BC
		else
		{
			dat->gI.dirichletBCFill(0, (dat->lambda_I*dat->MI[0]), dat->uo);
			dat->gE.dirichletBCFill(0, (dat->lambda_E*dat->ME[0]), dat->uo);
			
			dat->NI[0] = dat->NI[0] - (dat->lambda_I * dat->MI[0]);
			dat->NE[0] = dat->NE[0] - (dat->lambda_E * dat->ME[0]);
			dat->OI[0] = dat->OI[0] + dat->MI[0];
			dat->OE[0] = dat->OE[0] + dat->ME[0];
		}
		dat->MI[dat->LN-1] = dat->MI[dat->LN-1] + dat->OI[dat->LN-1];
		dat->ME[dat->LN-1] = dat->ME[dat->LN-1] + dat->OE[dat->LN-1];
	}
	//Negative Velocity
	else
	{
		//Fill in Input Dirichlet BC
		if (dat->Dirichlet == true)
		{
			dat->gI.dirichletBCFill(dat->LN-1, dat->OI[dat->LN-1], dat->uo);
			dat->gE.dirichletBCFill(dat->LN-1, dat->OE[dat->LN-1], dat->uo);
		}
		//Fill in Input Neumann BC
		else
		{
			dat->gI.dirichletBCFill(dat->LN-1, (-dat->lambda_I*dat->OI[dat->LN-1]), dat->uo);
			dat->gE.dirichletBCFill(dat->LN-1, (-dat->lambda_E*dat->OE[dat->LN-1]), dat->uo);
			
			dat->NI[dat->LN-1] = dat->NI[dat->LN-1] + (dat->lambda_I*dat->OI[dat->LN-1]);
			dat->NE[dat->LN-1] = dat->NE[dat->LN-1] + (dat->lambda_E*dat->OE[dat->LN-1]);
			dat->MI[dat->LN-1] = dat->MI[dat->LN-1] + dat->OI[dat->LN-1];
			dat->ME[dat->LN-1] = dat->ME[dat->LN-1] + dat->OE[dat->LN-1];
		}
		
		//Fill in Exit Neumann BC (first row, 2nd column)
		if (dat->d == 0)
		{
			dat->OI[0] = dat->MI[0] + dat->OI[0];
			dat->OE[0] = dat->ME[0] + dat->OE[0];
		}
		else
		{
			dat->NI[0] = 1.0;
			dat->NE[0] = 1.0;
			dat->OI[0] = -1.0;
			dat->OE[0] = -1.0;
		}
	}
	
	return success;
}

//Function to evalute the residual (res) of the current iterate (x) using user data
int default_res(const Matrix<double> &x, Matrix<double> &res, const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	//Set Up all parameters
	success = (*dat->setparams) (dat->param_data);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Set Up all Coefficients and slope limiters
	success = (*dat->discretize) (dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Set up the BCs for the system
	success = (*dat->setbcs) (dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Form the residual line by line
	for (int l=0; l<dat->LN; l++)
	{
		if (l==0)
		{
			res.edit(l, 0, (dat->NE[l]*dat->un(l,0) + dat->OE[l]*dat->un(l+1,0) + dat->gE(l,0) - dat->gI(l,0) + dat->Fn(l,0) + dat->Sn(l,0)) - (dat->NI[l]*x(l,0) + dat->OI[l]*x(l+1,0) + dat->Fnp1(l,0)) + dat->Snp1(l,0));
		}
		else if (l==dat->LN-1)
		{
			res.edit(l, 0, (dat->ME[l]*dat->un(l-1,0) + dat->NE[l]*dat->un(l,0) + dat->gE(l,0) - dat->gI(l,0) + dat->Fn(l,0) + dat->Sn(l,0)) - (dat->MI[l]*x(l-1,0) + dat->NI[l]*x(l,0) + dat->Fnp1(l,0)) + dat->Snp1(l,0));
		}
		else
		{
			res.edit(l, 0, (dat->ME[l]*dat->un(l-1,0) + dat->NE[l]*dat->un(l,0) + dat->OE[l]*dat->un(l+1,0) + dat->gE(l,0) - dat->gI(l,0) + dat->Fn(l,0) + dat->Sn(l,0)) - (dat->MI[l]*x(l-1,0) + dat->NI[l]*x(l,0) + dat->OI[l]*x(l+1,0) + dat->Fnp1(l,0)) + dat->Snp1(l,0));
		}
	}
	
	return success;
}

//Function to solve the Tridiagonal Matrix<double> Directly (Can be used as preconditioner)
int default_precon(const Matrix<double> &b, Matrix<double> &p, const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	double d[dat->LN], a[dat->LN], c[dat->LN];
	double dp[dat->LN], ap[dat->LN];
	
	//Forward Sweep
	for (int l=0; l<dat->LN; l++)
	{
		if (dat->NI[l] == 0.0)
		{
			mError(singular_matrix);
			return -1;
		}
		
		if (l==0)
		{
			a[l] = 0.0;
			c[l] = dat->OI[l]/dat->NI[l];
		}
		else if (l==dat->LN-1)
		{
			a[l] = dat->MI[l]/dat->NI[l];
			c[l] = 0.0;
		}
		else
		{
			a[l] = dat->MI[l]/dat->NI[l];
			c[l] = dat->OI[l]/dat->NI[l];
		}
		
		d[l] = b(l,0)/dat->NI[l];
	}
	
	//Reverse Sweep
	for (int l=dat->LN-1; l>=0; l--)
	{
		if (l==dat->LN-1)
		{
			dp[l] = d[l];
			ap[l] = a[l];
		}
		else if (l==0)
		{
			dp[l] = ( d[l] - (c[l]*dp[l+1]) ) / ( 1.0 - (c[l]*ap[l+1]) );
			ap[l] = 0.0;
		}
		else
		{
			dp[l] = ( d[l] - (c[l]*dp[l+1]) ) / ( 1.0 - (c[l]*ap[l+1]) );
			ap[l] = ( a[l] ) / ( 1.0 - (c[l]*ap[l+1]) );
		}
	}
	
	//Forward Sweep
	for (int l=0; l<dat->LN; l++)
	{
		if (l==0)
			p.edit(l,0,dp[l]);
		else
			p.edit(l,0,dp[l] - (ap[l]*p(l-1,0)) );
	}
	
	return success;
}

//Default Post-Process Function
int default_postprocess(const void *user_data)
{
	std::cout << "Postprocess Actions..." << std::endl;
	return 0;
}

int default_reset(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	dat->unm1 = dat->un;
	dat->un = dat->unp1;
	dat->vn = dat->vnp1;
	dat->Dn = dat->Dnp1;
	dat->kn = dat->knp1;
	dat->Rn = dat->Rnp1;
	dat->kfn = dat->kfnp1;
	dat->uT_old = dat->uT;
	dat->uAvg_old = dat->uAvg;
	dat->dt_old = dat->dt;
	dat->t_old = dat->t;
	
	return success;
}

//Test case initial conditions for non-linear problem
int buckley_leverett_ic(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	double zl;
	
	dat->uo = 1.0;
	dat->vo = dat->uo;
	for (int l=0; l<dat->LN; l++)
	{
		if (dat->SteadyState == true)
			dat->Rn.edit(l, 0, 0.0);
		else
			dat->Rn.edit(l, 0, dat->RIC);
		zl = (double)(l+1)*dat->dz;
		
		if (zl <= 0.33333)
		{
			dat->un.edit(l, 0, (1.0 - 3.0*zl));
			dat->unm1.edit(l, 0, dat->un(l,0));
			dat->unp1.edit(l, 0, dat->un(l,0));
		}
		else
		{
			dat->un.edit(l, 0, 0.0);
			dat->unm1.edit(l, 0, dat->un(l,0));
			dat->unp1.edit(l, 0, dat->un(l,0));
		}
		
		dat->Dn.edit(l, 0, (4.0*0.01*(1.0 - dat->un(l,0))) );
		dat->vn.edit(l, 0, (dat->un(l,0) / (pow(dat->un(l,0), 2.0) + pow((1.0-dat->un(l,0)), 2.0)) ) );
		
		dat->kn.edit(l, 0, 0.0);
		dat->Sn.edit(l, 0, 0.0);
		
	}
	
	if (dat->SteadyState == true)
		dat->unp1.ConstantICFill(0.0);
	
	success = uTotal(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	success = uAverage(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	dat->uT_old = dat->uT;
	dat->uAvg_old = dat->uAvg;
	
	return success;
}

//Test case parameters for non-linear problem
int buckley_leverett_params(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	double zl;
	
	dat->vo = dat->uo;
	dat->Do = 0.0;
	for (int l=0; l<dat->LN; l++)
	{
		if (dat->SteadyState == true)
			dat->Rnp1.edit(l, 0, 0.0);
		else
			dat->Rnp1.edit(l, 0, dat->RIC);
		zl = (double)(l+1)*dat->dz;
		
		dat->Dnp1.edit(l, 0, (4.0*0.01*(1.0 - dat->unp1(l,0))) );
		dat->vnp1.edit(l, 0, (dat->unp1(l,0) / (pow(dat->unp1(l,0), 2.0) + pow((1.0-dat->unp1(l,0)), 2.0)) ) );
		
		dat->knp1.edit(l, 0, 0.0);
		dat->Snp1.edit(l, 0, 0.0);
		
	}
	
	return success;
}

//Function for the custom initial conditions for the inviscous Burger's Problem
int burgers_ic(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	double zl;
	
	dat->uo = 0.5;
	dat->vo = 0.5*dat->uo;
	for (int l=0; l<dat->LN; l++)
	{
		if (dat->SteadyState == true)
			dat->Rn.edit(l, 0, 0.0);
		else
			dat->Rn.edit(l, 0, 1.0);
		zl = (double)(l+1)*dat->dz;
		
		dat->Dn.edit(l, 0, 0.0);
		dat->kn.edit(l, 0, 0.0);
		dat->Sn.edit(l, 0, 0.0);
		
		dat->un.edit(l, 0, 0.5 + sin(zl));
		dat->unm1.edit(l, 0, dat->un(l,0));
		dat->unp1.edit(l, 0, dat->un(l,0));
		
		dat->vn.edit(l, 0, 0.5*dat->un(l,0));
	}
	
	success = uTotal(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	success = uAverage(dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	dat->uT_old = dat->uT;
	dat->uAvg_old = dat->uAvg;
	
	return success;
}

//Function for the custom parameters for the inviscous Burger's Problem
int burgers_params(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	dat->uo = dat->unp1(dat->LN-1,0);
	dat->vo = 0.5*dat->uo;
	for (int l=0; l<dat->LN; l++)
	{
		if (dat->SteadyState == true)
			dat->Rnp1.edit(l, 0, 0.0);
		else
			dat->Rnp1.edit(l, 0, 1.0);
		
		dat->Dnp1.edit(l, 0, 0.0);
		dat->knp1.edit(l, 0, 0.0);
		dat->Snp1.edit(l, 0, 0.0);
		dat->vnp1.edit(l, 0, 0.5*dat->unp1(l,0));
	}
	
	return success;
}

//Function for the custom boundary conditions for the inviscous Burger's Problem
int burgers_bcs(const void *user_data)
{
	int success = 0;
	FINCH_DATA *dat = (FINCH_DATA *) user_data;
	
	dat->uo = dat->unp1(dat->LN-1,0);
	dat->vo = 0.5*dat->uo;
	
	if (dat->Dirichlet == false)
	{
		mError(invalid_boolean);
		return -1;
	}
	else
	{
		dat->gI.dirichletBCFill(0, dat->MI[0], dat->unp1(dat->LN-1,0));
		dat->gE.dirichletBCFill(0, dat->ME[0], dat->un(dat->LN-1,0));
		
		dat->gI.edit(dat->LN-1, 0, dat->OI[dat->LN-1]*dat->unp1(0,0));
		dat->gE.edit(dat->LN-1, 0, dat->OE[dat->LN-1]*dat->un(0,0));
	}
	
	/*
		NOTE: This demonstrate how one would handle a Periodic BC in FINCH. For this
	 case, the Dirichlet BC would be used and we specify that the value at
	 the input node is equal to the value at the output node. However, the
	 input node is not expilictly solved for, whereas the output node is. Thus,
	 we solve for the output node in the mesh and use it's value as the input
	 node. Best to use an iterative method for this type of BC.
	 */
	
	
	return success;
}

//Function runs the FINCH tests for convergence, accuracy, and stability when tasked
int FINCH_TESTS()
{
	int success = 0;
	
	/* 					Testing of the scheme				*/
	
	//Declarations
	FINCH_DATA dat;
	double time;
	FILE *Output;
	double exponent = 0.0;
	double truncErr;
	int time_steps = 0;
	
	//Initializations
	time = clock();
	Output = fopen("output/FINCH_TEST_Output.txt","w+");
	if (Output == nullptr)
	{
		system("mkdir output");
		Output = fopen("output/FINCH_TEST_Output.txt","w+");
	}
	
	//Change Parameters for Testing
	dat.uIC = 0.0;
	//dat.uIC = 1.0;
	dat.uo = 1.0;
	//dat.uo = 0.0;
	//dat.vIC = 992.2941164;
	//dat.vo = 992.2941164;
	dat.vIC = 2.0;
	dat.vo = 2.0;
	//dat.vIC = 1.0;
	//dat.vo = 1.0;
	//dat.DIC = 0.546244074;
	//dat.Do = 0.546244074;
	dat.DIC = 0.01;
	dat.Do = 0.01;
	//dat.DIC = 0.0;
	//dat.Do = 0.0;
	//dat.kIC = 0.0;
	//dat.ko = 0.0;
	dat.kIC = 5.0;
	dat.ko = 5.0;
	dat.RIC = 1.0;
	dat.Ro = 1.0;
	//dat.RIC = 285991.8319;
	//dat.Ro = 285991.8319;
	dat.kfn = 0.0;
	dat.kfnp1 = 0.0;
	//dat.L = 0.127;
	dat.L = 1.0;
	//dat.L = 2.0*M_PI;
	dat.s = 1.0;
	dat.T = 0.2;
	//dat.T = 60.0;
	dat.LN = 40;
	dat.t_old = 0.0;
	dat.dt_old = 0.0;
	dat.d = 0;
	
	//Boolean Statments
	dat.Dirichlet = true;
	dat.CheckMass = false;
	dat.Iterative = true;
	dat.SteadyState = false;
	dat.NormTrack = true;
	
	//Iterative Methods
	dat.nl_method = LARK_PJFNK; //0 = FINCH_Picard, 1 = LARK_Picard, 2 = LARK_PJFNK
	dat.pjfnk_dat.nl_tol_rel = 1e-6;
	dat.pjfnk_dat.nl_tol_abs = 1e-6;
	dat.pjfnk_dat.linear_solver = QR;
	//dat.pjfnk_dat.L_Output = true;
	//dat.pjfnk_dat.lin_tol = 1e-10;
	dat.pjfnk_dat.LineSearch = true;
	dat.pjfnk_dat.Bounce = true;
	
	/*
	 After extensive testing, we can show that our Picard iteration is the most
		efficient solution method. However, PJFNK is still good and will be useful
	 for solving more complex, non-linear systems.
	 */
	
	//Used in determining truncation error
	if (dat.CN == true)
		exponent = exponent + 1.0;
	else
		exponent = exponent + 0.5;
	if (dat.ExplicitFlux == false)
		exponent = exponent + 1.0;
	else
		exponent = exponent + 0.5;
	
	//Set up the FINCH_DATA
	
	//Buckley-Leverett Non-Linear Tests with Default Dirichlet BCs
	success = setup_FINCH_DATA(default_execution,buckley_leverett_ic,default_timestep,default_preprocess,default_solve,buckley_leverett_params,vanAlbada_discretization,default_bcs,default_res,default_precon,default_postprocess,default_reset,&dat,(void *)&dat);
	
	//Inviscous Burger's Non-Linear Tests with Periodic BCs
	//success = setup_FINCH_DATA(default_execution,burgers_ic,default_timestep,default_preprocess,default_solve,burgers_params,minmod_discretization,burgers_bcs,default_res,default_precon,default_postprocess,default_reset,&dat,(void *)&dat);
	
	//Below uses minmod discretization (least dispersive, least oscillatory, worst convergence)
	//success = setup_FINCH_DATA(default_execution,default_ic,default_timestep,default_preprocess,default_solve,default_params,minmod_discretization,default_bcs,default_res,default_precon,default_postprocess,default_reset,&dat,(void *)&dat);
	
	//Below uses Ospre discretization (less dispersive, less oscillatory, better convergence)
	//success = setup_FINCH_DATA(default_execution,default_ic,default_timestep,default_preprocess,default_solve,default_params,ospre_discretization,default_bcs,default_res,default_precon,default_postprocess,default_reset,&dat,(void *)&dat);
	
	//Below uses van Albada discretization (most dispersive, most oscillations, best convergence)
	//success = setup_FINCH_DATA(default_execution,default_ic,default_timestep,default_preprocess,default_solve,default_params,vanAlbada_discretization,default_bcs,default_res,default_precon,default_postprocess,default_reset,&dat,(void *)&dat);
	
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Make header file for output
	print2file_dim_header(Output, &dat);
	print2file_newline(Output, &dat);
	print2file_time_header(Output, &dat);
	print2file_newline(Output, &dat);
	
	//Set Initial Conditions
	success = (*dat.setic) (&dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Print out ICs
	print2file_result_old(Output, &dat);
	print2file_newline(Output, &dat);
	
	//Loop to solve for each time step in simulation
	do
	{
		//Check to see if system needs updating
		if (dat.Update == true)
		{
			success = (*dat.resettime) ((void *)&dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		
		//Step size based of off CFL condition
		success = (*dat.settime) ((void *)&dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		//dat.dt = 0.1;
		if (dat.SteadyState == false)
			dat.t = dat.t_old + dat.dt;
		else
			dat.t = INFINITY;
		
		//Call the routine
		std::cout << "Evaluating Time: " << dat.t << std::endl;
		success = (*dat.callroutine) ((void *)&dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			dat.Update = true; //Be sure to set update = true if simulation successful!
		}
		else {mError(simulation_fail); dat.Update = false; return -1;}
		
		//Print out simulation results
		print2file_result_new(Output, &dat);
		print2file_newline(Output, &dat);
		
		time_steps++;
		
	} while (dat.t < (dat.T) && dat.SteadyState == false);
	
	//END PROGRAM
	fclose(Output);
	time = clock() - time;
	
	if (dat.SteadyState == false)
		truncErr = pow(dat.dz,2.0) + pow(dat.dt,exponent);
	else
		truncErr = pow(dat.dz,2.0);
	
	//Display performance metrics for tests
	std::cout << "Runtime Time (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "Truncation Error:\t" << truncErr << std::endl;
	std::cout << "Total Iterations:\t" << dat.total_iter << std::endl;
	std::cout << "Total Time Steps:\t" << time_steps+1 << std::endl;
	std::cout << "Average Iterations:\t" << (double)dat.total_iter/(time_steps+1) << std::endl;
	std::cout << "Complexity (ms):\t" << (time / CLOCKS_PER_SEC)/(double)(dat.total_iter+time_steps)*1000.0 << std::endl;
	return success;
}
