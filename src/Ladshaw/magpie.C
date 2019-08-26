/*!
 *  \file magpie.h magpie.cpp
 *	\brief Multicomponent Adsorption Generalized Procedure for Isothermal Equilibria
 *  \author Austin Ladshaw
 *	\date 12/17/2013
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "magpie.h"

//Function to evaluate the adsorption capacity of the GSTA model given a pressure and the component index
double qo(double po, const void *data, int i)
{
	double qo = 0, top = 0, bot = 1;
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	
	for (int n=0; n<dat->gsta_dat[i].m; n++)
	{
		double tempKo = exp( lnKo(dat->gsta_dat[i].dHo[n], dat->gsta_dat[i].dSo[n], dat->sys_dat.T) );
		top = top + ( (n+1) * tempKo * pow( (po / Po) , (n+ 1) ) );
		bot = bot + (tempKo * pow( (po / Po) ,(n+1)));
	}
	qo = (dat->gsta_dat[i].qmax / dat->gsta_dat[i].m) * (top / bot);
	
	return qo;
}

//Function calculates the gradient of the GSTA isotherm
double dq_dp(double p, const void *data, int i)
{
	double gradient = 0;
	double sum1 = 0.0, sum2 = 1.0, sum3 = 0.0, sum4 = 0.0;
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	
	for (int n=0; n<dat->gsta_dat[i].m; n++)
	{
		double tempKo = exp( lnKo(dat->gsta_dat[i].dHo[n], dat->gsta_dat[i].dSo[n], dat->sys_dat.T) );
		sum1 = sum1 + ( (n+1) * (n+1) * tempKo * pow( (p / Po) , (n) ) );
		sum2 = sum2 + (tempKo * pow( (p / Po) ,(n+1)));
		sum3 = sum3 + ( (n+1) * tempKo * pow( (p / Po) , (n+1) ) );
		sum4 = sum4 + ( (n+1) * tempKo * pow( (p / Po) , (n) ) );
	}
	gradient = (dat->gsta_dat[i].qmax / dat->gsta_dat[i].m) *
	( (sum1*sum2) - (sum3*sum4) ) / pow(sum2,2.0) / Po;
	
	return gradient;
}

//Calculate the value of q/p for a given species
double q_p(double p, const void *data, int i)
{
	double ratio = 0.0;
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	
	double tempK1 = exp( lnKo(dat->gsta_dat[i].dHo[0], dat->gsta_dat[i].dSo[0], dat->sys_dat.T) );
	
	if (p <= 0.0)
	{
		ratio = He(dat->gsta_dat[i].qmax,tempK1,dat->gsta_dat[i].m);
	}
	else
	{
		ratio = qo(p,data,i)/p;
	}
	
	return ratio;
}

//Function to evaluate the spreading pressure term PI from the AST integral based on a give po and component index
double PI(double po, const void *data, int i)
{
	double PI = 0, bot = 1;
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	
	for(int n=0; n<dat->gsta_dat[i].m; n++)
	{
		double tempKo = exp( lnKo(dat->gsta_dat[i].dHo[n], dat->gsta_dat[i].dSo[n], dat->sys_dat.T) );
		bot = bot + (tempKo * pow( (po / Po) , (n+1) ));
	}
	PI = (dat->gsta_dat[i].qmax / dat->gsta_dat[i].m) * log(bot);
	
	return PI;
}

//Function to calculate the maximum lateral energies
double eMax(const void *data, int i)
{
	double mu = 0;
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	
	if (dat->gsta_dat[i].m == 1)
	{
		mu = 1e-6;
	}
	else
	{
		mu = (2 * (-dat->gsta_dat[i].dHo[(dat->gsta_dat[i].m-1)] +
				   dat->gsta_dat[i].dHo[(dat->gsta_dat[i].m-2)] + dat->gsta_dat[i].dHo[0]) ) /
		(Z * dat->mspd_dat[i].s);
	}
	return (2.0*mu);
}

//Function to evaluate the isosteric heat of adsorption based on the pure component isotherm parameters given the index
double Qst(double po, const void *data, int i)
{
	double Qst = 0, top = 0, bot = 0;
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	
	double phi = qo(po, dat, i) / dat->gsta_dat[i].qmax;
	
	for(int n=0; n<dat->gsta_dat[i].m; n++)
	{
		double tempKo = exp( lnKo(dat->gsta_dat[i].dHo[n], dat->gsta_dat[i].dSo[n], dat->sys_dat.T) );
		top = top + ( ( (dat->gsta_dat[i].m*phi) - (n+1) ) * tempKo * pow((po/Po),(n+1)) * (-1*dat->gsta_dat[i].dHo[n]) );
		bot = bot + ( ( (dat->gsta_dat[i].m*phi) - (n+1) ) * (n+1) * tempKo * pow((po/Po),(n+1)) );
	}
	if (po <= 0.0)
		Qst = -1*dat->gsta_dat[i].dHo[0];
	else
		Qst = top / bot;
	
	return Qst;
}

//Function to evaluate the activities from the SPD model using the correct geometric mean
double lnact_mSPD(const double *par, const void *data, int i, volatile double PI)
{
	/*
	 * Parameter Convention:
	 * ---------------------
	 * par[0] = PI			par[(i+1)] = x[i]
	 * fabs					fabs
	 *
	 * n_par = 1 + dat.N
	 * m_dat = 1 + dat.N
	 *
	 * Make sure that all parameters are set accordingly before this function is called
	 * otherwise, errors will occur...
	 */
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	double LNact = 0, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
	int N = dat->sys_dat.N;
	double theta[N];
	double T[N][N], e[N][N];
	double alpha[N][N];
	
	//Forward: GPAST evaluations; Reverse: All Evaluations
	if (par != NULL)
	{
		//Loop to fill out theta array & e[j][j] & T[j][j]
		for(int j=0; j<N; j++)
		{
			for(int l=0; l<N; l++)
			{
				if (dat->sys_dat.Recover == false)
					sum4 = sum4 + (dat->mspd_dat[l].s * (fabs(par[(l+1)])) );
				else
					sum4 = sum4 + (dat->mspd_dat[l].s * dat->gpast_dat[l].x );
				if(j==l)
				{
					double po[1];
					dat->sys_dat.J = j;
					
					if (dat->sys_dat.Recover == false)
					{
						po[0] = (dat->sys_dat.PT * dat->gpast_dat[j].y) / (fabs(par[(j+1)]) );
						dat->sys_dat.PI = fabs(par[0]);
					}
					else
					{
						po[0] = (dat->sys_dat.PT * fabs(par[(j+1)])) / ( dat->gpast_dat[j].x );
						dat->sys_dat.PI = PI;
					}
					lm_status_struct status;
					lm_control_struct control = lm_control_double;
					lmmin(1, po, 1, dat, eval_po_PI, &control, &status, lm_printout_std);
					dat->sys_dat.total_eval = dat->sys_dat.total_eval + status.nfev;
					if (status.info > 5) {std::cout << "\nWarning! Non-convergent intermediate steps!" << std::endl; return 0.0;}
					
					e[j][j] = (2 * (Qst(fabs(po[0]), dat, j) + dat->gsta_dat[j].dHo[0]) ) / (Z * dat->mspd_dat[j].s);
					
					T[j][j] = 1.0;
				}
			}
			if (dat->sys_dat.Recover == false)
				theta[j] = (dat->mspd_dat[j].s * (fabs(par[(j+1)]))) / sum4;
			else
				theta[j] = (dat->mspd_dat[j].s * dat->gpast_dat[j].x) / sum4;
			sum4 = 0;
		}
		//Loop to fill out e[i][j] & T[i][j]
		for (int j=0; j<N; j++)
		{
			for(int l=0; l<N; l++)
			{
				if(j!=l)
				{
					//Calculate the shift factor for geometric mean
					double shift = sqrt( fabs(dat->mspd_dat[l].eMax * dat->mspd_dat[j].eMax) );
					
					if (dat->mspd_dat[j].eta[l] == 0 && dat->mspd_dat[l].eta[j] == 0)
					{
						alpha[l][j] = 0;
					}
					else
					{
						if (dat->sys_dat.Recover == false)
						{
							alpha[l][j] = ( dat->mspd_dat[j].eta[l] - dat->mspd_dat[l].eta[j]) *
							( fabs(par[(l+1)]) /(fabs(par[(l+1)]) + fabs(par[(j+1)])) ) +
							dat->mspd_dat[l].eta[j];
						}
						else
						{
							alpha[l][j] = ( dat->mspd_dat[j].eta[l] - dat->mspd_dat[l].eta[j]) *
							( dat->gpast_dat[l].x /(dat->gpast_dat[l].x + dat->gpast_dat[j].x) ) +
							dat->mspd_dat[l].eta[j];
						}
					}
					
					e[l][j] =  sqrt( ( fabs(dat->mspd_dat[l].eMax) + e[l][l] ) * ( fabs(dat->mspd_dat[j].eMax) + e[j][j] ) ) - (shift * alpha[l][j] );
					
					e[j][l] = e[l][j];
					
					T[l][j] = exp(-(Z * (e[l][j] - e[j][j])) / (2 * R * dat->sys_dat.T));
					T[j][l] = exp(-(Z * (e[j][l] - e[l][l])) / (2 * R * dat->sys_dat.T));
				}
			}
		}
	}
	//Forward: Gradient Estimations; Reverse: NOT USED
	else
	{
		//Loop to fill out theta array & e[j][j] & T[j][j]
		for(int j=0; j<N; j++)
		{
			for(int l=0; l<N; l++)
			{
				sum4 = sum4 + (dat->mspd_dat[l].s * dat->gpast_dat[l].x );
				if(j==l)
				{
					//Set the index for eval_po_PI
					double po[1];
					dat->sys_dat.J = j;
					if (dat->sys_dat.Recover == false)
					{
						po[0] = (dat->sys_dat.PT * dat->gpast_dat[j].y) / (dat->gpast_dat[j].x);
						dat->sys_dat.PI = PI;
						lm_status_struct status;
						lm_control_struct control = lm_control_double;
						control.printflags = 0;
						lmmin(1, po, 1, dat, eval_po_PI, &control, &status, lm_printout_std);
						dat->sys_dat.total_eval = dat->sys_dat.total_eval + status.nfev;
						if (status.info > 5) {std::cout << "\nWarning! Gradient Estimation Failed!" << std::endl; return 0.0;}
						e[j][j] = (2 * (Qst(fabs(po[0]), dat, j) + dat->gsta_dat[j].dHo[0]) ) / (Z * dat->mspd_dat[j].s);
					}
					else
					{
						mError(magpie_reverse_error); return 0.0; //Note: This return forces activity to be unity
					}
					T[j][j] = 1.0;
				}
			}
			theta[j] = (dat->mspd_dat[j].s * dat->gpast_dat[j].x ) / sum4;
			sum4 = 0;
		}
		//Loop to fill out e[i][j] & T[i][j]
		for (int j=0; j<N; j++)
		{
			for(int l=0; l<N; l++)
			{
				if(j!=l)
				{
					//Calculate the shift factor for geometric mean
					double shift = sqrt( fabs(dat->mspd_dat[l].eMax * dat->mspd_dat[j].eMax) );
					
					if (dat->mspd_dat[j].eta[l] == 0 && dat->mspd_dat[l].eta[j] == 0)
					{
						alpha[l][j] = 0;
					}
					else
					{
						alpha[l][j] = ( dat->mspd_dat[j].eta[l] - dat->mspd_dat[l].eta[j]) *
						( dat->gpast_dat[l].x /(dat->gpast_dat[l].x  + dat->gpast_dat[j].x) ) +
						dat->mspd_dat[l].eta[j];
					}
					
					e[l][j] =  sqrt( ( fabs(dat->mspd_dat[l].eMax) + e[l][l] ) * ( fabs(dat->mspd_dat[j].eMax) + e[j][j] ) ) - (shift * alpha[l][j] );
					
					e[j][l] = e[l][j];
					
					T[l][j] = exp(-(Z * (e[l][j] - e[j][j])) / (2 * R * dat->sys_dat.T));
					T[j][l] = exp(-(Z * (e[j][l] - e[l][l])) / (2 * R * dat->sys_dat.T));
				}
			}
		}
	}
	//Loop for solving other summations
	for(int j=0; j<N; j++)
	{
		sum1 = sum1 + (theta[j] * T[j][i]);
		for(int k=0; k<N; k++)
		{
			sum2 = sum2 + (theta[k] * T[k][j]);
		}
		sum3 = sum3 + ( (theta[j] * T[i][j]) / sum2 );
		sum2 = 0;
	}
	LNact = dat->mspd_dat[i].s * (1 - log(sum1) - sum3);
	return LNact;
}

//Function to estimate numerically the gradient for AST_SPD
double grad_mSPD(const double *par, const void *data, int i)
{
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	volatile double h, xph, xmh, dx, lnact_ph, lnact_mh;
	double grad=0;
	
	if (dat->sys_dat.Recover == false)
	{
		h = sqrt(DBL_EPSILON) * dat->sys_dat.PI;
		xph = dat->sys_dat.PI + h;
		xmh = dat->sys_dat.PI - h;
		lnact_ph = lnact_mSPD(NULL, dat, i, xph);
		lnact_mh = lnact_mSPD(NULL, dat, i, xmh);
	}
	else
	{
		h = sqrt(DBL_EPSILON) * fabs(par[0]);
		xph = fabs(par[0]) + h;
		xmh = fabs(par[0]) - h;
		lnact_ph = lnact_mSPD(par, dat, i, xph);
		lnact_mh = lnact_mSPD(par, dat, i, xmh);
	}
	
	dx = xph - xmh;
	if (dx == 0.0)
		grad = 0.0;
	else
		grad = (lnact_ph - lnact_mh) / dx;
	return grad;
}

//Function to evaluate the total adsorbed amount from AST_SPD numerically
double qT(const double *par, const void *data)
{
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	double qT = 0, sum1 = 0, sum2 = 0;
	double po[dat->sys_dat.N], act[dat->sys_dat.N];
	for(int i=0; i<dat->sys_dat.N; i++)
	{
		if (dat->sys_dat.Recover == false)
		{
			po[i] = (dat->sys_dat.PT * dat->gpast_dat[i].y) / (dat->gpast_dat[i].x * dat->mspd_dat[i].gama);
			sum1 = sum1 + (dat->gpast_dat[i].x / qo(po[i], dat, i));
			if (dat->sys_dat.Ideal == false)
				sum2 = sum2 + (dat->gpast_dat[i].x * grad_mSPD(NULL,dat, i));
			else
				sum2 = 0.0;
		}
		else
		{
			if (dat->sys_dat.Ideal == false)
				act[i] = exp(lnact_mSPD(par,dat,i,fabs(par[0])));
			else
				act[i] = 1.0;
			po[i] = (dat->sys_dat.PT * fabs(par[(i+1)])) / (dat->gpast_dat[i].x * act[i]);
			sum1 = sum1 + (dat->gpast_dat[i].x / qo(po[i], dat, i));
			if (dat->sys_dat.Ideal == false)
				sum2 = sum2 + (dat->gpast_dat[i].x * grad_mSPD(par,dat, i));
			else
				sum2 = 0.0;
		}
	}
	qT = pow((sum1+sum2), -1);
	return qT;
}

//Educated Guessing routine for initial values used in eval_AST_SPD
void initialGuess_mSPD(double *par, const void *data)
{
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	double sum = 0;
	double temp;
	for(int i=0; i<dat->sys_dat.N; i++)
	{
		if (dat->sys_dat.Recover == false)
		{
			par[(i+1)] = (1.0 / dat->sys_dat.N);
		}
		else
		{
			if (dat->sys_dat.Carrier == false)
				par[(i+1)] = (1.0 / dat->sys_dat.N);
			else
			{
				par[(i+1)] = (1.0 / (dat->sys_dat.N+1));
				par[(dat->sys_dat.N+1)] = (1.0 / (dat->sys_dat.N+1));
			}
		}
	}
	
	for(int i=0; i<dat->sys_dat.N; i++)
	{
		if (dat->sys_dat.Recover == false)
			temp = (dat->sys_dat.PT * dat->gpast_dat[i].y) / ((par[(i+1)]));
		else
		{
			temp = (dat->sys_dat.PT * par[(i+1)]) / (dat->gpast_dat[i].x);
			
		}
		sum = sum + PI( temp, dat, i );
	}
	par[0] = sum / dat->sys_dat.N;
}

//Function to evaluate po based on PI used to solve for activities in SPD
void eval_po_PI(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	//NOTE: sys_dat.J must be set appropriately before this function can be called
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	fvec[0] = dat->sys_dat.PI - PI(fabs(par[0]), dat, dat->sys_dat.J);
}

//Function to evaluate the reference state pressure given the adsorbed amount
void eval_po_qo(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	fvec[0] = dat->gpast_dat[dat->sys_dat.I].qo - qo(fabs(par[0]),dat,dat->sys_dat.I);
}

//Function to be used with lmmin to evaluate the value of po based on subsystem
void eval_po(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	fvec[0] = dat->gpast_dat[dat->sys_dat.J].PIo - PI(fabs(par[0]),dat,dat->sys_dat.I);
}

//Evaluation of the Beta binary parameters of the mSPD model
void eval_eta(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	/*
	 * Parameter Convention
	 * --------------------
	 * Forward: i,j -> par[0] = eta[i][j]
	 * Reverse: j,i -> par[1] = eta[j][i]
	 *
	 * ln(gama_inf[i][j]) = s[i] * (1 - ln(tau[j][i]) - tau[i][j])
	 * ln(gama_inf[j][i]) = s[j] * (1 - ln(tau[i][j]) - tau[j][i])
	 *
	 * tau,e,alpha[0][1][k] -> dat->sys_dat.I,dat->sys_dat.J, k index is for which system
	 *
	 * 0 maps to I, 1 maps to J
	 *
	 * NOTE: for gama_inf[i][j], x[i] = 0 and x[j] = 1
	 * 		 for gama_inf[j][i], x[i] = 1 and x[j] = 0
	 */
	double ejj[2], eii[2], eij[2];
	double Tij[2], Tji[2];
	double res[2];
	int i = dat->sys_dat.I;
	int j = dat->sys_dat.J;
	double shift = sqrt( fabs(dat->mspd_dat[i].eMax * dat->mspd_dat[j].eMax) );
	
	//Infinite Dilution of i will be indexed 0
	ejj[0] = ( 2 * (Qst(dat->gpast_dat[j].po[j], dat, j) - (-1*dat->gsta_dat[j].dHo[0]) ) ) / ( Z * dat->mspd_dat[j].s);
	eii[0] = ( 2 * (Qst(dat->gpast_dat[i].po[j], dat, i) - (-1*dat->gsta_dat[i].dHo[0]) ) ) / ( Z * dat->mspd_dat[i].s);
	eij[0] = sqrt( (fabs(dat->mspd_dat[i].eMax) + eii[0]) * (fabs(dat->mspd_dat[j].eMax) + ejj[0]) ) - ( fabs(par[0]) * shift);
	
	Tij[0] = exp( - (Z * (eij[0]-ejj[0])) / (2*R*dat->sys_dat.T) );
	Tji[0] = exp( - (Z * (eij[0]-eii[0])) / (2*R*dat->sys_dat.T) );
	res[0] = dat->mspd_dat[i].s * (1 - log(Tji[0]) - Tij[0]);
	
	//Infinite Dilution of j will be indexed 1
	ejj[1] = ( 2 * (Qst(dat->gpast_dat[j].po[i], dat, j) - (-1*dat->gsta_dat[j].dHo[0]) ) ) / ( Z * dat->mspd_dat[j].s);
	eii[1] = ( 2 * (Qst(dat->gpast_dat[i].po[i], dat, i) - (-1*dat->gsta_dat[i].dHo[0]) ) ) / ( Z * dat->mspd_dat[i].s);
	eij[1] = sqrt( (fabs(dat->mspd_dat[i].eMax) + eii[1]) * (fabs(dat->mspd_dat[j].eMax) + ejj[1])) - ( fabs(par[1]) * shift);
	
	Tij[1] = exp( - (Z * (eij[1]-ejj[1])) / (2*R*dat->sys_dat.T) );
	Tji[1] = exp( - (Z * (eij[1]-eii[1])) / (2*R*dat->sys_dat.T) );
	res[1] = dat->mspd_dat[j].s * (1 - log(Tij[1]) - Tji[1]);
	
	//Residual Evaluations for each infinite dilution
	fvec[0] = log(dat->gpast_dat[i].gama_inf[j]) - res[0];
	fvec[1] = log(dat->gpast_dat[j].gama_inf[i]) - res[1];
}

//Function for evaluation of AST using SPD model for activities based on shifted geo mean
void eval_GPAST(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	/*
	 * Parameter Convention:
	 * ---------------------
	 * par[0] = PI			par[(i+1)] = x[i]
	 *
	 * n_par = 1 + dat->sys_dat.N; m_dat = n_par
	 */
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	
	if (dat->sys_dat.Recover == false)
	{
		double sum1[dat->sys_dat.N], sum2 = 0;
		double act[dat->sys_dat.N];
		//Loop for number of components
		for (int i=0; i<dat->sys_dat.N; i++)
		{
			if (dat->sys_dat.Ideal == false)
				act[i] = exp( (lnact_mSPD(par, dat, i, -1)) );
			else
				act[i] = 1.0;
			sum2 = sum2 + fabs(par[(i+1)]);
			sum1[i] = 1;
			//Loop for number of parameters for ith component
			for(int n=0; n<dat->gsta_dat[i].m; n++)
			{
				double tempKo = exp( lnKo(dat->gsta_dat[i].dHo[n], dat->gsta_dat[i].dSo[n], dat->sys_dat.T) );
				sum1[i] = sum1[i] + (tempKo * pow( ( ( dat->sys_dat.PT * dat->gpast_dat[i].y ) /
													((fabs(par[(i+1)])) * act[i] * Po) ), (n+1) ) );
			}
			fvec[(i+1)] =  fabs(par[0]) - ( (dat->gsta_dat[i].qmax * log(sum1[i])) / (dat->gsta_dat[i].m) );
		}
		fvec[0] = (1.0 - (sum2) );
	}
	else
	{
		dat->sys_dat.PI = fabs(par[0]);
		double act[dat->sys_dat.N], po[dat->sys_dat.N];
		double ysum = 0;
		for (int i=0; i<dat->sys_dat.N; i++)
		{
			ysum = ysum + fabs(par[(i+1)]);
			if (dat->sys_dat.Ideal == false)
				act[i] = exp(lnact_mSPD(par,dat,i,fabs(par[0])));
			else
				act[i] = 1.0;
			po[i] = (dat->sys_dat.PT * fabs(par[(i+1)])) / (dat->gpast_dat[i].x * act[i]);
			fvec[(i+1)] = fabs(par[0]) - PI(po[i],dat,i);
		}
		if (dat->sys_dat.Carrier == true)
			ysum = ysum + fabs(par[(dat->sys_dat.N+1)]);
		fvec[0] = dat->sys_dat.qT - qT(par,dat);
		fvec[(dat->sys_dat.N+1)] = 1 - ysum;
	}
}

//Function to solve GPAST given a set of inputs that have been initialized
int MAGPIE(const void *data)
{
	//NOTE: This function assumes that all relevant data in the data structure as been initialized
	//Boundary Checking for variables must be completed prior to calling this routine
	int success = 0;
	MAGPIE_DATA *dat = (MAGPIE_DATA *) data;
	dat->sys_dat.total_eval = 0;
	lm_status_struct status;
	status.nfev = 0;
	lm_control_struct control = lm_control_double;
	control.epsilon = sqrt(DBL_EPSILON);
	control.ftol = 1e-6;
	control.gtol = 1e-6;
	control.xtol = 1e-6;
	control.stepbound = 100.0;
	control.maxcall = 1000;
	control.printflags = 0;
	
	//STEP 0: Check and correct for components whose mole fractions are zero at this sample point
	std::vector<double> qmax0;
	std::vector<int> m0;
	std::vector<std::vector<double> > dHo0;
	std::vector<std::vector<double> > dSo0;
	std::vector<double> v0;
	std::vector<double> xy0;
	std::vector<double> present0;
	int N0 = dat->sys_dat.N;
	
	//Temp objects to hold solutions
	std::vector<double> xy1;
	std::vector<double> q1;
	std::vector<double> gama1;
	
	//Copy all static info and locate index of non-existent species
	qmax0.resize(N0);
	m0.resize(N0);
	v0.resize(N0);
	present0.resize(N0);
	xy0.resize(N0);
	
	//Copy the necessary information into Temp working space
	for (int i=0; i<N0; i++)
	{
		dHo0.push_back(std::vector<double> ());
		dSo0.push_back(std::vector<double> ());
		
		qmax0[i] = dat->gsta_dat[i].qmax;
		m0[i] = dat->gsta_dat[i].m;
		v0[i] = dat->mspd_dat[i].v;
		
		if (dat->sys_dat.Recover == false)
		{
			xy0[i] = dat->gpast_dat[i].y;
			if (dat->sys_dat.PT == 0.0)
				xy0[i] = 0.0;
		}
		else
		{
			xy0[i] = dat->gpast_dat[i].x;
			if (dat->sys_dat.qT == 0.0)
				xy0[i] = 0.0;
		}
		
		if (xy0[i] == 0.0)
		{
			dat->gpast_dat[i].present = false;
			present0[i] = false;
			dat->sys_dat.N--;
		}
		else
		{
			dat->gpast_dat[i].present = true;
			present0[i] = true;
		}
		
		for (int n=0; n<m0[i]; n++)
		{
			dHo0[i].push_back(dat->gsta_dat[i].dHo[n]);
			dSo0[i].push_back(dat->gsta_dat[i].dSo[n]);
		}
	}//END Copy
	
	//STEP 0.5: Reallocation of working space
	int k = 0;
	if (N0 > dat->sys_dat.N)
	{
		//Reallocation
		if (dat->sys_dat.Output == true)
			std::cout << "\nNot all components present. Reallocating working space.\n" << std::endl;
		
		dat->gsta_dat.resize(dat->sys_dat.N);
		dat->gpast_dat.resize(dat->sys_dat.N);
		dat->mspd_dat.resize(dat->sys_dat.N);
		xy1.resize(dat->sys_dat.N);
		q1.resize(dat->sys_dat.N);
		gama1.resize(dat->sys_dat.N);
		
		for (int i=0; i<N0; i++)
		{
			if (dat->gpast_dat[i].present == true)
			{
				//Allocations
				dat->mspd_dat[k].eta.resize( dat->sys_dat.N );
				dat->gpast_dat[k].gama_inf.resize( dat->sys_dat.N );
				dat->gpast_dat[k].po.resize( dat->sys_dat.N );
				dat->gsta_dat[k].dHo.resize( m0[i]);
				dat->gsta_dat[k].dSo.resize( m0[i]);
				
				//Initializations
				dat->gsta_dat[k].qmax = qmax0[i];
				dat->gsta_dat[k].m = m0[i];
				dat->mspd_dat[k].v = v0[i];
				for (int n=0; n<dat->gsta_dat[k].m; n++)
				{
					dat->gsta_dat[k].dHo[n] = dHo0[i][n];
					dat->gsta_dat[k].dSo[n] = dSo0[i][n];
				}
				if (dat->sys_dat.Recover == false)
					dat->gpast_dat[k].y = xy0[i];
				else
					dat->gpast_dat[k].x = xy0[i];
				k++;
			}
			else
			{
				//Component is not present and should be left out of the evaluations
			}
		}
	}
	else
	{
		//Do nothing
	}//END Reallocation
	
	//Check the number of components
	if (dat->sys_dat.N > 1)
	{
		
		//Step 1: Calculate Henry's Coeffs., Ref.Capacities, Spreading Pressures, Shape Factors, and Maximum Energies
		for (int i=0; i<dat->sys_dat.N; i++)
		{
			//Both
			double tempKo = exp( lnKo(dat->gsta_dat[i].dHo[0], dat->gsta_dat[i].dSo[0], dat->sys_dat.T) );
			dat->gpast_dat[i].He = He(dat->gsta_dat[i].qmax,tempKo,dat->gsta_dat[i].m);
			dat->mspd_dat[i].s = shapeFactor( dat->mspd_dat[i].v );
			dat->mspd_dat[i].eMax = eMax(dat,i);
			
			//Forward
			if (dat->sys_dat.Recover == false)
			{
				dat->gpast_dat[i].qo = qo((dat->sys_dat.PT),dat,i);
				dat->gpast_dat[i].PIo = PI((dat->sys_dat.PT),dat,i);
			}
			//Reverse
			else
			{
				dat->gpast_dat[i].qo = qo((dat->sys_dat.PT),dat,i);
				dat->gpast_dat[i].PIo = PI((dat->sys_dat.PT),dat,i);
			}
		}
		
		//Step 2: Evaluate the Reference state pressures and infinite dilution activities
		for (int i=0; i<dat->sys_dat.N; i++)
		{
			for (int j=0; j<dat->sys_dat.N; j++)
			{
				if (i==j)
				{
					dat->gpast_dat[i].gama_inf[j] = 1.0;
					dat->gpast_dat[i].po[j] = (dat->sys_dat.PT);
				}
				else
				{
					//Must be solved for using iterative method
					dat->sys_dat.I = i; dat->sys_dat.J = j;
					double par_po[1];
					if (dat->sys_dat.Recover == false)
						par_po[0] = (dat->sys_dat.PT*dat->gpast_dat[i].y);
					else
						par_po[0] = (dat->sys_dat.PT*dat->gpast_dat[i].x);
					if (dat->sys_dat.Output == true)
						std::cout << "\nEvaluating Reference State Pressures...\n";
					lmmin(1,par_po,1,dat,eval_po,&control,&status,lm_printout_std);
					dat->sys_dat.total_eval = dat->sys_dat.total_eval + status.nfev;
					dat->sys_dat.avg_norm = dat->sys_dat.avg_norm + status.fnorm;
					if (dat->sys_dat.max_norm < fabs(status.fnorm))
						dat->sys_dat.max_norm = status.fnorm;
					if (dat->sys_dat.Output == true)
					{
						std::cout << lm_infmsg[status.info] << std::endl;
						std::cout << "E.Norm: "<< status.fnorm << std::endl;
					}
					//if (status.info > 5) {mError(simulation_fail); return status.info;}
					dat->gpast_dat[i].po[j] = fabs(par_po[0]);
					dat->gpast_dat[i].gama_inf[j] = (dat->gpast_dat[j].qo / (dat->gpast_dat[i].He * dat->gpast_dat[i].po[j]));
				}
			}
		}
		
		//Step 2.5: Check for ideality
		for (int i=0; i<dat->sys_dat.N; i++)
		{
			if(dat->gsta_dat[i].m != 1)
			{
				dat->sys_dat.Ideal = false;
				break;
			}
			else if (dat->gsta_dat[i].m == 1)
			{
				dat->sys_dat.Ideal = true;
			}
			else
			{mError(indexing_error); return -1;}
		}
		
		if (dat->sys_dat.Ideal == true)
			if (dat->sys_dat.Output == true)
				std::cout << "\nBased on single component isotherms, system is expected to behave ideally..." << std::endl;
		
		//Step 3: Evaluate the eta's using the infinite dilution activities
		for (int i=0; i<dat->sys_dat.N; i++)
		{
			for (int j=0; j<dat->sys_dat.N; j++)
			{
				if (i==j)
				{
					dat->mspd_dat[i].eta[j] = 0.0;
				}
				else if (i<j)
				{
					//Must be solved for using iterative method
					dat->sys_dat.I = i; dat->sys_dat.J = j;
					
					//Check to see if solutions will behave ideally
					if (dat->gsta_dat[i].m == 1 && dat->gsta_dat[j].m == 1)
					{
						if (dat->sys_dat.Output == true)
							std::cout << "\nBinary pair will behave ideally..." << std::endl;
						dat->mspd_dat[i].eta[j] = 0.0;
						dat->mspd_dat[j].eta[i] = 0.0;
					}
					else
					{
						double par_eta[2];
						par_eta[0] = 1.0; par_eta[1] = 1.0;
						if (dat->sys_dat.Output == true)
							std::cout << "\nEvaluating Binary Interaction Parameters...\n";
						lmmin(2,par_eta,2,dat,eval_eta,&control,&status,lm_printout_std);
						dat->sys_dat.total_eval = dat->sys_dat.total_eval + status.nfev;
						dat->sys_dat.avg_norm = dat->sys_dat.avg_norm + status.fnorm;
						if (dat->sys_dat.max_norm < fabs(status.fnorm))
							dat->sys_dat.max_norm = status.fnorm;
						if (dat->sys_dat.Output == true)
						{
							std::cout << lm_infmsg[status.info] << std::endl;
							std::cout << "E.Norm: "<< status.fnorm << std::endl;
						}
						//if (status.info > 5) {mError(simulation_fail); return status.info;}
						dat->mspd_dat[i].eta[j] = fabs(par_eta[0]);
						dat->mspd_dat[j].eta[i] = fabs(par_eta[1]);
					}
				}
				else
				{
					//No Calculation
				}
			}
		}
		
		//Step 4: Solve the resulting GPAST system of equations using all the information now gathered
		/*
		 * Parameter Convention: Forward
		 * -----------------------------
		 * PI = par[0]		x[i] = par[(i+1)]
		 * i goes from 0 to N+1, N  being the number of components
		 *
		 * Parameter Convention: Reverse
		 * -----------------------------
		 * PI = par[0]		y[i] = par[(i+1)]
		 *
		 */
		int n_par, m_dat;
		
		if (dat->sys_dat.Carrier == false)
			n_par = 1+dat->sys_dat.N;
		else
			n_par = 1+dat->sys_dat.N;
		
		if (dat->sys_dat.Carrier == true && dat->sys_dat.Recover == true)
			n_par++;
		
		if (dat->sys_dat.Recover == false)
			m_dat = 1+dat->sys_dat.N;
		else
			m_dat = 2+dat->sys_dat.N;
		
		if (dat->sys_dat.Output == true)
			std::cout << "\nAttemping the GPAST solution of the system..." << std::endl;
		double par_gpast[n_par];
		initialGuess_mSPD(par_gpast,dat);
		
		lmmin(n_par,par_gpast,m_dat,dat,eval_GPAST,&control,&status,lm_printout_std);
		dat->sys_dat.total_eval = dat->sys_dat.total_eval + status.nfev;
		dat->sys_dat.avg_norm = dat->sys_dat.avg_norm + status.fnorm;
		if (dat->sys_dat.max_norm < fabs(status.fnorm))
			dat->sys_dat.max_norm = status.fnorm;
		if (dat->sys_dat.Output == true)
		{
			std::cout << lm_infmsg[status.info] << std::endl;
			std::cout << "E.Norm: "<< status.fnorm << std::endl;
		}
		//if (status.info > 5) {mError(simulation_fail); return status.info;}
		//success = status.info;
		
		dat->sys_dat.PI = fabs(par_gpast[0]);
		xy1.resize(dat->sys_dat.N);
		gama1.resize(dat->sys_dat.N);
		q1.resize(dat->sys_dat.N);
		for (int i=0; i<dat->sys_dat.N; i++)
		{
			if (dat->sys_dat.Recover == false)
			{
				dat->gpast_dat[i].x = fabs(par_gpast[(i+1)]);
				if (dat->sys_dat.Ideal == false)
					dat->mspd_dat[i].gama = exp( lnact_mSPD(par_gpast,dat,i,-1) );
				else
					dat->mspd_dat[i].gama = 1.0;
				xy1[i] = dat->gpast_dat[i].x;
				gama1[i] = dat->mspd_dat[i].gama;
			}
			else
			{
				dat->gpast_dat[i].y = fabs(par_gpast[(i+1)]);
				if (dat->sys_dat.Ideal == false)
					dat->mspd_dat[i].gama = exp(lnact_mSPD(par_gpast,dat,i,fabs(par_gpast[0])));
				else
					dat->mspd_dat[i].gama = 1.0;
				xy1[i] = dat->gpast_dat[i].y;
				gama1[i] = dat->mspd_dat[i].gama;
			}
		}
		if (dat->sys_dat.Recover == false)
			dat->sys_dat.qT = qT(NULL,dat);
		for (int i=0; i<dat->sys_dat.N; i++)
		{
			dat->gpast_dat[i].q = dat->sys_dat.qT * dat->gpast_dat[i].x;
			q1[i] = dat->gpast_dat[i].q;
		}
		
		//Restore original data
		k = 0;
		if (N0 >= dat->sys_dat.N)
		{
			dat->sys_dat.N = N0;
			//Allocate back the original space
			dat->gsta_dat.resize(N0);
			dat->gpast_dat.resize(N0);
			dat->mspd_dat.resize(N0);
			
			for (int i=0; i<N0; i++)
			{
				//Allocations
				dat->mspd_dat[i].eta.resize( N0 );
				dat->gpast_dat[i].gama_inf.resize( N0 );
				dat->gpast_dat[i].po.resize( N0 );
				
				dat->gsta_dat[i].dHo.resize( m0[i]);
				dat->gsta_dat[i].dSo.resize( m0[i]);
				
				//Restore original information
				dat->gsta_dat[i].qmax = qmax0[i];
				dat->gsta_dat[i].m = m0[i];
				dat->mspd_dat[i].v = v0[i];
				for (int n=0; n<dat->gsta_dat[i].m; n++)
				{
					dat->gsta_dat[i].dHo[n] = dHo0[i][n];
					dat->gsta_dat[i].dSo[n] = dSo0[i][n];
				}
				
				//Restore Answers
				if (present0[i] == false)
				{
					if (dat->sys_dat.Recover == false)
					{
						dat->gpast_dat[i].x = 0.0;
						dat->gpast_dat[i].y = xy0[i];
						dat->gpast_dat[i].q = 0.0;
						dat->mspd_dat[i].gama = 1.0;
						dat->gpast_dat[i].poi = 0.0;
					}
					else
					{
						dat->gpast_dat[i].y = 0.0;
						dat->gpast_dat[i].x = xy0[i];
						dat->gpast_dat[i].q = 0.0;
						dat->mspd_dat[i].gama = 1.0;
						dat->gpast_dat[i].poi = 0.0;
					}
				}
				else
				{
					if (dat->sys_dat.Recover == false)
					{
						dat->gpast_dat[i].x = xy1[k];
						dat->gpast_dat[i].y = xy0[i];
						dat->gpast_dat[i].q = q1[k];
						dat->mspd_dat[i].gama = gama1[k];
						dat->gpast_dat[i].poi = (dat->sys_dat.PT * dat->gpast_dat[i].y) /
						(dat->gpast_dat[i].x * dat->mspd_dat[i].gama);
					}
					else
					{
						dat->gpast_dat[i].y = xy1[k];
						dat->gpast_dat[i].x = xy0[i];
						dat->gpast_dat[i].q = q1[k];
						dat->mspd_dat[i].gama = gama1[k];
						dat->gpast_dat[i].poi = (dat->sys_dat.PT * dat->gpast_dat[i].y) /
						(dat->gpast_dat[i].x * dat->mspd_dat[i].gama);
					}
					k++;
				}
			}
		}
		else
		{
			mError(simulation_fail); success = -1;
		}//END Restoration
		
	}//End Evaluation for multiple components
	
	//Perform evaluation for single component system
	else if (dat->sys_dat.N == 1)
	{
		xy1.resize(1);
		q1.resize(1);
		gama1.resize(1);
		if (dat->sys_dat.Output == true)
			std::cout << "Single Adsorbable Component! Using Standard GSTA Isotherm..." << std::endl;
		if (dat->sys_dat.Recover == false)
		{
			//Given PT, T, y
			dat->gpast_dat[0].poi = dat->sys_dat.PT * dat->gpast_dat[0].y;
			dat->gpast_dat[0].x = 1.0;
			dat->mspd_dat[0].gama = 1.0;
			dat->sys_dat.qT = qo(dat->gpast_dat[0].poi,dat,0);
			dat->gpast_dat[0].q = dat->sys_dat.qT;
			dat->gpast_dat[0].qo = dat->gpast_dat[0].q;
			dat->sys_dat.PI = PI(dat->gpast_dat[0].poi,dat,0);
			dat->sys_dat.total_eval++;
			success = 1;
			
			//Store Temp Ans
			xy1[0] = dat->gpast_dat[0].x;
			q1[0] = dat->gpast_dat[0].q;
			gama1[0] = dat->mspd_dat[0].gama;
			
		}
		else
		{
			//Given PT, T, x, qT
			dat->sys_dat.I = 0;
			dat->gpast_dat[0].qo = dat->sys_dat.qT;
			dat->gpast_dat[0].q = dat->sys_dat.qT;
			dat->gpast_dat[0].x = 1.0;
			dat->mspd_dat[0].gama = 1.0;
			
			double par_po0[1];
			par_po0[0] = dat->sys_dat.PT;
			lmmin(1,par_po0,1,dat,eval_po_qo,&control,&status,lm_printout_std);
			dat->sys_dat.total_eval = dat->sys_dat.total_eval + status.nfev;
			dat->sys_dat.avg_norm = dat->sys_dat.avg_norm + status.fnorm;
			if (dat->sys_dat.max_norm < fabs(status.fnorm))
				dat->sys_dat.max_norm = status.fnorm;
			if (dat->sys_dat.Output == true)
			{
				std::cout << lm_infmsg[status.info] << std::endl;
				std::cout << "E.Norm: "<< status.fnorm << std::endl;
			}
			
			dat->gpast_dat[0].poi = fabs(par_po0[0]);
			dat->sys_dat.PI = PI(dat->gpast_dat[0].poi,dat,0);
			dat->gpast_dat[0].y = dat->gpast_dat[0].poi / dat->sys_dat.PT;
			
			//Store Temp Ans
			xy1[0] = dat->gpast_dat[0].y;
			q1[0] = dat->gpast_dat[0].q;
			gama1[0] = dat->mspd_dat[0].gama;
			
			//if (status.info > 5) {mError(simulation_fail); return status.info;}
			//success = status.info;
		}
		
		//Restore original data
		k = 0;
		if (N0 >= dat->sys_dat.N)
		{
			dat->sys_dat.N = N0;
			//Allocate back the original space
			dat->gsta_dat.resize(N0);
			dat->gpast_dat.resize(N0);
			dat->mspd_dat.resize(N0);
			
			for (int i=0; i<N0; i++)
			{
				//Allocations
				dat->mspd_dat[i].eta.resize( N0 );
				dat->gpast_dat[i].gama_inf.resize( N0 );
				dat->gpast_dat[i].po.resize( N0 );
				
				dat->gsta_dat[i].dHo.resize( m0[i]);
				dat->gsta_dat[i].dSo.resize( m0[i]);
				
				//Restore original information
				dat->gsta_dat[i].qmax = qmax0[i];
				dat->gsta_dat[i].m = m0[i];
				dat->mspd_dat[i].v = v0[i];
				for (int n=0; n<dat->gsta_dat[i].m; n++)
				{
					dat->gsta_dat[i].dHo[n] = dHo0[i][n];
					dat->gsta_dat[i].dSo[n] = dSo0[i][n];
				}
				
				//Restore Answers
				if (present0[i] == false)
				{
					if (dat->sys_dat.Recover == false)
					{
						dat->gpast_dat[i].x = 0.0;
						dat->gpast_dat[i].y = xy0[i];
						dat->gpast_dat[i].q = 0.0;
						dat->mspd_dat[i].gama = 1.0;
						dat->gpast_dat[i].poi = 0.0;
					}
					else
					{
						dat->gpast_dat[i].y = 0.0;
						dat->gpast_dat[i].x = xy0[i];
						dat->gpast_dat[i].q = 0.0;
						dat->mspd_dat[i].gama = 1.0;
						dat->gpast_dat[i].poi = 0.0;
					}
				}
				else
				{
					if (dat->sys_dat.Recover == false)
					{
						dat->gpast_dat[i].x = xy1[k];
						dat->gpast_dat[i].y = xy0[i];
						dat->gpast_dat[i].q = q1[k];
						dat->mspd_dat[i].gama = gama1[k];
						dat->gpast_dat[i].poi = (dat->sys_dat.PT * dat->gpast_dat[i].y) /
						(dat->gpast_dat[i].x * dat->mspd_dat[i].gama);
					}
					else
					{
						dat->gpast_dat[i].y = xy1[k];
						dat->gpast_dat[i].x = xy0[i];
						dat->gpast_dat[i].q = q1[k];
						dat->mspd_dat[i].gama = gama1[k];
						dat->gpast_dat[i].poi = (dat->sys_dat.PT * dat->gpast_dat[i].y) /
						(dat->gpast_dat[i].x * dat->mspd_dat[i].gama);
					}
					k++;
				}
			}
		}
		else
		{
			mError(simulation_fail); return -1;
		}//END Restoration
	}
	
	else if (dat->sys_dat.N == 0)
	{
		//Loop for the original number of components
		dat->sys_dat.PI = 0.0;
		dat->sys_dat.qT = 0.0;
		success = 1;
		
		//Restore original data
		k = 0;
		if (N0 > dat->sys_dat.N)
		{
			dat->sys_dat.N = N0;
			//Allocate back the original space
			dat->gsta_dat.resize(N0);
			dat->gpast_dat.resize(N0);
			dat->mspd_dat.resize(N0);
			
			for (int i=0; i<N0; i++)
			{
				//Allocations
				dat->mspd_dat[i].eta.resize( N0 );
				dat->gpast_dat[i].gama_inf.resize( N0 );
				dat->gpast_dat[i].po.resize( N0 );
				dat->gsta_dat[i].dHo.resize( m0[i]);
				dat->gsta_dat[i].dSo.resize( m0[i]);
				
				//Initializations
				dat->gsta_dat[i].qmax = qmax0[i];
				dat->gsta_dat[i].m = m0[i];
				dat->mspd_dat[i].v = v0[i];
				for (int n=0; n<dat->gsta_dat[i].m; n++)
				{
					dat->gsta_dat[i].dHo[n] = dHo0[i][n];
					dat->gsta_dat[i].dSo[n] = dSo0[i][n];
				}
				
				dat->gpast_dat[i].y = 0.0;
				dat->gpast_dat[i].x = 0.0;
				
				dat->gpast_dat[i].q = 0.0;
				dat->mspd_dat[i].gama = 1.0;
				dat->gpast_dat[i].qo = 0.0;
				dat->gpast_dat[i].poi = 0.0;
			}
		}
		else
		{
			mError(simulation_fail); success = -1;
		}
	}
	//Report error if number of components is less than one
	else
	{
		mError(simulation_fail); success = -1;
	}
	
	//Clear out temporary memory
	qmax0.clear();
	m0.clear();
	v0.clear();
	dHo0.clear();
	dSo0.clear();
	xy1.clear();
	q1.clear();
	gama1.clear();
	dat->sys_dat.N = N0;
	
	return success;
}

int MAGPIE_SCENARIOS(const char *inputFileName, const char *sceneFileName)
{
	int success = 0;
	std::string inputName, sceneName;
	//Check to see if files are given
	if (inputFileName == NULL || sceneFileName == NULL)
	{
		std::cout << "Enter the name of the input file: ";
		std::cin >> inputName;
		std::cout << "Enter the name of the scenario file: ";
		std::cin >> sceneName;
		std::cout << "\n";
		
		inputFileName = inputName.c_str();
		sceneFileName = sceneName.c_str();
	}
	std::ifstream inputFile( inputFileName );
	std::ifstream sceneFile( sceneFileName );
	
	//Check to see if files exist
	if (inputFile.good()==false || sceneFile.good()==false)
	{
		mError(file_dne);
		return -1;
	}
	
	
	//Declarations
	MAGPIE_DATA dat;
	double d_read;
	int i_read;
	double time;
	int num_scene;
	int solnFlag = 0;
	FILE *sceneResults;
	
	//Initializations
	time = clock();
	dat.sys_dat.total_eval = 0;
	dat.sys_dat.avg_norm = 0;
	dat.sys_dat.max_norm = 0;
	dat.sys_dat.Recover = false;
	dat.sys_dat.Carrier = false;
	dat.sys_dat.Ideal = false;
	dat.sys_dat.Output = true;
	
	//Check to see if file exists
	if (inputFile.good()==false || sceneFile.good()==false)
	{
		mError(file_dne);
		std::cout << "Check file names, then re-run program..." << std::endl;
		return -1;
	}
	sceneResults = fopen("output/MAGPIE_Results.txt","w+");
	if (sceneResults == nullptr)
	{
		system("mkdir output");
		sceneResults = fopen("output/MAGPIE_Results.txt","w+");
	}
	
	//Read in data from Input File
	inputFile >> i_read; dat.sys_dat.N = i_read;
	if (dat.sys_dat.N < 0) {mError(invalid_components); return -1;}
	//resize working space for data structures
	dat.gsta_dat.resize( dat.sys_dat.N );
	dat.gpast_dat.resize( dat.sys_dat.N );
	dat.mspd_dat.resize( dat.sys_dat.N );
	for (int i=0; i<dat.sys_dat.N; i++)
	{
		//resize working space for parameters
		dat.mspd_dat[i].eta.resize( dat.sys_dat.N );
		dat.gpast_dat[i].gama_inf.resize( dat.sys_dat.N );
		dat.gpast_dat[i].po.resize( dat.sys_dat.N );
		inputFile >> d_read; dat.mspd_dat[i].v = d_read;
		inputFile >> d_read; dat.gsta_dat[i].qmax = d_read;
		inputFile >> i_read; dat.gsta_dat[i].m = i_read;
		dat.gsta_dat[i].dHo.resize( dat.gsta_dat[i].m );
		dat.gsta_dat[i].dSo.resize( dat.gsta_dat[i].m );
		for (int n=0; n<dat.gsta_dat[i].m; n++)
		{
			inputFile >> d_read; dat.gsta_dat[i].dHo[n] = d_read;
			inputFile >> d_read; dat.gsta_dat[i].dSo[n] = d_read;
		}
	}
	//END of Input Read
	inputFile.close();
	
	//Read in Scenario file and run all scenarios
	sceneFile >> i_read;
	if (i_read == 0)
	{
		dat.sys_dat.Recover = false;
		dat.sys_dat.Carrier = false;
	}
	else if (i_read == 1)
		dat.sys_dat.Recover = true;
	else
	{mError(invalid_boolean); return -1;}
	
	//Create a Header for the output file
	fprintf(sceneResults, "T(K)\tPT(kPa)\t");
	for (int i=0; i<dat.sys_dat.N; i++)
		fprintf(sceneResults, "y[%i]\t", (i+1));
	fprintf(sceneResults, "qT[mol/kg]\t");
	for (int i=0; i<dat.sys_dat.N; i++)
		fprintf(sceneResults, "x[%i]\t", (i+1));
	fprintf(sceneResults, "PI[mol/kg]\n");
	
	//Continue Reading in Scenario file
	sceneFile >> i_read; num_scene = i_read;
	dat.sys_dat.Par = (dat.sys_dat.N * (dat.sys_dat.N - 1));
	dat.sys_dat.Sys = dat.sys_dat.Par / 2;
	std::cout << "Total Number of Components: " << dat.sys_dat.N << std::endl;
	std::cout << "Total Number of Scenarios: " << num_scene << std::endl;
	std::cout << "Number of Interaction Parameters: " << dat.sys_dat.Par << " per Scenario" << std::endl;
	std::cout << "Total Number of Systems to Solve: " << (dat.sys_dat.Sys * num_scene) << std::endl;
	std::cout << "--------------------------------------------------\n" << std::endl;
	int sceneCount = 0;
	double total_iter = 0;
	do
	{
		sceneFile >> d_read; dat.sys_dat.PT = d_read;
		sceneFile >> d_read; dat.sys_dat.T = d_read;
		
		//Start Simulation Run
		std::cout << "Scenario #" << (sceneCount+1) << std::endl;
		std::cout << "PT(kPa): " << dat.sys_dat.PT << "\tT(K): " << dat.sys_dat.T;
		if (dat.sys_dat.Recover == true)
		{
			sceneFile >> d_read; dat.sys_dat.qT = d_read;
			std::cout << "\tqT[M/M]: " << dat.sys_dat.qT << std::endl;
			sceneFile >> i_read;
			if (i_read == 0)
				dat.sys_dat.Carrier = false;
			else if (i_read == 1)
				dat.sys_dat.Carrier = true;
			else
			{mError(invalid_boolean); return -1;}
		}
		else
			std::cout << std::endl;
		
		//Read in the Gas or Solid Mole Fractions
		double y_check = 0;
		double xsum = 0;
		for (int i=0; i<dat.sys_dat.N; i++)
		{
			sceneFile >> d_read;
			if (dat.sys_dat.Recover == false)
			{
				dat.gpast_dat[i].y = d_read;
				if (dat.gpast_dat[i].y > 1.0 || dat.gpast_dat[i].y < 0.0)
				{mError(invalid_molefraction); return -1;}
				y_check = y_check + dat.gpast_dat[i].y;
				xsum = 1.0;
			}
			else
			{
				dat.gpast_dat[i].x = d_read;
				if (dat.gpast_dat[i].x > 1.0 || dat.gpast_dat[i].x < 0.0)
				{mError(invalid_molefraction); return -1;}
				xsum = xsum + dat.gpast_dat[i].x;
			}
		}
		if (y_check > (1.0 + 1e-06))
		{mError(invalid_gas_sum); return -1;}
		if ((xsum > (1.0 + 1e-06) || xsum < (1.0 - 1e-06)) && dat.sys_dat.qT != 0.0)
		{mError(invalid_solid_sum); return -1;}
		
		//All Reads and initializations are now complete: Call MAGPIE Routine
		solnFlag = MAGPIE((void *)&dat);
		total_iter += dat.sys_dat.total_eval;
		if (solnFlag < 4) {std::cout << "\nScenario simulation successful!\n" << std::endl;}
		else {mError(scenario_fail);}
		
		//Retrieve data from MAGPIE simulation for file output
		fprintf(sceneResults,"%.6g\t",dat.sys_dat.T);
		fprintf(sceneResults,"%.6g\t",dat.sys_dat.PT);
		for (int i=0; i<dat.sys_dat.N; i++)
			fprintf(sceneResults,"%.6g\t",dat.gpast_dat[i].y);
		fprintf(sceneResults,"%.6g\t",dat.sys_dat.qT);
		for (int i=0; i<dat.sys_dat.N; i++)
			fprintf(sceneResults,"%.6g\t",dat.gpast_dat[i].x);
		fprintf(sceneResults,"%.6g\n",dat.sys_dat.PI);
		
		num_scene--;
		sceneCount++;
		std::cout << "--------------------------------------------------\n" << std::endl;
	} while(num_scene>0);
	//END of Scenario Simulations
	sceneFile.close();
	
	//END of Program
	fclose(sceneResults);
	
	//Clean Memory
	dat.gsta_dat.clear();
	dat.gpast_dat.clear();
	dat.mspd_dat.clear();
	
	//Display Stats
	time = clock() - time;
	if (dat.sys_dat.total_eval > 0)
	{
		std::cout << "\nTotal Non-Linear Evaluations: " << total_iter << std::endl;
		std::cout << "\nMaximum E. Norm Calculated: " << dat.sys_dat.max_norm << std::endl;
		std::cout << "\nAverage Euclidean Norm: " <<
		(dat.sys_dat.avg_norm / total_iter) << std::endl;
	}
	std::cout << "\nTotal Runtime: " << (time/ CLOCKS_PER_SEC) << " seconds" << std::endl;
	
	return success;
}
