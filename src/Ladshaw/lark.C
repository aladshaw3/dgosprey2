/*!
 *  \file lark.cpp lark.h
 *	\brief Linear Algebra Residual Kernels
 *  \author Austin Ladshaw
 *	\date 10/14/2014
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "lark.h"

/// \cond

//Example matrix vector product function
int matvec_ex01(const Matrix<double>& v, Matrix<double>& w, const void *data)
{
	int success = 0;
	EX01_DATA *dat = (EX01_DATA *) data;
	
	w = dat->M * v;
	
	return success;
}

//Example for preconditioning
int precon_ex01(const Matrix<double>& b, Matrix<double>& p, const void *data)
{
	int success = 0;
	EX01_DATA *dat = (EX01_DATA *) data;
	
	//Example 1.1: Symmetric Gauss-Siedel Preconditioning (Slowest, but improves convergence)
	//Matrix<double> interim;
	//p.lowerTriangularSolve(dat->M, interim.upperTriangularSolve(dat->M, b));
	
	//Example 1.2: Upper Triangular Preconditioning
	//p.upperTriangularSolve(dat->M, b);
	
	//Example 1.3: Lower Triangular Precondtioning
	//p.lowerTriangularSolve(dat->M, b);
	
	//Example 1.4: Jacobi Preconditioning (Fastest, but requires more iterations)
	p.diagonalSolve(dat->M, b);
	
	//Example 1.5: Tridiagonal Preconditioning
	//p.ladshawSolve(dat->M, b);
	
	return success;
}

//Example matrix vector product function
int matvec_ex02(const Matrix<double>& v, Matrix<double>& w, const void *data)
{
	int success = 0;
	EX02_DATA *dat = (EX02_DATA *) data;
	
	w = dat->M * v;
	
	return success;
}

//Example matrix vector product function
int matvec_ex04(const Matrix<double>& v, Matrix<double>& w, const void *data)
{
	int success = 0;
	EX04_DATA *dat = (EX04_DATA *) data;
	w = dat->M * v;
	return success;
}

//Tridiagonal preconditioning for a 3D Laplacian
int precon_ex04(const Matrix<double>& b, Matrix<double>& p, const void *data)
{
	int success = 0;
	EX04_DATA *dat = (EX04_DATA *) data;
	
	//Example 1.1: Symmetric Gauss-Siedel Preconditioning (Slowest, but improves convergence)
	Matrix<double> interim;
	p.lowerTriangularSolve(dat->M, interim.upperTriangularSolve(dat->M, b));
	
	//Example 1.2: Upper Triangular Preconditioning
	//p.upperTriangularSolve(dat->M, b);
	
	//Example 1.3: Lower Triangular Precondtioning
	//p.lowerTriangularSolve(dat->M, b);
	
	//Example 1.4: Jacobi Preconditioning (Fastest, but requires more iterations)
	//p.diagonalSolve(dat->M, b);
	
	//Example 1.5: Tridiagonal Preconditioning
	//p.ladshawSolve(dat->M, b);
	
	return success;
}

//Function to approximate the solution to x for the Picard Iteration
int evalx_ex09(const Matrix<double> &x, Matrix<double>& G, const void *data)
{
	int success = 0;
	EX09_DATA *dat = (EX09_DATA *) data;
	
	if (G.rows() != dat->N)
	{
		mError(dim_mis_match);
		return -1;
	}
	
	for (int i=0; i<dat->N; i++)
	{
		if (i == 0)
		{
			G(i,0) = 0.5*x(i+1,0) + 0.5*dat->h*dat->h*dat->k*exp(x(i,0));
		}
		else if (i==dat->N-1)
		{
			G(i,0) = 0.5*x(i-1,0) + 0.5*dat->h*dat->h*dat->k*exp(x(i,0));
		}
		else
		{
			G(i,0) = 0.5*x(i+1,0) + 0.5*x(i-1,0) + 0.5*dat->h*dat->h*dat->k*exp(x(i,0));
		}
	}
	
	return success;
}

//Function evaluation for Picard example 9
int funeval_ex09(const Matrix<double> &x, Matrix<double>& F, const void *data)
{
	int success = 0;
	EX09_DATA *dat = (EX09_DATA *) data;
	if (F.rows() != dat->N)
	{
		mError(dim_mis_match);
		return -1;
	}
	
	for (int i=0; i<dat->N; i++)
	{
		if (i == 0)
		{
			F(i,0) = -x(i+1,0) + 2*x(i,0) - dat->h*dat->h*dat->k*exp(x(i,0));
		}
		else if (i==dat->N-1)
		{
			F(i,0) = 2*x(i,0) - x(i-1,0) - dat->h*dat->h*dat->k*exp(x(i,0));
		}
		else
		{
			F(i,0) = -x(i+1,0) + 2*x(i,0) - x(i-1,0) - dat->h*dat->h*dat->k*exp(x(i,0));
		}
	}
	
	return success;
}

//Function evaluation for PJFNK method example 10
int funeval_ex10(const Matrix<double> &x, Matrix<double>& F, const void *data)
{
	int success = 0;
	EX09_DATA *dat = (EX09_DATA *) data;
	if (F.rows() != dat->N)
	{
		mError(dim_mis_match);
		return -1;
	}
	
	for (int i=0; i<dat->N; i++)
	{
		if (i == 0)
		{
			F(i,0) = -x(i+1,0) + 2*x(i,0) - dat->h*dat->h*dat->k*exp(x(i,0));
		}
		else if (i==dat->N-1)
		{
			F(i,0) = 2*x(i,0) - x(i-1,0) - dat->h*dat->h*dat->k*exp(x(i,0));
		}
		else
		{
			F(i,0) = -x(i+1,0) + 2*x(i,0) - x(i-1,0) - dat->h*dat->h*dat->k*exp(x(i,0));
		}
	}
	
	return success;
}

//Preconditioner for the inner linear iterates for PJFNK example 10
int precon_ex10(const Matrix<double> &r, Matrix<double>& p, const void *data)
{
	int success = 0;
	EX09_DATA *dat = (EX09_DATA *) data;
	p.ladshawSolve(dat->M, r);
	return success;
}

int matvec_ex15(const Matrix<double>& v, Matrix<double>& w, const void *data)
{
	int success = 0;
	EX15_DATA *dat = (EX15_DATA *) data;
	int r = dat->m;
	int r2 = dat->m*dat->m;
	int r3 = r2;
	int r4 = 0;
	//Perform action A*v and store in vector w
	if (dat->N != w.rows())
	{
		mError(matvec_mis_match);
		return -1;
	}
	for (int i=0; i<dat->N; i++)
	{
		//If statements for tridiagonal portion
		w.edit(i, 0, 6*v(i,0));
		if (i == 0)
		{
			w.edit(i, 0, w(i,0)-1*v(i+1,0));
		}
		else if (i == dat->N-1)
		{
			w.edit(i, 0, w(i,0)-1*v(i-1,0));
		}
		else if (i == r-1)
		{
			w.edit(i, 0, w(i,0)-1*v(i-1,0));
		}
		else if (i == r)
		{
			w.edit(i, 0, w(i,0)-1*v(i+1,0));
			r = r + dat->m;
		}
		else
		{
			w.edit(i, 0, w(i,0)-1*v(i-1,0));
			w.edit(i, 0, w(i,0)-1*v(i+1,0));
		}
		
		//If statements for 2nd diagonal bands
		if (i > dat->m-1)
		{
			if (i <= r3-1)
			{
				w.edit(i, 0, w(i,0)-1*v(i-dat->m,0));
			}
			else if (i > r3-1)
			{
				r4 = r4+1;
				if (r4 == dat->m-1)
				{
					r3 = r2;
					r4 = 0;
				}
			}
		}
		if (i <= dat->N-dat->m-1 && i <= r2-dat->m-1)
		{
			w.edit(i, 0, w(i,0)-1*v(i+dat->m,0));
		}
		if (i == r2-1)
		{
			r2 = r2+(dat->m*dat->m);
		}
		
		//If statements for 3rd diagonal bands
		if (i > (dat->m*dat->m)-1)
		{
			w.edit(i, 0, w(i,0)-1*v(i-(dat->m*dat->m),0));
		}
		if (i <= dat->N-(dat->m*dat->m)-1)
		{
			w.edit(i, 0, w(i,0)-1*v(i+(dat->m*dat->m),0));
		}
	}
	return success;
}

int precon_ex15(const Matrix<double>& w, Matrix<double>& p, const void *data)
{
	int success = 0;
	EX15_DATA *dat = (EX15_DATA *) data;
	if (dat->N != p.rows())
	{
		mError(matvec_mis_match);
		return -1;
	}
	for (int i=0; i<dat->N; i++)
	{
		p.edit(i, 0, w(i,0)/6.0);
	}
	return success;
}

/// \endcond

//Function to compute the updated solution given the matrix-vector arguments
int update_arnoldi_solution(Matrix<double>& x, Matrix<double>& x0, ARNOLDI_DATA *arnoldi_dat)
{
	int success = 0;
	
	//Check for wrong matrix sizes
	if (arnoldi_dat->Vk.size() == 0 || arnoldi_dat->Vk.size() != arnoldi_dat->yk.rows() || arnoldi_dat->Vk[0].rows() != x0.rows())
	{
		mError(dim_mis_match);
		return -1;
	}
	if (x.rows() != x0.rows())
	{
		x.set_size(x0.rows(), 1);
	}
	
	//perform transformations to get the search direction yk
	arnoldi_dat->Hkp1.upperHessenberg2Triangular(arnoldi_dat->e1);
	arnoldi_dat->Hkp1.rowShrink();
	arnoldi_dat->e1.rowShrink();
	arnoldi_dat->yk.upperTriangularSolve(arnoldi_dat->Hkp1, arnoldi_dat->e1);
	
	//Loop over all rows
	double sum = 0.0;
	for (int i=0; i<x.rows(); i++)
	{
		sum = 0.0;
		for (int k=0; k<arnoldi_dat->Vk.size(); k++)
		{
			if (i >= arnoldi_dat->Vk[k].rows())
				sum = sum + 0.0;
			else
				sum = sum + (arnoldi_dat->Vk[k](i,0) * arnoldi_dat->yk(k,0));
		}
		x.edit(i, 0, sum + x0(i,0));
	}
	
	return success;
}

//Function to construct the orthonormal basis with optional preconditioning
int arnoldi( int (*matvec) (const Matrix<double>& v, Matrix<double> &w, const void *data),
			int (*precon) (const Matrix<double>& b, Matrix<double> &p, const void *data),
			Matrix<double> &r0, ARNOLDI_DATA *arnoldi_dat, const void *matvec_data,
			const void *precon_data )
{
	int success = 0;
	double h = 0;
	
	//Initialize the krylov data
	if (arnoldi_dat->k < 2 || arnoldi_dat->k > r0.rows())
		arnoldi_dat->k = r0.rows();
	if (r0.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	if ( (*matvec) == NULL)
	{
		success = -1;
		mError(nullptr_func);
		return success;
	}
	if ( arnoldi_dat->Vk.size() != arnoldi_dat->k )
	{
		arnoldi_dat->Vk.resize(arnoldi_dat->k);
	}
	if ( arnoldi_dat->w.rows() != r0.rows() )
	{
		arnoldi_dat->w.set_size(r0.rows(), 1);
	}
	if ( arnoldi_dat->v.rows() != r0.rows() )
	{
		arnoldi_dat->v.set_size(r0.rows(), 1);
	}
	if ( arnoldi_dat->sum.rows() != r0.rows() )
	{
		arnoldi_dat->sum.set_size(r0.rows(), 1);
	}
	if ( arnoldi_dat->Hkp1.rows() != arnoldi_dat->k+1 )
	{
		arnoldi_dat->Hkp1.set_size(arnoldi_dat->k+1, arnoldi_dat->k);
	}
	if ( arnoldi_dat->e1.rows() != arnoldi_dat->k+1)
	{
		arnoldi_dat->e1.set_size(arnoldi_dat->k+1, 1);
	}
	if ( arnoldi_dat->yk.rows() != arnoldi_dat->k )
	{
		arnoldi_dat->yk.set_size(arnoldi_dat->k, 1);
	}
	
	//Check for and apply preconditioning
	if ( (*precon) != NULL)
	{
		success = (*precon) (r0, arnoldi_dat->v, precon_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
	}
	else
	{
		arnoldi_dat->v = r0;
	}
	
	//Initialize data before loop
	arnoldi_dat->beta = arnoldi_dat->v.norm();
	arnoldi_dat->e1.edit(0, 0, arnoldi_dat->beta);
	if (arnoldi_dat->Vk[0].rows() != r0.rows())
		arnoldi_dat->Vk[0].set_size(r0.rows(), 1);
	for (int n=0; n<r0.rows(); n++)
		arnoldi_dat->Vk[0].edit(n, 0, arnoldi_dat->v(n,0)/arnoldi_dat->beta);
	
	//Begin looping for the size of the subspace
	arnoldi_dat->iter = 0;
	for (int j=0; j<arnoldi_dat->k; j++)
	{
		if (arnoldi_dat->Output == true)
		{
			std::cout << "Arnoldi vector " << j+1 << " being built..." << std::endl;
		}
		arnoldi_dat->sum.zeros();
		success = (*matvec) (arnoldi_dat->Vk[j], arnoldi_dat->w, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
		if ( (*precon) != NULL)
		{
			arnoldi_dat->v = arnoldi_dat->w;
			success = (*precon) (arnoldi_dat->v, arnoldi_dat->w, precon_data);
			if (success != 0)
			{
				mError(simulation_fail);
				return success;
			}
		}
		else {/* No change to w */}
		
		//Inner loop for formation of Hessenberg matrix
		for (int i=0; i<=j; i++)
		{
			h = arnoldi_dat->w.inner_product(arnoldi_dat->Vk[i]);
			arnoldi_dat->Hkp1.edit(i, j, h);
			for (int n=0; n<r0.rows(); n++)
				arnoldi_dat->sum.edit(n, 0, arnoldi_dat->sum(n,0) + (arnoldi_dat->Vk[i](n,0)*h));
			
		} //END inner loop
		
		arnoldi_dat->w = arnoldi_dat->w - arnoldi_dat->sum;
		arnoldi_dat->hp1 = arnoldi_dat->w.norm();
		arnoldi_dat->Hkp1.edit(j+1, j, arnoldi_dat->hp1);
		if (j<arnoldi_dat->k-1 && arnoldi_dat->Vk[j+1].rows() != r0.rows())
			arnoldi_dat->Vk[j+1].set_size(r0.rows(), 1);
		
		if (arnoldi_dat->hp1 == 0.0)
		{
			if (j==0 && arnoldi_dat->Hkp1(j,j)==0.0)
				arnoldi_dat->Hkp1.edit(j, j, 1.0);
			for (int i=j+1; i<arnoldi_dat->k; i++)
			{
				arnoldi_dat->Hkp1.edit(i, i, 1.0);
				if (j<arnoldi_dat->k-1 && arnoldi_dat->Vk[j+1].rows() != r0.rows())
					arnoldi_dat->Vk[j+1].set_size(r0.rows(), 1);
			}
			break;
		}
		else
		{
			if (j<arnoldi_dat->k-1)
			{
				for (int n=0; n<r0.rows(); n++)
					arnoldi_dat->Vk[j+1].edit(n, 0, arnoldi_dat->w(n,0)/arnoldi_dat->hp1);
			}
			else
				arnoldi_dat->v = arnoldi_dat->w/arnoldi_dat->hp1;
		}
		
		arnoldi_dat->iter++;
		
	} //END Subspace loop
	
	//Solve the resulting least squares problem
	if (arnoldi_dat->Output == true)
	{
		std::cout << "Krylov subspace construction complete!\n" << std::endl;
	}
	
	return success;
}

//Function to perform the Restarted GMRES algorithm for iteratively solving a linear system
int gmresLeftPreconditioned( int (*matvec) (const Matrix<double>& v, Matrix<double> &w, const void *data),
							int (*precon) (const Matrix<double>& b, Matrix<double> &P, const void *data),
							Matrix<double> &b, GMRESLP_DATA *gmreslp_dat, const void *matvec_data,
							const void *precon_data )
{
	int success = 0;
	double res_old;
	
	//Check arguments
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		success = -1;
		return success;
	}
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	
	//Initialize Krylov Data and GMRES Data
	if (gmreslp_dat->restart < 2 || gmreslp_dat->restart > b.rows())
	{
		gmreslp_dat->restart = std::min(20,b.rows());
	}
	if (gmreslp_dat->maxit < 1 || gmreslp_dat->maxit > b.rows())
	{
		gmreslp_dat->maxit = std::min(1000,b.rows());
	}
	if (gmreslp_dat->x.rows() != b.rows())
	{
		gmreslp_dat->x.set_size(b.rows(), 1);
	}
	if (gmreslp_dat->tol_rel >= 1)
	{
		gmreslp_dat->tol_rel = 1e-6;
	}
	if (gmreslp_dat->tol_abs >= 1)
	{
		gmreslp_dat->tol_abs = 1e-6;
	}
	if (gmreslp_dat->tol_rel < MIN_TOL)
	{
		gmreslp_dat->tol_rel = MIN_TOL;
	}
	if (gmreslp_dat->tol_abs < MIN_TOL)
	{
		gmreslp_dat->tol_abs = MIN_TOL;
	}
	gmreslp_dat->arnoldi_dat.k = gmreslp_dat->restart;
	gmreslp_dat->arnoldi_dat.Output = gmreslp_dat->Output;
	
	//Check to see if the data has been properly initialized
	if (gmreslp_dat->arnoldi_dat.Vk.size() != gmreslp_dat->arnoldi_dat.k)
	{
		gmreslp_dat->arnoldi_dat.Vk.resize(gmreslp_dat->arnoldi_dat.k);
	}
	if ( gmreslp_dat->arnoldi_dat.w.rows() != b.rows() )
	{
		gmreslp_dat->arnoldi_dat.w.set_size(b.rows(), 1);
	}
	if ( gmreslp_dat->arnoldi_dat.v.rows() != b.rows() )
	{
		gmreslp_dat->arnoldi_dat.v.set_size(b.rows(), 1);
	}
	if ( gmreslp_dat->arnoldi_dat.sum.rows() != b.rows() )
	{
		gmreslp_dat->arnoldi_dat.sum.set_size(b.rows(), 1);
	}
	if ( gmreslp_dat->arnoldi_dat.Hkp1.rows() != gmreslp_dat->arnoldi_dat.k+1 )
	{
		gmreslp_dat->arnoldi_dat.Hkp1.set_size(gmreslp_dat->arnoldi_dat.k+1, gmreslp_dat->arnoldi_dat.k);
	}
	if ( gmreslp_dat->arnoldi_dat.e1.rows() != gmreslp_dat->arnoldi_dat.k+1)
	{
		gmreslp_dat->arnoldi_dat.e1.set_size(gmreslp_dat->arnoldi_dat.k+1, 1);
	}
	if ( gmreslp_dat->arnoldi_dat.yk.rows() != gmreslp_dat->arnoldi_dat.k )
	{
		gmreslp_dat->arnoldi_dat.yk.set_size(gmreslp_dat->arnoldi_dat.k, 1);
	}
	gmreslp_dat->iter = 0;
	gmreslp_dat->steps = 0;
	
	//Form first matrix vector product
	success = (*matvec) (gmreslp_dat->x, gmreslp_dat->arnoldi_dat.w, matvec_data);
	if (success != 0)
	{
		mError(simulation_fail);
		return success;
	}
	gmreslp_dat->r = b - gmreslp_dat->arnoldi_dat.w;
	gmreslp_dat->relres_base = gmreslp_dat->r.norm();
	gmreslp_dat->res = gmreslp_dat->relres_base;
	gmreslp_dat->bestres = gmreslp_dat->res;
	gmreslp_dat->bestx = gmreslp_dat->x;
	gmreslp_dat->relres = 1.0;
	res_old = gmreslp_dat->res;
	if (gmreslp_dat->Output == true)
		std::cout << "\nRelRes[" << 0 << "] =\t" << gmreslp_dat->relres << std::endl;
	
	//Check for imediate convergence
	if (gmreslp_dat->res < gmreslp_dat->tol_abs)
	{
		if (gmreslp_dat->Output == true)
		{
			std::cout << "\nSolution converged after 0 GMRES iterations witin residual tolerance!" << std::endl;
			std::cout << "AbsRes[0] =\t" << gmreslp_dat->res << std::endl << std::endl;
		}
		return success;
	}
	
	//If no solution found imediately, begin looping to search for solution
	for (int m=0; m<gmreslp_dat->maxit; m++)
	{
		//Call the Krylov Function
		success = arnoldi(matvec,precon,gmreslp_dat->r,&gmreslp_dat->arnoldi_dat,matvec_data,precon_data);
		if (success != 0)
		{
			mError(simulation_fail);
			std::cout << "Vk = \n";
			for (int n=0; n<b.rows(); n++)
			{
				for (int k=0; k<gmreslp_dat->arnoldi_dat.k; k++)
				{
					std::cout << "\t" << gmreslp_dat->arnoldi_dat.Vk[k](n,0);
				}
				std::cout << std::endl;
			}
			gmreslp_dat->arnoldi_dat.Hkp1.Display("H");
			gmreslp_dat->arnoldi_dat.yk.Display("yk");
			gmreslp_dat->x.Display("x0");
			success = update_arnoldi_solution(gmreslp_dat->x, gmreslp_dat->x, &gmreslp_dat->arnoldi_dat);
			gmreslp_dat->x.Display("x");
			success = (*matvec) (gmreslp_dat->x, gmreslp_dat->arnoldi_dat.w, matvec_data);
			gmreslp_dat->r = b - gmreslp_dat->arnoldi_dat.w;
			gmreslp_dat->res = gmreslp_dat->r.norm();
			std::cout << "RelRes at stop = \t" << gmreslp_dat->relres << std::endl;
			std::cout << "AbsRes at stop = \t" << gmreslp_dat->res << std::endl << std::endl;
			success = -1;
			return success;
		}
		
		//From the new solution and residual vector
		success = update_arnoldi_solution(gmreslp_dat->x, gmreslp_dat->x, &gmreslp_dat->arnoldi_dat);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
		success = (*matvec) (gmreslp_dat->x, gmreslp_dat->arnoldi_dat.w, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
		gmreslp_dat->r = b - gmreslp_dat->arnoldi_dat.w;
		gmreslp_dat->res = gmreslp_dat->r.norm();
		gmreslp_dat->relres = gmreslp_dat->res / gmreslp_dat->relres_base;
		gmreslp_dat->iter++;
		gmreslp_dat->steps = gmreslp_dat->steps + gmreslp_dat->iter + gmreslp_dat->arnoldi_dat.iter;
		if (gmreslp_dat->Output == true)
			std::cout << "RelRes[" << m+1 << "] =\t" << gmreslp_dat->relres << std::endl;
		
		//Check residual for convergence
		if (gmreslp_dat->res < gmreslp_dat->tol_abs)
		{
			if (gmreslp_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << m+1 << " GMRES restart(s) witin residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << m+1 << "] =\t" << gmreslp_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else if ( gmreslp_dat->relres < gmreslp_dat->tol_rel)
		{
			if (gmreslp_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << m+1 << " GMRES restart(s) witin relative residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << m+1 << "] =\t" << gmreslp_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else
		{
			if (gmreslp_dat->res < gmreslp_dat->bestres)
			{
				gmreslp_dat->bestres = gmreslp_dat->res;
				gmreslp_dat->bestx = gmreslp_dat->x;
			}
			res_old = gmreslp_dat->res;
		}
		
		//Check for failures
		if (isnan(gmreslp_dat->res) || isinf(gmreslp_dat->res))
		{
			if (gmreslp_dat->Output == true)
				std::cout << "\nGMRES broke down at the " << m+1 << " restart!" << std::endl;
			success = -1;
			break;
		}
	}
	
	//Post Loop messages
	if (gmreslp_dat->iter >= gmreslp_dat->maxit)
	{
		gmreslp_dat->res = gmreslp_dat->bestres;
		gmreslp_dat->x = gmreslp_dat->bestx;
		if (gmreslp_dat->Output == true)
		{
			std::cout << "\nGMRES did not find a solution within " << gmreslp_dat->maxit << " GMRES restarts!" << std::endl;
			std::cout << "Returning the best found solution..." << std::endl;
			std::cout << "Best absolute residual =\t" << gmreslp_dat->res << std::endl << std::endl;
		}
		success = 0;
	}
	if (success == -1)
	{
		gmreslp_dat->res = gmreslp_dat->bestres;
		gmreslp_dat->x = gmreslp_dat->bestx;
		if (gmreslp_dat->Output == true)
		{
			std::cout << "\nGMRES did not find a solution within " << gmreslp_dat->iter << " GMRES restarts!" << std::endl;
			std::cout << "Returning the best found solution..." << std::endl;
			std::cout << "Best absolute residual =\t" << gmreslp_dat->res << std::endl << std::endl;
		}
		success = 0;
	}
	return success;
}

//Function to perform the Unrestarted GMRES for directly solving a linear system using krylov function
int fom( int (*matvec) (const Matrix<double>& v, Matrix<double> &w, const void *data),
		int (*precon) (const Matrix<double>& b, Matrix<double> &P, const void *data),
		Matrix<double> &b, GMRESLP_DATA *gmreslp_dat, const void *matvec_data,
		const void *precon_data )
{
	/*
	 NOTE: This implementation is the Restarted GMRES without any restarting. Thus, this
	 method is gaurenteed to converge to the exact (or machine precision) solution
	 to the linear system as long as (i) the linear system is non-singular and (ii)
	 there is enough computer memory available to construct the full space.
	 
	 If you find that the GMRES implementation with restarting is unable to solve
	 the linear system, then try using this method. If this method cannot solve the
	 system, then the system is too large or is singular.
	 
	 This is named FOM for Full Othogonalization Method. However, GMRES and FOM are actually
	 different algorithms all together. FOM and GMRES are the same only if the full orthogonal
	 subspace is constructed, as is the case here.
	 */
	int success = 0;
	gmreslp_dat->restart = b.rows();
	gmreslp_dat->maxit = 1;
	success = gmresLeftPreconditioned(matvec,precon,b,gmreslp_dat,matvec_data,precon_data);
	if (success != 0) {mError(simulation_fail); return success;}
	return success;
}

//Implementation of the Compact GMRES algorithm
int gmresRightPreconditioned( int (*matvec) (const Matrix<double>& v, Matrix<double> &w, const void *data),
							 int (*precon) (const Matrix<double>& b, Matrix<double> &p, const void *data),
							 Matrix<double> &b, GMRESRP_DATA *gmresrp_dat, const void *matvec_data,
							 const void *precon_data )
{
	int success = 0;
	double res_old, beta, beta1, beta0 = 0.0, eta, h;
	
	//Check arguments
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		success = -1;
		return success;
	}
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	
	//Initialize Data
	if (gmresrp_dat->restart <= 1)
	{
		gmresrp_dat->restart = std::min(b.rows(),20);
	}
	else if (gmresrp_dat->restart > b.rows() && b.rows() <= 1000)
	{
		gmresrp_dat->restart = b.rows();
	}
	if (gmresrp_dat->restart == b.rows())
	{
		gmresrp_dat->maxit = 1;
	}
	else if (gmresrp_dat->maxit <= 0 || gmresrp_dat->maxit > b.rows())
	{
		gmresrp_dat->maxit = std::min(b.rows(),1000);
	}
	if (gmresrp_dat->tol_rel >= 1)
	{
		gmresrp_dat->tol_rel = 1e-6;
	}
	if (gmresrp_dat->tol_abs >= 1)
	{
		gmresrp_dat->tol_abs = 1e-6;
	}
	if (gmresrp_dat->tol_rel < MIN_TOL)
	{
		gmresrp_dat->tol_rel = MIN_TOL;
	}
	if (gmresrp_dat->tol_abs < MIN_TOL)
	{
		gmresrp_dat->tol_abs = MIN_TOL;
	}
	if (gmresrp_dat->x.rows() != b.rows())
	{
		gmresrp_dat->x.set_size(b.rows(), 1);
	}
	if (gmresrp_dat->w.rows() != b.rows())
	{
		gmresrp_dat->w.set_size(b.rows(), 1);
	}
	if (gmresrp_dat->v.rows() != b.rows())
	{
		gmresrp_dat->v.set_size(b.rows(), 1);
	}
	if (gmresrp_dat->sum.rows() != b.rows())
	{
		gmresrp_dat->sum.set_size(b.rows(), 1);
	}
	gmresrp_dat->iter_inner = 0;
	gmresrp_dat->iter_outer = 0;
	
	//Form first matrix vector product
	success = (*matvec) (gmresrp_dat->x, gmresrp_dat->w, matvec_data);
	if (success != 0)
	{
		mError(simulation_fail);
		return success;
	}
	gmresrp_dat->r = b - gmresrp_dat->w;
	gmresrp_dat->relres_base = gmresrp_dat->r.norm();
	gmresrp_dat->res = gmresrp_dat->relres_base;
	gmresrp_dat->bestres = gmresrp_dat->res;
	gmresrp_dat->bestx = gmresrp_dat->x;
	gmresrp_dat->relres = 1.0;
	res_old = gmresrp_dat->res;
	if (gmresrp_dat->Output == true)
		std::cout << "\nRelRes[" << 0 << "] =\t" << gmresrp_dat->relres << "\t---Start[0]---" << std::endl;
	
	//Check for imediate convergence
	if (gmresrp_dat->res < gmresrp_dat->tol_abs)
	{
		if (gmresrp_dat->Output == true)
		{
			std::cout << "\nSolution converged after 0 GMRES iterations witin residual tolerance!" << std::endl;
			std::cout << "AbsRes[0] =\t" << gmresrp_dat->res << std::endl << std::endl;
		}
		return success;
	}
	
	int step = 0;
	//Begin iterative procedure
	for (int m=0; m<gmresrp_dat->maxit; m++)
	{
		//Prepare the norms at restart
		beta = gmresrp_dat->r.norm();
		if (m == 0)
		{
			beta0 = beta;
			eta = ( beta0==0.0 ? 0.0 : 1.0 );
		}
		else
		{
			eta = beta/beta0;
		}
		beta1 = beta;
		
		if (gmresrp_dat->e0.size() < 1)
			gmresrp_dat->e0.push_back(beta);
		else
			gmresrp_dat->e0[0] = beta;
		
		if (m == 0)
		{
			std::vector<double> temp;
			temp.push_back(0.0);
			if (gmresrp_dat->H.size() < 1)
				gmresrp_dat->H.push_back(temp);
			else
				gmresrp_dat->H[0] = temp;
		}
		
		//Normalize and prepare vector space
		if (gmresrp_dat->Vk.size() < 1)
			gmresrp_dat->Vk.push_back(gmresrp_dat->r * (1.0/beta));
		else
			gmresrp_dat->Vk[0] = gmresrp_dat->r * (1.0/beta);
		
		//Begin inner column loops
		gmresrp_dat->iter_inner = 0;
		for (int j=0; j<gmresrp_dat->restart; j++)
		{
			gmresrp_dat->sum.zeros();
			
			//Orthogonalization and Right Preconditioning
			if ( (*precon) != NULL)
			{
				success = (*precon) (gmresrp_dat->Vk[j], gmresrp_dat->v, precon_data);
				if (success != 0)
				{
					mError(simulation_fail);
					return success;
				}
			}
			else
			{
				gmresrp_dat->v = gmresrp_dat->Vk[j];
			}
			//Store result in growing preconditioned vector space
			if (gmresrp_dat->Zk.size() < j+1)
				gmresrp_dat->Zk.push_back(gmresrp_dat->v);
			else
				gmresrp_dat->Zk[j] = gmresrp_dat->v;
			success = (*matvec) (gmresrp_dat->v, gmresrp_dat->w, matvec_data);
			if (success != 0)
			{
				mError(simulation_fail);
				return success;
			}
			
			//Row loops for formation of Hessenberg matrix
			for (int i=0; i<=j; i++)
			{
				h = gmresrp_dat->w.inner_product(gmresrp_dat->Vk[i]);
				
				if (gmresrp_dat->H[j].size() < i+1)
					gmresrp_dat->H[j].push_back(h);
				else
					gmresrp_dat->H[j][i] = h;
				
				for (int n=0; n<b.rows(); n++)
					gmresrp_dat->sum.edit(n, 0, gmresrp_dat->sum(n,0) + (gmresrp_dat->Vk[i](n,0)*h));
				
			} //END Hessenberg loop
			
			//Form the j+1 row of the j column of H
			gmresrp_dat->w = gmresrp_dat->w - gmresrp_dat->sum;
			h = gmresrp_dat->w.norm();
			if (gmresrp_dat->H[j].size() < j+2)
				gmresrp_dat->H[j].push_back(h);
			else
				gmresrp_dat->H[j][j+1] = h;
			
			if (gmresrp_dat->e0.size() < j+2)
				gmresrp_dat->e0.push_back(0.0);
			else
				gmresrp_dat->e0[j+1] = 0.0;
			
			//Copy constructors to retain the original data during rotations
			gmresrp_dat->e0_bar = gmresrp_dat->e0;
			gmresrp_dat->H_bar = gmresrp_dat->H;
			
			//Start the Givens Rotations of e0 and H and store in e0_bar and H_bar
			double s,c;
			double value_i, value_ip1;
			//Row Loop
			for (int i=0; i<gmresrp_dat->H[j].size()-1; i++)
			{
				s = gmresrp_dat->H_bar[i][i+1] / sqrt( pow(gmresrp_dat->H_bar[i][i],2.0) + pow(gmresrp_dat->H_bar[i][i+1],2.0) );
				c = gmresrp_dat->H_bar[i][i] / sqrt( pow(gmresrp_dat->H_bar[i][i],2.0) + pow(gmresrp_dat->H_bar[i][i+1],2.0) );
				
				if (isnan(s) || isnan(c) || isinf(s) || isinf(c))
				{
					mError(singular_matrix);
					return -1;
				}
				
				value_i = ((gmresrp_dat->e0_bar[i]*c) + (gmresrp_dat->e0_bar[i+1]*s));
				value_ip1 = (-(gmresrp_dat->e0_bar[i]*s) + (gmresrp_dat->e0_bar[i+1]*c));
				
				if (gmresrp_dat->e0_bar.size() < i+1)
					gmresrp_dat->e0_bar.push_back(value_i);
				else
					gmresrp_dat->e0_bar[i] = value_i;
				if (gmresrp_dat->e0_bar.size() < i+2)
					gmresrp_dat->e0_bar.push_back(value_ip1);
				else
					gmresrp_dat->e0_bar[i+1] = value_ip1;
				
				//Column Loop
				for (int k=i; k<gmresrp_dat->H.size(); k++)
				{
					value_i = ((gmresrp_dat->H_bar[k][i]*c) + (gmresrp_dat->H_bar[k][i+1]*s));
					value_ip1 = (-(gmresrp_dat->H_bar[k][i]*s) + (gmresrp_dat->H_bar[k][i+1]*c));
					gmresrp_dat->H_bar[k][i] = value_i;
					gmresrp_dat->H_bar[k][i+1] = value_ip1;
				}
				
			} //END Givens Rotations
			
			beta1 = fabs(gmresrp_dat->e0_bar[j+1]);
			eta = beta1/beta0;
			gmresrp_dat->relres = eta;
			if (gmresrp_dat->Output == true)
				std::cout << "RelRes[" << step+m+1 << "] =\t" << eta << std::endl;
			
			//Stop the routine if subspace finishes or if within tolerance
			if (h == 0.0  || eta <= gmresrp_dat->tol_rel )
			{
				step++;
				gmresrp_dat->iter_inner++;
				break;
			}
			//Form the j+1 orthonormal vector
			else
			{
				if (j < gmresrp_dat->restart-1)
				{
					if (gmresrp_dat->Vk.size() < j+2)
						gmresrp_dat->Vk.push_back(gmresrp_dat->w * (1.0/h));
					else
						gmresrp_dat->Vk[j+1] = gmresrp_dat->w * (1.0/h);
					
					if (gmresrp_dat->H.size() < j+2)
						gmresrp_dat->H.push_back(gmresrp_dat->H[j]);
					else
						gmresrp_dat->H[j+1] = gmresrp_dat->H[j];
					
				}
			}
			
			step++;
			gmresrp_dat->iter_inner++;
			
		}//END Inner Loop
		
		gmresrp_dat->iter_total = gmresrp_dat->iter_inner + gmresrp_dat->iter_total;
		
		//Solve Least Squares Problem
		if (gmresrp_dat->y.size() < gmresrp_dat->iter_inner)
			gmresrp_dat->y.resize(gmresrp_dat->iter_inner);
		double sum;
		for (int i=gmresrp_dat->iter_inner-1; i>=0; i--)
		{
			sum = 0.0;
			for (int j=gmresrp_dat->iter_inner-1; j>i; j--)
			{
				sum = sum + gmresrp_dat->H_bar[j][i] * gmresrp_dat->y[j];
			}
			if (gmresrp_dat->H_bar[i][i] == 0.0)
			{
				//Check for special case of exact solution when applying preconditioner
				if (i == 0 && gmresrp_dat->iter_inner == 1)
				{
					gmresrp_dat->y[i] = 1.0;
				}
				else
				{
					mError(singular_matrix)
					return -1;
				}
			}
			else
			{
				gmresrp_dat->y[i] = (gmresrp_dat->e0_bar[i] - sum) / gmresrp_dat->H_bar[i][i];
			}
		}
		
		//Form the new solution vector (x = x0 + Zk*y)
		for (int n=0; n<b.rows(); n++)
		{
			sum = 0.0;
			for (int i=0; i<gmresrp_dat->iter_inner; i++)
			{
				sum = sum + (gmresrp_dat->Zk[i](n,0) * gmresrp_dat->y[i]);
			}
			gmresrp_dat->w.edit(n, 0, sum);
		}
		gmresrp_dat->x = gmresrp_dat->x + gmresrp_dat->w;
		
		//Check the solution
		success = (*matvec) (gmresrp_dat->x, gmresrp_dat->w, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
		gmresrp_dat->r = b - gmresrp_dat->w;
		gmresrp_dat->res = gmresrp_dat->r.norm();
		gmresrp_dat->relres = gmresrp_dat->res / gmresrp_dat->relres_base;
		gmresrp_dat->iter_outer++;
		gmresrp_dat->iter_total = step+m+1;
		if (gmresrp_dat->Output == true)
			std::cout << "RelRes[" << step+m+1 << "] =\t" << gmresrp_dat->relres << "\t---Restart[" << m+1 << "]---" << std::endl;
		
		//Check residual for convergence
		if (gmresrp_dat->res < gmresrp_dat->tol_abs)
		{
			if (gmresrp_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << step+m+1 << " GMRES iteration(s) witin residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << step+m+1 << "] =\t" << gmresrp_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else if ( gmresrp_dat->relres < gmresrp_dat->tol_rel)
		{
			if (gmresrp_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << step+m+1 << " GMRES iteration(s) witin relative residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << step+m+1 << "] =\t" << gmresrp_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else
		{
			if (gmresrp_dat->res < gmresrp_dat->bestres)
			{
				gmresrp_dat->bestres = gmresrp_dat->res;
				gmresrp_dat->bestx = gmresrp_dat->x;
			}
			res_old = gmresrp_dat->res;
		}
		
		//Check for failures
		if (isnan(gmresrp_dat->res) || isinf(gmresrp_dat->res))
		{
			if (gmresrp_dat->Output == true)
				std::cout << "\nGMRES broke down at the " << step+m+1 << " iterate!" << std::endl;
			success = -1;
			break;
		}
		
	}//END Outer Loop
	
	//Post Loop messages
	if (gmresrp_dat->iter_outer >= gmresrp_dat->maxit)
	{
		gmresrp_dat->res = gmresrp_dat->bestres;
		gmresrp_dat->x = gmresrp_dat->bestx;
		if (gmresrp_dat->Output == true)
		{
			std::cout << "\nGMRES did not find a solution within " << gmresrp_dat->maxit+step << " GMRES iterations!" << std::endl;
			std::cout << "Returning the best found solution..." << std::endl;
			std::cout << "Best absolute residual =\t" << gmresrp_dat->res << std::endl << std::endl;
		}
		success = 0;
	}
	if (success == -1)
	{
		gmresrp_dat->res = gmresrp_dat->bestres;
		gmresrp_dat->x = gmresrp_dat->bestx;
		if (gmresrp_dat->Output == true)
		{
			std::cout << "\nGMRES did not find a solution within " << gmresrp_dat->iter_outer+step << " GMRES iterations!" << std::endl;
			std::cout << "Returning the best found solution..." << std::endl;
			std::cout << "Best absolute residual =\t" << gmresrp_dat->res << std::endl << std::endl;
		}
		success = 0;
	}
	
	return success;
}

//Implementation of the PCG algorithm for iteratively solving symmetric linear systems
int pcg( int (*matvec) (const Matrix<double>& p, Matrix<double> &Ap, const void *data),
		int (*precon) (const Matrix<double>& r, Matrix<double> &z, const void *data),
		Matrix<double> &b, PCG_DATA *pcg_dat, const void *matvec_data,
		const void *precon_data )
{
	int success = 0;
	double res_old;
	
	//Check input arguments
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		success = -1;
		return success;
	}
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	
	//Setup the working space and start the method
	if (pcg_dat->maxit <= 0 || pcg_dat->maxit > b.rows())
	{
		pcg_dat->maxit = std::min(1000,b.rows());
	}
	if (pcg_dat->tol_rel >= 1)
	{
		pcg_dat->tol_rel = 1e-6;
	}
	if (pcg_dat->tol_abs >= 1)
	{
		pcg_dat->tol_abs = 1e-6;
	}
	if (pcg_dat->tol_rel < MIN_TOL)
	{
		pcg_dat->tol_rel = MIN_TOL;
	}
	if (pcg_dat->tol_abs < MIN_TOL)
	{
		pcg_dat->tol_abs = MIN_TOL;
	}
	if (pcg_dat->x.rows() != b.rows())
	{
		pcg_dat->x.set_size(b.rows(), 1);
	}
	if (pcg_dat->r.rows() != b.rows())
	{
		pcg_dat->r.set_size(b.rows(), 1);
	}
	if (pcg_dat->r_old.rows() != b.rows())
	{
		pcg_dat->r_old.set_size(b.rows(), 1);
	}
	if (pcg_dat->z.rows() != b.rows())
	{
		pcg_dat->z.set_size(b.rows(), 1);
	}
	if (pcg_dat->z_old.rows() != b.rows())
	{
		pcg_dat->z_old.set_size(b.rows(), 1);
	}
	if (pcg_dat->p.rows() != b.rows())
	{
		pcg_dat->p.set_size(b.rows(), 1);
	}
	if (pcg_dat->Ap.rows() != b.rows())
	{
		pcg_dat->Ap.set_size(b.rows(), 1);
	}
	
	//Check for an arbitrary solution at this stage
	success = (*matvec) (pcg_dat->x, pcg_dat->Ap, matvec_data);
	if (success != 0)
	{
		mError(simulation_fail);
		return success;
	}
	pcg_dat->r_old = b - pcg_dat->Ap;
	pcg_dat->res = pcg_dat->r_old.norm();
	pcg_dat->bestres = pcg_dat->res;
	pcg_dat->bestx = pcg_dat->x;
	res_old = pcg_dat->res;
	pcg_dat->relres_base = pcg_dat->res;
	pcg_dat->relres = 1.0;
	if (pcg_dat->Output == true)
		std::cout << "\nRelRes[" << 0 << "] =\t" << pcg_dat->relres << std::endl;
	if (pcg_dat->res <= pcg_dat->tol_abs)
	{
		if (pcg_dat->Output == true)
		{
			std::cout << "\nSolution converged after 0 PCG iterations within residual tolerance!\n";
			std::cout << "AbsRes[0] =\t" << pcg_dat->res << std::endl << std::endl;
		}
		return success;
	}
	
	//Otherwise, begin the PCG algorithm
	if ( (*precon) != NULL)
	{
		success = (*precon) (pcg_dat->r_old, pcg_dat->z_old, precon_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
	}
	else
	{
		pcg_dat->z_old = pcg_dat->r_old;
	}
	pcg_dat->p = pcg_dat->z_old;
	pcg_dat->iter = 0;
	
	//Begin looping till convergence
	for (int k=0; k<pcg_dat->maxit; k++)
	{
		success = (*matvec) (pcg_dat->p, pcg_dat->Ap, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
		pcg_dat->alpha = pcg_dat->r_old.inner_product(pcg_dat->z_old) / pcg_dat->p.inner_product(pcg_dat->Ap);
		pcg_dat->x = pcg_dat->x + (pcg_dat->p * pcg_dat->alpha);
		pcg_dat->r = pcg_dat->r_old - (pcg_dat->Ap * pcg_dat->alpha);
		pcg_dat->iter++;
		
		//Check for convergence
		pcg_dat->res = pcg_dat->r.norm();
		pcg_dat->relres = pcg_dat->res / pcg_dat->relres_base;
		if (pcg_dat->Output == true)
			std::cout << "RelRes[" << k+1 << "] =\t" << pcg_dat->relres << std::endl;
		if (pcg_dat->res <= pcg_dat->tol_abs)
		{
			if (pcg_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << k+1 << " PCG iterations within residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << k+1 << "] =\t" << pcg_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else if (pcg_dat->relres <= pcg_dat->tol_rel)
		{
			if (pcg_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << k+1 << " PCG iterations within relative residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << k+1 << "] =\t" << pcg_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else
		{
			//Check for best residual
			if (pcg_dat->res < pcg_dat->bestres)
			{
				pcg_dat->bestres = pcg_dat->res;
				pcg_dat->bestx = pcg_dat->x;
			}
			res_old = pcg_dat->res;
			
			//Update and prepare for top of loop
			if ( (*precon) != NULL)
			{
				success = (*precon) (pcg_dat->r, pcg_dat->z, precon_data);
				if (success != 0)
				{
					mError(simulation_fail);
					return success;
				}
			}
			else
			{
				pcg_dat->z = pcg_dat->r;
			}
			pcg_dat->beta = pcg_dat->z.inner_product(pcg_dat->r) / pcg_dat->z_old.inner_product(pcg_dat->r_old);
			pcg_dat->p = pcg_dat->z + (pcg_dat->p * pcg_dat->beta);
			
			//Set all old values to current values
			pcg_dat->r_old = pcg_dat->r;
			pcg_dat->z_old = pcg_dat->z;
		}
		
		//Check for failures
		if (isnan(pcg_dat->res) || isinf(pcg_dat->res))
		{
			if (pcg_dat->Output == true)
				std::cout << "\nPCG broke down at the " << k+1 << " iterate!" << std::endl;
			success = -1;
			break;
		}
	}
	
	//Post loop messages
	if (pcg_dat->iter >= pcg_dat->maxit)
	{
		pcg_dat->res = pcg_dat->bestres;
		pcg_dat->x = pcg_dat->bestx;
		if (pcg_dat->Output == true)
		{
			std::cout << "\nPCG did not find a solution within " << pcg_dat->maxit << " iterations!" << std::endl;
			std::cout << "Best absolute residual =\t" << pcg_dat->res << std::endl << std::endl;
		}
		success = 0;
	}
	if (success == -1)
	{
		pcg_dat->res = pcg_dat->bestres;
		pcg_dat->x = pcg_dat->bestx;
		std::cout << "\nPCG did not find a solution within " << pcg_dat->iter << " iterations!" << std::endl;
		std::cout << "Best absolute residual =\t" << pcg_dat->res << std::endl << std::endl;
		success = 0;
	}
	
	return success;
}

//Implementation of the BiCGSTAB iterative method with preconditioning
int bicgstab( int (*matvec) (const Matrix<double>& p, Matrix<double> &Ap, const void *data),
			 int (*precon) (const Matrix<double>& r, Matrix<double> &z, const void *data),
			 Matrix<double> &b, BiCGSTAB_DATA *bicg_dat, const void *matvec_data,
			 const void *precon_data )
{
	int success = 0;
	
	//Check input arguments
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		success = -1;
		return success;
	}
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	
	//Setup the working space and start the method
	if (bicg_dat->maxit <= 0)
	{
		bicg_dat->maxit = std::min(1000,2*b.rows());
	}
	if (bicg_dat->tol_rel >= 1)
	{
		bicg_dat->tol_rel = 1e-6;
	}
	if (bicg_dat->tol_abs >= 1)
	{
		bicg_dat->tol_abs = 1e-6;
	}
	if (bicg_dat->tol_rel < MIN_TOL)
	{
		bicg_dat->tol_rel = MIN_TOL;
	}
	if (bicg_dat->tol_abs < MIN_TOL)
	{
		bicg_dat->tol_abs = MIN_TOL;
	}
	if (bicg_dat->x.rows() != b.rows())
	{
		bicg_dat->x.set_size(b.rows(), 1);
	}
	if (bicg_dat->r.rows() != b.rows())
	{
		bicg_dat->r.set_size(b.rows(), 1);
	}
	if (bicg_dat->r0.rows() != b.rows())
	{
		bicg_dat->r0.set_size(b.rows(), 1);
	}
	if (bicg_dat->v.rows() != b.rows())
	{
		bicg_dat->v.set_size(b.rows(), 1);
	}
	if (bicg_dat->p.rows() != b.rows())
	{
		bicg_dat->p.set_size(b.rows(), 1);
	}
	if (bicg_dat->y.rows() != b.rows())
	{
		bicg_dat->y.set_size(b.rows(), 1);
	}
	if (bicg_dat->s.rows() != b.rows())
	{
		bicg_dat->s.set_size(b.rows(), 1);
	}
	if (bicg_dat->z.rows() != b.rows())
	{
		bicg_dat->z.set_size(b.rows(), 1);
	}
	if (bicg_dat->t.rows() != b.rows())
	{
		bicg_dat->t.set_size(b.rows(), 1);
	}
	
	//Check for the arbitray solution
	success = (*matvec) (bicg_dat->x, bicg_dat->t, matvec_data);
	if (success != 0)
	{
		mError(simulation_fail);
		return success;
	}
	bicg_dat->r = b - bicg_dat->t;
	bicg_dat->r0 = bicg_dat->r;
	bicg_dat->res = bicg_dat->r.norm();
	bicg_dat->bestres = bicg_dat->res;
	bicg_dat->bestx = bicg_dat->x;
	bicg_dat->relres_base = bicg_dat->res;
	bicg_dat->relres = 1.0;
	if (bicg_dat->Output == true)
		std::cout << "\nRelRes[" << 0 << "] =\t" << bicg_dat->relres << std::endl;
	if (bicg_dat->res <= bicg_dat->tol_abs)
	{
		if (bicg_dat->Output == true)
		{
			std::cout << "\nSolution converged after 0 BiCGSTAB iterations within residual tolerance!\n";
			std::cout << "AbsRes[" << 0 << "] =\t" << bicg_dat->res << std::endl << std::endl;
		}
		return success;
	}
	
	//Otherwise, begin the BiCGSTAB algorithm
	bicg_dat->breakdown = false;
	bicg_dat->rho_old = 1.0;
	bicg_dat->alpha = 1.0;
	bicg_dat->omega_old = 1.0;
	bicg_dat->p.set_size(b.rows(), 1);
	bicg_dat->v.set_size(b.rows(), 1);
	bicg_dat->iter = 0;
	
	//Begin looping till convergence
	for (int k=0; k<bicg_dat->maxit; k++)
	{
		bicg_dat->iter++;
		bicg_dat->rho = bicg_dat->r0.inner_product(bicg_dat->r);
		bicg_dat->beta = (bicg_dat->rho / bicg_dat->rho_old) * (bicg_dat->alpha / bicg_dat->omega_old);
		if (isinf(bicg_dat->beta) || isnan(bicg_dat->beta))
		{
			bicg_dat->breakdown = true;
			if (bicg_dat->Output == true)
				std::cout << "\nBiCGSTAB broke down at the " << k+1 << " iterate...\n" << std::endl;
			break;
		}
		bicg_dat->p = bicg_dat->r + (bicg_dat->p * bicg_dat->beta) - (bicg_dat->v * (bicg_dat->beta * bicg_dat->omega_old));
		if ( (*precon) != NULL )
		{
			success = (*precon) (bicg_dat->p, bicg_dat->y, precon_data);
			if (success != 0)
			{
				mError(simulation_fail);
				return -1;
			}
		}
		else
		{
			bicg_dat->y = bicg_dat->p;
		}
		success = (*matvec) (bicg_dat->y, bicg_dat->v, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return -1;
		}
		bicg_dat->alpha = bicg_dat->rho / bicg_dat->r0.inner_product(bicg_dat->v);
		if (isinf(bicg_dat->alpha) || isnan(bicg_dat->alpha))
		{
			bicg_dat->breakdown = true;
			if (bicg_dat->Output == true)
				std::cout << "\nBiCGSTAB broke down at the " << k+1 << " iterate...\n" << std::endl;
			break;
		}
		bicg_dat->s = bicg_dat->r - (bicg_dat->v * bicg_dat->alpha);
		if ( (*precon) != NULL )
		{
			success = (*precon) (bicg_dat->s, bicg_dat->z, precon_data);
			if (success != 0)
			{
				mError(simulation_fail);
				return -1;
			}
		}
		else
		{
			bicg_dat->z = bicg_dat->s;
		}
		success = (*matvec) (bicg_dat->z, bicg_dat->t, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return -1;
		}
		bicg_dat->omega = bicg_dat->t.inner_product(bicg_dat->s)/bicg_dat->t.inner_product(bicg_dat->t);
		if (isinf(bicg_dat->omega) || isnan(bicg_dat->omega))
		{
			bicg_dat->breakdown = true;
			if (bicg_dat->Output == true)
				std::cout << "\nBiCGSTAB broke down at the " << k+1 << " iterate...\n" << std::endl;
			break;
		}
		bicg_dat->x = bicg_dat->x + (bicg_dat->y * bicg_dat->alpha) + (bicg_dat->z * bicg_dat->omega);
		bicg_dat->r = bicg_dat->s - (bicg_dat->t * bicg_dat->omega);
		
		//Check for convergence
		bicg_dat->res = bicg_dat->r.norm();
		bicg_dat->relres = bicg_dat->res / bicg_dat->relres_base;
		if (bicg_dat->Output == true)
			std::cout << "RelRes[" << k+1 << "] =\t" << bicg_dat->relres << std::endl;
		if (bicg_dat->res <= bicg_dat->tol_abs)
		{
			if (bicg_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << k+1 << " BiCGSTAB iterations within residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << k+1 << "] =\t" << bicg_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else if (bicg_dat->relres <= bicg_dat->tol_rel)
		{
			if (bicg_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << k+1 << " BiCGSTAB iterations within relative residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << k+1 << "] =\t" << bicg_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else
		{
			bicg_dat->rho_old = bicg_dat->rho;
			bicg_dat->omega_old = bicg_dat->omega;
			if (bicg_dat->res < bicg_dat->bestres)
			{
				bicg_dat->bestres = bicg_dat->res;
				bicg_dat->bestx = bicg_dat->x;
			}
		}
		
		//Check for failures
		if (isnan(bicg_dat->res) || isinf(bicg_dat->res))
		{
			if (bicg_dat->Output == true)
				std::cout << "\nBiCGSTAB broke down at the " << k+1 << " iterate!" << std::endl;
			success = -1;
			bicg_dat->breakdown = true;
			break;
		}
	}
	
	//Post loop messages
	if (bicg_dat->iter >= bicg_dat->maxit)
	{
		if (bicg_dat->Output == true)
		{
			std::cout << "\nBiCGSTAB did not find a solution within " << bicg_dat->maxit << " iterations!" << std::endl;
			std::cout << "Best absolute residual =\t" << bicg_dat->bestres << std::endl << std::endl;
		}
		bicg_dat->x = bicg_dat->bestx;
		bicg_dat->res = bicg_dat->bestres;
		success = 0;
	}
	if (bicg_dat->breakdown == true)
	{
		if (bicg_dat->Output == true)
			std::cout << "\nReturing the best found solution with norm = \t" << bicg_dat->bestres << std::endl << std::endl;
		bicg_dat->x = bicg_dat->bestx;
		bicg_dat->res = bicg_dat->bestres;
		success = 0;
	}
	
	return success;
}

//Implementation of the CGS method for solving non-symmetric linear systems
int cgs( int (*matvec) (const Matrix<double>& p, Matrix<double> &Ap, const void *data),
		int (*precon) (const Matrix<double>& r, Matrix<double> &z, const void *data),
		Matrix<double> &b, CGS_DATA *cgs_dat, const void *matvec_data,
		const void *precon_data )
{
	int success = 0;
	
	//Check input arguments
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		success = -1;
		return success;
	}
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	
	//Setup the working space and start the method
	if (cgs_dat->maxit <= 0)
	{
		cgs_dat->maxit = std::min(1000,2*b.rows());
	}
	if (cgs_dat->tol_rel >= 1)
	{
		cgs_dat->tol_rel = 1e-6;
	}
	if (cgs_dat->tol_abs >= 1)
	{
		cgs_dat->tol_abs = 1e-6;
	}
	if (cgs_dat->tol_rel < MIN_TOL)
	{
		cgs_dat->tol_rel = MIN_TOL;
	}
	if (cgs_dat->tol_abs < MIN_TOL)
	{
		cgs_dat->tol_abs = MIN_TOL;
	}
	if (cgs_dat->x.rows() != b.rows())
	{
		cgs_dat->x.set_size(b.rows(), 1);
	}
	if (cgs_dat->r.rows() != b.rows())
	{
		cgs_dat->r.set_size(b.rows(), 1);
	}
	if (cgs_dat->r0.rows() != b.rows())
	{
		cgs_dat->r0.set_size(b.rows(), 1);
	}
	if (cgs_dat->u.rows() != b.rows())
	{
		cgs_dat->u.set_size(b.rows(), 1);
	}
	if (cgs_dat->w.rows() != b.rows())
	{
		cgs_dat->w.set_size(b.rows(), 1);
	}
	if (cgs_dat->v.rows() != b.rows())
	{
		cgs_dat->v.set_size(b.rows(), 1);
	}
	if (cgs_dat->p.rows() != b.rows())
	{
		cgs_dat->p.set_size(b.rows(), 1);
	}
	if (cgs_dat->c.rows() != b.rows())
	{
		cgs_dat->c.set_size(b.rows(), 1);
	}
	if (cgs_dat->z.rows() != b.rows())
	{
		cgs_dat->z.set_size(b.rows(), 1);
	}
	
	//Check for the arbitray solution
	success = (*matvec) (cgs_dat->x, cgs_dat->p, matvec_data);
	if (success != 0)
	{
		mError(simulation_fail);
		return success;
	}
	cgs_dat->r = b - cgs_dat->p;
	cgs_dat->r0 = cgs_dat->r;
	cgs_dat->res = cgs_dat->r.norm();
	cgs_dat->bestres = cgs_dat->res;
	cgs_dat->bestx = cgs_dat->x;
	cgs_dat->relres_base = cgs_dat->res;
	cgs_dat->relres = 1.0;
	if (cgs_dat->Output == true)
		std::cout << "\nRelRes[" << 0 << "] =\t" << cgs_dat->relres << std::endl;
	if (cgs_dat->res <= cgs_dat->tol_abs)
	{
		if (cgs_dat->Output == true)
		{
			std::cout << "\nSolution converged after 0 CGS iterations within residual tolerance!\n";
			std::cout << "AbsRes[" << 0 << "] =\t" << cgs_dat->res << std::endl << std::endl;
		}
		return success;
	}
	
	//Otherwise, begin the CGS algorithm
	cgs_dat->u.set_size(b.rows(), 1);
	cgs_dat->w.set_size(b.rows(), 1);
	cgs_dat->alpha = 1.0;
	cgs_dat->sigma = 1.0;
	cgs_dat->breakdown = false;
	cgs_dat->iter = 0;
	
	//Begin looping till convergence
	for (int k=0; k<cgs_dat->maxit; k++)
	{
		cgs_dat->iter++;
		cgs_dat->rho = cgs_dat->r.inner_product(cgs_dat->r0);
		cgs_dat->beta = (1.0/cgs_dat->alpha) * (cgs_dat->rho/cgs_dat->sigma);
		if (isinf(cgs_dat->beta) || isnan(cgs_dat->beta))
		{
			cgs_dat->breakdown = true;
			if (cgs_dat->Output == true)
				std::cout << "\nCGS broke down at the " << k+1 << " iterate...\n" << std::endl;
			break;
		}
		cgs_dat->v = cgs_dat->r + (cgs_dat->u * cgs_dat->beta);
		cgs_dat->w = cgs_dat->v + (cgs_dat->u * cgs_dat->beta) + (cgs_dat->w * (cgs_dat->beta * cgs_dat->beta));
		if ( (*precon) != NULL )
		{
			success = (*precon) (cgs_dat->w, cgs_dat->p, precon_data);
			if (success != 0)
			{
				mError(simulation_fail);
				return -1;
			}
		}
		else
		{
			cgs_dat->p = cgs_dat->w;
		}
		success = (*matvec) (cgs_dat->p, cgs_dat->c, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return -1;
		}
		cgs_dat->sigma = cgs_dat->c.inner_product(cgs_dat->r0);
		cgs_dat->alpha = cgs_dat->rho/cgs_dat->sigma;
		if (isinf(cgs_dat->alpha) || isnan(cgs_dat->alpha))
		{
			cgs_dat->breakdown = true;
			if (cgs_dat->Output == true)
				std::cout << "\nCGS broke down at the " << k+1 << " iterate...\n" << std::endl;
			break;
		}
		cgs_dat->u = cgs_dat->v - (cgs_dat->c * cgs_dat->alpha);
		cgs_dat->z = cgs_dat->u + cgs_dat->v;
		if ( (*precon) != NULL )
		{
			success = (*precon) (cgs_dat->z, cgs_dat->p, precon_data);
			if (success != 0)
			{
				mError(simulation_fail);
				return -1;
			}
		}
		else
		{
			cgs_dat->p = cgs_dat->z;
		}
		cgs_dat->x = cgs_dat->x + (cgs_dat->p * cgs_dat->alpha);
		success = (*matvec) (cgs_dat->x, cgs_dat->p, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return -1;
		}
		cgs_dat->r = b - cgs_dat->p;
		
		//Check for convergence
		cgs_dat->res = cgs_dat->r.norm();
		cgs_dat->relres = cgs_dat->res / cgs_dat->relres_base;
		if (cgs_dat->Output == true)
			std::cout << "RelRes[" << k+1 << "] =\t" << cgs_dat->relres << std::endl;
		if (cgs_dat->res <= cgs_dat->tol_abs)
		{
			if (cgs_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << k+1 << " CGS iterations within residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << k+1 << "] =\t" << cgs_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else if (cgs_dat->relres <= cgs_dat->tol_rel)
		{
			if (cgs_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << k+1 << " CGS iterations within relative residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << k+1 << "] =\t" << cgs_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else
		{
			if (cgs_dat->res < cgs_dat->bestres)
			{
				cgs_dat->bestres = cgs_dat->res;
				cgs_dat->bestx = cgs_dat->x;
			}
		}
		
		//Check for failures
		if (isnan(cgs_dat->res) || isinf(cgs_dat->res))
		{
			if (cgs_dat->Output == true)
				std::cout << "\nCGS broke down at the " << k+1 << " iterate!" << std::endl;
			success = -1;
			cgs_dat->breakdown = true;
			break;
		}
	}
	
	//Post loop messages
	if (cgs_dat->iter >= cgs_dat->maxit)
	{
		if (cgs_dat->Output == true)
		{
			std::cout << "\nCGS did not find a solution within " << cgs_dat->maxit << " iterations!" << std::endl;
			std::cout << "Best absolute residual =\t" << cgs_dat->bestres << std::endl << std::endl;
		}
		cgs_dat->x = cgs_dat->bestx;
		cgs_dat->res = cgs_dat->bestres;
		success = 0;
	}
	if (cgs_dat->breakdown == true)
	{
		if (cgs_dat->Output == true)
		{
			std::cout <<"\nReturing the best found solution with norm = \t" << cgs_dat->bestres << std::endl << std::endl;
		}
		cgs_dat->x = cgs_dat->bestx;
		cgs_dat->res = cgs_dat->bestres;
		success = 0;
	}
	
	return success;
}

//Function for forming the tranpose multiplication operation for a given linear operator
int operatorTranspose(int(*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data),
					  Matrix<double> &r, Matrix<double> &u, OPTRANS_DATA *transpose_dat,
					  const void *matvec_data)
{
	int success = 0;
	
	//Check arguments
	if ((*matvec) == NULL)
	{
		mError(nullptr_error);
		return -1;
	}
	if (r.rows() < 2)
	{
		mError(matrix_too_small);
		return -1;
	}
	if (u.rows() != r.rows())
	{
		u.set_size(r.rows(), 1);
	}
	if (transpose_dat->Ai.rows() != u.rows())
	{
		transpose_dat->Ai.set_size(u.rows(), 1);
	}
	if (transpose_dat->Ii.rows() != u.rows())
	{
		transpose_dat->Ii.set_size(u.rows(), 1);
	}
	transpose_dat->Ii.zeros();
	
	//Form the matrix elements of u as r^T*A*Ii == A^T*r @ i
	transpose_dat->Ii.edit(0, 0, 1.0);
	for (int i = 0; i < u.rows(); i++)
	{
		//Form the matrix-vector product of A with the identity vector
		success = (*matvec) (transpose_dat->Ii, transpose_dat->Ai, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
		u.edit(i, 0, r.inner_product(transpose_dat->Ai));
		
		//Reset the identity vector for next loop
		transpose_dat->Ii.edit(i, 0, 0.0);
		if (i < (u.rows()-1)) transpose_dat->Ii.edit(i+1, 0, 1.0);
		
	}
	
	return success;
}

//Function for solving a non-symmetric linear system using GCR
int gcr( int (*matvec) (const Matrix<double>& x, Matrix<double> &Ax, const void *data),
		int (*precon) (const Matrix<double>& r, Matrix<double> &Mr, const void *data),
		Matrix<double> &b, GCR_DATA *gcr_dat, const void *matvec_data,
		const void *precon_data )
{
	int success = 0;
	double c_norm = 0.0;
	
	//Check arguments
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		success = -1;
		return success;
	}
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	
	//Initialize Data
	if (gcr_dat->restart <= 0)
	{
		gcr_dat->restart = std::min(b.rows(),20);
	}
	else if (gcr_dat->restart > b.rows() && b.rows() <= 1000)
	{
		gcr_dat->restart = b.rows();
	}
	if (gcr_dat->restart == b.rows())
	{
		gcr_dat->maxit = 1;
	}
	else if (gcr_dat->maxit <= 0 || gcr_dat->maxit > b.rows())
	{
		gcr_dat->maxit = std::min(b.rows(),1000);
	}
	if (gcr_dat->tol_rel >= 1)
	{
		gcr_dat->tol_rel = 1e-6;
	}
	if (gcr_dat->tol_abs >= 1)
	{
		gcr_dat->tol_abs = 1e-6;
	}
	if (gcr_dat->tol_rel < MIN_TOL)
	{
		gcr_dat->tol_rel = MIN_TOL;
	}
	if (gcr_dat->tol_abs < MIN_TOL)
	{
		gcr_dat->tol_abs = MIN_TOL;
	}
	if (gcr_dat->x.rows() != b.rows())
	{
		gcr_dat->x.set_size(b.rows(), 1);
	}
	if (gcr_dat->c_temp.rows() != b.rows())
	{
		gcr_dat->c_temp.set_size(b.rows(), 1);
	}
	if (gcr_dat->u_temp.rows() != b.rows())
	{
		gcr_dat->u_temp.set_size(b.rows(), 1);
	}
	gcr_dat->iter_inner = 0;
	gcr_dat->iter_outer = 0;
	gcr_dat->breakdown = false;
	
	//Form the first matrix-vector product and check for arbitray solution
	success = (*matvec) (gcr_dat->x,gcr_dat->c_temp,matvec_data);
	if (success != 0)
	{
		mError(simulation_fail);
		return success;
	}
	gcr_dat->r = b - gcr_dat->c_temp;
	gcr_dat->relres_base = gcr_dat->r.norm();
	gcr_dat->res = gcr_dat->relres_base;
	gcr_dat->bestres = gcr_dat->res;
	gcr_dat->bestx = gcr_dat->x;
	gcr_dat->relres = 1.0;
	if (gcr_dat->Output == true)
		std::cout << "\nRelRes[" << 0 << "] =\t" << gcr_dat->relres << "\t---Start[0]---" << std::endl;
	
	//Check for imediate convergence
	if (gcr_dat->res < gcr_dat->tol_abs)
	{
		if (gcr_dat->Output == true)
		{
			std::cout << "\nSolution converged after 0 GCR iterations witin residual tolerance!" << std::endl;
			std::cout << "AbsRes[0] =\t" << gcr_dat->res << std::endl << std::endl;
		}
		return success;
	}
	
	//Begin outer iterations to find solution
	int step = 0;
	for (int m=0; m<gcr_dat->maxit; m++)
	{
		//Begin inner iterations to build the subspace
		for (int k=0; k<gcr_dat->restart; k++)
		{
			//Apply preconditioner  if available
			if ( (*precon) != NULL)
			{
				success = (*precon) (gcr_dat->r, gcr_dat->u_temp, precon_data);
				if (success != 0)
				{
					mError(simulation_fail);
					return success;
				}
			}
			else
			{
				gcr_dat->u_temp = gcr_dat->r;
			}
			
			//Check for a breakdown of the algorithm
			if (gcr_dat->u_temp.norm() == 0.0)
			{
				success = operatorTranspose(matvec, gcr_dat->r, gcr_dat->u_temp, &gcr_dat->transpose_dat, matvec_data);
				if (success != 0)
				{
					mError(simulation_fail);
					return success;
				}
			}
			
			//Form initial c vector
			success = (*matvec) (gcr_dat->u_temp,gcr_dat->c_temp,matvec_data);
			if (success != 0)
			{
				mError(simulation_fail);
				return success;
			}
			
			//Begin subspace orthogonal vector construction
			for (int i=0; i<k; i++)
			{
				gcr_dat->alpha = gcr_dat->c[i].inner_product(gcr_dat->c_temp);
				for (int n=0; n<b.rows(); n++)
				{
					gcr_dat->c_temp.edit(n, 0, gcr_dat->c_temp(n,0) - (gcr_dat->c[i](n,0) * gcr_dat->alpha));
					gcr_dat->u_temp.edit(n, 0, gcr_dat->u_temp(n,0) - (gcr_dat->u[i](n,0) * gcr_dat->alpha));
				}
				
			}//End Subspace Construction
			
			//Check for a breakdown of the algorithm
			c_norm = gcr_dat->c_temp.norm();
			if (c_norm == 0.0)
			{
				gcr_dat->breakdown = true;
				if (gcr_dat->Output == true)
					std::cout << "\nGCR may have broken down at the " << m+step+1 << " iterate...\n" << std::endl;
				break;
			}
			
			//Form and store the new vectors
			if (gcr_dat->c.size() < k+1)
			{
				gcr_dat->c.push_back( gcr_dat->c_temp * (1.0/c_norm) );
			}
			else
			{
				gcr_dat->c[k] = gcr_dat->c_temp * (1.0/c_norm);
			}
			if (gcr_dat->u.size() < k+1)
			{
				gcr_dat->u.push_back( gcr_dat->u_temp * (1.0/c_norm) );
			}
			else
			{
				gcr_dat->u[k] = gcr_dat->u_temp * (1.0/c_norm);
			}
			
			//Update solution and residuals
			gcr_dat->beta = gcr_dat->c[k].inner_product(gcr_dat->r);
			for (int n=0; n<b.rows(); n++)
			{
				gcr_dat->x.edit(n, 0, gcr_dat->x(n,0) + (gcr_dat->u[k](n,0) * gcr_dat->beta));
				gcr_dat->r.edit(n, 0, gcr_dat->r(n,0) - (gcr_dat->c[k](n,0) * gcr_dat->beta));
			}
			gcr_dat->res = gcr_dat->r.norm();
			gcr_dat->relres = gcr_dat->res / gcr_dat->relres_base;
			gcr_dat->total_iter = m+step+1;
			gcr_dat->iter_inner++;
			
			//Check residuals for convergence
			if (gcr_dat->Output == true)
				std::cout << "RelRes[" << m+step+1 << "] =\t" << gcr_dat->relres << std::endl;
			if (gcr_dat->res <= gcr_dat->tol_abs)
			{
				if (gcr_dat->Output == true)
				{
					std::cout << "\nSolution converged after " << m+step+1 << " GCR iterations within residual tolerance!" << std::endl;
					std::cout << "AbsRes[" << m+step+1 << "] =\t" << gcr_dat->res << std::endl << std::endl;
				}
				return success;
			}
			else if (gcr_dat->relres <= gcr_dat->tol_rel)
			{
				if (gcr_dat->Output == true)
				{
					std::cout << "\nSolution converged after " << m+step+1 << " GCR iterations within relative residual tolerance!" << std::endl;
					std::cout << "AbsRes[" << m+step+1 << "] =\t" << gcr_dat->res << std::endl << std::endl;
				}
				return success;
			}
			else
			{
				if (gcr_dat->res < gcr_dat->bestres)
				{
					gcr_dat->bestres = gcr_dat->res;
					gcr_dat->bestx = gcr_dat->x;
				}
			}
			
			//Check for failures
			if (isnan(gcr_dat->res) || isinf(gcr_dat->res))
			{
				if (gcr_dat->Output == true)
					std::cout << "\nGCR broke down at the " << step+m+1 << " iterate!" << std::endl;
				gcr_dat->breakdown = true;
				break;
			}
			
			step++;
		}// END inner iterations
		
		//Setup for restart if not yet converged and check for actual convergence
		success = (*matvec) (gcr_dat->x,gcr_dat->c_temp,matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
		gcr_dat->r = b - gcr_dat->c_temp;
		gcr_dat->res = gcr_dat->r.norm();
		gcr_dat->relres = gcr_dat->res / gcr_dat->relres_base;
		gcr_dat->iter_outer++;
		gcr_dat->total_iter = gcr_dat->iter_outer + gcr_dat->iter_inner;
		
		//Check residuals for convergence
		if (gcr_dat->Output == true)
			std::cout << "RelRes[" << step+m+1 << "] =\t" << gcr_dat->relres << "\t---Restart[" << m+1 << "]---" << std::endl;
		if (gcr_dat->res <= gcr_dat->tol_abs)
		{
			if (gcr_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << m+step+1 << " GCR iterations within residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << m+step+1 << "] =\t" << gcr_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else if (gcr_dat->relres <= gcr_dat->tol_rel)
		{
			if (gcr_dat->Output == true)
			{
				std::cout << "\nSolution converged after " << m+step+1 << " GCR iterations within relative residual tolerance!" << std::endl;
				std::cout << "AbsRes[" << m+step+1 << "] =\t" << gcr_dat->res << std::endl << std::endl;
			}
			return success;
		}
		else
		{
			if (gcr_dat->res < gcr_dat->bestres)
			{
				gcr_dat->bestres = gcr_dat->res;
				gcr_dat->bestx = gcr_dat->x;
			}
		}
		
		//Check for failures
		if (isnan(gcr_dat->res) || isinf(gcr_dat->res))
		{
			if (gcr_dat->Output == true)
				std::cout << "\nGCR broke down at the " << step+m+1 << " iterate!" << std::endl;
			gcr_dat->breakdown = true;
			break;
		}
		
	}// END outer iterations
	
	//Post loop messages
	if (gcr_dat->iter_outer >= gcr_dat->maxit)
	{
		if (gcr_dat->Output == true)
		{
			std::cout << "\nGCR did not find a solution within " << step+gcr_dat->maxit << " iterations!" << std::endl;
			std::cout << "Best absolute residual =\t" << gcr_dat->bestres << std::endl << std::endl;
		}
		gcr_dat->x = gcr_dat->bestx;
		gcr_dat->res = gcr_dat->bestres;
		success = 0;
	}
	if (gcr_dat->breakdown == true)
	{
		if (gcr_dat->Output == true)
			std::cout << "\nReturing the best found solution with norm = \t" << gcr_dat->bestres << std::endl << std::endl;
		gcr_dat->x = gcr_dat->bestx;
		gcr_dat->res = gcr_dat->bestres;
		success = 0;
	}
	
	return success;
}

//Function for preconditioning GCR in the GMRESR application
int gmresrPreconditioner( const Matrix<double>& r, Matrix<double> &Mr, const void *data)
{
	int success = 0;
	GMRESR_DATA *dat = (GMRESR_DATA *) data;
	
	//Check for errors
	if ( dat->matvec == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	dat->gmres_dat.restart = dat->gmres_restart;
	dat->gmres_dat.maxit = dat->gmres_maxit;
	dat->gmres_dat.tol_abs = dat->gmres_tol * dat->gcr_dat.res;
	dat->gmres_dat.tol_rel = dat->gmres_tol;
	
	//Enfore a zero initial guess
	if (dat->gmres_dat.x.rows() != dat->N)
	{
		dat->gmres_dat.x.set_size(dat->N, 1);
	}
	else
	{
		dat->gmres_dat.x.zeros();
	}
	dat->arg = r;
	//Call GMRES to precondition
	success = gmresRightPreconditioned(dat->matvec, dat->terminal_precon, dat->arg, &dat->gmres_dat, dat->matvec_data, dat->term_precon);
	if (success != 0) {mError(simulation_fail); return -1;}
	dat->iter_inner+=dat->gmres_dat.iter_total;
	
	//Set the vector to the updated value
	Mr = dat->gmres_dat.x;
	
	return success;
}

//Function for the implementation of the GMRESR algorithm
int gmresr( int (*matvec) (const Matrix<double>& x, Matrix<double> &Ax, const void *data),
		   int (*terminal_precon) (const Matrix<double>& r, Matrix<double> &Mr, const void *data),
		   Matrix<double> &b, GMRESR_DATA *gmresr_dat, const void *matvec_data,
		   const void *term_precon_data )
{
	int success = 0;
	
	//Check input args for errors
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	else
	{
		gmresr_dat->matvec = (*matvec);
	}
	gmresr_dat->matvec_data = matvec_data;
	gmresr_dat->terminal_precon = terminal_precon;
	gmresr_dat->term_precon = term_precon_data;
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	else
	{
		gmresr_dat->N = b.rows();
	}
	gmresr_dat->gmres_dat.Output = gmresr_dat->GMRES_Output;
	gmresr_dat->gcr_dat.Output = gmresr_dat->GCR_Output;
	
	//GMRESR Precondtioner Tolerance
	if (gmresr_dat->gmres_tol >= 1)
	{
		gmresr_dat->gmres_tol = 0.1;
	}
	if (gmresr_dat->gmres_tol < 1e-4)
	{
		gmresr_dat->gmres_tol = 1e-4;
	}
	
	//GCR Tolerances
	if (gmresr_dat->gcr_rel_tol >= 1)
	{
		gmresr_dat->gcr_rel_tol = 1e-6;
	}
	if (gmresr_dat->gcr_abs_tol >= 1)
	{
		gmresr_dat->gcr_abs_tol = 1e-6;
	}
	if (gmresr_dat->gcr_abs_tol < MIN_TOL)
	{
		gmresr_dat->gcr_abs_tol = MIN_TOL;
	}
	if (gmresr_dat->gcr_rel_tol < MIN_TOL)
	{
		gmresr_dat->gcr_rel_tol = MIN_TOL;
	}
	gmresr_dat->gcr_dat.tol_abs = gmresr_dat->gcr_abs_tol;
	gmresr_dat->gcr_dat.tol_rel = gmresr_dat->gcr_rel_tol;
	if (gmresr_dat->gcr_restart <= 0)
	{
		gmresr_dat->gcr_restart = std::min(b.rows(),20);
	}
	else if (gmresr_dat->gcr_restart > b.rows() && b.rows() <= 1000)
	{
		gmresr_dat->gcr_restart = b.rows();
	}
	if (gmresr_dat->gcr_restart == b.rows())
	{
		gmresr_dat->gcr_maxit = 1;
	}
	else if (gmresr_dat->gcr_maxit <= 0 || gmresr_dat->gcr_maxit > b.rows())
	{
		gmresr_dat->gcr_maxit = std::min(b.rows(),1000);
	}
	gmresr_dat->gcr_dat.restart = gmresr_dat->gcr_restart;
	gmresr_dat->gcr_dat.maxit = gmresr_dat->gcr_maxit;
	if (gmresr_dat->gmres_restart <= 1)
	{
		gmresr_dat->gmres_restart = std::min(b.rows(),20);
	}
	else if (gmresr_dat->gmres_restart > b.rows())
	{
		gmresr_dat->gmres_restart = std::min(b.rows(),20);
	}
	if (gmresr_dat->gmres_restart == b.rows())
	{
		gmresr_dat->gmres_maxit = 1;
	}
	else if (gmresr_dat->gmres_maxit <= 0 || gmresr_dat->gmres_maxit > b.rows())
	{
		gmresr_dat->gmres_maxit = std::min(b.rows(),5);
	}
	gmresr_dat->gmres_dat.restart = gmresr_dat->gmres_restart;
	gmresr_dat->gmres_dat.maxit = gmresr_dat->gmres_maxit;
	gmresr_dat->total_iter = 0;
	gmresr_dat->iter_inner = 0;
	gmresr_dat->iter_outer = 0;
	
	//Call GCR with the GMRES preconditioner
	success = gcr(gmresr_dat->matvec, gmresrPreconditioner, b, &gmresr_dat->gcr_dat, gmresr_dat->matvec_data, (void *)gmresr_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	gmresr_dat->iter_outer = gmresr_dat->gcr_dat.total_iter;
	gmresr_dat->total_iter = gmresr_dat->iter_outer + gmresr_dat->iter_inner;
	
	return success;
}

// Preconditioner function for the Krylov Multi-Space
int kmsPreconditioner( const Matrix<double>& r, Matrix<double> &Mr, const void *data)
{
	int success = 0;
	KMS_DATA * kms_dat = (KMS_DATA *) data;
	
	//Check for errors
	if ( kms_dat->matvec == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	kms_dat->level++;
	kms_dat->gmres_in[kms_dat->level-1].Output = kms_dat->Output_inner;
	
	//Formulate the best number of inner and outer preconditioning iterates based on current status
	double inner_par = 0;
	
	if (kms_dat->level == 1)
		inner_par += 3;
	else if (kms_dat->level > 1 && kms_dat->level <= 2)
		inner_par += 5;
	else if (kms_dat->level > 2 && kms_dat->level <= 3)
		inner_par += 8;
	else if (kms_dat->level > 3 && kms_dat->level <= 4)
		inner_par += 12;
	else
		inner_par += 20;
	
	if (kms_dat->max_level == 1)
		inner_par += 20;
	else if (kms_dat->max_level > 1 && kms_dat->max_level <= 2)
		inner_par += 12;
	else if (kms_dat->max_level > 2 && kms_dat->max_level <= 3)
		inner_par += 8;
	else if (kms_dat->max_level > 3 && kms_dat->max_level <= 4)
		inner_par += 5;
	else
		inner_par += 3;
	
	kms_dat->gmres_in[kms_dat->level - 1].restart = (int) (inner_par / 2.0);
	
	if (kms_dat->gmres_in[kms_dat->level - 1].restart <= 20 && kms_dat->gmres_in[kms_dat->level - 1].restart > 10)
		kms_dat->gmres_in[kms_dat->level - 1].maxit = 1;
	else
		kms_dat->gmres_in[kms_dat->level - 1].maxit = 2;
	
	if (kms_dat->level > 1)
	{
		kms_dat->gmres_in[kms_dat->level-1].tol_abs = kms_dat->gmres_in[kms_dat->level-2].res * kms_dat->inner_reltol;
		kms_dat->gmres_in[kms_dat->level-1].tol_rel = kms_dat->inner_reltol;
	}
	else
	{
		kms_dat->gmres_in[kms_dat->level-1].tol_abs = kms_dat->gmres_out.res * kms_dat->inner_reltol;
		kms_dat->gmres_in[kms_dat->level-1].tol_rel = kms_dat->inner_reltol;
	}
	
	Matrix<double> temp = r;
	if (kms_dat->gmres_in[kms_dat->level - 1].x.rows() != Mr.rows())
		kms_dat->gmres_in[kms_dat->level - 1].x.set_size(Mr.rows(), 1);
	else
		kms_dat->gmres_in[kms_dat->level - 1].x.zeros();
	if (kms_dat->level == kms_dat->max_level)
		success = gmresRightPreconditioned(kms_dat->matvec, kms_dat->terminal_precon, temp, &kms_dat->gmres_in[kms_dat->level-1], kms_dat->matvec_data, kms_dat->term_precon);
	else
		success = gmresRightPreconditioned(kms_dat->matvec, kmsPreconditioner, temp, &kms_dat->gmres_in[kms_dat->level-1], kms_dat->matvec_data, (void *)kms_dat);
	kms_dat->inner_iter += kms_dat->gmres_in[kms_dat->level-1].iter_total;
	Mr = kms_dat->gmres_in[kms_dat->level-1].x;
	kms_dat->level--;
	
	return success;
}

// Function to iteratively solve a non-symmetric, indefinite linear system with KMS
int krylovMultiSpace( int (*matvec) (const Matrix<double>& x, Matrix<double> &Ax, const void *data),
					 int (*terminal_precon) (const Matrix<double>& r, Matrix<double> &Mr, const void *data),
					 Matrix<double> &b, KMS_DATA *kms_dat, const void *matvec_data,
					 const void *term_precon_data )
{
	int success = 0;
	
	//Check input args for errors
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	else
	{
		kms_dat->matvec = (*matvec);
	}
	kms_dat->matvec_data = matvec_data;
	kms_dat->terminal_precon = terminal_precon;
	kms_dat->term_precon = term_precon_data;
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	kms_dat->gmres_out.Output = kms_dat->Output_outer;
	
	//Preconditioner Tolerances
	if (kms_dat->inner_reltol >= 0.5)
	{
		kms_dat->inner_reltol = 0.5;
	}
	if (kms_dat->inner_reltol < 1e-4)
	{
		kms_dat->inner_reltol = 1e-4;
	}
	
	//Outer Step Tolerances
	if (kms_dat->outer_reltol >= 1)
	{
		kms_dat->outer_reltol = 1e-6;
	}
	if (kms_dat->outer_abstol >= 1)
	{
		kms_dat->outer_abstol = 1e-6;
	}
	if (kms_dat->outer_abstol < MIN_TOL)
	{
		kms_dat->outer_abstol = MIN_TOL;
	}
	if (kms_dat->outer_reltol < MIN_TOL)
	{
		kms_dat->outer_reltol = MIN_TOL;
	}
	kms_dat->gmres_out.tol_abs = kms_dat->outer_abstol;
	kms_dat->gmres_out.tol_rel = kms_dat->outer_reltol;
	
	//Method Parameters
	if (kms_dat->restart <= 0)
	{
		kms_dat->restart = std::min(b.rows(),20);
	}
	else if (kms_dat->restart > b.rows() && b.rows() <= 1000)
	{
		kms_dat->restart = b.rows();
	}
	if (kms_dat->restart == b.rows())
	{
		kms_dat->maxit = 1;
	}
	else if (kms_dat->maxit <= 0 || kms_dat->maxit > b.rows())
	{
		kms_dat->maxit = std::min(b.rows(),1000);
	}
	if (kms_dat->max_level < 0)
		kms_dat->max_level = 0;
	if (kms_dat->max_level > 5)
		kms_dat->max_level = 5;
	kms_dat->gmres_out.restart = kms_dat->restart;
	kms_dat->gmres_out.maxit = kms_dat->maxit;
	kms_dat->level = 0;
	kms_dat->total_iter = 0;
	kms_dat->inner_iter = 0;
	kms_dat->outer_iter = 0;
	
	if (kms_dat->gmres_in.size() != kms_dat->max_level)
		kms_dat->gmres_in.resize(kms_dat->max_level);
	
	if (kms_dat->max_level == 0)
	{
		success = gmresRightPreconditioned(kms_dat->matvec, kms_dat->terminal_precon, b, &kms_dat->gmres_out, kms_dat->matvec_data, kms_dat->term_precon);
		if (success != 0) { mError(simulation_fail); return -1; }
	}
	else
	{
		success = gmresRightPreconditioned(kms_dat->matvec, kmsPreconditioner, b, &kms_dat->gmres_out, kms_dat->matvec_data, (void *)kms_dat);
		if (success != 0) { mError(simulation_fail); return -1; }
	}
	kms_dat->outer_iter = kms_dat->gmres_out.iter_total;
	kms_dat->total_iter = kms_dat->outer_iter + kms_dat->inner_iter;
	
	return success;
}

//Function to perform the QR factorization of an invertable linear operator
int QRsolve( int (*matvec) (const Matrix<double>& x, Matrix<double> &Ax, const void *data),
			Matrix<double> &b, QR_DATA *qr_dat, const void *matvec_data)
{
	int success = 0;
	
	//Check for errors
	if ( (*matvec) == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	if (b.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	
	//Setup working environment
	qr_dat->ek.set_size(b.rows(), 1);
	qr_dat->x.set_size(b.rows(), 1);
	qr_dat->Ro.set_size(b.rows(), b.rows());
	qr_dat->ek.zeros();
	
	//Initialize the matrix R with linear operator entries
	for (int k=0; k<b.rows(); k++)
	{
		qr_dat->ek.edit(k, 0, 1.0);
		success = (*matvec) (qr_dat->ek, qr_dat->x, matvec_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return success;
		}
		qr_dat->Ro.columnReplace(k, qr_dat->x);
		qr_dat->ek.edit(k, 0, 0.0);
	}
	
	//Call the QR solve function from macaw
	qr_dat->x.qrSolve(qr_dat->Ro, b);
	
	return success;
}

//Function for solving a non-linear system using a Picard or Fixed-Point iteration
int picard( int (*res) (const Matrix<double>& x, Matrix<double> &r, const void *data),
		   int (*evalx) (const Matrix<double>& x0, Matrix<double> &x, const void *data),
		   Matrix<double> &x, PICARD_DATA *picard_dat, const void *res_data,
		   const void *evalx_data )
{
	int success = 0;
	double res_old;
	
	//Check for existance of user supplied function
	if ((*res) == NULL)
	{
		mError(nullptr_func);
		success = -1;
		return success;
	}
	if ((*evalx) == NULL)
	{
		mError(nullptr_func);
		success = -1;
		return success;
	}
	
	//Setup working space
	if ( picard_dat->maxit <= 0 )
	{
		picard_dat->maxit = std::min(3*x.rows(),1000);
	}
	if (picard_dat->tol_rel >= 1.0)
	{
		picard_dat->tol_rel = 1e-6;
	}
	if (picard_dat->tol_abs >= 1.0)
	{
		picard_dat->tol_abs = 1e-6;
	}
	if (picard_dat->tol_rel < MIN_TOL)
	{
		picard_dat->tol_rel = MIN_TOL;
	}
	if (picard_dat->tol_abs < MIN_TOL)
	{
		picard_dat->tol_abs = MIN_TOL;
	}
	if ( picard_dat->x0.rows() != x.rows() )
	{
		picard_dat->x0.set_size(x.rows(), 1);
	}
	if ( picard_dat->bestx.rows() != x.rows() )
	{
		picard_dat->bestx.set_size(x.rows(), 1);
	}
	if ( picard_dat->r.rows() != x.rows() )
	{
		picard_dat->r.set_size(x.rows(), 1);
	}
	picard_dat->iter = 0;
	
	//Start the method by forming first residual
	success = (*res) (x,picard_dat->r,res_data);
	if (success != 0)
	{
		mError(simulation_fail);
		return -1;
	}
	picard_dat->bestx = x;
	picard_dat->res = picard_dat->r.norm();
	res_old = picard_dat->res;
	picard_dat->bestres = picard_dat->res;
	picard_dat->relres = 1.0;
	picard_dat->relres_base = picard_dat->res;
	if (picard_dat->Output == true)
	{
		std::cout << "Begin Picard iterations...\n" << std::endl;
		std::cout << "Norm[" << 0 << "] =\t" << picard_dat->res << std::endl;
	}
	picard_dat->x0 = x;
	
	//Loop till convergence
	for (int k=0; k<picard_dat->maxit; k++)
	{
		picard_dat->iter++;
		
		//Form the approximation to the solution with user supplied function
		success = (*evalx) (picard_dat->x0,x,evalx_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return -1;
		}
		
		//Call user residual function
		success = (*res) (x,picard_dat->r,res_data);
		if (success != 0)
		{
			mError(simulation_fail);
			return -1;
		}
		picard_dat->res = picard_dat->r.norm();
		picard_dat->relres = picard_dat->res/picard_dat->relres_base;
		
		//Check for convergence
		if (picard_dat->Output == true)
			std::cout << "Norm[" << k+1 << "] =\t" << picard_dat->res << std::endl;
		if (picard_dat->res <= picard_dat->tol_abs)
		{
			if (picard_dat->Output == true)
				std::cout << "\nSolution converged after " << k+1 << " Picard iterations within residual tolerance!" << std::endl;
			return success;
		}
		else if (picard_dat->relres <= picard_dat->tol_rel)
		{
			if (picard_dat->Output == true)
				std::cout << "\nSolution converged after " << k+1 << " Picard iterations within relative residual tolerance!" << std::endl;
			return success;
		}
		else
		{
			if (picard_dat->res < picard_dat->bestres)
			{
				picard_dat->bestres = picard_dat->res;
				picard_dat->bestx = x;
			}
			picard_dat->x0 = x;
		}
		res_old = picard_dat->res;
		
		//Check for failures
		if (isnan(picard_dat->res) || isinf(picard_dat->res))
		{
			if (picard_dat->Output == true)
				std::cout << "\nPicard broke down at the " << k+1 << " iterate!" << std::endl;
			success = -1;
			break;
		}
	}
	
	//Post loop messages
	if (picard_dat->iter >= picard_dat->maxit)
	{
		if (picard_dat->Output == true)
		{
			std::cout << "\nPicard did not find exact solution within " << picard_dat->maxit << " iterations!" << std::endl;
			std::cout << "Best reported norm =\t" << picard_dat->bestres << std::endl;
		}
		x = picard_dat->bestx;
		picard_dat->res = picard_dat->bestres;
		success = 0;
	}
	if (success == -1)
	{
		if (picard_dat->Output == true)
		{
			std::cout << "\nPicard did not find exact solution within " << picard_dat->iter << " iterations!" << std::endl;
			std::cout << "Best reported norm =\t" << picard_dat->bestres << std::endl;
		}
		x = picard_dat->bestx;
		picard_dat->res = picard_dat->bestres;
		success = 0;
	}
	
	return success;
}

//Jacobian-Vector Product function used by the PJFNK implementation
int jacvec(const Matrix<double>& v, Matrix<double>& Jv, const void *data)
{
	int success = 0;
	PJFNK_DATA *dat = (PJFNK_DATA *) data;
	double Jv_i;
	//Forward Difference
	for (int i=0; i<Jv.rows(); i++)
	{
		dat->v.edit(i,0,dat->x(i,0) + (dat->eps * v(i,0)));
	}
	success = (*dat->funeval) ( dat->v, dat->Fv, dat->res_data);
	if (success != 0)
	{
		mError(simulation_fail);
		return -1;
	}
	for (int i=0; i<Jv.rows(); i++)
	{
		Jv_i = (dat->Fv(i,0) - dat->F(i,0)) / dat->eps;
		Jv.edit(i,0,Jv_i);
	}
	
	return success;
}

//Implementation of the Backtracking Linesearch Algorithm for Newton's Method
int backtrackLineSearch( int (*feval) (const Matrix<double>& x, Matrix<double> &F, const void *data),
						Matrix<double> &Fkp1, Matrix<double> &xkp1, Matrix<double> &pk, double normFk,
						BACKTRACK_DATA *backtrack_dat, const void *feval_data)
{
	int success = 0;
	double slope, lambda_old, lambda, norm_old, lambda_temp;
	double a, b, disc;
	
	//Check input arguments for problems
	if ( (*feval) == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	if (Fkp1.rows() < 2 || xkp1.rows() < 2 || pk.rows() < 2)
	{
		mError(matrix_too_small);
		return -1;
	}
	if (normFk < 0)
	{
		mError(invalid_norm);
		return -1;
	}
	if (Fkp1.rows() != xkp1.rows() && Fkp1.rows() != pk.rows())
	{
		mError(dim_mis_match);
		return -1;
	}
	
	//Check an initialize all Backtracking Data
	if (backtrack_dat->constRho == true)
	{
		if (backtrack_dat->rho < 0.1 || backtrack_dat->rho > 0.5)
		{
			backtrack_dat->rho = 0.1;
		}
	}
	if (backtrack_dat->alpha <= 0.0 || backtrack_dat->alpha >= 0.5)
	{
		backtrack_dat->alpha = 1e-4;
	}
	if (backtrack_dat->lambdaMin > 1.0 || backtrack_dat->lambdaMin < DBL_EPSILON)
	{
		backtrack_dat->lambdaMin = DBL_EPSILON;
	}
	backtrack_dat->fun_call = 0;
	
	//Start the backtracking algorithm
	backtrack_dat->xk = xkp1;
	backtrack_dat->Fk = Fkp1;
	slope = backtrack_dat->Fk.inner_product(pk*(-1.0));
	if (isnan(slope) || isinf(slope))
	{
		return 1;
	}
	lambda = 1.0;
	lambda_old = lambda;
	
	/*
	 
	 Important Note: The pk search vector being passed to this algorithm is the POSITIVE
	 gradient and therefore not the decent direction. To make the result into a decent direction
	 the scalar multiple, lambda, must be a negative number between -1.0 and -lambdaMin OR we
	 need to subtract pk from xk instead of adding it.
	 
	 This is an artifact of the PJFNK routine, which uses a Krylov subspace to solve the
	 linear system (J*s = F) instead of (J*s = -F). The reason we solve for the positive gradient
	 with the linear system is that our residual vector F cannot be passed to the subspace methods
	 as -F. Therefore, we solve for the positive gradient to maintain efficiency, then just subtract
	 the positive gradient as opposed to adding the negative gradient. Mathematically, it is exactly
	 the same.
	 
	 */
	
	do
	{
		xkp1 = backtrack_dat->xk - (pk*lambda);
		success = (*feval) (xkp1,Fkp1,feval_data);
		backtrack_dat->fun_call++;
		if (success != 0) {mError(simulation_fail); return -1;}
		backtrack_dat->normFkp1 = Fkp1.norm();
		norm_old = backtrack_dat->normFkp1;
		
		if ( isnan(backtrack_dat->normFkp1) || isinf(backtrack_dat->normFkp1) )
		{
			lambda = lambda_old * 0.1;
			lambda_old = lambda;
		}
		
	} while ((isnan(backtrack_dat->normFkp1) && lambda > backtrack_dat->lambdaMin) ||
			 (isinf(backtrack_dat->normFkp1) && lambda > backtrack_dat->lambdaMin));
	
	//Check for imediate success
	if (backtrack_dat->normFkp1 <= normFk + (backtrack_dat->alpha * lambda * slope))
	{
		return 0;
	}
	//Begin to reduce the step size
	else
	{
		//Form first new lambda and new solution and norm based on that lambda
		if (backtrack_dat->constRho == false)
			lambda_temp = -slope / (2.0 * (backtrack_dat->normFkp1 - normFk - slope));
		else
			lambda_temp = backtrack_dat->rho * lambda;
		if (lambda_temp < 0.1*lambda)
			lambda = 0.1 * lambda;
		else
			lambda = lambda_temp;
		if (isnan(lambda) || isinf(lambda))
		{
			return 1;
		}
		
		do
		{
			xkp1 = backtrack_dat->xk - (pk*lambda);
			success = (*feval) (xkp1,Fkp1,feval_data);
			backtrack_dat->fun_call++;
			if (success != 0) {mError(simulation_fail); return -1;}
			backtrack_dat->normFkp1 = Fkp1.norm();
			
			if ( isnan(backtrack_dat->normFkp1) || isinf(backtrack_dat->normFkp1) )
			{
				lambda = lambda_old * 0.1;
				lambda_old = lambda;
			}
			
		} while ((isnan(backtrack_dat->normFkp1) && lambda > backtrack_dat->lambdaMin) ||
				 (isinf(backtrack_dat->normFkp1) && lambda > backtrack_dat->lambdaMin));
		
		//Loop until we have a suitable step size
		while (backtrack_dat->normFkp1 > normFk + (backtrack_dat->alpha * lambda * slope))
		{
			//Form new lambda
			if (backtrack_dat->constRho == false)
			{
				a = (1.0 / (lambda - lambda_old)) * ( (backtrack_dat->normFkp1/pow(lambda,2.0)) - (normFk/pow(lambda,2.0)) - (slope/lambda) - (norm_old/pow(lambda_old,2.0)) + (normFk/pow(lambda_old,2.0)) + (slope/lambda_old) );
				
				b = (1.0 / (lambda - lambda_old)) * ( -((backtrack_dat->normFkp1*lambda_old)/pow(lambda,2.0)) + ((normFk*lambda_old)/pow(lambda,2.0)) + ((slope*lambda_old)/lambda) + ((norm_old*lambda)/pow(lambda_old,2.0)) - ((normFk*lambda)/pow(lambda_old,2.0)) - ((slope*lambda)/lambda_old) );
				
				disc = pow(b,2.0) - (3.0*a*slope);
				if (disc < 0.0)
					disc = 0.0;
				
				//Prepare for next iteration
				lambda_old = lambda;
				norm_old = backtrack_dat->normFkp1;
				
				//Check to see if form is cubic or quadratic
				if (a == 0.0)
					lambda_temp = -slope / (2.0*b);
				else
					lambda_temp = (-b + sqrt(disc)) / (3.0*a);
				
				//Check lambda_temp for violations in size
				if (lambda_temp > 0.5 * lambda)
					lambda_temp	= 0.5 * lambda;
			}
			else
			{
				lambda_temp = backtrack_dat->rho * lambda;
			}
			if (lambda_temp < 0.1*lambda || isnan(lambda_temp))
				lambda = 0.1 * lambda;
			else
				lambda = lambda_temp;
			
			//Check the value of lambda against the minimum
			if (lambda <= backtrack_dat->lambdaMin)
			{
				lambda = backtrack_dat->lambdaMin;
				xkp1 = backtrack_dat->xk - (pk*lambda);
				success = (*feval) (xkp1,Fkp1,feval_data);
				backtrack_dat->fun_call++;
				if (success != 0) {mError(simulation_fail); return -1;}
				backtrack_dat->normFkp1 = Fkp1.norm();
				success = 1;
				break;
			}
			
			//Form new solution and norm
			xkp1 = backtrack_dat->xk - (pk*lambda);
			success = (*feval) (xkp1,Fkp1,feval_data);
			backtrack_dat->fun_call++;
			if (success != 0) {mError(simulation_fail); return -1;}
			backtrack_dat->normFkp1 = Fkp1.norm();
		}
	}
	
	return success;
}

//Implementation of the PJFNK method for soling non-linear systems of equations
int pjfnk( int (*res) (const Matrix<double>& x, Matrix<double> &F, const void *data),
		  int (*precon) (const Matrix<double>& r, Matrix<double>& p, const void *data),
		  Matrix<double> &x, PJFNK_DATA *pjfnk_dat, const void *res_data,
		  const void *precon_data )
{
	int success = 0;
	double res_old;
	bool LS_Flag = pjfnk_dat->LineSearch;
	bool LS_User = pjfnk_dat->LineSearch;
	int bounce = 0;
	
	//Check input arguments for problems
	if ( (*res) == NULL )
	{
		mError(nullptr_func);
		return -1;
	}
	else
	{
		pjfnk_dat->funeval = (*res);
	}
	pjfnk_dat->precon = (*precon);
	if (x.rows() < 2)
	{
		success = -1;
		mError(matrix_too_small);
		return success;
	}
	pjfnk_dat->res_data = res_data;
	pjfnk_dat->precon_data = precon_data;
	pjfnk_dat->nl_iter = 0;
	pjfnk_dat->l_iter = 0;
	pjfnk_dat->fun_call = 0;
	
	//Check and initialize all PJFNK data
	if (pjfnk_dat->F.rows() != x.rows())
		pjfnk_dat->F.set_size(x.rows(), 1);
	if (pjfnk_dat->Fv.rows() != x.rows())
		pjfnk_dat->Fv.set_size(x.rows(), 1);
	if (pjfnk_dat->v.rows() != x.rows())
		pjfnk_dat->v.set_size(x.rows(), 1);
	if (pjfnk_dat->x.rows() != x.rows())
		pjfnk_dat->x.set_size(x.rows(), 1);
	if (pjfnk_dat->bestx.rows() != x.rows())
		pjfnk_dat->bestx.set_size(x.rows(), 1);
	for (int i=0; i<x.rows(); i++)
	{
		pjfnk_dat->x.edit(i,0,x(i,0));
		pjfnk_dat->bestx.edit(i,0,x(i,0));
	}
	if (pjfnk_dat->nl_maxit <= 0)
		pjfnk_dat->nl_maxit = std::min(3*x.rows(),1000);
	if (pjfnk_dat->linear_solver < GMRESLP || pjfnk_dat->linear_solver > QR)
	{
		//Choose the best linear solver based on problem size and availability of preconditioning
		if (x.rows() >= 100 && (*precon) == NULL)
			pjfnk_dat->linear_solver = GMRESRP; //Restarted GMRES with Right Preconditioning
		else if (x.rows() >= 100 && (*precon) != NULL)
			pjfnk_dat->linear_solver = BiCGSTAB; //BiCGSTAB
		else if (x.rows() < 100 && (*precon) != NULL)
			pjfnk_dat->linear_solver = CGS; //CGS
		else
			pjfnk_dat->linear_solver = GMRESLP; //Restarted GMRES with Left Preconditioning
		
		//NOTE: you would only use PCG if the linear system were symmetric
		//Symmetry is not something we can check for, so it never defaults to PCG
		
		//NOTE: you would specify linear_solver = FOM if you wanted to take an exact
		//Newton step. Otherwise, all steps are only approximate.
	}
	if (pjfnk_dat->nl_tol_abs >= 1.0)
		pjfnk_dat->nl_tol_abs = 1e-6;
	if (pjfnk_dat->nl_tol_abs < MIN_TOL)
		pjfnk_dat->nl_tol_abs = MIN_TOL;
	if (pjfnk_dat->nl_tol_rel >= 1.0)
		pjfnk_dat->nl_tol_rel = 1e-6;
	if (pjfnk_dat->nl_tol_rel < MIN_TOL)
		pjfnk_dat->nl_tol_rel = MIN_TOL;
	if (pjfnk_dat->lin_tol_rel >= 1.0)
		pjfnk_dat->lin_tol_rel = 0.1;
	if (pjfnk_dat->lin_tol_rel < MIN_TOL)
		pjfnk_dat->lin_tol_rel = MIN_TOL;
	if (pjfnk_dat->lin_tol_abs >= 1.0)
		pjfnk_dat->lin_tol_abs = 0.1;
	if (pjfnk_dat->lin_tol_abs < MIN_TOL)
		pjfnk_dat->lin_tol_abs = MIN_TOL;
	if (pjfnk_dat->eps >= 1.0 || pjfnk_dat->eps < sqrt(DBL_EPSILON))
		pjfnk_dat->eps = sqrt(DBL_EPSILON);
	if (pjfnk_dat->Bounce == true)
		pjfnk_dat->LineSearch = true;
	
	//Start the method by calling the users F(x) function to form an initial residual
	success = (*pjfnk_dat->funeval) (pjfnk_dat->x,pjfnk_dat->F,pjfnk_dat->res_data);
	pjfnk_dat->fun_call++;
	if (success != 0) {mError(simulation_fail); return -1;}
	pjfnk_dat->nl_res = pjfnk_dat->F.norm();
	pjfnk_dat->nl_bestres = pjfnk_dat->nl_res;
	pjfnk_dat->nl_res_base = pjfnk_dat->nl_res;
	pjfnk_dat->nl_relres = 1.0;
	res_old = pjfnk_dat->nl_res;
	
	if (pjfnk_dat->NL_Output == true)
	{
		std::cout << "Begin PJFNK iterations with ";
		if (pjfnk_dat->linear_solver == GMRESLP)
			std::cout << "Restarted Left-Preconditioned GMRES linear solver";
		else if (pjfnk_dat->linear_solver == PCG)
			std::cout << "PCG linear solver";
		else if (pjfnk_dat->linear_solver == BiCGSTAB)
			std::cout << "BiCGSTAB linear solver";
		else if (pjfnk_dat->linear_solver == CGS)
			std::cout << "CGS linear solver";
		else if (pjfnk_dat->linear_solver == FOM)
			std::cout << "Full Orthogonalization Method linear solver";
		else if (pjfnk_dat->linear_solver == GMRESRP)
			std::cout << "Restarted Right-Preconditioned GMRES linear solver";
		else if (pjfnk_dat->linear_solver == GCR)
			std::cout << "GCR linear solver";
		else if (pjfnk_dat->linear_solver == GMRESR)
			std::cout << "GMRESR linear solver";
		else if (pjfnk_dat->linear_solver == KMS)
			std::cout << "KMS linear solver";
		else if (pjfnk_dat->linear_solver == QR)
			std::cout << "QR linear solver";
		else
			return -1;
		if (pjfnk_dat->LineSearch == true)
			std::cout << " and Backtracking Linesearch...\n" << std::endl;
		else
			std::cout << "...\n" << std::endl;
		std::cout << "NL_Norm[" << 0 << "] =\t" << pjfnk_dat->nl_res << std::endl;
	}
	//Start looping for convergence
	for (int k=0; k<pjfnk_dat->nl_maxit; k++)
	{
		//Solve: J(x_k)*s_k = F(x_k) for s_k using a Krylov Subspace
		//Next iterate: x_kp1 = x_k - s_k
		//Check New Residual: F(x_kp1) then repeat
		pjfnk_dat->nl_iter++;
		
		//Restarted GMRES with Left Preconditioning - Best Convergence for any non-singluar matrix
		if (pjfnk_dat->linear_solver == GMRESLP)
		{
			pjfnk_dat->gmreslp_dat.Output = pjfnk_dat->L_Output;
			pjfnk_dat->gmreslp_dat.tol_abs = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->gmreslp_dat.tol_rel = pjfnk_dat->lin_tol_rel;
			pjfnk_dat->gmreslp_dat.maxit = x.rows();
			success = gmresLeftPreconditioned(jacvec, pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->gmreslp_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->gmreslp_dat.steps;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->gmreslp_dat.steps;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			//Form new solution
			if (pjfnk_dat->LineSearch == false)
			{
				if (success != 0) {mError(simulation_fail); success = -1; break;}
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->gmreslp_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->gmreslp_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//PCG - Convergence Guaranteed for symmetric, posivitive definite matrix
		else if (pjfnk_dat->linear_solver == PCG)
		{
			pjfnk_dat->pcg_dat.Output = pjfnk_dat->L_Output;
			pjfnk_dat->pcg_dat.tol_abs = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->pcg_dat.tol_rel = pjfnk_dat->lin_tol_rel;
			pjfnk_dat->pcg_dat.maxit = 2*x.rows();
			success = pcg(jacvec, pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->pcg_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->pcg_dat.iter;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->pcg_dat.iter;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			//Form new solution
			if (pjfnk_dat->LineSearch == false)
			{
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->pcg_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->pcg_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//BiCGSTAB - Good Convergence for non-symmetric, positive definite matrix
		else if (pjfnk_dat->linear_solver == BiCGSTAB)
		{
			pjfnk_dat->bicgstab_dat.Output = pjfnk_dat->L_Output;
			pjfnk_dat->bicgstab_dat.tol_abs = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->bicgstab_dat.tol_rel = pjfnk_dat->lin_tol_rel;
			pjfnk_dat->bicgstab_dat.maxit = 2*x.rows();
			success = bicgstab(jacvec, pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->bicgstab_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->bicgstab_dat.iter;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->bicgstab_dat.iter;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			if (pjfnk_dat->LineSearch == false)
			{
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->bicgstab_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->bicgstab_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//CGS - Good Convergence, if preconditioned, for non-symmetric, positive definite matix
		else if (pjfnk_dat->linear_solver == CGS)
		{
			pjfnk_dat->cgs_dat.Output = pjfnk_dat->L_Output;
			pjfnk_dat->cgs_dat.tol_abs = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->cgs_dat.tol_rel = pjfnk_dat->lin_tol_rel;
			pjfnk_dat->cgs_dat.maxit = 2*x.rows();
			success = cgs(jacvec, pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->cgs_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->cgs_dat.iter;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->cgs_dat.iter;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			if (pjfnk_dat->LineSearch == false)
			{
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->cgs_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->cgs_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//Unrestarted GMRES with Left Preconditioning (i.e. Full Orthogonalization) - Basically Equivalent to an Exact Newton Step
		else if (pjfnk_dat->linear_solver == FOM)
		{
			pjfnk_dat->gmreslp_dat.Output = pjfnk_dat->L_Output;
			pjfnk_dat->gmreslp_dat.tol_abs = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->gmreslp_dat.tol_rel = pjfnk_dat->lin_tol_rel;
			success = fom(jacvec, pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->gmreslp_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->gmreslp_dat.steps;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->gmreslp_dat.steps;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			if (pjfnk_dat->LineSearch == false)
			{
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->gmreslp_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->gmreslp_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//Restarted GMRES with Right Preconditioning - May be more efficient than Left Preconditioned GMRES for larger problems
		else if (pjfnk_dat->linear_solver == GMRESRP)
		{
			pjfnk_dat->gmresrp_dat.Output = pjfnk_dat->L_Output;
			pjfnk_dat->gmresrp_dat.tol_abs = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->gmresrp_dat.tol_rel = pjfnk_dat->lin_tol_rel;
			pjfnk_dat->gmresrp_dat.maxit = x.rows();
			success=gmresRightPreconditioned(jacvec, pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->gmresrp_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->gmresrp_dat.iter_total;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->gmresrp_dat.iter_total;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			if (pjfnk_dat->LineSearch == false)
			{
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->gmresrp_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->gmresrp_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//GCR iterative solver with restarting
		else if (pjfnk_dat->linear_solver == GCR)
		{
			pjfnk_dat->gcr_dat.Output = pjfnk_dat->L_Output;
			pjfnk_dat->gcr_dat.tol_abs = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->gcr_dat.tol_rel = pjfnk_dat->lin_tol_rel;
			pjfnk_dat->gcr_dat.maxit = x.rows();
			success = gcr(jacvec, pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->gcr_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->gcr_dat.total_iter;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->gcr_dat.total_iter;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			if (pjfnk_dat->LineSearch == false)
			{
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->gcr_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->gcr_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//GMRESR iterative solver with restarting
		else if (pjfnk_dat->linear_solver == GMRESR)
		{
			pjfnk_dat->gmresr_dat.GCR_Output = pjfnk_dat->L_Output;
			pjfnk_dat->gmresr_dat.GMRES_Output = false;
			pjfnk_dat->gmresr_dat.gcr_abs_tol = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->gmresr_dat.gcr_rel_tol = pjfnk_dat->lin_tol_rel;
			pjfnk_dat->gmresr_dat.gcr_maxit = x.rows();
			success = gmresr(jacvec, pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->gmresr_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->gmresr_dat.total_iter;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->gmresr_dat.total_iter;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			if (pjfnk_dat->LineSearch == false)
			{
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->gmresr_dat.gcr_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->gmresr_dat.gcr_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//KMS iterative solver with restarting
		else if (pjfnk_dat->linear_solver == KMS)
		{
			pjfnk_dat->kms_dat.Output_outer = pjfnk_dat->L_Output;
			pjfnk_dat->kms_dat.Output_inner = false;
			pjfnk_dat->kms_dat.outer_abstol = pjfnk_dat->lin_tol_abs;
			pjfnk_dat->kms_dat.outer_reltol = pjfnk_dat->lin_tol_rel;
			pjfnk_dat->kms_dat.maxit = x.rows();
			success = krylovMultiSpace(jacvec,pjfnk_dat->precon, pjfnk_dat->F, &pjfnk_dat->kms_dat, pjfnk_dat, pjfnk_dat->precon_data);
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + pjfnk_dat->kms_dat.total_iter;
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->kms_dat.total_iter;
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			if (pjfnk_dat->LineSearch == false)
			{
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->kms_dat.gmres_out.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->kms_dat.gmres_out.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
		}
		
		//QR direct solver
		else if (pjfnk_dat->linear_solver == QR)
		{
			success = QRsolve(jacvec, pjfnk_dat->F, &pjfnk_dat->qr_dat, pjfnk_dat);
			pjfnk_dat->fun_call = pjfnk_dat->fun_call + x.rows();
			pjfnk_dat->l_iter = pjfnk_dat->l_iter + x.rows();
			if (success != 0) {mError(simulation_fail); success = -1; break;}
			
			if (pjfnk_dat->LineSearch == false)
			{
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, pjfnk_dat->x(i,0) - pjfnk_dat->qr_dat.x(i,0));
					x.edit(i, 0, pjfnk_dat->x(i,0));
				}
			}
			else
			{
				success	= backtrackLineSearch(pjfnk_dat->funeval, pjfnk_dat->F, x, pjfnk_dat->qr_dat.x, pjfnk_dat->nl_res, &pjfnk_dat->backtrack_dat, pjfnk_dat->res_data);
				pjfnk_dat->fun_call = pjfnk_dat->fun_call + pjfnk_dat->backtrack_dat.fun_call;
				if (success < 0) {mError(simulation_fail); success = -1; break;}
				//Form new solution
				for (int i=0; i<x.rows(); i++)
				{
					pjfnk_dat->x.edit(i, 0, x(i,0));
				}
				if (pjfnk_dat->Bounce == true)
					LS_Flag = false;
			}
			
		}
		
		//No other Krylov Method currently available
		else
		{
			mError(invalid_boolean);
			return -1;
		}
		
		//Form the next non-linear residual
		if (pjfnk_dat->LineSearch == false)
		{
			success = (*pjfnk_dat->funeval) (pjfnk_dat->x,pjfnk_dat->F,pjfnk_dat->res_data);
			pjfnk_dat->fun_call++;
			if (success != 0) {mError(simulation_fail); return -1;}
			pjfnk_dat->nl_res = pjfnk_dat->F.norm();
			
			if (LS_User == true)
			{
				pjfnk_dat->LineSearch = true;
				LS_Flag = pjfnk_dat->LineSearch;
				bounce++;
				
				if (res_old <= pjfnk_dat->nl_res + pjfnk_dat->nl_tol_abs )
				{
					pjfnk_dat->LineSearch = false;
					LS_Flag = pjfnk_dat->LineSearch;
					bounce--;
				}
			}
			
		}
		else
		{
			pjfnk_dat->nl_res = pjfnk_dat->backtrack_dat.normFkp1;
			
			//Check for saddle point
			if (success == 1 && res_old <= pjfnk_dat->nl_res + pjfnk_dat->nl_tol_abs)
			{
				if (pjfnk_dat->nl_res <= pjfnk_dat->nl_tol_abs)
				{
					success = 0;
				}
				else
				{
					if ( (bounce >= 5 || pjfnk_dat->Bounce == false) && res_old <= pjfnk_dat->nl_res + pjfnk_dat->nl_tol_abs)
					{
						if (pjfnk_dat->NL_Output == true)
							std::cout << "\nPotential saddle point found! Ending algorithm early!" << std::endl;
						success = -1;
						break;
					}
					else
					{
						bounce++;
						pjfnk_dat->LineSearch = LS_Flag;
					}
				}
			}
		}
		pjfnk_dat->nl_relres = pjfnk_dat->nl_res / pjfnk_dat->nl_res_base;
		if (pjfnk_dat->NL_Output == true)
			std::cout << "NL_Norm[" << k+1 << "] =\t" << pjfnk_dat->nl_res << std::endl;
		
		//Check for convergence
		if (pjfnk_dat->nl_res <= pjfnk_dat->nl_tol_abs)
		{
			if (pjfnk_dat->NL_Output == true)
				std::cout << "\nSolution converged after " << k+1 << " PJFNK iterations within residual tolerance!" << std::endl;
			return success;
		}
		else if (pjfnk_dat->nl_relres <= pjfnk_dat->nl_tol_rel)
		{
			if (pjfnk_dat->NL_Output == true)
				std::cout << "\nSolution converged after " << k+1 << " PJFNK iterations within relative residual tolerance!" << std::endl;
			return success;
		}
		else
		{
			if (pjfnk_dat->nl_res < pjfnk_dat->nl_bestres)
			{
				pjfnk_dat->nl_bestres = pjfnk_dat->nl_res;
				pjfnk_dat->bestx = x;
			}
		}
		res_old = pjfnk_dat->nl_res;
		
		//Check for failures
		if (isnan(pjfnk_dat->nl_res) || isinf(pjfnk_dat->nl_res))
		{
			if (pjfnk_dat->NL_Output == true)
				std::cout << "\nPJFNK broke down at the " << k+1 << " iterate!" << std::endl;
			success = -1;
			break;
		}
	}
	
	//Post loop messages
	if (pjfnk_dat->nl_iter >= pjfnk_dat->nl_maxit)
	{
		if (pjfnk_dat->NL_Output == true)
		{
			std::cout << "\nPJFNK did not find a solution within " << pjfnk_dat->nl_maxit << " iterations!" << std::endl;
			std::cout << "Best reported norm =\t" << pjfnk_dat->nl_bestres << std::endl;
		}
		x = pjfnk_dat->bestx;
		pjfnk_dat->nl_res = pjfnk_dat->nl_bestres;
		success = 0;
	}
	if (success == -1)
	{
		if (pjfnk_dat->NL_Output == true)
		{
			std::cout << "\nPJFNK did not find exact solution within " << pjfnk_dat->nl_iter << " iterations!" << std::endl;
			std::cout << "Best reported norm =\t" << pjfnk_dat->nl_bestres << std::endl;
		}
		x = pjfnk_dat->bestx;
		pjfnk_dat->nl_res = pjfnk_dat->nl_bestres;
		success = 0;
	}
	
	return success;
}

//Form a numerical jacobian matrix
int NumericalJacobian( int (*Func) (const Matrix<double> &x, Matrix<double> &F, const void *user_data),
					  const Matrix<double> &x, Matrix<double> &J, int Nx, int Nf, NUM_JAC_DATA *jac_dat,
					  const void *user_data)
{
	int success = 0;
	
	//Check input arguments for problems
	if ( (*Func) == NULL )
	{
		mError(nullptr_func);
		return -1;
	}
	if ( jac_dat == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	if (jac_dat->dxj.rows() != Nx)
	{
		jac_dat->dxj.set_size(Nx, 1);
	}
	if (J.rows() != Nf && J.columns() != Nx)
	{
		J.set_size(Nf, Nx);
	}
	if (jac_dat->Fx.rows() != Nf)
	{
		jac_dat->Fx.set_size(Nf, 1);
	}
	if (jac_dat->Fxp.rows() != Nf)
	{
		jac_dat->Fxp.set_size(Nf, 1);
	}
	if (jac_dat->eps < sqrt(DBL_EPSILON) || jac_dat->eps >= 1.0)
	{
		jac_dat->eps = sqrt(DBL_EPSILON);
	}
	
	//Form the first fuction evaluation
	success = (*Func) (x, jac_dat->Fx, user_data);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Create a copy of x to change piecewise
	jac_dat->dxj = x;
	for (int j=0; j<Nx; j++)
	{
		//Change the jth variable
		jac_dat->dxj(j,0) = x(j,0) + jac_dat->eps;
		success = (*Func) (jac_dat->dxj, jac_dat->Fxp, user_data);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Approximate each row of the first column of the jacobian
		for (int i=0; i<Nf; i++)
		{
			J(i,j) = (jac_dat->Fxp(i,0) - jac_dat->Fx(i,0)) / jac_dat->eps;
		}
		
		//Recover the jth variable before continuing
		jac_dat->dxj(j,0) = x(j,0);
	}
	
	return success;
}

//Testing grounds for the LARK Functions
int LARK_TESTS()
{
	int success = 0;
	
	//Example 1:------------------------------------------------------------------------
	//Building Kyrlov Subspace from linear system
	std::cout << "------------------Begin Example 1------------------" << std::endl;
	EX01_DATA ex01_dat;
	int rows = 100;
	int cols = 100;
	int maxk = 100;
	double bound = -1.0;
	double time;
	ex01_dat.M.set_size(rows, cols);
	ex01_dat.b.set_size(rows, 1);
	ex01_dat.M.tridiagonalFill(-1, 2, -1, false);
	ex01_dat.b.dirichletBCFill(0, 1, bound);
	
	Matrix<double> solution(rows,1), guess;
	guess.set_size(rows,1);
	time = clock();
	solution.ladshawSolve(ex01_dat.M,ex01_dat.b);
	time = clock() - time;
	std::cout << "Direct Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "Norm =\t" << (ex01_dat.b - ex01_dat.M*solution).norm() << std::endl; //(PASS)
	
	Matrix<double> r0(rows,1);
	r0 = ex01_dat.b - ex01_dat.M*guess;
	ARNOLDI_DATA arnoldi_dat;
	arnoldi_dat.k = maxk;
	
	time = clock();
	success = arnoldi(matvec_ex01,NULL,r0,&arnoldi_dat,(void *)&ex01_dat, (void *)&ex01_dat);
	
	Matrix<double> x;
	success = update_arnoldi_solution(x, guess, &arnoldi_dat);
	time = clock() - time;
	std::cout << "Krylov Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;;
	std::cout << "Norm =\t" << (ex01_dat.b - ex01_dat.M*x).norm() << std::endl; //(PASS)
	GMRESLP_DATA fom_dat;
	time = clock();
	success = fom(matvec_ex01,precon_ex01,ex01_dat.b,&fom_dat,(void *)&ex01_dat,(void *)&ex01_dat);
	time = clock() - time;
	std::cout << "FULL GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "FULL GMRES Norm: " << (ex01_dat.b - ex01_dat.M*fom_dat.x).norm() << std::endl;
	std::cout << "------------------END Example 1------------------\n" << std::endl;
	//END Example 1:-------------------------------------------------------------------
	
	//Example 2: ----------------------------------------------------------------------
	/*
	 Solve 2x2 system defined in Saad and Schultz (1986) paper
	 
	 Ax = b -> A = |0  1|		b = |1|		x0 = |0|
	 |-1 0|		    |1|			 |0|
	 
	 Solution: x = |-1|
	 | 1|
	 
	 Any other solution technique will break down for this system because the diagonals
	 of A are all zero. However, the system is not singular, it is just order strangely.
	 This demonstrates that the Krylov method can be used to solve linear systems directly,
	 where other methods may break down and fail to find a solution.
	 
	 */
	std::cout << "------------------Begin Example 2------------------" << std::endl;
	EX02_DATA ex02_dat;
	ex02_dat.M.set_size(2, 2);
	ex02_dat.b.set_size(2, 1);
	Matrix<double> x0(2,1);
	ex02_dat.M.edit(0, 0, 0); ex02_dat.M.edit(0, 1, 1);
	ex02_dat.M.edit(1, 0, -1); ex02_dat.M.edit(1, 1, 0);
	ex02_dat.b.ConstantICFill(1);
	x0.ConstantICFill(0);
	r0 = ex02_dat.b - ex02_dat.M*x0;
	arnoldi_dat.k = 2;
	arnoldi_dat.iter = 0;
	time = clock();
	success = arnoldi(matvec_ex02,NULL,r0,&arnoldi_dat,(void *)&ex02_dat,(void *)&ex02_dat);
	std::cout << "Krylov Exit Code: " << success << std::endl;
	success = update_arnoldi_solution(x, x0, &arnoldi_dat);
	time = clock() - time;
	std::cout << "Krylov Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	x.Display("x");
	std::cout << "Norm =\t" << (ex02_dat.b - ex02_dat.M*x).norm() << std::endl; //(PASS)
	std::cout << "------------------END Example 2------------------\n" << std::endl;
	//END Example 2:-------------------------------------------------------------------
	
	//Example 3:-------------------------------------------------------------------------
	/*
	 
	 Solve Example 1 iteratively using the default settings for GMRES
	 
	 */
	
	std::cout << "------------------Begin Example 3------------------" << std::endl;
	std::cout << "Solving Example 1 with default GMRES arguments" << std::endl;
	GMRESLP_DATA gmres_dat;
	gmres_dat.restart = 20;
	time = clock();
	success = gmresLeftPreconditioned(matvec_ex01,precon_ex01,ex01_dat.b,&gmres_dat,(void *)&ex01_dat,(void *)&ex01_dat); //(PASS)
	time = clock() - time;
	std::cout << "GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GMRES Norm: " << (ex01_dat.b - ex01_dat.M*gmres_dat.x).norm() << std::endl;
	std::cout << "GMRES Total Evaluations: " << gmres_dat.steps << std::endl;
	std::cout << "------------------END Example 3------------------\n" << std::endl;
	//END Example 3:-------------------------------------------------------------------------
	
	
	//Example 4: -------------------------------------------------------------------------------
	/*
	 We want to see how well the GMRES implementation does at solving a more complex, sparse
	 system of equations. One in which none of our current methods, aside from FOM, would be
	 able to solve directly.
	 
	 NOTE: if I am going to use DENSE matrices with these codes, I should not use a size
	 greater than 1000x1000. For sizes larger than this, you must use either Sparse matrices
	 or some kind of matrix free format for evaluating the functions.
	 */
	std::cout << "------------------Begin Example 4------------------" << std::endl;
	std::cout << "Solving a 3D Laplacian function with default GMRES" << std::endl;
	EX04_DATA ex04_dat;
	int m = 5;
	ex04_dat.M.naturalLaplacian3D(m); //Creates a mxmxm matrix for a 3D Laplacian function
	ex04_dat.b.set_size(m*m*m, 1);
	ex04_dat.b.edit(0, 0, -1);
	ex04_dat.b.edit(m*m-1, 0, -1);
	ex04_dat.b.edit(m*m*m-1, 0, -1);
	GMRESLP_DATA gmres_dat4;
	gmres_dat4.restart = 20;
	
	time = clock();
	success = gmresLeftPreconditioned(matvec_ex04,precon_ex04,ex04_dat.b,&gmres_dat4,(void *)&ex04_dat,(void *)&ex04_dat);
	time = clock() - time;
	x.ladshawSolve(ex04_dat.M, ex04_dat.b);
	std::cout << "Precon Norm: " << (ex04_dat.b - ex04_dat.M*x).norm() << std::endl;
	std::cout << "GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GMRES Norm: " << (ex04_dat.b - ex04_dat.M*gmres_dat4.x).norm() << std::endl;
	std::cout << "GMRES Total Evaluations: " << gmres_dat4.steps << std::endl;
	std::cout << "------------------END Example 4------------------\n" << std::endl;
	
	//END Example 4:----------------------------------------------------------------------------
	
	
	//Example 5: --------------------------------------------------------------------------------
	/*
	 Solve example 4 using PCG algorithm and compare computation time to that of GMRES.
	 Note that example 4 was solving a symmetric system with GMRES and performed very
	 well. Here, we expect PCG to perform better because it requires much less memory
	 and fewer matrix-vector multiplications
	 
	 This method did indeed outperform GMRES for the same symmetric problem!
	 */
	std::cout << "------------------Begin Example 5------------------" << std::endl;
	std::cout << "Solving a 3D Laplacian function with default PCG" << std::endl;
	PCG_DATA pcg_dat5;
	time = clock();
	success = pcg(matvec_ex04,precon_ex04,ex04_dat.b,&pcg_dat5,(void *)&ex04_dat,(void *)&ex04_dat);
	//success = pcg(matvec_ex04,NULL,ex04_dat.b,&pcg_dat5,(void *)&ex04_dat,(void *)&ex04_dat);
	time = clock() - time;
	std::cout << "PCG Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "PCG Norm: " << (ex04_dat.b - ex04_dat.M*pcg_dat5.x).norm() << std::endl;
	std::cout << "PCG Iterations: " << pcg_dat5.iter << std::endl;
	std::cout << "------------------END Example 5------------------\n" << std::endl;
	
	//END Example 5:----------------------------------------------------------------------------
	
	//Example 6: --------------------------------------------------------------------------------
	/*
	 Attempt to solve example 1 with PCG. Example 1 is a non-symmetric matrix, so we do
	 not anticipate this to converge to a reasonable solution. However, we want to see
	 what happens in the method when the matrix is non-symmetric. Does the residual
	 improve? What is the behavior like?
	 
	 The norm was reduced, but it never converged! The behavior of the norm was not
	 very good. Was noisy and unstable, which may contribute to non-convergence. For
	 examples 7 and 8, we will try to solve this same problem, but using BiCGSTAB.
	 */
	std::cout << "------------------Begin Example 6------------------" << std::endl;
	std::cout << "Solving a 1D Non-symmetric Laplacian function with default PCG" << std::endl;
	PCG_DATA pcg_dat6;
	time = clock();
	success = pcg(matvec_ex01,precon_ex01,ex01_dat.b,&pcg_dat6,(void *)&ex01_dat,(void *)&ex01_dat);
	time = clock() - time;
	std::cout << "PCG Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "PCG Norm: " << (ex01_dat.b - ex01_dat.M*pcg_dat6.x).norm() << std::endl;
	std::cout << "PCG Iterations: " << pcg_dat6.iter << std::endl;
	std::cout << "------------------END Example 6------------------\n" << std::endl;
	
	//END Example 6:----------------------------------------------------------------------------
	
	
	//Example 7: --------------------------------------------------------------------------------
	/*
	 Solve the non-symmetric 1-D Laplacian with BiCGSTAB algorithm. Previous example tried
	 and failed to solve this system with PCG. Example 3 demonstrated the the GMRES
	 implementation was able to solve this problem and reduce the residual monotonically.
	 Here, we expect BiCGSTAB to be able to solve the problem, but expect the residuals
	 to not decrease monotonically. Whether or not this is more computationally efficient
	 than GMRES will really depend on the problem we are solving.
	 
	 For this particular problem, BiCGSTAB too more iterations to converge than GMRES, but
	 solved the system faster than GMRES. This is because BiCGSTAB is a much more
	 computationally efficient algorithm that GMRES due to the recurrences in the code.
	 */
	std::cout << "------------------Begin Example 7------------------" << std::endl;
	std::cout << "Solving a 1D Non-symmetric Laplacian function with default BiCGSTAB" << std::endl;
	BiCGSTAB_DATA bicg_dat7;
	time = clock();
	success = bicgstab(matvec_ex01,precon_ex01,ex01_dat.b,&bicg_dat7,(void *)&ex01_dat,(void *)&ex01_dat);
	time = clock() - time;
	std::cout << "BiCGSTAB Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "BiCGSTAB Norm: " << (ex01_dat.b - ex01_dat.M*bicg_dat7.x).norm() << std::endl;
	std::cout << "BiCGSTAB Iterations: " << bicg_dat7.iter << std::endl;
	std::cout << "------------------END Example 7------------------\n" << std::endl;
	
	//END Example 7:----------------------------------------------------------------------------
	
	
	//Example 8: --------------------------------------------------------------------------------
	/*
	 Solve the non-symmetric 1-D Laplacian with the CGS algorithm. This is another
	 iterative method for solving non-symmetric linear systems that has often been
	 implemented by engineers and scientists. In theory, BiCGSTAB should be more
	 stable than CGS, but for some problems CGS will converge faster than BiCGSTAB.
	 
	 Case in point, for this example, CGS converged in 100 iterations whereas
	 BiCGSTAB converged in 110. However, the convergence behavior of BiCGSTAB was
	 much, much smoother than that of this CGS implementation.
	 
	 That is of course if we consider using preconditioning. Without preconditioning,
	 CGS for this problem diverges and will not converge in 200 iterations, while
	 BiCGSTAB still converges, but takes 159 iterations.
	 
	 The take away from this is that if a preconditioner is available, use CGS over
	 BiCGSTAB. However, if no preconditioner is available, then use BiCGSTAB. If the
	 linear system is not well conditioned, or if CGS and BiCGSTAB fail, then use
	 GMRES. If the system is small enough, just use FOM. For any symmetric matrix,
	 always use PCG. This will be a convention we adopt moving forward to developing
	 our implementation of PJFNK for non-linear systems.
	 */
	std::cout << "------------------Begin Example 8------------------" << std::endl;
	
	//NOTE: This also demonstrates use of pointers to data structures embeded in other structures
	PJFNK_DATA test;
	test.res_data = (void *)&ex01_dat;
	
	std::cout << "Solving a 1D Non-symmetric Laplacian function with default CGS" << std::endl;
	CGS_DATA cgs_dat8;
	time = clock();
	success = cgs(matvec_ex01,precon_ex01,ex01_dat.b,&cgs_dat8,test.res_data,(void *)&ex01_dat);
	time = clock() - time;
	std::cout << "CGS Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "CGS Norm: " << (ex01_dat.b - ex01_dat.M*cgs_dat8.x).norm() << std::endl;
	std::cout << "CGS Iterations: " << cgs_dat8.iter << std::endl;
	std::cout << "------------------END Example 8------------------\n" << std::endl;
	
	//END Example 8:----------------------------------------------------------------------------
	
	
	//Example 9 and 10: ------------------------------------------------------------------------
	
	/*
	 This example shows the use of both Picard and PJFNK methods for solving non-linear
	 systems of equations. In this example, the non-linear system arises from a non-linear
	 PDE of a 1-D Laplacian with a non-linear term k*exp(u);
	 
	 Solve: 		-d^2u/dx^x - k*exp(u) = 0
	 
	 F(u) = 0 = -u(i+1) + 2*u(i) - u(i-1) - h^2*k*exp(u(i))
	 BCs:	u(0) = u(1) = 0
	 
	 Split into linear and non-linear parts for both the Picard iteration and as a
	 preconditioner for the PJFNK method and compare the results.
	 
	 Note: Picard could not converge for this problem. May be too strongly non-linear...
	 
	 Was able to get Picard to converge after 291 Non-linear iterations, but problem size
	 is of order 10. Very poor convergence for this problem.
	 
	 */
	
	std::cout << "------------------Begin Example 9------------------" << std::endl;
	
	EX09_DATA ex09_dat;
	ex09_dat.k = 1.0;
	ex09_dat.N = 12;
	ex09_dat.h = 1.0/ex09_dat.N;
	ex09_dat.N = ex09_dat.N-2;
	ex09_dat.x.set_size(ex09_dat.N, 1);
	ex09_dat.s.set_size(ex09_dat.N, 1);
	ex09_dat.p.set_size(ex09_dat.N, 1);
	ex09_dat.M.set_size(ex09_dat.N, ex09_dat.N);
	ex09_dat.M.tridiagonalFill(-1, 2, -1, false);
	
	//First try to solve with Picard's Method
	PICARD_DATA picard_dat09;
	picard_dat09.maxit = 300;
	time = clock();
	success = picard(funeval_ex09, evalx_ex09, ex09_dat.x, &picard_dat09, (void *)&ex09_dat, (void *)&ex09_dat);
	time = clock() - time;
	success = funeval_ex09(ex09_dat.x, ex09_dat.p, (void *)&ex09_dat);
	std::cout << "PICARD Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "PICARD Norm: " << ex09_dat.p.norm() << std::endl;
	std::cout << "PICARD Iterations: " << picard_dat09.iter << std::endl;
	
	//std::cout << std::endl;
	//ex09_dat.x.Display("x");
	
	std::cout << "------------------END Example 9------------------\n" << std::endl;
	
	std::cout << "------------------Begin Example 10------------------" << std::endl;
	
	ex09_dat.k = 1.0;
	ex09_dat.N = 102;
	ex09_dat.h = 1.0/ex09_dat.N;
	ex09_dat.N = ex09_dat.N-2;
	ex09_dat.x.set_size(ex09_dat.N, 1);
	ex09_dat.s.set_size(ex09_dat.N, 1);
	ex09_dat.p.set_size(ex09_dat.N, 1);
	ex09_dat.M.set_size(ex09_dat.N, ex09_dat.N);
	ex09_dat.M.tridiagonalFill(-1, 2, -1, false);
	
	//Next try to solve with PJFNK method
	ex09_dat.x.ConstantICFill(0.0);
	PJFNK_DATA pjfnk_dat10;
	pjfnk_dat10.linear_solver = GMRESRP;
	//pjfnk_dat10.L_Output = true;
	//pjfnk_dat10.gmres_dat.Expand = true;
	
	//NOTE: You do not have to set these, they will be automatically set when PJFNK is called
	//pjfnk_dat10.linear_solver = GMRESLP;
	//pjfnk_dat10.lin_tol = 1e-4;
	//pjfnk_dat10.res_data = (void *)&ex09_dat;		//Required if your functions need your data
	//pjfnk_dat10.precon_data = (void *)&ex09_dat;	//Required if your preconditioner needs data
	//pjfnk_dat10.funeval = (*funeval_ex10);			//Always required
	//pjfnk_dat10.precon = (*precon_ex10);			//Optional, but recommended
	
	time = clock();
	//success = pjfnk(funeval_ex10, precon_ex10, ex09_dat.x, &pjfnk_dat10, (void *)&ex09_dat, (void *)&ex09_dat);
	success = pjfnk(funeval_ex10, NULL, ex09_dat.x, &pjfnk_dat10, (void *)&ex09_dat, NULL);
	time = clock() - time;
	success = funeval_ex10(ex09_dat.x, ex09_dat.p, (void *)&ex09_dat);
	std::cout << "PJFNK Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "PJFNK Norm: " << ex09_dat.p.norm() << std::endl;
	std::cout << "PJFNK Iterations: " << pjfnk_dat10.nl_iter << std::endl;
	
	//std::cout << std::endl;
	//ex09_dat.x.Display("x");
	
	std::cout << "------------------END Example 10------------------\n" << std::endl;
	
	//END Example 9 and 10: --------------------------------------------------------------------
	
	//Example 11: ----------------------------------------------------------------------
	/*
	 Solve 2x2 system defined in Saad and Schultz (1986) paper
	 
	 Ax = b -> A = |0  1|		b = |1|		x0 = |0|
	 |-1 0|		    |1|			 |0|
	 
	 Solution: x = |-1|
	 | 1|
	 
	 Any other solution technique will break down for this system because the diagonals
	 of A are all zero. However, the system is not singular, it is just order strangely.
	 This demonstrates that the Krylov method can be used to solve linear systems directly,
	 where other methods may break down and fail to find a solution.
	 
	 */
	std::cout << "------------------Begin Example 11------------------" << std::endl;
	std::cout << "Solving a 2x2 indefinite system with default Right-Preconditioned GMRES" << std::endl;
	ex02_dat.M.set_size(2, 2);
	ex02_dat.b.set_size(2, 1);
	ex02_dat.M.edit(0, 0, 0); ex02_dat.M.edit(0, 1, 1);
	ex02_dat.M.edit(1, 0, -1); ex02_dat.M.edit(1, 1, 0);
	ex02_dat.b.ConstantICFill(1);
	x0.ConstantICFill(0);
	r0 = ex02_dat.b - ex02_dat.M*x0;
	
	GMRESRP_DATA cgmres_dat;
	//cgmres_dat.restart = 2;
	//cgmres_dat.maxit = 1;
	time = clock();
	
	success = gmresRightPreconditioned(matvec_ex02,NULL,ex02_dat.b,&cgmres_dat,(void *)&ex02_dat,(void *)&ex02_dat);
	std::cout << "GMRES Exit Code: " << success << std::endl;
	
	time = clock() - time;
	std::cout << "GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	cgmres_dat.x.Display("gmres_x");
	std::cout << "Norm =\t" << (ex02_dat.b - ex02_dat.M*cgmres_dat.x).norm() << std::endl; //(PASS)
	
	std::cout << "\nSolving a 2x2 indefinite system with default GCR" << std::endl;
	GCR_DATA gcr11;
	
	time = clock();
	
	success = gcr(matvec_ex02,NULL,ex02_dat.b,&gcr11,(void *)&ex02_dat,(void *)&ex02_dat);
	std::cout << "GCR Exit Code: " << success << std::endl;
	
	time = clock() - time;
	std::cout << "GCR Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GCR inner iterations: " << gcr11.iter_inner << std::endl;
	std::cout << "GCR outer iterations: " << gcr11.iter_outer << std::endl;
	std::cout << "GCR total iterations: " << gcr11.total_iter << std::endl;
	gcr11.x.Display("gcr_x");
	std::cout << "Norm =\t" << (ex02_dat.b - ex02_dat.M*gcr11.x).norm() << std::endl; //FAIL!!!
	
	
	std::cout << "\nSolving a 2x2 indefinite system with default GMRESR" << std::endl;
	
	GMRESR_DATA gmresr11;
	
	time = clock();
	
	success = gmresr(matvec_ex02,NULL,ex02_dat.b,&gmresr11,(void *)&ex02_dat,(void *)&ex02_dat);
	std::cout << "GCR Exit Code: " << success << std::endl;
	
	time = clock() - time;
	std::cout << "GCR Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GCR inner iterations: " << gmresr11.gcr_dat.iter_inner << std::endl;
	std::cout << "GCR outer iterations: " << gmresr11.gcr_dat.iter_outer << std::endl;
	std::cout << "GCR total iterations: " << gmresr11.gcr_dat.total_iter << std::endl;
	gmresr11.gcr_dat.x.Display("gmresr_x");
	std::cout << "Norm =\t" << (ex02_dat.b - ex02_dat.M*gmresr11.gcr_dat.x).norm() << std::endl; //FAIL!!!
	
	std::cout << "------------------END Example 11------------------\n" << std::endl;
	//END Example 11:-------------------------------------------------------------------
	
	//Example 12: -------------------------------------------------------------------------------
	/*
	 We want to see how well the Compact GMRES implementation does at solving a more complex, sparse
	 system of equations. One in which none of our current methods, aside from FOM, would be
	 able to solve directly.
	 
	 NOTE: if I am going to use DENSE matrices with these codes, I should not use a size
	 greater than 1000x1000. For sizes larger than this, you must use either Sparse matrices
	 or some kind of matrix free format for evaluating the functions.
	 */
	std::cout << "------------------Begin Example 12------------------" << std::endl;
	std::cout << "Solving a 3D Laplacian function with default Right-Preconditioned GMRES" << std::endl;
	m = 5;
	ex04_dat.M.naturalLaplacian3D(m); //Creates a mxmxm matrix for a 3D Laplacian function
	ex04_dat.b.set_size(m*m*m, 1);
	ex04_dat.b.edit(0, 0, -1);
	ex04_dat.b.edit(m*m-1, 0, -1);
	ex04_dat.b.edit(m*m*m-1, 0, -1);
	GMRESRP_DATA cgmres_dat4;
	//cgmres_dat4.restart = 5;
	
	time = clock();
	success = gmresRightPreconditioned(matvec_ex04,precon_ex04,ex04_dat.b,&cgmres_dat4,(void *)&ex04_dat,(void *)&ex04_dat);
	//success = gmresRightPreconditioned(matvec_ex04,NULL,ex04_dat.b,&cgmres_dat4,(void *)&ex04_dat,(void *)&ex04_dat);
	time = clock() - time;
	x.ladshawSolve(ex04_dat.M, ex04_dat.b);
	std::cout << "Precon Norm: " << (ex04_dat.b - ex04_dat.M*x).norm() << std::endl;
	std::cout << "GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GMRES Norm: " << (ex04_dat.b - ex04_dat.M*cgmres_dat4.x).norm() << std::endl;
	std::cout << "GMRES Total Evaluations: " << cgmres_dat4.iter_total << std::endl;
	std::cout << "GMRES Restarts: " << cgmres_dat4.iter_outer << std::endl;
	
	std::cout << "\nSolving a 3D Laplacian function with default GCR" << std::endl;
	GCR_DATA gcr12;
	//gcr12.restart = 5;
	time = clock();
	success = gcr(matvec_ex04,precon_ex04,ex04_dat.b,&gcr12,(void *)&ex04_dat,(void *)&ex04_dat);
	//success = gcr(matvec_ex04,NULL,ex04_dat.b,&gcr12,(void *)&ex04_dat,(void *)&ex04_dat);
	time = clock() - time;
	std::cout << "GCR Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GCR Norm: " << (ex04_dat.b - ex04_dat.M*gcr12.x).norm() << std::endl; //PASS
	std::cout << "GCR inner iterations: " << gcr12.iter_inner << std::endl;
	std::cout << "GCR outer iterations: " << gcr12.iter_outer << std::endl;
	std::cout << "GCR total iterations: " << gcr12.total_iter << std::endl;
	
	std::cout << "------------------END Example 12------------------\n" << std::endl;
	
	//END Example 12:----------------------------------------------------------------------------
	
	//Example 13: --------------------------------------------------------------------------------
	/*
	 Solve the non-symmetric 1-D Laplacian with CGMRES algorithm. Previous example tried
	 and failed to solve this system with PCG. Example 3 demonstrated the the GMRES
	 implementation was able to solve this problem and reduce the residual monotonically.
	 Here, we expect BiCGSTAB to be able to solve the problem, but expect the residuals
	 to not decrease monotonically. Whether or not this is more computationally efficient
	 than GMRES will really depend on the problem we are solving.
	 
	 For this particular problem, BiCGSTAB too more iterations to converge than GMRES, but
	 solved the system faster than GMRES. This is because BiCGSTAB is a much more
	 computationally efficient algorithm that GMRES due to the recurrences in the code.
	 */
	std::cout << "------------------Begin Example 13------------------" << std::endl;
	std::cout << "Solving a 1D Non-symmetric Laplacian function with default Right-Preconditioned GMRES" << std::endl;
	GMRESRP_DATA cgmres_dat13;
	cgmres_dat13.restart = ex01_dat.b.rows();
	cgmres_dat13.x.set_size(ex01_dat.b.rows(), 1);
	cgmres_dat13.restart = 100;
	time = clock();
	success = gmresRightPreconditioned(matvec_ex01,precon_ex01,ex01_dat.b,&cgmres_dat13,(void *)&ex01_dat,(void *)&ex01_dat);
	//success = gmresRightPreconditioned(matvec_ex01,NULL,ex01_dat.b,&cgmres_dat13,(void *)&ex01_dat,(void *)&ex01_dat);
	time = clock() - time;
	std::cout << "GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GMRES Norm: " << (ex01_dat.b - ex01_dat.M*cgmres_dat13.x).norm() << std::endl;
	std::cout << "GMRES Iterations: " << cgmres_dat13.iter_total << std::endl;
	std::cout << "GMRES Restarts: " << cgmres_dat13.iter_outer << std::endl;
	std::cout << "------------------END Example 13------------------\n" << std::endl;
	
	std::cout << "\nSolving a 1D Non-symmetric Laplacian function with default GCR" << std::endl;
	
	GCR_DATA gcr13;
	gcr13.restart = 100;
	gcr13.x.set_size(ex01_dat.b.rows(), 1);
	time = clock();
	success = gcr(matvec_ex01,precon_ex01,ex01_dat.b,&gcr13,(void *)&ex01_dat,(void *)&ex01_dat);
	//success = gcr(matvec_ex01,NULL,ex01_dat.b,&gcr13,(void *)&ex01_dat,(void *)&ex01_dat);
	time = clock() - time;
	std::cout << "GCR Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GCR Norm: " << (ex01_dat.b - ex01_dat.M*gcr13.x).norm() << std::endl; //PASS
	std::cout << "GCR inner iterations: " << gcr13.iter_inner << std::endl;
	std::cout << "GCR outer iterations: " << gcr13.iter_outer << std::endl;
	std::cout << "GCR total iterations: " << gcr13.total_iter << std::endl;
	
	std::cout << "\nSolving a 1D Non-symmetric Laplacian function with default GMRESR" << std::endl;
	
	GMRESR_DATA gmresr13;
	//gmresr13.gmres_restart = 10;
	time = clock();
	success = gmresr(matvec_ex01,precon_ex01,ex01_dat.b,&gmresr13,(void *)&ex01_dat,(void *)&ex01_dat);
	time = clock() - time;
	std::cout << "GMRESR Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "GMRESR Norm: " << (ex01_dat.b - ex01_dat.M*gmresr13.gcr_dat.x).norm() << std::endl; //PASS
	std::cout << "GMRESR inner iterations: " << gmresr13.gcr_dat.iter_inner << std::endl;
	std::cout << "GMRESR outer iterations: " << gmresr13.gcr_dat.iter_outer << std::endl;
	std::cout << "GMRESR total iterations: " << gmresr13.gcr_dat.total_iter << std::endl;
	
	//END Example 13:----------------------------------------------------------------------------
	
	std::cout << "------------------Begin Example 14------------------" << std::endl;
	
	ex09_dat.k = 2.0;
	ex09_dat.N = 102;
	ex09_dat.h = 1.0/ex09_dat.N;
	ex09_dat.N = ex09_dat.N-2;
	ex09_dat.x.set_size(ex09_dat.N, 1);
	ex09_dat.s.set_size(ex09_dat.N, 1);
	ex09_dat.p.set_size(ex09_dat.N, 1);
	ex09_dat.M.set_size(ex09_dat.N, ex09_dat.N);
	ex09_dat.M.tridiagonalFill(-1, 2, -1, false);
	
	//Next try to solve with PJFNK method
	ex09_dat.x.ConstantICFill(0.0);
	PJFNK_DATA pjfnk_dat14;
	
	//NOTE: You do not have to set these, they will be automatically set when PJFNK is called
	pjfnk_dat14.L_Output = true;
	pjfnk_dat14.linear_solver = KMS; //Choice 1 (PCG) and 4 (GMRES_FULL) are not choosen by PJFNK
	pjfnk_dat14.LineSearch = true;
	//pjfnk_dat14.Bounce = true;
	//pjfnk_dat14.backtrack_dat.constRho = true;
	
	//KMS Options
	pjfnk_dat14.kms_dat.max_level = 2;
	pjfnk_dat14.kms_dat.inner_reltol = 0.001;
	
	time = clock();
	success = pjfnk(funeval_ex10, precon_ex10, ex09_dat.x, &pjfnk_dat14, (void *)&ex09_dat, (void *)&ex09_dat);
	//success = pjfnk(funeval_ex10, NULL, ex09_dat.x, &pjfnk_dat14, (void *)&ex09_dat, (void *)&ex09_dat);
	time = clock() - time;
	success = funeval_ex10(ex09_dat.x, ex09_dat.p, (void *)&ex09_dat);
	std::cout << "PJFNK Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	std::cout << "PJFNK Norm: " << ex09_dat.p.norm() << std::endl;
	std::cout << "PJFNK Iterations: " << pjfnk_dat14.nl_iter << std::endl;
	
	//std::cout << std::endl;
	//ex09_dat.x.Display("x");
	
	std::cout << "------------------END Example 14------------------\n" << std::endl;
	
	
	
	
	
	std::cout << "------------------Begin Example 15------------------" << std::endl;
	
	/*
	 *			SOLVE a 3D Laplacian without using matrices... Sparse System Handling
	 */
	
	std::cout << "Solving a 125,000 Linear Equations with GMRESRP, PCG, BiCGSTAB, CGS, GMRESLP, GMRESR, and KMS" << std::endl;
	
	//Setup Data Structure for Sparse Laplacian
	EX15_DATA ex15_dat;
	ex15_dat.m =50;
	ex15_dat.N = ex15_dat.m*ex15_dat.m*ex15_dat.m;
	ex15_dat.b.set_size(ex15_dat.N, 1);
	ex15_dat.b.edit(0, 0, -1);
	ex15_dat.b.edit(ex15_dat.m*ex15_dat.m-1, 0, -1);
	ex15_dat.b.edit(ex15_dat.m*ex15_dat.m*ex15_dat.m-1, 0, -1);
	x.set_size(ex15_dat.N, 1);
	
	std::cout << "\n------------------START GMRES------------------" << std::endl;
	
	//Setup GMRES
	GMRESRP_DATA gmres_dat15;
	
	time = clock();
	success = gmresRightPreconditioned(matvec_ex15,NULL,ex15_dat.b,&gmres_dat15,(void *)&ex15_dat,(void *)&ex15_dat);
	time = clock() - time;
	std::cout << "GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	success = matvec_ex15(gmres_dat15.x, x, (void *)&ex15_dat);
	std::cout << "GMRES Norm: " << (ex15_dat.b - x).norm() << std::endl;
	std::cout << "GMRES Total Evaluations: " << gmres_dat15.iter_total << std::endl;
	std::cout << "GMRES Restarts: " << gmres_dat15.iter_outer << std::endl;
	
	//Setup PCG
	PCG_DATA pcg_dat15;
	
	std::cout << "\n------------------START PCG------------------" << std::endl;
	
	time = clock();
	success = pcg(matvec_ex15, NULL, ex15_dat.b, &pcg_dat15, (void *)&ex15_dat, (void *)&ex15_dat);
	time = clock() - time;
	std::cout << "PCG Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	success = matvec_ex15(pcg_dat15.x, x, (void *)&ex15_dat);
	std::cout << "PCG Norm: " << (ex15_dat.b - x).norm() << std::endl;
	std::cout << "PCG Total Evaluations: " << pcg_dat15.iter << std::endl;
	
	//Setup BiCGSTAB
	BiCGSTAB_DATA bicg_dat15;
	
	std::cout << "\n------------------START BiCGSTAB------------------" << std::endl;
	
	time = clock();
	success = bicgstab(matvec_ex15, NULL, ex15_dat.b, &bicg_dat15, (void *)&ex15_dat, (void *)&ex15_dat);
	time = clock() - time;
	std::cout << "BiCGSTAB Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	success = matvec_ex15(bicg_dat15.x, x, (void *)&ex15_dat);
	std::cout << "BiCGSTAB Norm: " << (ex15_dat.b - x).norm() << std::endl;
	std::cout << "BiCGSTAB Total Evaluations: " << bicg_dat15.iter << std::endl;
	
	//Setup CGS
	CGS_DATA cgs_dat15;
	
	std::cout << "\n------------------START CGS------------------" << std::endl;
	
	time = clock();
	success = cgs(matvec_ex15, NULL, ex15_dat.b, &cgs_dat15, (void *)&ex15_dat, (void *)&ex15_dat);
	time = clock() - time;
	std::cout << "CGS Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	success = matvec_ex15(cgs_dat15.x, x, (void *)&ex15_dat);
	std::cout << "CGS Norm: " << (ex15_dat.b - x).norm() << std::endl;
	std::cout << "CGS Total Evaluations: " << cgs_dat15.iter << std::endl;
	
	
	std::cout << "\n------------------START Left Preconditioned GMRES-------------------" << std::endl;
	
	GMRESLP_DATA gmres_rec15;
	
	time = clock();
	success = gmresLeftPreconditioned(matvec_ex15,NULL,ex15_dat.b,&gmres_rec15,(void *)&ex15_dat,NULL);
	time = clock() - time;
	std::cout << "GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	success = matvec_ex15(gmres_rec15.x, x, (void *)&ex15_dat);
	std::cout << "GMRES Norm: " << (ex15_dat.b - x).norm() << std::endl;
	std::cout << "GMRES Total Evaluations: " << gmres_rec15.steps << std::endl;
	std::cout << "GMRES Restarts: " << gmres_rec15.iter << std::endl;
	
	std::cout << "\n------------------START GCR-------------------" << std::endl;
	
	GCR_DATA gcr15;
	
	time = clock();
	success = gcr(matvec_ex15,NULL,ex15_dat.b,&gcr15,(void *)&ex15_dat,(void *)&ex15_dat);
	time = clock() - time;
	std::cout << "GCR Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	success = matvec_ex15(gcr15.x, x, (void *)&ex15_dat);
	std::cout << "GCR Norm: " << (ex15_dat.b - x).norm() << std::endl; //PASS
	std::cout << "GCR inner iterations: " << gcr15.iter_inner << std::endl;
	std::cout << "GCR outer iterations: " << gcr15.iter_outer << std::endl;
	std::cout << "GCR total iterations: " << gcr15.total_iter << std::endl;
	
	
	std::cout << "\n------------------START GMRESR-------------------" << std::endl;
	
	GMRESR_DATA gmresr15;
	//gmresr15.GMRES_Output = true;
	gmresr15.gmres_maxit = 1;
	
	time = clock();
	success = gmresr(matvec_ex15,NULL,ex15_dat.b,&gmresr15,(void *)&ex15_dat,(void *)&ex15_dat);
	time = clock() - time;
	std::cout << "GMRESR Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	success = matvec_ex15(gmresr15.gcr_dat.x, x, (void *)&ex15_dat);
	std::cout << "GMRESR Norm: " << (ex15_dat.b - x).norm() << std::endl; //PASS
	std::cout << "GMRESR inner iterations: " << gmresr15.iter_inner << std::endl;
	std::cout << "GMRESR outer iterations: " << gmresr15.iter_outer	<< std::endl;
	std::cout << "GMRESR total iterations: " << gmresr15.total_iter << std::endl;
	
	std::cout << "\n------------------START KMS-------------------" << std::endl;
	
	KMS_DATA kms15;
	kms15.Output_inner = false;
	kms15.max_level = 2;
	kms15.inner_reltol = 0.001;
	
	time = clock();
	success = krylovMultiSpace(matvec_ex15,NULL,ex15_dat.b,&kms15,(void *)&ex15_dat,(void *)&ex15_dat);
	time = clock() - time;
	std::cout << "KMS Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
	success = matvec_ex15(kms15.gmres_out.x, x, (void *)&ex15_dat);
	std::cout << "KMS Norm: " << (ex15_dat.b - x).norm() << std::endl; //PASS
	std::cout << "KMS inner iterations: " << kms15.inner_iter << std::endl;
	std::cout << "KMS outer iterations: " << kms15.outer_iter << std::endl;
	std::cout << "KMS total iterations: " << kms15.total_iter << std::endl;
	
	
	std::cout << "------------------END Example 15------------------\n" << std::endl;
	
	
	return success;
}
