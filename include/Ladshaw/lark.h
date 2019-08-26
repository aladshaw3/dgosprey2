/*!
 *  \file lark.h lark.cpp
 *	\brief Linear Algebra Residual Kernels
 *	\details	The functions contained within are designed to solve generic linear and
 non-linear square systems of equations given a function argument and
 data from the user. Optionally, the user can also provide a function to
 return a preconditioning result that will be applied to the system.
 
 Having the user define how the preconditioning is carried out provides
 two major advantages: (1) we do not need to store and large, sparse
 preconditioning matrices and instead only store the preconditioned
 vector result and (2) this allows the user to use any kind of preconditioner
 they see fit for their problem.
 
 The Arnoldi function is typically not called by the user, but can be if
 desired. It accepts the function arguments and a residual vector to form
 an orthonormal basis of the Krylov subspace using the Modified Gram-Schmidt
 process (aka Arnoldi Iteration). This function is called by GMRES to iteratively
 solve a linear system of equations. Note that you can use this function to
 directly solve the linear system as long as that system is not too large.
 Construction of the basis is expensive, which is why this is used as a sub-function
 of an iterative method.
 
 The Restarted GMRES function will accept function arguments for a linear system
 and attempt to solve said system iteratively by constructing an orthonormal
 basis from the Krylov function. Note that this GMRES function does support
 restarting and will use restarting by default if the linear system is too
 large.
 
 Also included is a GMRES algorithm without restarting. This will directly solve
 the linear system within residual tolerance using a Full Orthogonal basis set
 of that system. It is equivalent to calling the Krylov method with the k parameter
 equal to N (i.e. the number of equations). This method is nick-named the Full
 Othogonalization Method (FOM), although the true FOM algorithm in literature is
 slightly different.
 
 The PJFNK function will accept function arguments for a square, non-linear
 system of equations and attempt to solve it iteratively using both the
 GMRES and Krylov functions with Newton's method to convert the non-linear
 system into a linear system.
 
 Also built here is a PCG implementation for solving symmetric linear systems.
 Can also be called by PJFNK if we know that the linear system (i.e. the
 Jacobian) is symmetric. This algorithm is significantly more efficient
 than GMRES, but is only valid if the system of equations is symmetric.
 
 Other linear solvers implemented in this work are the BiCGSTAB and CGS algorithms
 for non-symmetric, positive definite matrices. These algorithms are significantly
 more computationally efficient than GMRES or FOM. However, they can both break down
 if the linear system is poorly conditioned. In general, you only want to use these
 methods if you have preconditioning available and your linear system is very, very
 large. Otherwise, you will be better suited to using GMRES or FOM.
 
 There is also an implementation of the Generalized Conjugate Residual (GCR) method
 with and without restarting. This is a GMRES-like method that should give the
 exact solution within N iterations, where N is the original size of the matrix. Built
 ontop of the GCR method is a GMRESR (or GMRES Recursive) algorithm that uses GCR as
 the base method and performs GMRESRP iterations as a preconditioner at each iteration
 of GCR. This is the only linear solver that has built-in preconditioning. As a result,
 it may be slower than other algorithms for simple problems, but generally will have
 much better convergence behavior and will almost always give better residual reduction,
 even for hard to solve problems.
 
 We have also developed a novel/experimental iterative method based on the idea of
 recursively preconditioning a Krylov Subspace with more Krylov Subspaces. We have
 called with algorithm the Krylov Multi-Space (KMS) method. This algorithm is based on
 publications from Vorst and Vuik (1991) and Saad (1993). The idea is too use the FGMRES
 algorithm developed by Saad (1993) and precondition it with more FGMRES steps, i.e.,
 nesting the iterations as Vorst and Vuik (1991) had proposed. In this way, we have
 created a generalized Krylov Subspace method that has it's own variable preconditioner
 that can be adjusted depending on the user's desired complexity and convergence rate.
 If the levels of recursion requested is zero, then this algorithm is exactly equal to
 GMRES with right preconditioning. If the level is one, then it is FGMRES with a GMRES
 preconditioner. However, we allow the levels of recursion to reach up to 5, thus allowing
 us to precondition the preconitioners with more GMRES steps. This can result is significantly
 faster convergence rates, but is typically only necessary for very large or difficult
 to solve problems.
 
 NOTE: There are three GMRES implementations: (i) gmresLP, (ii) fom, and
 (iii) gmresRP. GMRESLP is a restarted GMRES implementation that is left
 preconditioned and only checks the residual on the outer loops. This may
 be less efficient than GMRESRP, which can check both outer and inner loop
 residuals. However, GMRESRP has to use right preconditioning, which also
 slightly changes the convergence behavior of the linear system. GMRES with
 left preconditioning and without restarting will just build the full
 subspace by default, thus solving the system exactly, but may require too
 much memory. You can do a GMRESRP unrestarted by specifying that the
 restart parameter be equal to the size of the problem.
 
 Basic Implementation Details:
 -----------------------------
 Linear Solvers -> Solve Ax=b for x \n
 Non-Linear Solvers -> Solve F(x)=0 for x \n
 All implementations require system size to be 2 or greater
 
 *  \author Austin Ladshaw
 *	\date 10/14/2014
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#pragma once

#ifndef LARK_HPP_
#define LARK_HPP_

#include "macaw.h"
#include <float.h>

#ifndef MIN_TOL
#define MIN_TOL 1e-15	///< Minimum Allowable Tolerance for linear and non-linear problems
#endif

/// Data structure for the construction of the Krylov subspaces for a linear system
/** C-style object used in conjunction with the Arnoldi algorithm to construct an orthonormal
	basis and upper Hessenberg representation of a given linear operator. This is used to solve
	a linear system both iteratively (i.e., in conjunction with GMRESLP) and directly (i.e., in
	conjunction with FOM). Alternatively, you can just store the factorized components for later
	use in another routine. */
typedef struct
{
	int k;              ///< Desired size of the Krylov subspace
	int iter;           ///< Actual size of the Krylov subspace
	
	double beta;        ///< Normalization parameter
	double hp1;         ///< Additional row element of H (separate storage for holding)
	
	bool Output = true; ///< True = print messages to console
	
	std::vector< Matrix<double> > Vk; ///< (N) x (k) orthonormal vector basis stored as a vector of column matrices
	Matrix<double> Hkp1;			///< (k+1) x (k) upper Hessenberg matrix
	Matrix<double> yk;          	///< (k) x (1) vector search direction
	Matrix<double>e1;				///< (k) x (1) orthonormal vector with 1 in first position
	Matrix<double> w;           	///< (N) x (1) interim result of the matrix_vector multiplication
	Matrix<double> v;				///< (N) x (1) holding cell for the column entries of Vk and other interims
	Matrix<double> sum;				///< (N) x (1) running sum of subspace vectors for use in altering w
	
}ARNOLDI_DATA;

/// Data structure for implementation of the Restarted GMRES algorithm with Left Preconditioning
/** C-style object used in conjunction with Generalized Minimum RESidual Left-Precondtioned
	(GMRESLP) and Full Orthogonalization Method (FOM) algorithms to iteratively or directly
	solve a linear system of equations. When using with GMRESLP, you can only check/observe the linear
	residuals before a restart or after the Arnoldi space is constructed. This is because this object
	uses Left-side Preconditioning. A faster routine may be GMRESRP, which is able to construct residuals
	after each Arnoldi iteration. */
typedef struct
{
	int restart = -1;    		///< Restart parameter - default = min(vector_size,20)
	int maxit = 0;          	///< Maximum allowable iterations - default = min(vector_size,1000)
	int iter = 0;           	///< Number of iterations needed for convergence
	int steps = 0;				///< Total number of gmres iterations and krylov iterations
	
	double tol_rel = 1e-6;		///< Relative tolerance for convergence - default = 1e-6
	double tol_abs = 1e-6;		///< Absolution tolerance for convergence - default = 1e-6
	double res;					///< Absolution redisual norm of the linear system
	double relres;      		///< Relative residual norm of the linear system
	double relres_base; 		///< Initial residual norm of the linear system
	double bestres;				///< Best found residual norm of the linear system
	
	bool Output = true;			///< True = print messages to console
	
	Matrix<double> x;           ///< Current solution to the linear system
	Matrix<double> bestx;		///< Best found solution to the linear system
	Matrix<double> r;			///< Residual vector for the linear system
	
	ARNOLDI_DATA arnoldi_dat;	///< Data structure for the kyrlov subspace
	
}GMRESLP_DATA;

/// Data structure for the Restarted GMRES algorithm with Right Preconditioning
/** C-style object used in conjunction with Generalized Minimum RESidual Right Preconditioned
	(GMRESRP) algorithm to iteratively solve a linear system of equations. Unlike GMRESLP, the
	GMRESRP method is capable of checking linear residuals at both the inner and outer steps. As
	a result, this algorithm may terminate earlier than GMRESLP if it has found a suitable solution
	during one of the inner steps. */
typedef struct
{
	int restart = -1;			///< Restart parameter - default = min(20,vector_size)
	int maxit = 0;				///< Maximum allowable outer iterations
	int iter_outer = 0;			///< Total number of outer iterations
	int iter_inner = 0;			///< Total number of inner iterations
	int iter_total = 0;			///< Total number of overall iterations
	
	double tol_rel = 1e-6;		///< Relative tolerance for convergence - default = 1e-6
	double tol_abs = 1e-6;		///< Absolute tolerance for convergence - default = 1e-6
	double res;					///< Absolute residual norm for linear system
	double relres;				///< Relative residual norm for linear system
	double relres_base;			///< Initial residual norm of the linear system
	double bestres;				///< Best found residual norm of the linear system
	
	bool Output = true;			///< True = print messages to console
	
	Matrix<double> x;					///< Current solution to the linear system
	Matrix<double> bestx;				///< Best found solution to the linear system
	Matrix<double> r;					///< Residual vector for the linear system
	
	std::vector< Matrix<double> > Vk;			///< (N x k) orthonormal vector basis
	std::vector< Matrix<double> > Zk;			///< (N x k) preconditioned vector set
	std::vector< std::vector< double > > H;		///< (k+1 x k) upper Hessenberg storage matrix
	std::vector< std::vector< double > > H_bar;	///< (k+1 x k) Factorized matrix
	std::vector< double > y;					///< (k x 1) Vector search direction
	std::vector< double > e0;					///< (k+1 x 1) Normalized vector with residual info
	std::vector< double > e0_bar;				///< (k+1 x 1) Factorized normal vector
	
	Matrix<double> w;           	///< (N) x (1) interim result of the matrix_vector multiplication
	Matrix<double> v;				///< (N) x (1) holding cell for the column entries of Vk and other interims
	Matrix<double> sum;				///< (N) x (1) running sum of subspace vectors for use in altering w
	
}GMRESRP_DATA;

/// Data structure for implementation of the PCG algorithms for symmetric linear systems
/** C-style object used in conjunction with the Preconditioned Conjugate Gradient (PCG)
	algorithm to iteratively solve a symmetric linear system of equations. This algorithm
	is optimal if your linear system is symmetric, but will not work at all if your system is
	asymmetric. For asymmetric systems, use one of the other linear methods. */
typedef struct
{
	int maxit = 0;				///< Maximum allowable iterations - default = min(vector_size,1000)
	int iter = 0;				///< Actual number of iterations taken
	
	double alpha;				///< Step size for new solution
	double beta;				///< Step size for new search direction
	double tol_rel = 1e-6;		///< Relative tolerance for convergence - default = 1e-6
	double tol_abs = 1e-6;		///< Absolution tolerance for convergence - default = 1e-6
	double res;					///< Absolute residual norm
	double relres;				///< Relative residual norm
	double relres_base;			///< Initial residual norm
	double bestres;				///< Best found residual norm
	
	bool Output = true;			///< True = print messages to console
	
	Matrix<double> x;						///< Current solution to the linear system
	Matrix<double> bestx;					///< Best found solution to the linear system
	Matrix<double> r;						///< Residual vector for the linear system
	Matrix<double> r_old;					///< Previous residual vector
	Matrix<double> z;						///< Preconditioned residual vector (result of precon function)
	Matrix<double> z_old;					///< Previous preconditioned residual vector
	Matrix<double> p;						///< Search direction
	Matrix<double> Ap;						///< Result of matrix-vector multiplication
	
}PCG_DATA;

/// Data structure for the implementation of the BiCGSTAB algorithm for non-symmetric linear systems
/** C-style object used in conjunction with the Bi-Conjugate Gradient STABalized (BiCGSTAB) algorithm
	to solve a linear system of equations. This algorithm is generally more efficient than any GMRES
	or GCR variant, but may not always reduce the residual at each step. However, if used with preconditioning,
	then this algorithm is very efficient, especially when used for solving grid-based linear systems. */
typedef struct
{
	int maxit = 0;				///< Maximum allowable iterations - default = min(2*vector_size,1000)
	int iter = 0;				///< Actual number of iterations
	bool breakdown;				///< Boolean to determine if the method broke down
	
	double alpha;				///< Step size parameter for next solution
	double beta;				///< Step size parameter for search direction
	double rho;					///< Scaling parameter for alpha and beta
	double rho_old;				///< Previous scaling parameter for alpha and beta
	double omega;				///< Scaling parameter and additional step length
	double omega_old;			///< Previous scaling parameter and step length
	double tol_rel = 1e-6;		///< Relative tolerance for convergence - default = 1e-6
	double tol_abs = 1e-6;		///< Absolution tolerance for convergence - default = 1e-6
	double res;					///< Absolute residual norm
	double relres;				///< Relative residual norm
	double relres_base;			///< Initial residual norm
	double bestres;				///< Best found residual norm
	
	bool Output = true;			///< True = print messages to console
	
	Matrix<double> x;					///< Current solution to the linear system
	Matrix<double> bestx;				///< Best found solution to the linear system
	Matrix<double> r;					///< Residual vector for the linear system
	Matrix<double> r0;					///< Initial residual vector
	Matrix<double> v;					///< Search direction for p
	Matrix<double> p;					///< Search direction for updating
	Matrix<double> y;					///< Preconditioned search direction
	Matrix<double> s;					///< Residual updating vector
	Matrix<double> z;					///< Preconditioned residual updating vector
	Matrix<double> t;					///< Search direction for resdidual updates
	
}BiCGSTAB_DATA;

/// Data structure for the implementation of the CGS algorithm for non-symmetric linear systems
/** C-style object to be used in conjunction with the Conjugate Gradient Squared (CGS) algorithm
	to solve linear systems of equations. This algorithm is slightly less computational work than
	BiCGSTAB, but is much less stable. As a result, I do not recommend using this algorithm unless
	you also use some form of preconditioning. */
typedef struct
{
	int maxit = 0;				///< Maximum allowable iterations - default = min(2*vector_size,1000)
	int iter = 0;				///< Actual number of iterations
	bool breakdown;				///< Boolean to determine if the method broke down
	
	double alpha;				///< Step size parameter for next solution
	double beta;				///< Step size parameter for search direction
	double rho;					///< Scaling parameter for alpha and beta
	double sigma;				///< Scaling parameter and additional step length
	double tol_rel = 1e-6;		///< Relative tolerance for convergence - default = 1e-6
	double tol_abs = 1e-6;		///< Absolution tolerance for convergence - default = 1e-6
	double res;					///< Absolute residual norm
	double relres;				///< Relative residual norm
	double relres_base;			///< Initial residual norm
	double bestres;				///< Best found residual norm
	
	bool Output = true;			///< True = print messages to console
	
	Matrix<double> x;					///< Current solution to the linear system
	Matrix<double> bestx;				///< Best found solution to the linear system
	Matrix<double> r;					///< Residual vector for the linear system
	Matrix<double> r0;					///< Initial residual vector
	Matrix<double> u;					///< Search direction for v
	Matrix<double> w;					///< Updates sigma and u
	Matrix<double> v;					///< Search direction for x
	Matrix<double> p;					///< Preconditioning result for w, z, and matvec for Ax
	Matrix<double> c;					///< Holds the matvec result between A and p
	Matrix<double> z;					///< Full search direction for x
	
}CGS_DATA;

/// Data structure for implementation of linear operator transposition
/** C-style object used in conjunction with the Operator Transpose algorithm to form
	an action of A^T*r when A is only available as a linear operator and not a matrix.
	This is a sub-routine required by GCR and GMRESR to stabilize the outer iterations. */
typedef struct
{
	Matrix<double> Ii;					///< The ith column vector of the identity operator
	Matrix<double> Ai;					///< The ith column vector of the user's linear operator
	
}OPTRANS_DATA;

/// Data structure for the implementation of the GCR algorithm for non-symmetric linear systems
/** C-style object used in conjunction with the Generalized Conjugate Residual (GCR) algorithm for
	solving a non-symmetric linear system of equations. When the linear system in question has a
	positive-definite-symmetric component to it, then this algorithm is equivalent to GMRESRP. However,
	it is generally less efficient than GMRESRP and can suffer breakdowns. */
typedef struct
{
	int restart = -1;			///< Restart parameter for outer iterations - default = 20
	int maxit = 0;				///< Maximum allowable outer iterations
	int iter_outer = 0;			///< Number of outer iterations taken
	int iter_inner = 0;			///< Number of inner iterations taken
	int total_iter = 0;			///< Total number of iterations taken
	bool breakdown = false;		///< Boolean to determine if a step has failed
	
	double alpha;				///< Inner iteration step size
	double beta;				///< Outer iteration step size
	double tol_rel = 1e-6;		///< Relative tolerance for convergence - default = 1e-6
	double tol_abs = 1e-6;		///< Absolute tolerance for convergence - default = 1e-6
	double res;					///< Absolute residual norm for linear system
	double relres;				///< Relative residual norm for linear system
	double relres_base;			///< Initial residual norm of the linear system
	double bestres;				///< Best found residual norm of the linear system
	
	bool Output = true;			///< True = print messages to the console
	
	Matrix<double> x;					///< Current solution to the linear system
	Matrix<double> bestx;				///< Best found solution to the linear system
	Matrix<double> r;					///< Residual Vector
	Matrix<double> c_temp;				///< Temporary c vector to be updated
	Matrix<double> u_temp;				///< Temporary u vector to be updated
	std::vector<Matrix<double> > u;		///< Vector span for updating x
	std::vector<Matrix<double> > c;		///< Vector span for updating r
	
	OPTRANS_DATA transpose_dat;			///< Data structure for Operator Transposition
	
}GCR_DATA;

/// Data structure for the implementation of GCR with Nested GMRES preconditioning (i.e., GMRESR)
/** C-style object to be used in conjunction with the Generalized Minimum RESidual Recurive (GMRESR)
	algorithm. Although the name suggests that this method used GMRES recursively, what it is actually
	doing is nesting GMRESRP iterations inside the GCR method to form a preconditioner for GCR. The name
	GMRESR came from literature (Vorst and Vuik, "GMRESR: A family of nested GMRES methods", 1991).*/
typedef struct
{
	int gcr_restart = -1;		///< Number of GCR restarts (default = 20, max = N)
	int gcr_maxit = 0;			///< Number of GCR iterations
	int gmres_restart = -1;		///< Number of GMRES restarts (max = 20)
	int gmres_maxit = 1;		///< Number of GMRES iterations (max = 5, default = 1)
	int N;						///< Dimension of the linear system
	int total_iter;				///< Total GMRES and GCR iterations
	int iter_outer;				///< Total GCR iterations
	int iter_inner;				///< Total GMRES iterations
	
	bool GCR_Output = true;		///< True = print GCR messages
	bool GMRES_Output = false;	///< True = print GMRES messages
	
	double gmres_tol = 0.1;		///< Tolerance relative to GCR iterations
	double gcr_rel_tol = 1e-6;	///< Relative outer residual tolerance
	double gcr_abs_tol = 1e-6;	///< Absolute outer residual tolerance
	
	Matrix<double> arg;			///< Argument matrix passed between preconditioner and iterator
	
	GCR_DATA gcr_dat;			///< Data structure for the outer GCR steps
	GMRESRP_DATA gmres_dat;		///< Data structure for the inner GMRES steps
	
	/// User supplied matrix-vector product function
	int (*matvec) (const Matrix<double> &x, Matrix<double> &Ax, const void *matvec_data);
	/// Optional user supplied terminal preconditioner
	int (*terminal_precon) (const Matrix<double> &r, Matrix<double> &p, const void *precon_data);
	
	const void *matvec_data;	///< Data structure for the user's matvec function
	const void *term_precon;	///< Data structure for the user's terminal preconditioner
	
}GMRESR_DATA;

/// Data structure for the implemenation of the Krylov Multi-Space (KMS) Method
/** C-style object to be used in conjunction with the Krylov Multi-Space (KMS) Algorithm to iteratively
	solve non-symmetric, indefinite linear systems. This method was inspired by the Flexible GMRES (FGMRES)
	and Recursive GMRES (GMRESR) methods proposed by Saad (1993) and Vorst and Vuik (1991), respectively.
	The idea behind this method is to recursively call FGMRES to solve a linear system with pregressively
	smaller Krylov Subspaces built by a Right-Preconditioned GMRES algorithm. Thus creating a "V-cycle"
	of iteration similar to that seen in Multi-Grid algorithms. */
typedef struct
{
	int level = 0;			///< Current level in the recursion
	int max_level = 0;		///< Maximum allowable recursion levels (Default = 0 -> GMRES, Max = 5)
	int restart = -1;		///< Restart parameter for the outer iterates (Default = 20, Max = N)
	int maxit = 0;			///< Maximum allowable iterations for the outer steps
	int inner_iter = 0;		///< Number of inner steps taken
	int outer_iter = 0;		///< Number of outer steps taken
	int total_iter = 0;		///< Total number of iterations in all steps
	
	double outer_reltol = 1e-6;		///< Relative residual tolerance for outer steps (Default = 1e-6)
	double outer_abstol = 1e-6;		///< Absolute residual tolerance for outer steps (Default = 1e-6)
	double inner_reltol = 0.1;		///< Residual tolerance for inner steps made relative to outer steps (Default = 0.1)
	
	bool Output_outer = true;			///< True = Print the outer steps residuals
	bool Output_inner = false;			///< True = Print the inner steps residuals
	
	GMRESRP_DATA gmres_out;					///< Data structure for the outer steps
	std::vector<GMRESRP_DATA> gmres_in;		///< Data structures for each recursion level
	
	/// User supplied matrix-vector product function
	int (*matvec) (const Matrix<double> &x, Matrix<double> &Ax, const void *matvec_data);
	/// Optional user supplied terminal preconditioner
	int (*terminal_precon) (const Matrix<double> &r, Matrix<double> &p, const void *precon_data);
	
	const void *matvec_data;	///< Data structure for the user's matvec function
	const void *term_precon;	///< Data structure for the user's terminal preconditioner
	
}KMS_DATA;

///Data structure for the implementation of a QR solver given some invertable linear operator
/** C-style object to be used in conjuction with a QR solver for invertable linear operators.
	This method will extract columns from the linear operator and use Householder Reflections
	to factor the operator into an upper triangular matrix and a unitary reflection matrix.
	It is generally less efficient to use this method for sparse systems, but is more stable
	and occassionally more efficient for dense systems. */
typedef struct
{
	Matrix<double> ek;				///< Unit vector used to extract columns from the linear operator
	Matrix<double> Ro;				///< Upper triangular matrix formed from factoring the linear operator
	Matrix<double> x;				///< Solution to the linear system
	
}QR_DATA;

/// Data structure for the implementation of a Picard or Fixed-Point iteration for non-linear systems
/** C-style object used in conjunction with the Picard algorithm for solving a non-linear system of equations.
	This is an extradorinarily simple iterative method by which a weak or loose form of the non-linear system
	is solved based on an initial guess. User must supplied a residual function for the non-linear system and
	a function representing the weak solution. Generally, this method is less efficient than Newton methods,
	but is significantly cheaper. */
typedef struct
{
	int maxit = 0;					///< Maximum allowable iterations - default = min(3*vec_size,1000)
	int iter = 0;					///< Actual number of iterations
	
	double tol_rel = 1e-6;			///< Relative tolerance for convergence - default = 1e-6
	double tol_abs = 1e-6;			///< Absolution tolerance for convergence - default = 1e-6
	double res;						///< Residual norm of the iterate
	double relres;					///< Relative residual norm of the iterate
	double relres_base;				///< Initial residual norm
	double bestres;					///< Best found residual norm
	
	bool Output = true;			///< True = print messages to console
	
	Matrix<double> x0;						///< Previous iterate solution vector
	Matrix<double> bestx;					///< Best found solution vector
	Matrix<double> r;						///< Residual of the non-linear system
	
}PICARD_DATA;

/// Data structure for the implementation of Backtracking Linesearch
/** C-style object used in conjunction with the Backtracking Linesearch algorithm to smooth out convergence
	of Netwon based iterative methods for non-linear systems of equations. The actual algorithm has been separated
	from the interior of the Newton method so that it can be included in any future Newton based iterative methods
	being developed. */
typedef struct
{
	int fun_call = 0;					///< Number of function calls made during line search
	double alpha = 1e-4;				///< Scaling parameter for determination of search step size
	double rho = 0.1;					///< Scaling parameter for to change step size by
	double lambdaMin=DBL_EPSILON;		///< Smallest allowable step length
	double normFkp1;					///< New residual norm of the Newton step
	
	bool constRho = false;				///< True = use a constant value for rho
	
	Matrix<double> Fk;					///< Old residual vector of the Newton step
	Matrix<double> xk;					///< Old solution vector of the Newton step
	
}BACKTRACK_DATA;

/// Enum of definitions for linear solver types in PJFNK
/** Enum delineates the available Krylov Subspace methods that can be used to
	solve the linear sub-problem at each non-linear iteration in a Newton method. */
typedef enum
{
	GMRESLP,
	PCG,
	BiCGSTAB,
	CGS,
	FOM,
	GMRESRP,
	GCR,
	GMRESR,
	KMS,
	QR
} krylov_method;

/// Data structure for the implementation of the PJFNK algorithm for non-linear systems
/** C-style object to be used in conjunction with the Preconditioned Jacobian-Free Newton-Krylov
	(PJFNK) method for solving a non-linear system of equations. You can use any of the Krylov methods
	listed in the krylov_method enum to solve the linear sub-problem. When FOM is specified as the
	Krylov method, this algorithm becomes equivalent to an exact Newton method. If no Krylov method
	is specified, then the algorithm will try to pick a method based on the problem size and availability
	of preconditioning. */
typedef struct
{
	int nl_iter = 0;         	///< Number of non-linear iterations
	int l_iter = 0;				///< Number of linear iterations
	int fun_call = 0;			///< Actual number of function calls made
	int nl_maxit = 0;			///< Maximum allowable non-linear steps
	int linear_solver = -1;		///< Flag to denote which linear solver to use - default = PJFNK Chooses
	
	double nl_tol_abs = 1e-6;   ///< Absolute Convergence tolerance for non-linear system - default = 1e-6
	double nl_tol_rel = 1e-6;	///< Relative Convergence tol for the non-linear system - default = 1e-6
	double lin_tol_rel = 1e-6;	///< Relative tolerance of the linear solver - default = 1e-6
	double lin_tol_abs = 1e-6;	///< Absolute tolerance of the linear solver - default = 1e-6
	double nl_res;				///< Absolute redidual norm for the non-linear system
	double nl_relres;       	///< Relative residual for the non-linear system
	double nl_res_base;     	///< Initial residual norm for the non-linear system
	double nl_bestres;			///< Best found residual norm
	double eps=sqrt(DBL_EPSILON);///< Value of epsilon used jacvec - default = sqrt(DBL_EPSILON)
	
	bool NL_Output = true;		///< True = print PJFNK messages to console
	bool L_Output = false;		///< True = print Linear messages to console
	bool LineSearch = false;	///< True = use Backtracking Linesearch for global convergence
	bool Bounce = false;		///< True = allow Linesearch to go outside local well, False = Strict local convergence
	
	Matrix<double> F;					///< Stored fuction evaluation at x (also the residual)
	Matrix<double> Fv;					///< Stored function evaluation at x+eps*v
	Matrix<double> v;					///< Stored vector of x+eps*v
	Matrix<double> x;					///< Current solution vector for the non-linear system
	Matrix<double> bestx;				///< Best found solution vector to the non-linear system
	
	//The PJFNK implementation will choose the linear method best suited for your problem
	GMRESLP_DATA gmreslp_dat;   	///< Data structure for the GMRESLP method
	PCG_DATA pcg_dat;				///< Data structure for the PCG method
	BiCGSTAB_DATA bicgstab_dat;		///< Data structure for the BiCGSTAB method
	CGS_DATA cgs_dat;				///< Data structure for the CGS method
	GMRESRP_DATA gmresrp_dat;		///< Data structure for the GMRESRP method
	GCR_DATA gcr_dat;				///< Data structure for the GCR method
	GMRESR_DATA gmresr_dat;			///< Data structure for the GMRESR method
	KMS_DATA kms_dat;				///< Data structure for the KMS method
	QR_DATA qr_dat;					///< Data structure for the QR solve method
	
	
	BACKTRACK_DATA backtrack_dat;	///< Data structure for the Backtracking Linesearch algorithm
	
	//------------Below is a list of pointers the user must provide---------------
	
	/// Data structure pointer for user's residual data
	const void *res_data;
	/// Data structure pointer for user's preconditioning data
	const void *precon_data;
	/// Function pointer for the user's function F(x) using there data
	int (*funeval) (const Matrix<double> &x, Matrix<double> &F, const void *res_data);
	/// Function pointer for the user's preconditioning function for the linear system
	int (*precon) (const Matrix<double> &r, Matrix<double> &p, const void *precon_data);
	
}PJFNK_DATA;

/// Data structure to form a numerical jacobian matrix with finite differences
/** C-style object to be used in conjunction with the Numerical Jacobian algorithm. This algorithm
	will used double-precision finite-differences to formulate an approximate Jacobian matrix at the
	given variable state for the given residual/non-linear function. */
typedef struct
{
	double eps = sqrt(DBL_EPSILON);		///< Perturbation value
	Matrix<double> Fx;					///< Vector of function evaluations at x
	Matrix<double> Fxp;					///< Vector of function evaluations at x+eps
	Matrix<double> dxj;					///< Vector of perturbed x values
	
}NUM_JAC_DATA;

/// \cond

//Data structures used by examples to test the algorithms of LARK
typedef struct
{
	Matrix<double> M;               //(N) x (N) coefficient matrix
	Matrix<double> b;               //(N) x (1) right side vector solution
	
}EX01_DATA;

typedef struct
{
	Matrix<double> M;               //(N) x (N) coefficient matrix
	Matrix<double> b;               //(N) x (1) right side vector solution
	
}EX02_DATA;

typedef struct
{
	Matrix<double> M;               //(N) x (N) coefficient matrix
	Matrix<double> b;               //(N) x (1) right side vector solution
	
}EX04_DATA;

typedef struct
{
	double h;				//Mesh size for non-linear PDE
	double k;				//Scaling factor for non-linear portion of PDE
	int N;					//Number of nodes in mesh
	Matrix<double> x;				//Nx1 solution vector
	Matrix<double> M;				//NxN preconditioning matrix
	Matrix<double> s;				//Search direction for Picard
	Matrix<double> p;				//Conditioned search direction
}EX09_DATA;

typedef struct
{
	int N;					//Number of nodes in the mesh
	int m;					//Dimensional length for natural laplacian
	Matrix<double> b;				//Boundary Conditions for the problem
}EX15_DATA;

//These functions are used to run LARK Tests
int matvec_ex01(const Matrix<double>& v, Matrix<double>& w, const void *data);

int precon_ex01(const Matrix<double>& b, Matrix<double>& p, const void *data);

int matvec_ex02(const Matrix<double>& v, Matrix<double>& w, const void *data);

int matvec_ex04(const Matrix<double>& v, Matrix<double>& w, const void *data);

int precon_ex04(const Matrix<double>& b, Matrix<double>& p, const void *data);

int evalx_ex09(const Matrix<double> &x, Matrix<double>& G, const void *data);

int funeval_ex09(const Matrix<double> &x, Matrix<double>& F, const void *data);

int funeval_ex10(const Matrix<double> &x, Matrix<double>& F, const void *data);

int precon_ex10(const Matrix<double> &r, Matrix<double>& p, const void *data);

int matvec_ex15(const Matrix<double>& v, Matrix<double>& w, const void *data);

int precon_ex15(const Matrix<double>& w, Matrix<double>& p, const void *data);

/// \endcond

//---------------Below are the actual functions available in LARK---------------------

/// Function to update the linear vector x based on the Arnoldi Krylov subspace
/** This function will update a solution vector x based on the previous solution x0 given
	the orthonormal basis and upper Hessenberg matrix formed in the Arnoldi algorithm. Updating
	is automatically called by the GMRESLP function. It is expected that the Arnoldi algorithm
	has already been called prior to calling this function.
 
	\param x matrix that will hold the new updated solution to the linear system
	\param x0 matrix that holds the previous solution to the linear system
	\param arnoldi_dat pointer to the ARNOLDI_DATA data structure*/
int update_arnoldi_solution(Matrix<double>& x, Matrix<double>& x0, ARNOLDI_DATA *arnoldi_dat);

/// Function to factor a linear operator into an orthonormal basis and upper Hessenberg matrix
/** This function performs the Arnoldi algorithm to factor a linear operator into an orthonormal
	basis and upper Hessenberg matrix. Each orthonormal vector is formed using a Modified Gram-Schmidt
	procedure. When used in conjunction with GMRESLP, user may supply a preconditioning operator to
	improve convergence of the linear system. However, this function can be used by itself to
	factor the user's linear operator.
 
	\param matvec user supplied linear operator given as an int function
	\param precon user supplied preconditioning operator given as an int function
	\param r0 user supplied vector to serve as the first basis vector in the orthonormal basis
	\param arnoldi_dat pointer to the ARNOLDI_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int arnoldi( int (*matvec) (const Matrix<double>& v, Matrix<double> &w, const void *data),
			int (*precon) (const Matrix<double>& b, Matrix<double> &p, const void *data),
			Matrix<double> &r0, ARNOLDI_DATA *arnoldi_dat, const void *matvec_data,
			const void *precon_data );

/// Function to iteratively solve a non-symmetric, indefinite linear system with GMRESLP
/** This function iteratively solves a non-symmetric, indefinite linear system using the
	Generalized Minimum RESidual method with Left Preconditioning (GMRESLP). It calls the
	Arnoldi algorithm to factor a linear operator into an orthonormal basis and upper
	Hessenberg matrix, then uses that factorization to form an approximation to the linear
	system. Because this algorithm uses left-side preconditioning, it can only check the
	linear residuals at the outer iterations.
 
	\param matvec user supplied linear operator given as an int function
	\param precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param gmreslp_dat pointer to the GMRESLP_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int gmresLeftPreconditioned( int (*matvec) (const Matrix<double>& v, Matrix<double> &w, const void *data),
							int (*precon) (const Matrix<double>& b, Matrix<double> &p, const void *data),
							Matrix<double> &b, GMRESLP_DATA *gmreslp_dat, const void *matvec_data,
							const void *precon_data );

/// Function to directly solve a non-symmetric, indefinite linear system with FOM
/** This function directly solves a non-symmetric, indefinite linear system using the
	Full Orthogonalization Method (FOM). This algorithm is exactly equivalent to GMRESLP
	without restarting. Therefore, it uses the GMRESLP_DATA structure and calls the
	GMRESLP algorithm without using restarts. As a result, it never checks linear residuals.
	However, this should give the exact solution upon completion, assuming the linear
	operator is not singular.
 
	\param matvec user supplied linear operator given as an int function
	\param precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param gmreslp_dat pointer to the GMRESLP_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int fom( int (*matvec) (const Matrix<double>& v, Matrix<double> &w, const void *data),
		int (*precon) (const Matrix<double>& b, Matrix<double> &p, const void *data),
		Matrix<double> &b, GMRESLP_DATA *gmreslp_dat, const void *matvec_data,
		const void *precon_data );

/// Function to iteratively solve a non-symmetric, indefinite linear system with GMRESRP
/** This function iteratively solves a non-symmetric, indefinite linear system using the
	Generalized Minimum RESidual method with Right Preconditioning (GMRESRP). Because this
	algorithm uses right preconditioning, it is able to check the linear residuals at both
	the outer and inner iterations. This may be much for efficient compared to GMRESLP. In
	order to check inner residuals, this algorithm has to perform it's own internal Modified
	Gram-Schmidt procedure and will not call the Arnoldi algorithm.
 
	\param matvec user supplied linear operator given as an int function
	\param precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param gmresrp_dat pointer to the GMRESRP_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int gmresRightPreconditioned( int (*matvec) (const Matrix<double>& v, Matrix<double> &w, const void *data),
							 int (*precon) (const Matrix<double>& b, Matrix<double> &p, const void *data),
							 Matrix<double> &b, GMRESRP_DATA *gmresrp_dat, const void *matvec_data,
							 const void *precon_data );

/// Function to iteratively solve a symmetric, definite linear system with PCG
/** This function iteratively solves a symmetric, definite linear system using the
	Preconditioned Conjugate Gradient (PCG) method. The PCG algorithm is optimal in
	terms of efficiency and residual reduction, but only if the linear system is
	symmetric. PCG will fail if the linear operator is non-symmetric!
 
	\param matvec user supplied linear operator given as an int function
	\param precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param pcg_dat pointer to the PCG_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int pcg( int (*matvec) (const Matrix<double>& p, Matrix<double> &Ap, const void *data),
		int (*precon) (const Matrix<double>& r, Matrix<double> &z, const void *data),
		Matrix<double> &b, PCG_DATA *pcg_dat, const void *matvec_data,
		const void *precon_data );

/// Function to iteratively solve a non-symmetric, definite linear system with BiCGSTAB
/** This function iteratively solves a non-symmetric, definite linear system using the
	Bi-Conjugate Gradient STABilized (BiCGSTAB) method. This is a highly efficient algorithm
	for solving non-symmetric problems, but will occassionally breakdown and fail. Most
	common failures are caused by poor preconditioning. Works very well for grid-based
	linear systems.
 
	\param matvec user supplied linear operator given as an int function
	\param precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param bicg_dat pointer to the BiCGSTAB_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int bicgstab( int (*matvec) (const Matrix<double>& p, Matrix<double> &Ap, const void *data),
			 int (*precon) (const Matrix<double>& r, Matrix<double> &z, const void *data),
			 Matrix<double> &b, BiCGSTAB_DATA *bicg_dat, const void *matvec_data,
			 const void *precon_data );

/// Function to iteratively solve a non-symmetric, definite linear system with CGS
/** This function iteratively solves a non-symmetric, definite linear system using the
	Conjugate Gradient Squared (CGS) method. This is an extremely efficient algorithm
	for solving non-symmetric problems, but will often breakdown and fail. Most common
	failures are caused by poor or no preconditioning. Works very will for grid-based
	linear systems.
 
	\param matvec user supplied linear operator given as an int function
	\param precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param cgs_dat pointer to the CGS_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int cgs( int (*matvec) (const Matrix<double>& p, Matrix<double> &Ap, const void *data),
		int (*precon) (const Matrix<double>& r, Matrix<double> &z, const void *data),
		Matrix<double> &b, CGS_DATA *cgs_dat, const void *matvec_data,
		const void *precon_data );

/// Function that is used to perform transposition of a linear operator and results in a new vector A^T*r=u
/** This function takes a user supplied linear operator and forms the result of that operator transposed
	and multiplied by a given vector r (A^T*r=u). Transposition is accomplised by reordering the transpose
	operator and multiplying the non-transposed operator by a complete set of orthonormal vectors. The end
	result gives the ith component of the vector u for each operation (u_i = r^T*A*i). Here, i is a vector
	made from the ith column of the identity matrix. If the linear system if sufficiently large, then this
	operation may take some time.
 
	\param matvec user supplied linear operator given as an int function
	\param r vector to be multiplied by the transpose of the operator
	\param u vector to store the result of the operator transposition (u=A^T*r)
	\param transpose_dat pointer to the OPTRANS_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/
int operatorTranspose(int(*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data),
					  Matrix<double> &r, Matrix<double> &u, OPTRANS_DATA *transpose_dat,
					  const void *matvec_data);

/// Function to iteratively solve a non-symmetric, definite linear system with GCR
/** This function iteratively solves a non-symmetric, definite linear system using the
	Generalized Conjugate Residual (GCR) method. Similar to GMRESRP, this algorithm
	will construct a growing orthonormal basis set that will eventually form the exact
	solution to the linear system. However, this algorithm is less efficient than GMRESRP
	and can suffer breakdowns if the linear system is indefinite.
 
	\param matvec user supplied linear operator given as an int function
	\param precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param gcr_dat pointer to the GCR_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int gcr( int (*matvec) (const Matrix<double>& x, Matrix<double> &Ax, const void *data),
		int (*precon) (const Matrix<double>& r, Matrix<double> &Mr, const void *data),
		Matrix<double> &b, GCR_DATA *gcr_dat, const void *matvec_data,
		const void *precon_data );

/// Function used in conjunction with GMRESR to apply GMRESRP iterations as a preconditioner
/** This function is required to take the form of the user supplied preconditioning functions
	for other iterative methods. However, it cannot be used in conjunction with any other
	Krylov method. It is only called by the GMRESR function when the preconditioner needs to
	be applied.
 
	\param r vector supplied to the preconditioner to operate on
	\param Mr vector to hold the result of the preconditioning operation
	\param data void pointer to the GMRESR_DATA data structure*/
int gmresrPreconditioner( const Matrix<double>& r, Matrix<double> &Mr, const void *data);

/// Function to iteratively solve a non-symmetric, indefinite linear system with GMRESR
/** This function iteratively solves a non-symmetric, indefinite linear system using the
	Generalized Minimum RESidual Recursive (GMRESR) method. This algorithm actually uses
	GCR at the outer iterations, but stabilizes GCR with GMRESRP inner iterations to
	implicitly form a variable preconditioner to the linear system. As such, this is one
	of only two methods that inherently includes preconditioning (the other is KMS), without
	any user supplied preconditioning operator. However, this algorithms is signficantly more
	computationally expensive than GCR or GMRESRP separately. It should only be used for solving
	very large or very hard to solve linear systems.
 
	\param matvec user supplied linear operator given as an int function
	\param terminal_precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param gmresr_dat pointer to the GMRESR_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param term_precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*terminal_precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int gmresr( int (*matvec) (const Matrix<double>& x, Matrix<double> &Ax, const void *data),
		   int (*terminal_precon) (const Matrix<double>& r, Matrix<double> &Mr, const void *data),
		   Matrix<double> &b, GMRESR_DATA *gmresr_dat, const void *matvec_data,
		   const void *term_precon_data );

/// Preconditioner function for the Krylov Multi-Space
/** This function is required to take the form of the user supplied preconditioning functions
	for other iterative methods. However, it cannot be used in conjunction with any other
	Krylov method. It is only called by the KMS function when the preconditioner needs to
	be applied.
 
	\param r vector supplied to the preconditioner to operate on
	\param Mr vector to hold the result of the preconditioning operation
	\param data void pointer to the KMS_DATA data structure*/
int kmsPreconditioner( const Matrix<double>& r, Matrix<double> &Mr, const void *data);

/// Function to iteratively solve a non-symmetric, indefinite linear system with KMS
/** This function iteratively solves a non-symmetric, indefinite linear system using the
	Krylov Multi-Space (KMS) method. This algorithm uses GMRESRP at both outer and inner
	iterations to implicitly form a variable preconditioner to the linear system. As such,
	this is one of only two methods that inherently includes preconditioning, without any user
	supplied preconditioning operator (the other being GMRESR). The advantage to this method
	over GMRESR is that this method is GMRES at its core, and will therefore never breakdown
	or need to be stabilized. Additionally, you can call this method and set it's max_level
	parameter (see KMS_DATA) to 0, which will make this algorithm exactly equal to GMRESRP.
	If the max_level is set to 1, then this algorithm is exactly FGMRES (Saad, 1993) with the
	GMRES algorithm as a preconditioner. However, you can set max_level higher to precondition
	the preconditioners with more preconditioners. Thus creating a method with any desired
	complexity or rate of convergence.
 
	\param matvec user supplied linear operator given as an int function
	\param terminal_precon user supplied preconditioning operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param kms_dat pointer to the KMS_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator
	\param term_precon_data user supplied void pointer to a data structure needed for the precondtioning operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*terminal_precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the above linear operator function and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the original
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int krylovMultiSpace( int (*matvec) (const Matrix<double>& x, Matrix<double> &Ax, const void *data),
					 int (*terminal_precon) (const Matrix<double>& r, Matrix<double> &Mr, const void *data),
					 Matrix<double> &b, KMS_DATA *kms_dat, const void *matvec_data,
					 const void *term_precon_data );

/// Function to solve a dense linear operator system using QR factorization
/** This function is used to solve a dense linear system using QR factorization. It should only be used
	if iterative methods are unstable or if the linear system is very dense. There will likely be memory
	limitations to using this method, since it is assumed that the matrix/operator is dense. This method
	may also be less efficient because it has to extract the matrix elements from the linear operator. So
	if the linear operator is large, then the setup cost for this method is high.
 
	Factorization is carried out using Householder Reflections. Each reflection matrix is iteratively
	applied to the operator and the vector b to convert the linear system to upper triangular. Then, the
	system is solved using backwards substitution.
 
	\param matvec user supplied linear operator given as an int function
	\param b matrix of boundary conditions in the linear system Ax=b
	\param qr_dat pointer to the QR_DATA data structure
	\param matvec_data user supplied void pointer to a data structure needed for the linear operator*/

/**	\note int (*matvec) (const Matrix<double>& v, Matrix<double> &Av, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a linear operator. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 v that will act on the linear operator a modified the matrix entries of Av to form the
 result of a matrix-vector product. Void pointer data is used to pass any user data structure
 that the function may need in order to perform the linear operation.
	\n --------------------------------------------------------------------------------*/
int QRsolve( int (*matvec) (const Matrix<double>& x, Matrix<double> &Ax, const void *data),
			Matrix<double> &b, QR_DATA *qr_dat, const void *matvec_data);

/// Function to iteratively solve a non-linear system using the Picard or Fixed-Point method
/** This function iteratively solves a non-linear system using the Picard method. User supplies
	a residual function and a weak solution form function. The weak form function is used to
	approximate the next solution vector for the non-linear system and the residual function is
	used to determine convergence. User also supplies an initial guess to the non-linear system
	as a matix x, which will also be used to store the solution. This algorithm is very simple
	and may not be sufficient to solve complex non-linear systems.
 
	\param res user supplied function for the non-linear residuals of the system
	\param evalx user supplied function for the weak form to estimate the next solution
	\param x user supplied matrix holding the initial guess to the non-linear system
	\param picard_dat pointer to the PICARD_DATA data structure
	\param res_data user supplied void pointer to a data structure used for residual evaluations
	\param evalx_data user supplied void pointer to a data structure used for evaluation of weak form*/

/**	\note int (*res) (const Matrix<double>& x, Matrix<double> &F, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for the non-linear residuals. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 x representing the current non-linear variables. Those variables are used to evaluate the
 users functions and return the residuals in the matrix F. The void pointer data is a data
 structure provided by the user to hold information the function may need in order to form
 the residuals.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*evalx) (const Matrix<double>& x0, Matrix<double> &x, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function to approximate the next solution vector x based on
 the previous solution vector x0. The x0 matrix is passed to this function and must
 be used to edit the entries of x based on the weak form of the problem. The user
 is free to define any weak form approximation. Void pointer data is the users
 data structure that may be used to pass additional information into this function
 in order to evaluate the weak form. \n
 
 Example Residual: F(x) = x^2 + x - 1    Goal is to make this function equal zero \n
 Example Weak Form: x = 1 - x0^2         Rearrage residual to form a weak solution
	\n --------------------------------------------------------------------------------*/
int picard( int (*res) (const Matrix<double>& x, Matrix<double> &r, const void *data),
		   int (*evalx) (const Matrix<double>& x0, Matrix<double> &x, const void *data),
		   Matrix<double> &x, PICARD_DATA *picard_dat, const void *res_data,
		   const void *evalx_data);

/// Function to form a linear operator of a Jacobian matrix used along with the PJFNK method
/** This function is used in conjunction with the PJFNK routine to form a linear operator that
	a Krylov method can operate on. This linear operator is formed from the current residual
	vector of the non-linear iteration in PJFNK using a finite difference approximation. \n
	
	Jacobian Linear Operator: J*v = ( F(x_k + eps*v) - F(x_k) ) / eps
 
	\param v vector to be multiplied by the Jacobian matrix
	\param Jv storage vector for the result of the Jacobi-vector product
	\param data void pointer to the PJFNK_DATA data structure holding solver information*/
int jacvec( const Matrix<double>& v, Matrix<double>& Jv, const void *data);

/// Function to perform a Backtracking Line Search operation to smooth out convergence of PJFNK
/** This function performs a simple backtracking line search operation on the residuals from the
	PJFNK method. The step size of the non-linear iteration is checked against a level of tolerance
	for residual reduction, then adjusted down if necessary. This method always starts out with the
	maximum allowable step size. If the largest step size is fine, then the algorithm does nothing.
	Otherwise, it iteratively adjusts the step size down, until a suitable step is found. In the
	case that no suitable step is found, this algorithm will report failure to the PJFNK method and
	PJFNK will decide whether to continue trying to find a global minimum or report that it is stuck
	in a local minimum.
 
	\param feval user supplied residual function for the non-linear system
	\param Fkp1 vector holding the residuals for the next non-linear step
	\param xkp1 vector holding the solution for the next non-linear step
	\param pk vector holding the current non-linear search direction
	\param normFk value of the current non-linear residual
	\param backtrack_dat pointer to the BACKTRACK_DATA data structure
	\param feval_data user supplied void pointer to the data structure needed for residual evaluation*/

/**	\note int (*feval) (const Matrix<double>& x, Matrix<double> &F, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for the non-linear residuals. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 x representing the current non-linear variables. Those variables are used to evaluate the
 users functions and return the residuals in the matrix F. The void pointer data is a data
 structure provided by the user to hold information the function may need in order to form
 the residuals.
	\n --------------------------------------------------------------------------------*/
int backtrackLineSearch( int (*feval) (const Matrix<double>& x, Matrix<double> &F, const void *data),
						Matrix<double> &Fkp1, Matrix<double> &xkp1, Matrix<double> &pk, double normFk,
						BACKTRACK_DATA *backtrack_dat, const void *feval_data);

/// Function to perform the PJFNK algorithm to solve a non-linear system of equations
/** This function solves a non-linear system of equations using the Preconditioned Jacobian-
	Free Newton-Krylov (PJFNK) algorithm. Each non-linear step of this method results in a
	linear sub-problem that is solved iteratively with one of the Krylov methods in the
	krylov_method enum. User must supplied a residual function that computes the non-linear
	residuals of the system given the current state of the variables x. Additionally, the
	user must also supplied an initial guess to the non-linear system. Optionally, the user
	may supply a preconditioning function for the linear sub-problem. \n
 
	Basic Steps: (i) Calc F(x_k), (ii) Solve J(x_k)*s_k=-F(x_k) for s_k, (iii) Form x_kp1 = x_k + s_k
	
	\param res user supplied residual function for the non-linear system
	\param precon user supplied preconditioning function for the linear sub-problems
	\param x user supplied initial guess and storage location of the solution
	\param pjfnk_dat pointer to the PJFNK_DATA data structure
	\param res_data user supplied void pointer to data structure used in residual function
	\param precon_data user supplied void pointer to data structure used in preconditioning function*/

/**	\note int (*res) (const Matrix<double>& x, Matrix<double> &F, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for the non-linear residuals. User's function must return an
 int of 0 upon success and anything else denotes a failure. The function accepts a matrix
 x representing the current non-linear variables. Those variables are used to evaluate the
 users functions and return the residuals in the matrix F. The void pointer data is a data
 structure provided by the user to hold information the function may need in order to form
 the residuals.
	\n --------------------------------------------------------------------------------*/

/**	\note int (*precon) (const Matrix<double>& b, Matrix<double> &Mb, const void *data)
	\n --------------------------------------------------------------------------------
	\n This is a user supplied function for a preconditioning operator. It has the same form as
 the linear operators from the Krylov methods and should have all the same properties. The only
 difference is that this function must form an approximate matrix inversion on the jacvec
 linear operator and modify the entries of Mb to represent the result of that approximate
 matrix inversion. The matrix b is given as the vector that this operator is acting on and
 the void pointer data is for any user data structure that the operator may need.
	\n --------------------------------------------------------------------------------*/
int pjfnk( int (*res) (const Matrix<double>& x, Matrix<double> &F, const void *data),
		  int (*precon) (const Matrix<double>& r, Matrix<double>& p, const void *data),
		  Matrix<double> &x, PJFNK_DATA *pjfnk_dat, const void *res_data,
		  const void *precon_data );

/// Function to form a full numerical Jacobian matrix from a given non-linear function
/** This function uses finite differences to form a full rank Jacobian matrix for a user
	supplied non-linear function. The Jacobian matrix will be formed at the current state
	of the non-linear variables x and stored in a full matrix J. Integers Nx and Nf are
	used to determine the size of the Jacobian matrix.
 
	\param Func user supplied function for evaluation of the non-linear system
	\param x matrix holding the current value of the non-linear variables
	\param J matrix that will store the numerical Jacobian result
	\param Nx number of non-linear variables in the system
	\param Nf number of non-linear functions in the system
	\param jac_dat pointer to the NUM_JAC_DATA data structure
	\param user_data user supplied void pointer to a data structure used in the non-linear function*/
int NumericalJacobian( int (*Func) (const Matrix<double> &x, Matrix<double> &F, const void *user_data),
					  const Matrix<double> &x, Matrix<double> &J, int Nx, int Nf, NUM_JAC_DATA *jac_dat,
					  const void *user_data);

/// Function that runs a variety of tests on all the functions in LARK
/** This function runs a variety of tests on the linear and non-linear methods
	developed in LARK. It can be called from the UI.*/
int LARK_TESTS();

#endif
