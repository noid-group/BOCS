/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file solv_lin_eqns.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef SOLV_LIN
#define SOLV_LIN


/*************************************************************************************************************************************************************************************************
 Functions for solving the system with LU decomposition and in general matrix format such that preconditioning can be applied 
*************************************************************************************************************************************************************************************************/

int solv_lin_eqns(FILE * fp_log, tW_system * sys);

/* Solve with LU decomposition */
extern void dgesv_(int *N, int *N_rhs, double AT[], int *LDA, int pivot[],
		   double x[], int *N_row, int *info);

/* Solve with LU decomposition and get error estimates */
extern void dgesvx_(char *FACT, char *TRANS, int *N, int *N_rhs,
		    double A[], int *LDA, double AF[], int *LDAF,
		    int pivot[], char *EQUED, double R[], double C[],
		    double b[], int *LDB, double phi[], int *LDX,
		    double *RCOND, double FERR[], double BERR[],
		    double work[], int iwork[], int *info);

/* Solve with LU decomposition and get extensive error estimates */
//extern void dgesvxx_( char *FACT, char *TRANS, int *N, int *N_rhs, double MT[], int *LDA, double AF[], int *LDAF, int pivot[], char *EQUED, double R[], double C[], double b[], int *LDB, double phi[], int *LDX, double *RCOND, double *RPVGRW, double BERR[], int *N_ERR_BNDS, double ERR_BNDS_NORM[], double ERR_BNDS_COMP[], int *N_PARAMS, double PARAMS[], double work[], int iwork[], int *info );


/*************************************************************************************************************************************************************************************************
 Functions for solving the system with either UU decomposition or Cholesky decomposition with the matrix in packed storage form 
*************************************************************************************************************************************************************************************************/

int solv_lin_eqns_symm(FILE * fp_log, tW_system * sys);

/* Solve with UU decomposition */
extern void dspsv_(char *UPLO, int *N, int *NRHS, double A[], int IPIV[],
		   double B[], int *LDB, int *INFO);

/* Solve with UU decomposition and get extensive error estimates */
extern void dspsvx_(char *FACT, char *UPLO, int *N, int *N_rhs,
		    double MT[], double AF[], int pivot[], double b[],
		    int *LDB, double phi[], int *LDX, double *RCOND,
		    double FERR[], double BERR[], double work[],
		    int iwork[], int *info);

/* Solve with Cholesky decomposition */
extern void dppsv_(char *UPLO, int *N, int *NRHS, double A[], double B[],
		   int *LDB, int *INFO);

/* Solve with Cholesky decomposition and get extensive error estimates */
extern void dppsvx_(char *FACT, char *UPLO, int *N, int *NRHS, double A[],
		    double AF[], char *EQUED, double S[], double B[],
		    int *LDB, double X[], int *LDX, double *RCOND,
		    double FERR[], double BERR[], double work[],
		    int iwork[], int *INFO);

/* Solve with Cholesky decomposition and get extensive error estimates */
//extern void dposvxx_( char *FACT, char *UPLO, int *N, int *NRHS, double A[], int *LDA, double AF[], int *LDAF, char *EQUED, double S[], double B[], int *LDB, double X[], int *LDX, double *RCOND, double *RPVGRW, double BERR[], int *N_ERR_BNDS, double ERR_BNDS_NORM[], double ERR_BNDS_COMP[], int *NPARAMS, double PARAMS[], double WORK[], int IWORK[], int *INFO );

/*************************************************************************************************************************************************************************************************
 Functions for solving the system with using the Bayesian inference method
*************************************************************************************************************************************************************************************************/

int solv_lin_eqns_Bayesian(FILE * fp_log, tW_system * sys);

int solv_lin_eqns_Bayesian_saveparams(FILE * fp_log, tW_system * sys,
				      double *alpha_save,
				      double *beta_save);

int solv_lin_eqns_constrain_tabdih_Bayesian(FILE * fp_log,
					    tW_system * sys);

extern void dgetri_(int *N, double A[], int *LDA, int pivot[],
		    double work[], int *lwork, int *INFO);

extern void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
		   double *ALPHA, double A[], int *LDA, double B[],
		   int *LDB, double *BETA, double C[], int *LDC);

double calc_Chi2(FILE * fp_log, tW_system * sys, double phi[], double b[]);

/*************************************************************************************************************************************************************************************************
 Functions for solving the system with preconditioning by rescaling phi
*************************************************************************************************************************************************************************************************/

int solv_lin_eqns_rescalephi(FILE * fp_log, tW_system * sys);


/*************************************************************************************************************************************************************************************************
 Functions for solving the system with constraints on the dihedral angle
*************************************************************************************************************************************************************************************************/

/* JFR - 12.20.13: This function is for any regularization and/or constraints, rendering several of the other functions obsolete */
/* See the bottom of solv_lin_eqns.c for obsolete functions and details */
int solv_lin_eqns_constrain_tabdih_regularize(FILE * fp_log,
					      tW_system * sys,
					      double *alpha, double beta);

int solv_lin_eqns_constrain_tabdih(FILE * fp_log, tW_system * sys);

int solv_lin_eqns_constrain_tabdih_Bayesian_params(FILE * fp_log,
						   tW_system * sys,
						   double *alpha,
						   double beta);

int solv_lin_eqns_constrain_tabdih_2(FILE * fp_log, tW_system * sys);

int solv_lin_eqns_constrain_tabdih_2_Bayesian_params(FILE * fp_log,
						     tW_system * sys,
						     double *alpha,
						     double beta);

extern void dgglse_(int *M, int *N, int *P, double A[], int *LDA,
		    double B[], int *LDB, double C[], double D[],
		    double X[], double WORK[], int *LWORK, int *INFO);

/*************************************************************************************************************************************************************************************************
 Functions for solving the system with SVD or just calculating the SVD's  
*************************************************************************************************************************************************************************************************/
int SVD(FILE * fp_log, tW_system * sys);

extern void dgesdd_(char *JOBZ, int *M, int *N, double A[], int *LDA,
		    double S[], double U[], int *LDU, double VT[],
		    int *LDVT, double work[], int *lwork, int iwork[],
		    int *info);

int solveSVD(FILE * fp_log, tW_system * sys, double rcond, tW_word b_info);

extern void dgelss_(int *M, int *N, int *NRHS, double A[], int *LDA,
		    double B[], int *LDB, double S[], double *RCOND,
		    int *RANK, double WORK[], int *LWORK, int *INFO);

/*************************************************************************************************************************************************************************************************
 Functions for iterative inversion calculation 
*************************************************************************************************************************************************************************************************/

int solv_lin_eqns_PT(FILE * fp_log, tW_system * sys);

int get_b_soln_PT(FILE * fp_log, tW_system * sys, int order,
		  double phi_struct[]);

int PT_eigen(int k, int N, int N_eigen, double **A);

/*************************************************************************************************************************************************************************************************
 Functions for eigenspectrum analysis
*************************************************************************************************************************************************************************************************/

int eigen(FILE * fp_log, tW_system * sys, int flag_Gbar);

//DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
extern void dsytrd_(char *UPLO, int *N, double A[], int *LDA, double D[],
		    double E[], double Tau[], double work[], int *lwork,
		    int *info);

//DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK,    INFO )
extern void dorgtr_(char *UPLO, int *N, double A[], int *LDA, double Tau[],
		    double work[], int *lwork, int *info);

//DSTEQR( COMPZ, N, D, E, Z,    LDZ, WORK, INFO )
extern void dsteqr_(char *COMPZ, int *N, double D[], double E[],
		    double Z[], int *LDZ, double work[], int *info);

/*************************************************************************************************************************************************************************************************
  Functions for uniqeness analysis
*************************************************************************************************************************************************************************************************/

int calc_Calpha(int N, double f[], double MT[], double D[],
		double Calpha[]);

int calc_fn_bn(int N_M, int N_L, double Calpha[], double MT[], double D[],
	       tW_system * sys);

int print_Mn(int N_M, int L, double D[], double MT[], tW_system * sys);

int calc_IPR(int N, double MT[]);

/*************************************************************************************************************************************************************************************************
  Preconditioning
*************************************************************************************************************************************************************************************************/

int Precondition(int N, double *M, double *b, double *phi,
		 double *norm_col, double *row_max, int flag_restore,
		 int flag_solve, tW_system * sys);

int Right_PC(int N, double *M, double *b, double *phi, double *norm_col,
	     double *row_max, int flag_restore, int flag_solve);

int Left_PC(int N, double *M, double *b, double *phi, double *norm_col,
	    double *row_max, int flag_restore, int flag_solve);

#endif
