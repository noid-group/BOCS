/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file solv_lin_eqns.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Functions related to solving sets of linear equations with lapack libraries
*/

//c library includes
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//local includes
#include "gmx-interface.h"
#include "cgff_types.h"
#include "solv_lin_eqns.h"
#include "safe_mem.h"
#include "io_output.h"
#include "wnoid_math.h"

/*****************************************************************************************
solv_lin_eqns(): 
*****************************************************************************************/
int solv_lin_eqns(FILE * fp_log, tW_system * sys)
{
    int i, j;
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    int LDB = N;
    int LDX = N;
    int pivot[N];		// pivot indices for permutation matrix
    int info = 0;			// reports on success of the calculation
    double *MT;			// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    char FACT = sys->ERR_var.FACT;
    char EQUED;
    /* variables for preconditioning */
    double row_max[N], norm_col[N];
    /* variables for the timing things */
    time_t start, end;
    double diff;

    int inv_flag = 0;
    if (sys->ERR_var.flag_ERR == TRUE) {
	inv_flag = 1;
    }

    fprintf(fp_log, "In solv_lin_eqns.\n");

    MT = (double *) ecalloc(N * N, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index];
	}
    }

    /*  Precondition the matrix */
    Precondition(N, MT, b, phi, norm_col, row_max,
		 FALSE /* => don't restore */ ,
		 TRUE /* => operate on b and phi */ , sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	phi[i] = b[i];
    }

    if (inv_flag == 0) {


	if (strcmp(sys->SOLN_var.SOLN_METH, "LU") == 0) {	// if solution method is LU, solve by LU decomposition
	    /*  Solve the system using LU Decomposition  */
	    fprintf(fp_log, " Inverting the matrix using LU decomposition in default solver \n");
	    time(&start);
	    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
	    time(&end);
	    diff = difftime(end, start);
	    if (info != 0) {
		fprintf(fp_log, "WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
		fprintf(stderr, "WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
	    }

	    fprintf(fp_log, " dgesv_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));
	} else {  // otherwise, solve by SVD
	    /* JFR - 02.06.13: variable for solving the system with SVD */
	    int N2 = N;
	    int N3 = N;
	    int rank;
	    double *S;
	    int lwork = 10 * N;
	    double work[lwork];
	    double rcond = sys->SVD_var.rcond;

	    S = (double *) ecalloc(N, sizeof(double));
	    if (S == NULL) {
		printf("STOP. S ecalloc failed. N: %d\n", N);
		exit(0);
	    }

	    /* Solve the overdetermined system using Singular Value Decomposition, setting all singular values under rcond to zero  */
	    fprintf(fp_log, " Inverting the matrix using Singular Value Decomposition in default solver \n");
	    time(&start);
	    dgelss_(&N, &N2, &N_rhs, MT, &LDA, phi, &N3, S, &rcond, &rank, work, &lwork, &info);
	    time(&end);
	    diff = difftime(end, start);
	    if (info != 0) {
		fprintf(fp_log, "WARNING: info = %d after dgelss_().  The inverse condition number may be less than machine precision \n", info);
		fprintf(stderr, "WARNING: info = %d after dgelss_().  The inverse condition number may be less than machine precision \n", info);
	    }
	    fprintf(fp_log, " dgelss_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));
	    fprintf(fp_log, " Rank(G) = %d \n", rank);
	    fprintf(fp_log, "\n");

	    free(S);
	} 
    } else if (inv_flag == 1) {
	/*  Solve the system using LU Decomposition and compute error estimates  */

	/* variables for matrix inversion with error estimates */
	char TRANS = 'N';
	int LDAF = N;
	int iwork[N];
	double work[4 * N];
	double AF[N * N];
	double R[N];
	double C[N];
	double RCOND;
	double FERR[1];
	double BERR[1];

	fprintf(fp_log, "Inverting the matrix using LU decomposition with error estimates \n");
	time(&start);
	dgesvx_(&FACT, &TRANS, &N, &N_rhs, MT, &LDA, AF, &LDAF, pivot, &EQUED, R, C, b, &LDB, phi, &LDX, &RCOND, FERR, BERR, work, iwork, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log, "WARNING: info = %d after dgesvx_().  The inverse condition number may be less than machine precision \n", info);
	    fprintf(stderr, "WARNING: info = %d after dgesvx_().  The inverse condition number may be less than machine precision \n", info);
	}
	fprintf(fp_log,"dgesvx_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));
	fprintf(fp_log, " RCOND = %.8e \n", RCOND);
	fprintf(fp_log, " FERR = %.8e \n", FERR[0]);
	fprintf(fp_log, " BERR = %.8e \n", BERR[0]);
	fprintf(fp_log, "\n");
    }

    else 
    {
	fprintf(fp_log, "ERROR: None of the matrix inversion options were valid!  The matrix problem will not be solved \n");
    }

    /* Undo the Preconditioning */
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE /* => restore */ ,
		 TRUE /* => operate on b and phi */ , sys);

    free(MT);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_symm(): 
*****************************************************************************************/
int solv_lin_eqns_symm(FILE * fp_log, tW_system * sys)
{
    int i;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDB = N;
    int LDX = N;
    int pivot[N];		// pivot indices for permutation matrix
    int info = 0;			// reports on success of the calculation
    double *MT;			// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    char UPLO = 'L';
    char FACT = sys->ERR_var.FACT;
    char EQUED;
    /* variables for the timing things */
    time_t start, end;
    double diff;

    int inv_flag = -1;

    if (strcmp(sys->SOLN_var.SOLN_METH, UU) == 0) {
	if (sys->ERR_var.flag_ERR == FALSE) {
	    inv_flag = 0;
	} else {
	    inv_flag = 1;
	}
    } else if (strcmp(sys->SOLN_var.SOLN_METH, CHOLESKY) == 0) {
	if (sys->ERR_var.flag_ERR == FALSE) {
	    inv_flag = 2;
	} else {
	    inv_flag = 3;
	}
    }

    fprintf(fp_log, "In solv_lin_eqns.\n");

    if (sys->MEM_var.flag_LOWMEM == FALSE) {
	MT = (double *) ecalloc(sys->N_pack, sizeof(double));
	if (MT == NULL) {
	    printf("STOP. MT ecalloc failed. N_pack: %d\n", sys->N_pack);
	    exit(0);
	}
    } else {
	MT = M;
    }

    /*  Precondition the matrix */
    /* Currently, there are no preconditiong options for the symmetric matrix */

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	phi[i] = b[i];
    }

    for (i = 0; i < sys->N_pack; i++) {
	MT[i] = M[i];
    }

    if (inv_flag == 0) {
	/*  Solve the system using UU Decomposition  */
	fprintf(fp_log, " Inverting the matrix using UU decomposition \n");
	time(&start);
	dspsv_(&UPLO, &N, &N_rhs, MT, pivot, phi, &LDB, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log, "WARNING: info = %d after dsysv_().  The inverse condition number may be less than machine precision \n", info);
	    fprintf(stderr, "WARNING: info = %d after dsysv_().  The inverse condition number may be less than machine precision \n", info);
	}
	fprintf(fp_log," dspsv_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));

    } else if (inv_flag == 1) {
	/*  Solve the system using UU Decomposition and compute error estimates  */

	/* variables for matrix inversion with error estimates */
	int iwork[N];
	double work[3 * N];
	double AF[(N * N + N) / 2];
	double RCOND;
	double FERR[1];
	double BERR[1];

	fprintf(fp_log,
		" Inverting the matrix using UU decomposition with error estimates \n");
	time(&start);
	dspsvx_(&FACT, &UPLO, &N, &N_rhs, MT, AF, pivot, b, &LDB, phi,
		&LDX, &RCOND, FERR, BERR, work, iwork, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log, "WARNING: info = %d after dsysvx_().  The inverse condition number may be less than machine precision \n", info);
	    fprintf(stderr, "WARNING: info = %d after dsysvx_().  The inverse condition number may be less than machine precision \n", info);
	}
	fprintf(fp_log,	" dspsvx_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));
	fprintf(fp_log, " RCOND = %.8e \n", RCOND);
	fprintf(fp_log, " FERR = %.8e \n", FERR[0]);
	fprintf(fp_log, " BERR = %.8e \n", BERR[0]);
	fprintf(fp_log, "\n");
    }

    else if (inv_flag == 2) {
	/*  Solve the system using Cholesky Decomposition  */
	fprintf(fp_log, " Inverting the matrix using Cholesky decomposition \n");
	time(&start);
	dppsv_(&UPLO, &N, &N_rhs, MT, phi, &LDB, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log, "WARNING: info = %d after dposv_().  The inverse condition number may be less than machine precision \n", info);
	    fprintf(stderr, "WARNING: info = %d after dposv_().  The inverse condition number may be less than machine precision \n", info);
	}
	fprintf(fp_log, " dppsv_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));

    } else if (inv_flag == 3) {
	/*  Solve the system using Cholesky Decomposition and compute error estimates  */

	/* variables for matrix inversion with error estimates */
	int LDAF = N;
	int iwork[N];
	double work[3 * N];
	double AF[LDAF * N];
	double RCOND;
	double FERR[1];
	double BERR[1];
	double S[N];

	fprintf(fp_log, " Inverting the matrix using Cholesky decomposition with error estimates \n");
	time(&start);
	dppsvx_(&FACT, &UPLO, &N, &N_rhs, MT, AF, &EQUED, S, b, &LDB, phi, &LDX, &RCOND, FERR, BERR, work, iwork, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log, "WARNING: info = %d after dposvx_().  The inverse condition number may be less than machine precision \n", info);
	    fprintf(stderr, "WARNING: info = %d after dposvx_().  The inverse condition number may be less than machine precision \n", info);
	}
	fprintf(fp_log,	" dppsvx_() took %lf seconds or %lf minutes or %lf hours \n",diff, diff / (60.0), diff / (3600.00));
	fprintf(fp_log, " RCOND = %.8e \n", RCOND);
	fprintf(fp_log, " FERR = %.8e \n", FERR[0]);
	fprintf(fp_log, " BERR = %.8e \n", BERR[0]);
	fprintf(fp_log, "\n");

    } else {
	fprintf(fp_log, "ERROR: None of the matrix inversion options were valid!  The matrix problem will not be solved \n");
    }

    free(MT);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_Bayesian_saveparams(): same as normal but the regularization parameters
are stored
*****************************************************************************************/
int solv_lin_eqns_Bayesian_saveparams(FILE * fp_log, tW_system * sys,
				      double *alpha_save,
				      double *beta_save)
{
    int i, j;
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N2 = N;
    int N3 = N;
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    int LDB = N;
    int LDC = N;
    int pivot[N];		// pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT, *MT2, *MT3;	// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double row_max[N], norm_col[N];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    int lwork = 10 * N;
    double work[lwork];
    /* variables for regularization */
    double alpha[N * N];
    double alpha_old[N * N];
    double beta = 0;
    double beta_old = 0;
    double tau_alpha = sys->REG_var.tau_alpha;
    double tau_beta = sys->REG_var.tau_beta;
    double f2 = 0;
    int nt = sys->REG_var.Nframes;
    int NT = 3 * nt * N;
    int Nmax = sys->REG_var.Nmax;
    int ctr;
    double trace;
    int flag_conv;
    int relerr;
    char TRANSA = 'N';
    char TRANSB = 'N';
    double aconst = 1.0;
    double bconst = 0.0;

    fprintf(fp_log, "In solv_lin_eqns.\n");

    MT = (double *) ecalloc(N * N, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }
    MT2 = (double *) ecalloc(N * N, sizeof(double));
    if (MT2 == NULL) {
	printf("STOP. MT2 ecalloc failed. N: %d\n", N);
	exit(0);
    }
    MT3 = (double *) ecalloc(N * N, sizeof(double));
    if (MT3 == NULL) {
	printf("STOP. MT3 ecalloc failed. N: %d\n", N);
	exit(0);
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index];
	}
    }

    /* Do the normal calculation */

    /*  Precondition the matrix */
    Precondition(N, MT, b, phi, norm_col, row_max,
		 FALSE /* => don't restore */ ,
		 TRUE /* => operate on b and phi */ , sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	phi[i] = b[i];
    }

    /*  Solve the system using LU Decomposition  */
    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
    time(&start);
    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log, "WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
	fprintf(stderr, "WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
    }
    fprintf(fp_log, " dgesv_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);
    /* Now calculate the input values of alpha and beta */

    for (i = 0; i < N; i++) {
	f2 += phi[i] * phi[i];
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    if (i == j) {
		alpha[i + N * j] = N / f2;
	    } else {
		alpha[i + N * j] = 0.0;
	    }
	}
    }

    beta = NT / calc_Chi2(fp_log, sys, phi, b);

    /* print out the initial regularization parameters */
    fprintf(fp_log, " Initial Regularization Parameters: \n");
    for (i = 0; i < N; i++) {
	fprintf(fp_log, " alpha_0[%d] = %lf \n", i, alpha[i + N * i]);
    }
    fprintf(fp_log, " beta_0 = %lf \n", beta);

    ctr = 0;
    do {

	/* Prepare to calculate mN */
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		index = index_Lpacked(i, j, N);
		MT[i + N * j] = M[index] + alpha[i + N * j] / beta;
	    }
	}

	/*  Precondition the matrix */
	//Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

	for (i = 0; i < N; i++) {
	    /* Copy b -> x which will be overwritten with soln. */
	    phi[i] = b[i];
	}

	/*  Solve the regularized system using LU Decomposition  */
	fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
	time(&start);
	dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log, "WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
	    fprintf(stderr, "WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
	}

	fprintf(fp_log,
		" dgesv__() took %lf seconds or %lf minutes or %lf hours \n",
		diff, diff / (60.0), diff / (3600.00));

	/* Undo the Preconditioning */
	//Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

	/* Calculate the inverse of the regularized matrix */
	dgetri_(&N, MT, &LDA, pivot, work, &lwork, &info);

	/* Calculate the new alpha and beta quantities */
	for (i = 0; i < N; i++) {
	    alpha_old[i + N * i] = alpha[i + N * i];
	    alpha[i + N * i] =
		(1.0 / (phi[i] * phi[i] + (MT[i + N * i] / beta)));
	}

	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		index = index_Lpacked(i, j, N);
		MT2[i + N * j] = M[index];
		MT3[i + N * j] = 0.0;
	    }
	}
	/* Do matrix multiplication */
	dgemm_(&TRANSA, &TRANSB, &N, &N2, &N3, &aconst, MT, &LDA, MT2,
	       &LDB, &bconst, MT3, &LDC);

	trace = 0.0;
	for (i = 0; i < N; i++) {
	    trace += beta * MT3[i + N * i];
	}

	beta_old = beta;
	beta = (NT - trace) / (calc_Chi2(fp_log, sys, phi, b));

	flag_conv = TRUE;
	for (i = 0; i < N; i++) {
	    relerr =
		fabs((alpha[i + N * i] -
		      alpha_old[i + N * i]) / alpha_old[i + N * i]);
	    if (relerr > tau_alpha) {
		flag_conv = FALSE;
	    }
	}

	relerr = fabs((beta - beta_old) / beta_old);
	if (relerr > tau_beta) {
	    flag_conv = FALSE;
	}

	ctr++;

    } while ((flag_conv == FALSE) && (ctr < Nmax));

    /* print out and store the final regularization parameters */
    fprintf(fp_log, " Final Regularization Parameters: \n");
    for (i = 0; i < N; i++) {
	alpha_save[i] = alpha[i + N * i];
	fprintf(fp_log, " alpha_f[%d] = %lf \n", i, alpha[i + N * i]);
    }
    fprintf(fp_log, " beta_f = %lf, nt = %d \n", beta, nt);
    *beta_save = beta;
    if (ctr < Nmax) {
	fprintf(fp_log,	" Optimization of regularization parameters converged after %d steps", ctr);
    } else {
	fprintf(fp_log,	" Optimization of regularization parameters did not converge after Nmax = %d steps", Nmax);
    }

    /* Don't need to solve the system a final time (below) but I will leave it for now */

    /* Prepare to calculate the FINAL mN */
    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index] + alpha[i + N * j] / beta;
	}
    }

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	phi[i] = b[i];
    }

    /*  Solve the regularized system using LU Decomposition  */
    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
    time(&start);
    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,	"WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
	fprintf(stderr,	"WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
    }
    fprintf(fp_log, " dgesv_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

    free(MT);
    free(MT2);
    free(MT3);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_constrain_tabdih_regularize(): 
*****************************************************************************************/
int solv_lin_eqns_constrain_tabdih_regularize(FILE * fp_log,
					      tW_system * sys,
					      double *alpha, double beta)
{
    int i, j;
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_row = N;
    int N_col = N;
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    int LDB = N;
    int pivot[N];		// pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT;			// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double row_max[N_row], norm_col[N_col];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    /* variables for constraining dihedrals */
    double lambda = 0.00;
    double dlambda = 0.00;
    double tau_dih = 0.00;
    int flag_conv = TRUE;
    int Nmax = 1;
    if (sys->CONSTRAIN_var.flag_CONSTRAIN == TRUE) {
	lambda = sys->CONSTRAIN_var.lambda;
	dlambda = sys->CONSTRAIN_var.dlambda;
	tau_dih = sys->CONSTRAIN_var.tau_dih;
	Nmax = sys->CONSTRAIN_var.Nmax;
    }
    int k, kk, ii;
    int ctr;
    int i_flag, row_i, row_k;
    /* variables for regularization */
    double *C;

    fprintf(fp_log, "In solv_lin_eqns.\n");

    fprintf(fp_log, "In solv_lin_eqns REG and CONST, tau_dih = %lf .\n", tau_dih);

    /* JFR - 02.06.13: variable for solving the system with SVD */
    int N2 = N;
    int N3 = N;
    int rank;
    double *S;
    int lwork = 10 * N;
    double work[lwork];
    double rcond = sys->SVD_var.rcond;
    FILE *fp_S;

    fp_S = fopen("S_CONST.dat", "w");

    S = (double *) ecalloc(N, sizeof(double));
    if (S == NULL) {
	printf("STOP. S ecalloc failed. N: %d\n", N);
	exit(0);
    }

    MT = (double *) ecalloc(N_row * N_col, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }
    C = (double *) ecalloc(N_row * N_col, sizeof(double));
    if (C == NULL) {
	printf("STOP. C ecalloc failed. N: %d\n", N);
	exit(0);
    }

    ctr = 0;
    do {
	if (sys->REG_var.flag_REG == TRUE) {
	    if (strcmp(sys->REG_var.type, "BAYES") == 0) {	/* Regularize with Bayesian inference method */
		for (j = 0; j < N; j++) {
		    for (i = 0; i < N; i++) {
			if (i == j) {
			    C[j + N * j] = alpha[i] / beta;
			}
		    }
		}
	    } else if (strcmp(sys->REG_var.type, "UNCERT") == 0) {	/* Regularize with simple scaling of MT uncertainties */
		for (j = 0; j < N; j++) {
		    for (i = 0; i < N; i++) {
			index = index_Lpacked(i, j, N);
			if (sys->M_cnt[index] > 1.0) {
			    C[j + N * j] +=
				sqrt(sys->d2M[index]) /
				sqrt(((sys->M_cnt[index] - 1.0)));
			}
		    }
		    C[j + N * j] *= beta;
		}

		/* Only regularize nb and intra-nb interactions, for now only do this for the (harsher) uncertainty regularization. */
		/* To extend to BAYES, just move this out of the if statements */
		for (i = 0; i < sys->N_Inter_Types; i++) {
		    i_flag =
			strcmp(B_DIHEDRAL, sys->Inter_Types[i].inter_type);
		    for (ii = 0; ii < sys->Inter_Types[i].N_pts; ii++) {
			row_i = sys->Inter_Types[i].i_0 + ii;
			if ((strcmp
			     (B_NB_PAIR_BOND,
			      sys->Inter_Types[i].inter_type) != 0)
			    &&
			    (strcmp
			     (NB_PAIR,
			      sys->Inter_Types[i].inter_type) != 0)) {
			    C[row_i + N * (row_i)] = 0.00;
			}
		    }
		}
	    } else {
		printf("ERROR: Regularization type incorrect");
		exit(0);
	    }
	}

	if (sys->CONSTRAIN_var.flag_CONSTRAIN == TRUE) {	/* constrain tabulated dihedrals */
	    for (i = 0; i < sys->N_Inter_Types; i++) {
		if (strcmp(B_DIHEDRAL, sys->Inter_Types[i].inter_type) !=
		    0) {
		    continue;
		}		/* check that i is a dihedral interaction */
		if (strcmp(TOY_DIHED_NAME, sys->Inter_Types[i].basis) == 0) {
		    continue;
		}		/* check that i is not using a TOY basis */
		if (strcmp
		    (RYCKAERT_BELLEMANS_NAME,
		     sys->Inter_Types[i].basis) == 0) {
		    continue;
		}		/* check that i is not using a RB basis */
		for (k = 0; k < sys->N_Inter_Types; k++) {
		    if (strcmp(B_DIHEDRAL, sys->Inter_Types[k].inter_type)
			!= 0) {
			continue;
		    }		/* check that k is a dihedral interaction */
		    if (strcmp(TOY_DIHED_NAME, sys->Inter_Types[k].basis)
			== 0) {
			continue;
		    }		/* check that k is not using a TOY basis */
		    if (strcmp
			(RYCKAERT_BELLEMANS_NAME,
			 sys->Inter_Types[k].basis) == 0) {
			continue;
		    }		/* check that k is not using a RB basis */
		    if (strcmp
			(sys->Inter_Types[i].inter_name,
			 sys->Inter_Types[k].inter_name) != 0) {
			continue;
		    }		/*  check that i and k are the same dihedral */
		    for (ii = 0; ii < sys->Inter_Types[i].N_pts; ii++) {
			row_i = sys->Inter_Types[i].i_0 + ii;
			for (kk = 0; kk < sys->Inter_Types[k].N_pts; kk++) {
			    row_k = sys->Inter_Types[k].i_0 + kk;
			    C[row_i + N * (row_k)] += lambda;
//fprintf( fp_log, "constraining row_i = %d and row_k = %d .\n", row_i, row_k );
			}
		    }
		}
	    }
	}

	/* Update the metric tensor with regularization and order constraint info */
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N_col; j++) {
		index = index_Lpacked(i, j, N);
		MT[i + N * j] = M[index] + C[i + N * j];
	    }
	}

	/*  Precondition the matrix */
	//Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

	for (i = 0; i < N; i++) {
	    /* Copy b -> x which will be overwritten with soln. */
	    phi[i] = b[i];
	}

	if (strcmp(sys->SOLN_var.SOLN_METH, "LU") == 0) {	
	    /*  Solve the system using LU Decomposition  */
	    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
	    time(&start);
	    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
	    time(&end);
	    diff = difftime(end, start);
	    if (info != 0) {
		fprintf(fp_log, "WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
		fprintf(stderr, "WARNING: info = %d after dgesv_().  The inverse condition number may be less than machine precision \n", info);
	    }
	    fprintf(fp_log, " dgesv_() took %lf seconds or %lf minutes or %lf hours \n", diff, diff / (60.0), diff / (3600.00));

	} else {
	    /* Solve the overdetermined system using Singular Value Decomposition, setting all singular values under rcond to zero  */
	    fprintf(fp_log, " Inverting the matrix using Singular Value Decomposition in tabdih constrain reg\n");
	    time(&start);
	    dgelss_(&N, &N2, &N_rhs, MT, &LDA, phi, &N3, S, &rcond, &rank, work, &lwork, &info);
	    time(&end);
	    diff = difftime(end, start);
	    if (info != 0) {
		fprintf(fp_log, "WARNING: info = %d after dgelss_().  The inverse condition number may be less than machine precision \n", info);
		fprintf(stderr, "WARNING: info = %d after dgelss_().  The inverse condition number may be less than machine precision \n", info);
	    }
	    fprintf(fp_log,
		    " dgelss_() took %lf seconds or %lf minutes or %lf hours \n",
		    diff, diff / (60.0), diff / (3600.00));
	    fprintf(fp_log, " Rank(G) = %d \n", rank);
	    fprintf(fp_log, "\n");

	    /* print out the singular values */
	    for (i = 0; i < N; i++) {
		fprintf(fp_S, " %.16lf \n", S[i]);
	    }
	} 
	/* Undo the Preconditioning */
	//Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

	/* check that the matrix was not rank deficient */
	if (info != 0) {
	    ctr++;
	    lambda = lambda * dlambda;
	    flag_conv = FALSE;
	    continue;
	}

	/* increase lambda by dlambda */
	lambda = lambda * dlambda;

	ctr++;

    } while ((flag_conv == FALSE) && (ctr < Nmax));

    /* print out the final lambda scaling terms */
    fprintf(fp_log, " lambda = %lf \n", lambda / dlambda);
    if (ctr < Nmax) {
	fprintf(fp_log, " Optimization of lambda parameter converged after %d steps", ctr);
    } else {
	fprintf(fp_log,	" Optimization of lambda parameter did not converge after Nmax = %d steps", Nmax);
    }

    fclose(fp_S);

    free(MT);
    free(S);

    return info;
}

/*****************************************************************************************
solveSVD(): 
*****************************************************************************************/
int solveSVD(FILE * fp_log, tW_system * sys, double rcond, tW_word b_info)
{

    int i, j;
    int count;
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = sys->N_coeff;	// leading dimension of A as stored in AT.  NB A is square
    int info, rank;
    double *MT;
    double row_max[N], norm_col[N];
    double *b = sys->b;
    double *phi;
    double *M = sys->M;
    double *S;
    int lwork = 10 * N;
    double work[lwork];
    /* variables for the timing things */
    time_t start, end;
    double diff;

    FILE *fp_phi;		/* file for printing out the soln */
    tW_word fname;

    int index;

    fprintf(fp_log, "In solv_lin_eqns.\n");

    MT = (double *) ecalloc(N * N, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }

    phi = (double *) ecalloc(N, sizeof(double));

    S = (double *) ecalloc(N, sizeof(double));

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, sys->N_coeff);
	    MT[i + N * j] = M[index];
	}
    }

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	phi[i] = b[i];
    }

    /* Solve the overdetermined system using Singular Value Decomposition, setting all singular values under rcond to zero  */
    fprintf(fp_log,
	    " Inverting the matrix using Singular Value Decomposition in solveSVD \n");
    time(&start);
    dgelss_(&N, &N, &N_rhs, MT, &LDA, phi, &N, S, &rcond, &rank, work,
	    &lwork, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgelss_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgelss_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));
    fprintf(fp_log, " Rank(G) = %d \n", rank);
    fprintf(fp_log, "\n");

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

    /* Print out phi for the SVD soln with given rcond */
    count = 0;
    for (i = 0; i < sys->N_Inter_Types; i++) {
	sprintf(fname, "f_%s_SVD.%s.rcond%2.1e.dat", b_info,
		sys->Inter_Types[i].inter_name, rcond);
	fp_phi = fopen(fname, "w");

	for (j = 0; j < sys->Inter_Types[i].N_pts; j++) {
	    fprintf(fp_phi, "%10.5lf    %15.5lf \n",
		    sys->Inter_Types[i].ptr_x[j], phi[count]);
	    count++;
	}

	fclose(fp_phi);
    }

    free(MT);
    free(phi);
    free(S);

    return info;

}

/*****************************************************************************************
SVD(): 
*****************************************************************************************/
int SVD(FILE * fp_log, tW_system * sys)
{
    int i, j;
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N2 = N;
    int LDA = N;		// leading dimension of A as stored in AT.  NB A is square
    int LDU = N;		// leading dimension of U as stored in U.  NB U is square
    int LDVT = N;		// leading dimension of V as stored in VT.  NB V is square
    int info = 0;
    double *MT;
    double *S, *U, *VT;
    double row_max[N];
    double norm_col[N];
    double *M = sys->M;
    int lwork;
    if (N < 7) {
	lwork = 4 * N * N + 7 * N;
    } else {
	lwork = 5 * N * N;
    }
    double work[lwork];
    int iwork[8 * N];
    char jobz = 'A';		//A - compute and return U and V

    double *b = sys->b;		/* I need this for preconditioning now, I won't change it though */
    double *phi = NULL;

    tW_word fname;
    FILE *fp_U, *fp_V;
    FILE *fp_S;

    int index;

    fp_S = fopen("S.dat", "w");	//the singular values of M

    fprintf(fp_log, "In SVD.\n");

    MT = (double *) ecalloc(N * N, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }

    S = (double *) ecalloc(N, sizeof(double));
    U = (double *) ecalloc(N * N, sizeof(double));
    VT = (double *) ecalloc(N * N, sizeof(double));

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, sys->N_coeff);
	    MT[i + N * j] = M[index];
	}
    }

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, FALSE /* => don't operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, FALSE, sys);

    dgesdd_(&jobz, &N2, &N, MT, &LDA, S, U, &LDU, VT, &LDVT, work, &lwork,
	    iwork, &info);

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, FALSE /* => don't operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, FALSE, sys);

    /* print out the SVD */
    for (i = 0; i < N; i++) {

	if (sys->SVD_var.flag_printevecs == TRUE) {
	    sprintf(fname, "U.%d.dat", i);
	    fp_U = fopen(fname, "w");	//M col of MxM orthogonal matrix

	    sprintf(fname, "V.%d.dat", i);
	    fp_V = fopen(fname, "w");	//M rows of MxM orthogonal matrix
	}

	fprintf(fp_S, " %.16lf \n", S[i]);

	for (j = 0; j < N; j++) {
	    //M[i][j] = row_max[i] * M[i][j];  

	    if (sys->SVD_var.flag_printevecs == TRUE) {
		fprintf(fp_U, " %lf \n", U[i + j * N]);
		fprintf(fp_V, " %lf \n", VT[i * N + j]);
	    }
	    //if( j == N - 1 )
	    //{
	    //  fprintf( fp_U, " \n" );
	    //  fprintf( fp_V, " \n" );
	    //}

	}
	if (sys->SVD_var.flag_printevecs == TRUE) {
	    fclose(fp_U);
	    fclose(fp_V);
	}
    }

    //free( MT ); /* This is destroyed by dgesdd_ */
    free(S);
    free(U);
    free(VT);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_PT(): 
*****************************************************************************************/
int solv_lin_eqns_PT(FILE * fp_log, tW_system * sys)
{
    int i, j, k, l;
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    double *phi;
    double *M = sys->M;
    double *b = sys->b;
    double *g = sys->g;
    double **H1, **HN, **HN_tmp;	/* H^N = (G2inv*(-1.0*G3))^N */
    double *M2inv;
    double tmp;
    tW_word fname;

    int index;

    int N_PT = sys->PT_var.N_PT;

    fprintf(fp_log, "In solv_lin_eqns.\n");

    fprintf(fp_log, "N_PT = %d \n", sys->PT_var.N_PT);

    phi = (double *) ecalloc(N, sizeof(double));
    M2inv = (double *) ecalloc(N, sizeof(double));
    H1 = (double **) ecalloc(N, sizeof(double *));
    HN = (double **) ecalloc(N, sizeof(double *));
    HN_tmp = (double **) ecalloc(N, sizeof(double *));
    for (i = 0; i < N; i++) {
	H1[i] = (double *) ecalloc(N, sizeof(double));
	HN[i] = (double *) ecalloc(N, sizeof(double));
	HN_tmp[i] = (double *) ecalloc(N, sizeof(double));
    }

    /* decompose the matrix into 2 and 3-body terms */
    for (i = 0; i < N; i++) {
	M2inv[i] = 1.0 / g[i];
	HN_tmp[i][i] = M2inv[i];
    }

    sprintf(fname, "Minv.%dorder", 0);
    print_sep_M(fname, HN_tmp, sys);

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {

	    index = index_Lpacked(i, j, sys->N_coeff);

	    if (i == j) {
		H1[i][j] = M2inv[i] * (M[index] - g[i]);	/* This only works for the delta basis!!! */
		HN[i][j] = H1[i][j];
	    } else {
		H1[i][j] = M2inv[i] * (M[index]);
		HN[i][j] = H1[i][j];
	    }

	}
    }

    sprintf(fname, "Minv.%dorder", 1);
    print_sep_M(fname, H1, sys);

    /* k = 0 case */
    if (sys->PT_var.flag_eigen == TRUE) {
	PT_eigen(0, N, sys->PT_var.N_eigen, HN_tmp);
    }

    for (i = 0; i < N; i++) {
	phi[i] = b[i] * M2inv[i];

	/* set diagonal of HN_tmp back to zero */
	HN_tmp[i][i] = 0.00;
    }

    if (sys->PT_var.flag_MMOTF_SEP == TRUE) {
	//JFR - get the contribution to b from this order and decompose into each interaction
	get_b_soln_PT(fp_log, sys, 0, phi);
    }
    //JFR - print out the forces for this order
    sprintf(fname, "%s.%dorder", "f_forces", 0);
    print_sep_forces(fname, phi, sys);

    /* calculate the kth order solution */
    for (k = 1; k < N_PT; k++) {

	if (k % sys->PT_var.dPT == 0)	// print out results every dPT 
	{
	    if (sys->PT_var.flag_eigen == TRUE) {
		/* calculate the eigenspectrum of the kth order matrix */
		PT_eigen(k, N, sys->PT_var.N_eigen, HN);
	    }
	}

	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		if ((k % 2) == 0) {	/* k is even */
		    phi[i] += HN[i][j] * M2inv[j] * b[j];
		} else {	/* k is odd */

		    phi[i] -= HN[i][j] * M2inv[j] * b[j];
		}
	    }

	}

	if (k % sys->PT_var.dPT == 0)	// print out results every dPT 
	{

	    if (sys->PT_var.flag_MMOTF_SEP == TRUE) {
		//JFR - get the contribution to b from this order and decompose into each interaction
		get_b_soln_PT(fp_log, sys, k, phi);
	    }
	    //JFR - print out the forces for this order
	    sprintf(fname, "%s.%dorder", "f_forces", k);
	    print_sep_forces(fname, phi, sys);

	}

	/* update HN for the next order */
	if (k != (N_PT - 1)) {
	    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
		    tmp = 0.00;
		    for (l = 0; l < N; l++) {
			tmp += HN[i][l] * H1[l][j];
		    }
		    HN_tmp[i][j] = tmp;
		}
	    }

	    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
		    HN[i][j] = HN_tmp[i][j];
		    HN_tmp[i][j] *= M2inv[j];	/* So I can print out the full matrix for this order */
		}
	    }
	    if ((k + 1) % sys->PT_var.dPT == 0)	// print out results every dPT 
	    {
		sprintf(fname, "Minv.%dorder", k + 1);
		print_sep_M(fname, HN_tmp, sys);
	    }

	}

    }

    free(M2inv);
    free(HN);
    free(HN_tmp);

    return 0;
}

/*****************************************************************************************
get_b_soln_PT(): multiply the obtained forces with the original matrix to obtain the solns, b
*****************************************************************************************/
int get_b_soln_PT(FILE * fp_log, tW_system * sys, int order,
		  double phi_struct[])
{
    int i, j, k, l;
    double b_soln_struct = 0.00;
    //double b_soln_forces = 0.00;
    //double * phi_struct = sys->phi_struct;
    //double * phi_forces = sys->phi_forces;
    double *M = sys->M;
    double *g = sys->g;
    FILE *fp_soln_struct, *fp_soln_struct_diag;
    //FILE * fp_soln_forces, fp_soln_forces_diag;
    tW_word fname;

    int index;

    fprintf(fp_log, "In get_b_soln.\n");

    //JFR - Separate the contributions to b
    for (i = 0; i < sys->N_Inter_Types; i++) {

	sprintf(fname, "%s.%s.%dorder.dat", "b_soln_struct_diag",
		sys->Inter_Types[i].inter_name, order);
	fp_soln_struct_diag = fopen(fname, "w");

	//sprintf( fname, "%s.%s.dat", "b_soln_forces_diag", sys->Inter_Types[i].inter_name );
	//fp_soln_forces_diag = fopen( fname, "w" );

	for (j = 0; j < sys->N_Inter_Types; j++) {

	    sprintf(fname, "%s.%s.%s.%dorder.dat", "b_soln_struct",
		    sys->Inter_Types[i].inter_name,
		    sys->Inter_Types[j].inter_name, order);
	    fp_soln_struct = fopen(fname, "w");

	    //sprintf( fname, "%s.%s.%s.dat", "b_soln_forces", sys->Inter_Types[i].inter_name, sys->Inter_Types[j].inter_name );
	    //fp_soln_forces = fopen( fname, "w" );

	    for (k = 0; k < sys->Inter_Types[i].N_pts; k++) {

		if ((i == j)) {

		    fprintf(fp_soln_struct_diag, "%10.5lf    %15.5lf \n",
			    sys->Inter_Types[i].ptr_x[k],
			    (g[sys->Inter_Types[i].i_0 + k] *
			     phi_struct[sys->Inter_Types[i].i_0 +
					k]) / sys->Inter_Types[i].dr);
		    //fprintf( fp_soln_forces_diag, "%10.5lf    %15.5lf \n", sys->Inter_Types[i].ptr_x[k], ( g[sys->Inter_Types[i].i_0 + k] * phi_forces[sys->Inter_Types[i].i_0 + k] ) / sys->Inter_Types[i].dr );

		}

		for (l = 0; l < sys->Inter_Types[j].N_pts; l++) {
		    index =
			index_Lpacked((sys->Inter_Types[i].i_0 + k),
				      (sys->Inter_Types[j].i_0 + l),
				      sys->N_coeff);

		    if ((i == j) && (k == l)) {
			b_soln_struct +=
			    (M[index] -
			     g[sys->Inter_Types[i].i_0 +
			       k]) * phi_struct[sys->Inter_Types[j].i_0 +
						l];
			//b_soln_forces += ( M[index] - g[sys->Inter_Types[i].i_0 + k] ) * phi_forces[sys->Inter_Types[j].i_0 + l];
		    } else {
			b_soln_struct +=
			    M[index] * phi_struct[sys->Inter_Types[j].i_0 +
						  l];
			//b_soln_forces += M[index] * phi_forces[sys->Inter_Types[j].i_0 + l];

		    }

		}

		fprintf(fp_soln_struct, "%10.5lf    %15.5lf \n",
			sys->Inter_Types[i].ptr_x[k],
			b_soln_struct / sys->Inter_Types[i].dr);

		//fprintf( fp_soln_forces, "%10.5lf    %15.5lf \n", sys->Inter_Types[i].ptr_x[k], b_soln_forces / sys->Inter_Types[i].dr );

		b_soln_struct = 0.00;
		//b_soln_forces = 0.00;  
	    }
	    fclose(fp_soln_struct);
	    //fclose( fp_soln_forces );

	}
	fclose(fp_soln_struct_diag);
	//fclose( fp_soln_forces_diag );

    }

    return 0;
}

/*****************************************************************************************
PT_eigen(): 
*****************************************************************************************/
int PT_eigen(int k, int N, int N_eigen, double **A)
{
    int i, j;
    int N_row = N;		// matrix is square
    int N_col = N;		// matrix is square     
    int LDA = N;		// leading dimension of A as stored in AT.  NB A is square
    int info;
    double MT[N * N];
    double *MTp = MT;
    double Tau[N];
    double *Taup = Tau;
    double E[N];
    double *Ep = E;
    double D[N];
    double *Dp = D;
    int lwork = N * N;
    double work[lwork];
    char uplo = 'U';
    char compz = 'V';

    tW_word fname;
    FILE *fp_U, *fp_S;

    sprintf(fname, "eigenvalues.%dorder.dat", k);
    fp_S = fopen(fname, "w");	//M col of MxM orthogonal matrix

    //MT = ( double * ) ecalloc( N_row*N_col, sizeof( double ) ); 
    //if ( MT == NULL ) 
    //{ printf( "STOP. MT ecalloc failed. N_row: %d \t N_col: %d\n", N_row, N_col );  exit(0); } 

    //Tau  = ( double * ) ecalloc( N_row-1, sizeof( double ) );
    //E    = ( double * ) ecalloc( N_row-1, sizeof( double ) );
    //D    = ( double * ) ecalloc( N_row-1, sizeof( double ) ); 

    for (j = 0; j < N_col; j++) {
	for (i = 0; i < N_row; i++) {
	    MTp[i + N_row * j] = A[i][j];
	}
    }

    dsytrd_(&uplo, &N_col, MTp, &LDA, Dp, Ep, Taup, work, &lwork, &info);

    dorgtr_(&uplo, &N_col, MTp, &LDA, Taup, work, &lwork, &info);

    dsteqr_(&compz, &N_col, Dp, Ep, MTp, &LDA, work, &info);

    /* Print out eigenspectrum */

    for (i = 0; i < N_row; i++) {
	fprintf(fp_S, " %lf \n", Dp[i]);

	if ((i < N_eigen) || (i >= (N - N_eigen))) {
	    sprintf(fname, "eigenvector.%d.%dorder.dat", i, k);
	    fp_U = fopen(fname, "w");	//M col of MxM orthogonal matrix

	    for (j = 0; j < N_col; j++) {
		fprintf(fp_U, " %lf \n", MTp[i * N + j]);
	    }
	    fclose(fp_U);
	}

    }
    fclose(fp_S);

    //free( MT ); /* no longer allocate memory for this */
    //free( D );  /* no longer allocate memory for this */
    //free( E ); /* E has been destroyed by dsteqr? */
    //free( Tau );  /* Tau has been destroyed by dorgtr? */

    return info;
}

/*****************************************************************************************
eigen(): 
*****************************************************************************************/
int eigen(FILE * fp_log, tW_system * sys, int flag_Gbar)
{
    int i, j;
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int LDA = N;		// leading dimension of A as stored in AT.  NB A is square
    int info;
    double *M = sys->M;
    double MT[N * N];
    double *MTp = MT;
    double Tau[N];
    double *Taup = Tau;
    double E[N];
    double *Ep = E;
    double D[N];
    double *Dp = D;

    double Calpha[N];
    double *phi = sys->phi;

    int N_M = N;
    int N_L = 0;

    int lwork = N * N;
    double work[lwork];
    char uplo = 'U';
    char compz = 'V';

    int N_eigen = sys->Eigen_var.N_Eigen;
    tW_word matrix;

    int flag_norm = sys->Eigen_var.flag_norm;
    int i_flag, j_flag;
    double ri, rj;
    double normalized;

    int index;

    FILE *fp_MT;

    int flag = 0;

    if (flag_Gbar == TRUE) {
	strcpy(matrix, "Gbar");
    } else {
	strcpy(matrix, "G");

	fp_MT = fopen("MT.total.dat", "w");
    }

    tW_word fname;
    FILE *fp_U, *fp_S;

    sprintf(fname, "eigenvalues_%s.dat", matrix);
    fp_S = fopen(fname, "w");	//the singular values of M

    fprintf(fp_log, "In SVD.\n");

    if (flag_Gbar == TRUE) {
	for (j = 0; j < N; j++) {
	    j_flag = strcmp(NB_PAIR, sys->Inter_Types[j].inter_type);
	    rj = sys->Inter_Types[j].R_min + j * sys->Inter_Types[j].dr;
	    for (i = 0; i < N; i++) {
		i_flag = strcmp(NB_PAIR, sys->Inter_Types[i].inter_type);
		ri = sys->Inter_Types[i].R_min +
		    i * sys->Inter_Types[i].dr;
		if (i_flag == 0) {
		    if (j_flag == 0) {
			normalized = ri * ri * rj * rj;
		    } else {
			normalized = ri * ri;
		    }
		} else if (j_flag == 0) {
		    normalized = rj * rj;
		} else {
		    normalized = 1.0;
		}
		if (flag_norm == FALSE) {
		    normalized = 1.0;
		}

		index = index_Lpacked(i, j, sys->N_coeff);

		if (i == j) {
		    MTp[i + N * j] = (M[index] - sys->g[i]) / normalized;	/* This only works for the delta basis!!! */
		} else {
		    MTp[i + N * j] = M[index] / normalized;
		}
	    }
	}
    } else {
	for (j = 0; j < N; j++) {
	    j_flag = strcmp(NB_PAIR, sys->Inter_Types[j].inter_type);
	    rj = sys->Inter_Types[j].R_min + j * sys->Inter_Types[j].dr;
	    for (i = 0; i < N; i++) {
		i_flag = strcmp(NB_PAIR, sys->Inter_Types[i].inter_type);
		ri = sys->Inter_Types[i].R_min +
		    i * sys->Inter_Types[i].dr;
		if (i_flag == 0) {
		    if (j_flag == 0) {
			normalized = ri * ri * rj * rj;
		    } else {
			normalized = ri * ri;
		    }
		} else if (j_flag == 0) {
		    normalized = rj * rj;
		} else {
		    normalized = 1.0;
		}
		if (flag_norm == FALSE) {
		    normalized = 1.0;
		}

		index = index_Lpacked(i, j, sys->N_coeff);

		MTp[i + N * j] = M[index] / normalized;

		fprintf(fp_MT, "%.6e \n", MT[i + N * j] / normalized);
	    }
	}
    }

    dsytrd_(&uplo, &N, MTp, &LDA, Dp, Ep, Taup, work, &lwork, &info);

    dorgtr_(&uplo, &N, MTp, &LDA, Taup, work, &lwork, &info);

    dsteqr_(&compz, &N, Dp, Ep, MTp, &LDA, work, &info);

    if (info != 0) { 
	fprintf(fp_S, "WARNING: info = %d, something wrong in eigen calc! \n", info);
	fprintf(stderr, "WARNING: info = %d, something wrong in eigen calc! \n", info);
	//exit(0);
    } else {
	/* print out the SVD */

	for (i = 0; i < N; i++) {

	    fprintf(fp_S, " %lf \n", D[i]);

	    if ((i < N_eigen) || (i >= (N - N_eigen))) {
		sprintf(fname, "eigenvector_%s.%d.dat", matrix, i);
		fp_U = fopen(fname, "w");	//M col of MxM orthogonal matrix

		for (j = 0; j < N; j++) {
		    //M[i][j] = row_max[i] * M[i][j];  

		    fprintf(fp_U, " %lf \n", MT[i * N + j]);

		    //if( j == N - 1 )
		    //{
		    //fprintf( fp_U, " \n" );
		    //fprintf( fp_V, " \n" );
		    //}

		}
		fclose(fp_U);
	    }

	}

	if ((sys->Eigen_var.flag_printn == TRUE) && (flag_Gbar == FALSE)) {

	    calc_IPR(N, MT);
	    calc_Calpha(N, phi, MT, D, Calpha);

	    do {

		do {

		    if (N_M > N_L) {
			calc_fn_bn(N_M, N_L, Calpha, MT, D, sys);
			print_Mn(N_M, N_L, D, MT, sys);
		    }

		    N_L = N_L + sys->Eigen_var.DL;

		} while (N_L < N);

		if (flag == 0) {
		    N_M = N_M - (N % sys->Eigen_var.DM);
		    flag = 1;
		} else {
		    N_M = N_M - sys->Eigen_var.DM;
		}

		N_L = 0;

	    } while (N_M > 0);
	}
//    free( MT );
//    free( D ); 
//    free( E ); /* E has been destroyed by dsteqr */
//    free( Tau );  /* Tau has been destroyed by dorgtr? */
    }
    fclose(fp_MT);
    fclose(fp_S);

    return info;
}

/*****************************************************************************************
print_Mn():
*****************************************************************************************/
int print_Mn(int N_M, int L, double D[], double MT[], tW_system * sys)
{
    int i, j, k, ki, kj;
    int row_i, col_j;
    int N = sys->N_coeff;
    int i_flag, j_flag;
    double Mn[N][N];
    double ri, rj;
    tW_Inter_Types *inter_i, *inter_j;
    tW_word fname;
    FILE *fp_M, *fp_M2, *fp_T, *fp_T2;
    int diag_flag = FALSE;

    double normalized;

    //initialize Mn
    for (i = 0; i < N; i++) {

	for (j = 0; j < N; j++) {
	    Mn[i][j] = 0.00;
	}

    }

    //calculate Mn from N_M eigenvectors and eigenvalues
    for (i = 0; i < N; i++) {

	for (j = 0; j < N; j++) {

	    for (k = L; k < N_M; k++) {
		Mn[i][j] += MT[k * N + i] * D[k] * MT[k * N + j];
	    }
	}

    }

    sprintf(fname, "T.%dto%devecs.dat", L, N_M);
    fp_T = fopen(fname, "w");	//Whole matrix without diagonals

    sprintf(fname, "T2.%dto%devecs.dat", L, N_M);
    fp_T2 = fopen(fname, "w");	//whole matrix with diagonal elements

    for (i = 0; i < sys->N_Inter_Types; i++) {
	inter_i = &(sys->Inter_Types[i]);

	i_flag = strcmp(NB_PAIR, inter_i->inter_type);

	for (j = 0; j < sys->N_Inter_Types; j++) {
	    inter_j = &(sys->Inter_Types[j]);

	    j_flag = strcmp(NB_PAIR, inter_j->inter_type);

	    if ((inter_i->i_basis == DELTA_BASIS_INDEX)
		&& (inter_j->i_basis == DELTA_BASIS_INDEX)) {

		sprintf(fname, "M.%dto%devecs.%s.%s.dat", L, N_M,
			inter_i->inter_name, inter_j->inter_name);
		fp_M = fopen(fname, "w");

		sprintf(fname, "M2.%dto%devecs.%s.%s.dat", L, N_M,
			inter_i->inter_name, inter_j->inter_name);
		fp_M2 = fopen(fname, "w");	//M2 will be the matrix with diagonal elements

//        ri = inter_i->R_min;
		ri = inter_i->ptr_x[0];

		for (ki = 0; ki < inter_i->N_coeff; ki++) {
		    row_i = inter_i->i_0 + ki;

//          rj = inter_j->R_min;
		    rj = inter_j->ptr_x[0];
		    for (kj = 0; kj < inter_j->N_coeff; kj++) {
			col_j = inter_j->i_0 + kj;

			if (i_flag == 0) {
			    if (j_flag == 0) {
				normalized = ri * ri * rj * rj;
			    } else {
				normalized = ri * ri;
			    }
			} else if (j_flag == 0) {
			    normalized = rj * rj;
			} else {
			    normalized = 1.0;
			}

//            if( fabs( sys->M[row_i][col_j] ) > FLOAT_EPS )
			if ((ri > FLOAT_EPS) && (rj > FLOAT_EPS)) {

			    if (row_i == col_j) {
				fprintf(fp_M2, "%.6e  %.6e  %.6e\n", ri,
					rj, Mn[row_i][col_j] / normalized);
				fprintf(fp_T2, "%d  %d  %.6e\n", row_i,
					col_j,
					Mn[row_i][col_j] / normalized);
			    } else {
				fprintf(fp_M, "%.6e  %.6e  %.6e\n", ri, rj,
					Mn[row_i][col_j] / normalized);
				fprintf(fp_T, "%d  %d  %.6e\n", row_i,
					col_j,
					Mn[row_i][col_j] / normalized);
				fprintf(fp_M2, "%.6e  %.6e  %.6e\n", ri,
					rj, Mn[row_i][col_j] / normalized);
				fprintf(fp_T2, "%d  %d  %.6e\n", row_i,
					col_j,
					Mn[row_i][col_j] / normalized);
			    }
			}

			rj += inter_j->dr;
		    }
		    ri += inter_i->dr;
		}

		fclose(fp_M);
		fclose(fp_M2);
	    } else if (inter_i->i_basis == LINEAR_BASIS_INDEX
		       && inter_j->i_basis == LINEAR_BASIS_INDEX) {
		sprintf(fname, "M.%dto%devecs.%s.%s.dat", L, N_M,
			inter_i->inter_name, inter_j->inter_name);
		fp_M = fopen(fname, "w");

		sprintf(fname, "M2.%dto%devecs.%s.%s.dat", L, N_M,
			inter_i->inter_name, inter_j->inter_name);
		fp_M2 = fopen(fname, "w");	//M2 will be the matrix with diagonal elements

//        ri = inter_i->R_min;
		ri = inter_i->ptr_x[0];
		for (ki = 0; ki < inter_i->N_coeff; ki++) {
		    row_i = inter_i->i_0 + ki;

//          rj = inter_j->R_min;
		    rj = inter_j->ptr_x[0];
		    for (kj = 0; kj < inter_j->N_coeff; kj++) {
			col_j = inter_j->i_0 + kj;

			if (i_flag == 0) {
			    if (j_flag == 0) {
				normalized = ri * ri * rj * rj;
			    } else {
				normalized = ri * ri;
			    }
			} else if (j_flag == 0) {
			    normalized = rj * rj;
			} else {
			    normalized = 1.0;
			}

			if ((row_i == col_j) || (row_i == col_j + 1)
			    || (row_i == col_j - 1)) {
			    diag_flag = TRUE;
			}
//            if( fabs( sys->M[row_i][col_j] ) > FLOAT_EPS )
			if ((ri > FLOAT_EPS) && (rj > FLOAT_EPS)) {

			    if (diag_flag == TRUE) {
				fprintf(fp_M2, "%.6e  %.6e  %.6e\n", ri,
					rj, Mn[row_i][col_j] / normalized);
				fprintf(fp_T2, "%d  %d  %.6e\n", row_i,
					col_j,
					Mn[row_i][col_j] / normalized);
			    } else {
				fprintf(fp_M, "%.6e  %.6e  %.6e\n", ri, rj,
					Mn[row_i][col_j] / normalized);
				fprintf(fp_T, "%d  %d  %.6e\n", row_i,
					col_j,
					Mn[row_i][col_j] / normalized);
				fprintf(fp_M2, "%.6e  %.6e  %.6e\n", ri,
					rj, Mn[row_i][col_j] / normalized);
				fprintf(fp_T2, "%d  %d  %.6e\n", row_i,
					col_j,
					Mn[row_i][col_j] / normalized);
			    }
			}

			diag_flag = FALSE;

			rj += inter_j->dr;
		    }
		    ri += inter_i->dr;
		}

		fclose(fp_M);
		fclose(fp_M2);
	    }

	}
    }

    fclose(fp_T);
    fclose(fp_T2);

    return 0;
}


/*******************************************************************************************
calc_Calpha():
********************************************************************************************/
int calc_Calpha(int N, double f[], double MT[], double D[],
		double Calpha[])
{
    int i, j;
    double dummy = 0.00;
    FILE *fp, *fp2;

    fp = fopen("Calpha.dat", "w");
    fp2 = fopen("Calpha_lambda.dat", "w");

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    dummy += MT[i * N + j] * f[j];
	}

	Calpha[i] = dummy;
	fprintf(fp, "%15.5lf \n", dummy);
	fprintf(fp2, "%15.5lf \n", dummy * D[i]);

	dummy = 0;

    }

    fclose(fp);
    fclose(fp2);

    return 0;

}

/*******************************************************************************************
********************************************************************************************/

int calc_fn_bn(int N_M, int N_L, double Calpha[], double MT[], double D[],
	       tW_system * sys)
{

    int N = sys->N_coeff;
    int i, j;
    FILE *fp_bn, *fp_fn;
    tW_word fname;
    double fn[N];
    double bn[N];

    for (i = 0; i < N; i++) {
	fn[i] = 0.00;
	bn[i] = 0.00;
    }

    for (i = 0; i < N; i++) {
	for (j = N_L; j < N_M; j++) {
	    fn[i] += MT[j * N + i] * Calpha[j];
	    bn[i] += MT[j * N + i] * Calpha[j] * D[j];
	}
    }

    //JFR - Separate the interactions for output
    for (i = 0; i < sys->N_Inter_Types; i++) {

	sprintf(fname, "%s.%dto%devecs.%s.%s.%s.dat", "b", N_L, N_M,
		sys->Inter_Types[i].inter_type, "total",
		sys->Inter_Types[i].inter_name);
	fp_bn = fopen(fname, "w");

	sprintf(fname, "%s.%dto%devecs.%s.%s.%s.dat", "f", N_L, N_M,
		sys->Inter_Types[i].inter_type, "total",
		sys->Inter_Types[i].inter_name);
	fp_fn = fopen(fname, "w");

	for (j = 0; j < sys->Inter_Types[i].N_pts; j++) {
	    fprintf(fp_fn, "%10.5lf    %15.5lf \n",
		    sys->Inter_Types[i].ptr_x[j], fn[j]);
	    fprintf(fp_bn, "%10.5lf    %15.5lf \n",
		    sys->Inter_Types[i].ptr_x[j], bn[j]);
	}

	fclose(fp_fn);
	fclose(fp_bn);
    }

    return 0;
}

/*******************************************************************************************
********************************************************************************************/

int calc_IPR(int N, double MT[])
{

    int i, j;
    FILE *fp_IPR;
    double tmp = 0.00;

    fp_IPR = fopen("IPR.dat", "w");

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    tmp += pow(MT[N * i + j], 4);
	}

	fprintf(fp_IPR, "%15.5lf \n", tmp);
	tmp = 0.00;
    }

    fclose(fp_IPR);

    return 0;

}

/*******************************************************************************************
********************************************************************************************/

int Precondition(int N, double *M, double *b, double *phi,
		 double *norm_col, double *row_max, int flag_restore,
		 int flag_solve, tW_system * sys)
{
    tW_word NO, dimless, colnorm, rowmax, MTvar, bvar;
    strcpy(NO, "NO");
    strcpy(dimless, "dimless");
    strcpy(colnorm, "colnorm");
    strcpy(rowmax, "rowmax");
    strcpy(MTvar, "MTvar");
    strcpy(bvar, "bvar");

    int i, j;
    int index;
    double col_norm, row_norm, b_norm, M_norm, M_norm_avg, col_norm_avg;

    /* Right Preconditioning */
    if (flag_restore == FALSE) {
	/* initialize the preconditioners */
	for (j = 0; j < N; j++) {
	    norm_col[j] = 0.0;
	    row_max[j] = 0.0;
	}

	if (strcmp(NO, sys->PC_var.RPC) == 0) {
	    for (j = 0; j < N; j++) {
		norm_col[j] = 1.0;
	    }
	} else if (strcmp(dimless, sys->PC_var.RPC) == 0) {
	    for (j = 0; j < N; j++) {
		norm_col[j] = sqrt(M[j + j * N]);
		if (fabs(norm_col[j] < FLOAT_EPS)) {
		    printf ("ERROR: trying to make the equations dimensionless in Dimless_PC(), but one of the diagonal elements of the metric tensor is 0 \n");
		    exit(0);
		}
	    }
	} else if (strcmp(colnorm, sys->PC_var.RPC) == 0) {
	    for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++) {
		    norm_col[j] += M[i + j * N] * M[i + j * N];
		}
		norm_col[j] = sqrt(norm_col[j]);
	    }
	} else if (strcmp(MTvar, sys->PC_var.RPC) == 0) {
	    for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++) {
		    index = index_Lpacked(i, j, N);
		    if (sys->M_cnt[index] > 1.0) {
			norm_col[j] +=
			    sqrt(sys->d2M[index]) /
			    sqrt(((sys->M_cnt[index] - 1.0)));
		    }
//printf("norm_col[%d] = %.15lf \n", j, norm_col[j]);
		}
		printf("norm_col[%d] = %.15lf \n", j, norm_col[j]);
	    }
	} else {
	    printf
		("ERROR: LPC_type not compatible in function Precondition() \n");
	    exit(0);
	}

	if (sys->PC_var.flag_normphi == TRUE) {	/* rescale so the avg norm of the columns of MT remain unchanged */
	    col_norm_avg = 0.00;
	    M_norm_avg = 0.00;
	    for (j = 0; j < N; j++) {
		col_norm = 0.00;
		M_norm = 0.00;
		for (i = 0; i < N; i++) {
		    M_norm += M[i + j * N] * M[i + j * N];
		    col_norm +=
			(M[i + j * N] / norm_col[j]) * (M[i + j * N] /
							norm_col[j]);
		}
		col_norm_avg += (1.0 / ((double) N)) * sqrt(col_norm);
		M_norm_avg += (1.0 / ((double) N)) * sqrt(M_norm);
	    }
	    for (j = 0; j < N; j++) {
		norm_col[j] /= (col_norm_avg / M_norm_avg);
	    }

	}

    }

    Right_PC(N, M, b, phi, norm_col, row_max, flag_restore, flag_solve);

    /* Left Preconditioning */
    if (flag_restore == FALSE) {
	if (strcmp(NO, sys->PC_var.LPC) == 0) {
	    for (i = 0; i < N; i++) {
		row_max[i] = 1.0;
	    }
	} else if (strcmp(dimless, sys->PC_var.LPC) == 0) {
	    for (i = 0; i < N; i++) {
		row_max[i] = sqrt(M[i + i * N]);
		if (fabs(row_max[i] < FLOAT_EPS)) {
		    printf
			("ERROR: trying to make the equations dimensionless in Dimless_PC(), but one of the diagonal elements of the metric tensor is 0 \n");
		    exit(0);
		}
	    }
	} else if (strcmp(rowmax, sys->PC_var.LPC) == 0) {
	    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
		    if (fabs(M[i + j * N]) > row_max[i]) {
			row_max[i] = fabs(M[i + j * N]);
		    }
		}
	    }
	} else if (strcmp(bvar, sys->PC_var.LPC) == 0) {
	    for (i = 0; i < N; i++) {
		if ((sys->rescale[i] * ((double) sys->REG_var.Nframes)) >
		    1.0) {
		    row_max[i] =
			sqrt(sys->d2b[i]) /
			sqrt((sys->rescale[i] *
			      ((double) sys->REG_var.Nframes)) - (1.0));
		    printf("row_max[%d] = %.15lf \n", i, row_max[i]);
		}
	    }
	} else {
	    printf
		("ERROR: LPC_type not compatible in function Precondition() \n");
	    exit(0);
	}

	if (sys->PC_var.flag_normb == TRUE) {	/* rescale so the norm of PCb is same as b */
	    row_norm = 0.00;
	    b_norm = 0.00;
	    for (i = 0; i < N; i++) {
		row_norm += (b[i] / row_max[i]) * (b[i] / row_max[i]);
		b_norm += b[i] * b[i];
	    }
	    for (i = 0; i < N; i++) {
		row_max[i] /= sqrt(row_norm / b_norm);
	    }
	}

    }

    Left_PC(N, M, b, phi, norm_col, row_max, flag_restore, flag_solve);

    return 0;
}

/*******************************************************************************************
********************************************************************************************/

int Right_PC(int N, double *M, double *b, double *phi, double *norm_col,
	     double *row_max, int flag_restore, int flag_solve)
{
    int i, j;

    if (flag_restore == FALSE) {
	for (j = 0; j < N; j++) {	/* loop over columns */
	    for (i = 0; i < N; i++) {
		M[i + j * N] = M[i + j * N] / norm_col[j];
	    }
	}
    } else {
	/* undo the right preconditioning */
	for (j = 0; j < N; j++) {
	    if (flag_solve == TRUE) {
		phi[j] = phi[j] / norm_col[j];
	    }
	    for (i = 0; i < N; i++) {
		M[i + j * N] = M[i + j * N] * norm_col[j];
	    }
	}
    }

    return 0;
}

/*******************************************************************************************
********************************************************************************************/

int Left_PC(int N, double *M, double *b, double *phi, double *norm_col,
	    double *row_max, int flag_restore, int flag_solve)
{

    int i, j;

    if (flag_restore == FALSE) {
	for (i = 0; i < N; i++) {	/* loop over rows */
	    if (flag_solve == TRUE) {
		b[i] = b[i] / row_max[i];
	    }

	    /* Divide through by max. element. */
	    for (j = 0; j < N; j++) {
		M[i + j * N] = M[i + j * N] / row_max[i];
	    }
	}
    } else {
	/* undo left preconditioning */
	for (i = 0; i < N; i++) {	/* loop over rows */
	    if (flag_solve == TRUE) {
		b[i] = b[i] * row_max[i];
	    }
	    for (j = 0; j < N; j++) {
		M[i + j * N] = M[i + j * N] * row_max[i];
	    }
	}
    }

    return 0;
}

/*****************************************************************************************
calc_Chi2(): Calc Chi2.  
JFR - added 01.29.13 
*****************************************************************************************/
double calc_Chi2(FILE * fp_log, tW_system * sys, double phi[], double b[])
{
    int i, j;
    double Gphi[sys->N_coeff];
    double Chi2 = sys->Chi2;
    int index;

    fprintf(fp_log, "ff = %lf \n", Chi2);

    for (i = 0; i < sys->N_coeff; i++) {
	Gphi[i] = 0.00;
    }

    for (i = 0; i < sys->N_coeff; i++) {
	Chi2 -= 2.0 * b[i] * phi[i];

	for (j = 0; j < sys->N_coeff; j++) {
	    index = index_Lpacked(i, j, sys->N_coeff);
	    Gphi[i] += sys->M[index] * phi[j];
	}
    }

    fprintf(fp_log, "ff-2bphi = %lf \n", Chi2);

    for (i = 0; i < sys->N_coeff; i++) {
	Chi2 += phi[i] * Gphi[i];
    }

    fprintf(fp_log, "ff-2bphi+phiGphi = %lf \n", Chi2);

    return Chi2;

}

/*****************************************************************************************
JFR - 12.20.13:  Below this point are functions that I am not currently using, but I 
didn't want to commit to deleting them yet.
*****************************************************************************************/

/*****************************************************************************************
solv_lin_eqns_constrain_tabdih(): 
JFR - 12.20.13: I consolidated this function and solv_lin_eqns_constrain_tabdih_Bayesian_params
() into one function which has the options to regularize and/or constrain.
*****************************************************************************************/
int solv_lin_eqns_constrain_tabdih(FILE * fp_log, tW_system * sys)
{
    int i, j;
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_row = N;
    int N_col = N;
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    int LDB = N;
    int pivot[N];		// pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT;			// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double row_max[N_row], norm_col[N_col];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    /* variables for constraining dihedrals */
    double lambda = 0.00;
    double dlambda = 0.00;
    double tau_dih = 0.00;
    int Nmax = 1;
    if (sys->CONSTRAIN_var.flag_CONSTRAIN == TRUE) {
	lambda = sys->CONSTRAIN_var.lambda;
	dlambda = sys->CONSTRAIN_var.dlambda;
	tau_dih = sys->CONSTRAIN_var.tau_dih;
	Nmax = sys->CONSTRAIN_var.Nmax;
    }
    int k, kk, ii;
    double err = 0.0;
    int ctr;
    int i_flag, row_i, row_k, flag_conv;
    /* variables for regularization */
    double *C;
    double beta = 0.00;
    if ((sys->REG_var.flag_REG == TRUE)) {
	beta = sys->REG_var.tau_beta;
    }

    fprintf(fp_log, "In solv_lin_eqns.\n");

    /* JFR - 02.06.13: variable for solving the system with SVD */
    int N2 = N;
    int N3 = N;
    int rank;
    double *S;
    int lwork = 10 * N;
    double work[lwork];
    double rcond = sys->SVD_var.rcond;
    FILE *fp_S;

    fp_S = fopen("S_CONST.dat", "w");

    S = (double *) ecalloc(N, sizeof(double));
    if (S == NULL) {
	printf("STOP. S ecalloc failed. N: %d\n", N);
	exit(0);
    }

    MT = (double *) ecalloc(N_row * N_col, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }
    C = (double *) ecalloc(N_row * N_col, sizeof(double));
    if (C == NULL) {
	printf("STOP. C ecalloc failed. N: %d\n", N);
	exit(0);
    }

    ctr = 0;
    do {

	/* regularize using the uncertainties of M and b */
	for (j = 0; j < N; j++) {
	    for (i = 0; i < N; i++) {
		index = index_Lpacked(i, j, N);
		if (sys->M_cnt[index] > 1.0) {
		    C[j + N * j] +=
			sqrt(sys->d2M[index]) /
			sqrt(((sys->M_cnt[index] - 1.0)));
		}
	    }
	    C[j + N * j] *= beta;
	}

	/* only regularize nb and intra-nb interactions  */
	for (i = 0; i < sys->N_Inter_Types; i++) {
	    i_flag = strcmp(B_DIHEDRAL, sys->Inter_Types[i].inter_type);
	    for (ii = 0; ii < sys->Inter_Types[i].N_pts; ii++) {
		row_i = sys->Inter_Types[i].i_0 + ii;
		if ((strcmp(B_NB_PAIR_BOND, sys->Inter_Types[i].inter_type)
		     != 0)
		    && (strcmp(NB_PAIR, sys->Inter_Types[i].inter_type) !=
			0)) {
		    C[row_i + N * (row_i)] = 0.00;
		}
	    }
	}

	/* constrain tabulated dihedrals */
	for (i = 0; i < sys->N_Inter_Types; i++) {
	    if (strcmp(B_DIHEDRAL, sys->Inter_Types[i].inter_type) != 0) {
		continue;
	    }			/* check that i is a dihedral interaction */
	    if (strcmp(TOY_DIHED_NAME, sys->Inter_Types[i].basis) == 0) {
		continue;
	    }			/* check that i is not using a TOY basis */
	    if (strcmp(RYCKAERT_BELLEMANS_NAME, sys->Inter_Types[i].basis)) {
		continue;
	    }			/* check that i is not using a RB basis */
	    for (k = 0; k < sys->N_Inter_Types; k++) {
		if (strcmp(B_DIHEDRAL, sys->Inter_Types[k].inter_type) !=
		    0) {
		    continue;
		}		/* check that k is a dihedral interaction */
		if (strcmp(TOY_DIHED_NAME, sys->Inter_Types[k].basis) == 0) {
		    continue;
		}		/* check that k is not using a TOY basis */
		if (strcmp
		    (RYCKAERT_BELLEMANS_NAME, sys->Inter_Types[k].basis)) {
		    continue;
		}		/* check that k is not using a RB basis */
		if (strcmp
		    (sys->Inter_Types[i].inter_name,
		     sys->Inter_Types[k].inter_name) != 0) {
		    continue;
		}		/*  check that i and k are the same dihedral */
		for (ii = 0; ii < sys->Inter_Types[i].N_pts; ii++) {
		    row_i = sys->Inter_Types[i].i_0 + ii;
		    for (kk = 0; kk < sys->Inter_Types[k].N_pts; kk++) {
			row_k = sys->Inter_Types[k].i_0 + kk;
			C[row_i + N * (row_k)] += lambda;
		    }
		}
	    }
	}

	for (i = 0; i < N; i++) {
	    for (j = 0; j < N_col; j++) {
		index = index_Lpacked(i, j, N);
		MT[i + N * j] = M[index] + C[i + N * j];
	    }
	}

	/*  Precondition the matrix */
	//Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

	for (i = 0; i < N; i++) {
	    /* Copy b -> x which will be overwritten with soln. */
	    phi[i] = b[i];
	}
	if (strcmp(sys->SOLN_var.SOLN_METH, "LU") == 0) {
	    /*  Solve the system using LU Decomposition  */
	    fprintf(fp_log,
		    " Inverting the matrix using LU decomposition \n");
	    time(&start);
	    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
	    time(&end);
	    diff = difftime(end, start);
	    if (info != 0) {
		fprintf(fp_log,
			" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
			info);
	    }
	    fprintf(fp_log,
		    " dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
		    diff, diff / (60.0), diff / (3600.00));
	} else {
	    /* Solve the overdetermined system using Singular Value Decomposition, setting all singular values under rcond to zero  */
	    fprintf(fp_log,
		    " Inverting the matrix using Singular Value Decomposition in contstrain tabdih \n");
	    time(&start);
	    dgelss_(&N, &N2, &N_rhs, MT, &LDA, phi, &N3, S, &rcond, &rank,
		    work, &lwork, &info);
	    time(&end);
	    diff = difftime(end, start);
	    if (info != 0) {
		fprintf(fp_log,
			" info = %d after dgelss_().  The inverse condition number may be less than machine precision \n",
			info);
	    }
	    fprintf(fp_log,
		    " dgelss_() took %lf seconds or %lf minutes or %lf hours \n",
		    diff, diff / (60.0), diff / (3600.00));
	    fprintf(fp_log, " Rank(G) = %d \n", rank);
	    fprintf(fp_log, "\n");

	    /* print out the singular values */
	    for (i = 0; i < N; i++) {
		fprintf(fp_S, " %.16lf \n", S[i]);
	    }
	} 

	/* Undo the Preconditioning */
	//Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

	/* check that the matrix was not rank deficient */
	if (info != 0) {
	    ctr++;
	    lambda = lambda * dlambda;
	    flag_conv = FALSE;
	    continue;
	}

	/* check that the solution is not identically zero */
//    err = 0.00;
//    for ( i=N; i<N_row; i++ )
//    {
//      err += phi[i]*phi[i];
//    }
//    if ( fabs(err) < FLOAT_EPS ) { ctr++; lambda = lambda*dlambda; flag_conv = FALSE; fprintf( fp_log, " The solution is exactly zero, continue optimization of lambda parameter" ); continue; }

	/* check the error in the constraint */
	flag_conv = TRUE;
	for (i = 0; i < N; i++) {
	    err = 0.00;
	    for (j = 0; j < N_col; j++) {
		err += C[i + N * j] * phi[j];
	    }
	    if (fabs(err) > tau_dih) {
		flag_conv = FALSE;
	    }
	}

	/* increase lambda by dlambda */
	lambda = lambda * dlambda;

	ctr++;

    } while ((flag_conv == FALSE) && (ctr < Nmax));

    /* print out the final lambda scaling terms */
    fprintf(fp_log, " lambda = %lf \n", lambda / dlambda);
    if (ctr < Nmax) {
	fprintf(fp_log,
		" Optimization of lambda parameter converged after %d steps",
		ctr);
    } else {
	fprintf(fp_log,
		" Optimization of lambda parameter did not converge after Nmax = %d steps",
		Nmax);
    }

    fclose(fp_S);

    free(MT);
    free(S);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_constrain_tabdih_Bayesian_params(): 
JFR - 12.20.13: I consolidated this function and solv_lin_eqns_constrain_tabdih() into one
function which has the options to regularize and/or constrain.
*****************************************************************************************/
int solv_lin_eqns_constrain_tabdih_Bayesian_params(FILE * fp_log,
						   tW_system * sys,
						   double *alpha,
						   double beta)
{
    int i, j;
    //int Ndih = 0; 
    //for ( i=0; i<sys->N_Inter_Types; i++ )
    //{
    //  if ( strcmp( B_DIHEDRAL, sys->Inter_Types[i].inter_type ) ) { Ndih++; }
    //}
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_row = N;
    int N_col = N;
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    int LDB = N;
    int pivot[N];		// pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT;			// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double row_max[N_row], norm_col[N_col];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    /* variables for constraining dihedrals */
    double lambda = sys->CONSTRAIN_var.lambda;
    double dlambda = sys->CONSTRAIN_var.dlambda;
    double tau_dih = sys->CONSTRAIN_var.tau_dih;
    int Nmax = sys->CONSTRAIN_var.Nmax;
    int k, kk, ii;
    double err = 0.0;
    int ctr;
    int i_flag, k_flag, row_i, row_k, flag_conv;
    double *C;

    fprintf(fp_log, "In solv_lin_eqns.\n");

    /* JFR - 02.06.13: variable for solving the system with SVD */
    int N2 = N;
    int N3 = N;
    int rank;
    double *S;
    int lwork = 10 * N;
    double work[lwork];
    double rcond = sys->SVD_var.rcond;
    FILE *fp_S;

    fp_S = fopen("S_CONST.dat", "w");

    S = (double *) ecalloc(N, sizeof(double));
    if (S == NULL) {
	printf("STOP. S ecalloc failed. N: %d\n", N);
	exit(0);
    }

    MT = (double *) ecalloc(N_row * N_col, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }
    C = (double *) ecalloc(N_row * N_col, sizeof(double));
    if (C == NULL) {
	printf("STOP. C ecalloc failed. N: %d\n", N);
	exit(0);
    }

    ctr = 0;
    do {

	/* regularize using the Bayesian regularization parameters */
	for (j = 0; j < N; j++) {
	    for (i = 0; i < N; i++) {
		if (i == j) {
		    C[j + N * j] = alpha[i] / beta;
		}
	    }
	}

	for (i = 0; i < sys->N_Inter_Types; i++) {
	    i_flag = strcmp(B_DIHEDRAL, sys->Inter_Types[i].inter_type);
	    for (k = 0; k < sys->N_Inter_Types; k++) {
		k_flag =
		    strcmp(B_DIHEDRAL, sys->Inter_Types[k].inter_type);
		for (ii = 0; ii < sys->Inter_Types[i].N_pts; ii++) {
		    row_i = sys->Inter_Types[i].i_0 + ii;
		    for (kk = 0; kk < sys->Inter_Types[k].N_pts; kk++) {
			row_k = sys->Inter_Types[k].i_0 + kk;
			if ((i_flag == 0) && (k_flag == 0)) {
			    C[row_i + N * (row_k)] += lambda;
			} else {
			    C[row_i + N * (row_k)] += 0.0;
			}
		    }
		}
	    }
	}

	for (i = 0; i < N; i++) {
	    for (j = 0; j < N_col; j++) {
		index = index_Lpacked(i, j, N);
		MT[i + N * j] = M[index] + C[i + N * j];
	    }
	}

	/*  Precondition the matrix */
	//Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

	for (i = 0; i < N; i++) {
	    /* Copy b -> x which will be overwritten with soln. */
	    phi[i] = b[i];
	}

	if (strcmp(sys->SOLN_var.SOLN_METH, "LU") == 0) {
	    /*  Solve the system using LU Decomposition  */
	    fprintf(fp_log,
		    " Inverting the matrix using LU decomposition in default tabdih constrain Bayesian \n");
	    time(&start);
	    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
	    time(&end);
	    diff = difftime(end, start);
	    if (info != 0) {
		fprintf(fp_log,
			" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
			info);
	    }
	    fprintf(fp_log,
		    " dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
		    diff, diff / (60.0), diff / (3600.00));
	} else {
	    /* Solve the overdetermined system using Singular Value Decomposition, setting all singular values under rcond to zero  */
	    fprintf(fp_log,
		    " Inverting the matrix using Singular Value Decomposition in tabdih constrain Bayesian \n");
	    time(&start);
	    dgelss_(&N, &N2, &N_rhs, MT, &LDA, phi, &N3, S, &rcond, &rank,
		    work, &lwork, &info);
	    time(&end);
	    diff = difftime(end, start);
	    if (info != 0) {
		fprintf(fp_log,
			" info = %d after dgelss_().  The inverse condition number may be less than machine precision \n",
			info);
	    }
	    fprintf(fp_log,
		    " dgelss_() took %lf seconds or %lf minutes or %lf hours \n",
		    diff, diff / (60.0), diff / (3600.00));
	    fprintf(fp_log, " Rank(G) = %d \n", rank);
	    fprintf(fp_log, "\n");

	    /* print out the singular values */
	    for (i = 0; i < N; i++) {
		fprintf(fp_S, " %.16lf \n", S[i]);
	    }
	} 

	/* Undo the Preconditioning */
	//Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

	/* check that the matrix was not rank deficient */
	if (info != 0) {
	    ctr++;
	    lambda = lambda * dlambda;
	    flag_conv = FALSE;
	    continue;
	}

	/* check that the solution is not identically zero */
//    err = 0.00;
//    for ( i=N; i<N_row; i++ )
//    {
//      err += phi[i]*phi[i];
//    }
//    if ( fabs(err) < FLOAT_EPS ) { ctr++; lambda = lambda*dlambda; flag_conv = FALSE; fprintf( fp_log, " The solution is exactly zero, continue optimization of lambda parameter" ); continue; }

	/* check the error in the constraint */
	flag_conv = TRUE;
	for (i = 0; i < N; i++) {
	    err = 0.00;
	    for (j = 0; j < N_col; j++) {
		err += C[i + N * j] * phi[j];
	    }
	    if (fabs(err) > tau_dih) {
		flag_conv = FALSE;
	    }
	}

	/* increase lambda by dlambda */
	lambda = lambda * dlambda;

	ctr++;

    } while ((flag_conv == FALSE) && (ctr < Nmax));

    /* print out the final lambda scaling terms */
    fprintf(fp_log, " lambda = %lf \n", lambda / dlambda);
    if (ctr < Nmax) {
	fprintf(fp_log,
		" Optimization of lambda parameter converged after %d steps",
		ctr);
    } else {
	fprintf(fp_log,
		" Optimization of lambda parameter did not converge after Nmax = %d steps",
		Nmax);
    }

    fclose(fp_S);

    free(MT);
    free(S);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_Bayesian(): 
JFR - 12.20.13: This is the original implementation of the Bayesian inference regularization.
solv_lin_eqns_Bayesian_params() was created from this function to store the optimized 
parameters for additional calculations.
Now I see why I separated the Bayes + contraint calculation, this way I only have one function
implementating the Bayesian inference method so that future adjustments will be as simple
as possible.
*****************************************************************************************/
int solv_lin_eqns_Bayesian(FILE * fp_log, tW_system * sys)
{
    int i, j;
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N2 = N;
    int N3 = N;
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    int LDB = N;
    int LDC = N;
    int pivot[N];		// pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT, *MT2, *MT3;	// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double row_max[N], norm_col[N];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    int lwork = 10 * N;
    double work[lwork];
    /* variables for regularization */
    double alpha[N * N];
    double alpha_old[N * N];
    double beta = 0;
    double beta_old = 0;
    double tau_alpha = sys->REG_var.tau_alpha;
    double tau_beta = sys->REG_var.tau_beta;
    double f2 = 0;
    int nt = sys->REG_var.Nframes;
    int NT = 3 * nt * N;
    int Nmax = sys->REG_var.Nmax;
    int ctr;
    double trace;
    int flag_conv;
    int relerr;
    char TRANSA = 'N';
    char TRANSB = 'N';
    double aconst = 1.0;
    double bconst = 0.0;

    fprintf(fp_log, "In solv_lin_eqns.\n");

    MT = (double *) ecalloc(N * N, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }
    MT2 = (double *) ecalloc(N * N, sizeof(double));
    if (MT2 == NULL) {
	printf("STOP. MT2 ecalloc failed. N: %d\n", N);
	exit(0);
    }
    MT3 = (double *) ecalloc(N * N, sizeof(double));
    if (MT3 == NULL) {
	printf("STOP. MT3 ecalloc failed. N: %d\n", N);
	exit(0);
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index];
	}
    }

    /* Do the normal calculation */

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	phi[i] = b[i];
    }

    /*  Solve the system using LU Decomposition  */
    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
    time(&start);
    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

    /* Now calculate the input values of alpha and beta */

    for (i = 0; i < N; i++) {
	f2 += phi[i] * phi[i];
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    if (i == j) {
		alpha[i + N * j] = N / f2;
	    } else {
		alpha[i + N * j] = 0.0;
	    }
	}
    }

    beta = NT / calc_Chi2(fp_log, sys, phi, b);

    /* print out the initial regularization parameters */
    fprintf(fp_log, " Initial Regularization Parameters: \n");
    for (i = 0; i < N; i++) {
	fprintf(fp_log, " alpha_0[%d] = %lf \n", i, alpha[i + N * i]);
    }
    fprintf(fp_log, " beta_0 = %lf \n", beta);

    ctr = 0;
    do {

	/* Prepare to calculate mN */
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		index = index_Lpacked(i, j, N);
		MT[i + N * j] = M[index] + alpha[i + N * j] / beta;
	    }
	}

	/*  Precondition the matrix */
	//Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

	for (i = 0; i < N; i++) {
	    /* Copy b -> x which will be overwritten with soln. */
	    phi[i] = b[i];
	}

	/*  Solve the regularized system using LU Decomposition  */
	fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
	time(&start);
	dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log,
		    " info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		    info);
	}
	fprintf(fp_log,
		" dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
		diff, diff / (60.0), diff / (3600.00));

	/* Undo the Preconditioning */
	//Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

	/* Calculate the inverse of the regularized matrix */
	dgetri_(&N, MT, &LDA, pivot, work, &lwork, &info);

	/* Calculate the new alpha and beta quantities */
	for (i = 0; i < N; i++) {
	    alpha_old[i + N * i] = alpha[i + N * i];
	    alpha[i + N * i] =
		(1.0 / (phi[i] * phi[i] + (MT[i + N * i] / beta)));
	}

	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		index = index_Lpacked(i, j, N);
		MT2[i + N * j] = M[index];
		MT3[i + N * j] = 0.0;
	    }
	}
	/* Do matrix multiplication */
	dgemm_(&TRANSA, &TRANSB, &N, &N2, &N3, &aconst, MT, &LDA, MT2,
	       &LDB, &bconst, MT3, &LDC);

	trace = 0.0;
	for (i = 0; i < N; i++) {
	    trace += beta * MT3[i + N * i];
	}

	beta_old = beta;
	beta = (NT - trace) / (calc_Chi2(fp_log, sys, phi, b));

	flag_conv = TRUE;
	for (i = 0; i < N; i++) {
	    relerr =
		fabs((alpha[i + N * i] -
		      alpha_old[i + N * i]) / alpha_old[i + N * i]);
	    if (relerr > tau_alpha) {
		flag_conv = FALSE;
	    }
	}

	relerr = fabs((beta - beta_old) / beta_old);
	if (relerr > tau_beta) {
	    flag_conv = FALSE;
	}

	ctr++;

    } while ((flag_conv == FALSE) && (ctr < Nmax));

    /* print out the final regularization parameters */
    fprintf(fp_log, " Final Regularization Parameters: \n");
    for (i = 0; i < N; i++) {
	fprintf(fp_log, " alpha_f[%d] = %lf \n", i, alpha[i + N * i]);
    }
    fprintf(fp_log, " beta_f = %lf, nt = %d \n", beta, nt);
    if (ctr < Nmax) {
	fprintf(fp_log,
		" Optimization of regularization parameters converged after %d steps",
		ctr);
    } else {
	fprintf(fp_log,
		" Optimization of regularization parameters did not converge after Nmax = %d steps",
		Nmax);
    }

    /* Prepare to calculate the FINAL mN */
    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index] + alpha[i + N * j] / beta;
	}
    }

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	phi[i] = b[i];
    }

    /*  Solve the regularized system using LU Decomposition  */
    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
    time(&start);
    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

    free(MT);
    free(MT2);
    free(MT3);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_constrain_tabdih_Bayesian(): 
JFR - 12.20.13:  This should be equivalent to solv_lin_eqns_constrain_tabdih_Bayesian_params(),
except the optimization of the regularization parameters happens inside the function.
I can't remember why I separated the two calculations (for simplicity? or so I could do
multiple calculations with the same optimized parameters?).
At the moment, I think the separated calculation is conceptually simpler and easier to 
manipulate so this function is obsolete, at least for now.
*****************************************************************************************/
int solv_lin_eqns_constrain_tabdih_Bayesian(FILE * fp_log, tW_system * sys)
{
    int i, j;
    int index;
    int Ndih = 0;
    for (i = 0; i < sys->N_Inter_Types; i++) {
	if (strcmp(B_DIHEDRAL, sys->Inter_Types[i].inter_type) == 0) {
	    Ndih++;
	}
    }
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N2 = N;
    int N3 = N;
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    int LDB = N;
    int LDC = N;
    int LDE = Ndih;
    int pivot[N];		// pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT, *MT2, *MT3;	// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double row_max[N], norm_col[N];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    int lwork = 10 * N;
    double work[lwork];
    /* variables for regularization */
    double alpha[N * N];
    double alpha_old[N * N];
    double beta = 0;
    double beta_old = 0;
    double tau_alpha = sys->REG_var.tau_alpha;
    double tau_beta = sys->REG_var.tau_beta;
    double f2 = 0;
    int nt = sys->REG_var.Nframes;
    int NT = 3 * nt * N;
    int Nmax = sys->REG_var.Nmax;
    int ctr;
    double trace;
    int flag_conv;
    int relerr;
    char TRANSA = 'N';
    char TRANSB = 'N';
    double aconst = 1.0;
    double bconst = 0.0;
    /* variables for constraining dihedrals */
    int k, kk;
    int k_flag, row_k;
    double *C;
    double *d;
    double *binp;

    fprintf(fp_log, "In solv_lin_eqns.\n");

    MT = (double *) ecalloc(N * N, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }
    MT2 = (double *) ecalloc(N * N, sizeof(double));
    if (MT2 == NULL) {
	printf("STOP. MT2 ecalloc failed. N: %d\n", N);
	exit(0);
    }
    MT3 = (double *) ecalloc(N * N, sizeof(double));
    if (MT3 == NULL) {
	printf("STOP. MT3 ecalloc failed. N: %d\n", N);
	exit(0);
    }
    C = (double *) ecalloc(Ndih * N, sizeof(double));
    if (C == NULL) {
	printf("STOP. C ecalloc failed. N: %d\n", N);
	exit(0);
    }
    d = (double *) ecalloc(Ndih, sizeof(double));
    if (d == NULL) {
	printf("STOP. d ecalloc failed. N: %d\n", N);
	exit(0);
    }
    binp = (double *) ecalloc(N, sizeof(double));
    if (binp == NULL) {
	printf("STOP. binp ecalloc failed. N: %d\n", N);
	exit(0);
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index];
	}
    }

    /* set the constraints */
    for (i = 0; i < Ndih; i++) {
	d[i] = 0.00;
	for (k = 0; k < sys->N_Inter_Types; k++) {
	    k_flag = strcmp(B_DIHEDRAL, sys->Inter_Types[k].inter_type);
	    for (kk = 0; kk < sys->Inter_Types[k].N_pts; kk++) {
		row_k = sys->Inter_Types[k].i_0 + kk;
		if (k_flag == 0) {
		    C[i * N + row_k] = 1.0;
		} else {
		    C[i * N + row_k] = 0.0;
		}
	    }
	}
    }

    /* Do the normal calculation */

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	binp[i] = b[i];
    }

    /*  Solve the system with equality constraints  */
    fprintf(fp_log, " Inverting the matrix with equality constraints \n");
    time(&start);
    dgglse_(&N, &N2, &Ndih, MT, &LDA, C, &LDE, binp, d, phi, work, &lwork,
	    &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
    // Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

    /* Now calculate the input values of alpha and beta */

    for (i = 0; i < N; i++) {
	f2 += phi[i] * phi[i];
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    if (i == j) {
		alpha[i + N * j] = N / f2;
	    } else {
		alpha[i + N * j] = 0.0;
	    }
	}
    }

    beta = NT / calc_Chi2(fp_log, sys, phi, b);

    /* print out the initial regularization parameters */
    fprintf(fp_log, " Initial Regularization Parameters: \n");
    for (i = 0; i < N; i++) {
	fprintf(fp_log, " alpha_0[%d] = %lf \n", i, alpha[i + N * i]);
    }
    fprintf(fp_log, " beta_0 = %lf \n", beta);

    ctr = 0;
    do {

	/* Prepare to calculate mN */
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		index = index_Lpacked(i, j, N);
		MT[i + N * j] = M[index] + alpha[i + N * j] / beta;
	    }
	}

	/* d was killed by dgglse */
	d = (double *) ecalloc(Ndih, sizeof(double));
	if (d == NULL) {
	    printf("STOP. d ecalloc failed. N: %d\n", N);
	    exit(0);
	}

	/* set the constraints */
	for (i = 0; i < Ndih; i++) {
	    d[i] = 0.00;
	    for (k = 0; k < sys->N_Inter_Types; k++) {
		k_flag =
		    strcmp(B_DIHEDRAL, sys->Inter_Types[k].inter_type);
		for (kk = 0; kk < sys->Inter_Types[k].N_pts; kk++) {
		    row_k = sys->Inter_Types[k].i_0 + kk;
		    if (k_flag == 0) {
			C[i * N + row_k] = 1.0;
		    } else {
			C[i * N + row_k] = 0.0;
		    }
		}
	    }
	}

	/*  Precondition the matrix */
	//Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
	Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

	for (i = 0; i < N; i++) {
	    /* Copy b -> x which will be overwritten with soln. */
	    binp[i] = b[i];
	}

	/*  Solve the system with equality constraints  */
	fprintf(fp_log,
		" Inverting the matrix with equality constraints \n");
	time(&start);
	dgglse_(&N, &N2, &Ndih, MT, &LDA, C, &LDE, binp, d, phi, work,
		&lwork, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log,
		    " info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		    info);
	}
	fprintf(fp_log,
		" dgglse_() took %lf seconds or %lf minutes or %lf hours \n",
		diff, diff / (60.0), diff / (3600.00));

	/* Need to calculate the LU decomposition of MT */
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		index = index_Lpacked(i, j, N);
		MT[i + N * j] = M[index] + alpha[i + N * j] / beta;
	    }
	}

	for (i = 0; i < N; i++) {
	    /* Copy b -> x which will be overwritten with soln. */
	    binp[i] = b[i];
	}

	/*  Solve the regularized system using LU Decomposition  */
	fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
	time(&start);
	dgesv_(&N, &N_rhs, MT, &LDA, pivot, binp, &LDB, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log,
		    " info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		    info);
	}
	fprintf(fp_log,
		" dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
		diff, diff / (60.0), diff / (3600.00));

	/* Calculate the inverse of the regularized matrix */
	dgetri_(&N, MT, &LDA, pivot, work, &lwork, &info);

	/* Calculate the new alpha and beta quantities */
	for (i = 0; i < N; i++) {
	    alpha_old[i + N * i] = alpha[i + N * i];
	    alpha[i + N * i] =
		(1.0 / (phi[i] * phi[i] + (MT[i + N * i] / beta)));
	}

	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		index = index_Lpacked(i, j, N);
		MT2[i + N * j] = M[index];
		MT3[i + N * j] = 0.0;
	    }
	}

	/* Do matrix multiplication */
	dgemm_(&TRANSA, &TRANSB, &N, &N2, &N3, &aconst, MT, &LDA, MT2,
	       &LDB, &bconst, MT3, &LDC);

	trace = 0.0;
	for (i = 0; i < N; i++) {
	    trace += beta * MT3[i + N * i];
	}

	beta_old = beta;
	beta = (NT - trace) / (calc_Chi2(fp_log, sys, phi, b));

	flag_conv = TRUE;
	for (i = 0; i < N; i++) {
	    relerr =
		fabs((alpha[i + N * i] -
		      alpha_old[i + N * i]) / alpha_old[i + N * i]);
	    if (relerr > tau_alpha) {
		flag_conv = FALSE;
	    }
	}

	relerr = fabs((beta - beta_old) / beta_old);
	if (relerr > tau_beta) {
	    flag_conv = FALSE;
	}

	ctr++;

    } while ((flag_conv == FALSE) && (ctr < Nmax));

    /* print out the final regularization parameters */
    fprintf(fp_log, " Final Regularization Parameters: \n");
    for (i = 0; i < N; i++) {
	fprintf(fp_log, " alpha_f[%d] = %lf \n", i, alpha[i + N * i]);
    }
    fprintf(fp_log, " beta_f = %lf, nt = %d \n", beta, nt);
    if (ctr < Nmax) {
	fprintf(fp_log,
		" Optimization of regularization parameters converged after %d steps",
		ctr);
    } else {
	fprintf(fp_log,
		" Optimization of regularization parameters did not converge after Nmax = %d steps",
		Nmax);
    }

    /* Prepare to calculate the FINAL mN */
    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index] + alpha[i + N * j] / beta;
	}
    }

    /* d was killed by dgglse */
    d = (double *) ecalloc(Ndih, sizeof(double));
    if (d == NULL) {
	printf("STOP. d ecalloc failed. N: %d\n", N);
	exit(0);
    }

    /* set the constraints */
    for (i = 0; i < Ndih; i++) {
	d[i] = 0.00;
	for (k = 0; k < sys->N_Inter_Types; k++) {
	    k_flag = strcmp(B_DIHEDRAL, sys->Inter_Types[k].inter_type);
	    for (kk = 0; kk < sys->Inter_Types[k].N_pts; kk++) {
		row_k = sys->Inter_Types[k].i_0 + kk;
		if (k_flag == 0) {
		    C[i * N + row_k] = 1.0;
		} else {
		    C[i * N + row_k] = 0.0;
		}
	    }
	}
    }

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	binp[i] = b[i];
    }

    /*  Solve the system with equality constraints  */
    fprintf(fp_log, " Inverting the matrix with equality constraints \n");
    time(&start);
    dgglse_(&N, &N2, &Ndih, MT, &LDA, C, &LDE, binp, d, phi, work, &lwork,
	    &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgglse_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /* Need to calculate the LU decomposition of MT */
    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index] + alpha[i + N * j] / beta;
	}
    }

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	binp[i] = b[i];
    }

    /*  Solve the system with equality constraints  */
    fprintf(fp_log, " Inverting the matrix with equality constraints \n");
    time(&start);
    dgglse_(&N, &N2, &Ndih, MT, &LDA, C, &LDE, binp, d, phi, work, &lwork,
	    &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgglse_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

    free(MT);
    free(MT2);
    free(MT3);
    free(C);
    //free( d ); /* killed by dgglse */
    free(binp);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_constrain_tabdih_2(): 2 => use direct substition algorithm for constraints.
JFR - 12.20.13:  I found that this does not work well in practice
*****************************************************************************************/
int solv_lin_eqns_constrain_tabdih_2(FILE * fp_log, tW_system * sys)
{
    int i, j;
    int Ndih = 0;
    for (i = 0; i < sys->N_Inter_Types; i++) {
	if (strcmp(B_DIHEDRAL, sys->Inter_Types[i].inter_type) == 0) {
	    Ndih++;
	}
    }
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_row = N;
    int N_col = N;
    //int N_rhs = 1;             // no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    //int LDB = N;
    int LDC = Ndih;
    //int pivot[N];              // pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT;			// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double row_max[N_row], norm_col[N_col];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    /* variables for constraining dihedrals */
    int k, kk;
    int k_flag, row_k;
    double *C;
    double *d;
    double *binp;
    int lwork = 10 * N;
    double work[lwork];

    fprintf(fp_log, "In solv_lin_eqns.\n");

    MT = (double *) ecalloc(N_row * N_col, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }
    C = (double *) ecalloc(Ndih * N, sizeof(double));
    if (C == NULL) {
	printf("STOP. C ecalloc failed. N: %d\n", N);
	exit(0);
    }
    d = (double *) ecalloc(Ndih, sizeof(double));
    if (d == NULL) {
	printf("STOP. d ecalloc failed. N: %d\n", N);
	exit(0);
    }
    binp = (double *) ecalloc(N, sizeof(double));
    if (binp == NULL) {
	printf("STOP. binp ecalloc failed. N: %d\n", N);
	exit(0);
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N_col; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index];
	}
    }

    for (i = 0; i < Ndih; i++) {
	d[i] = 0.00;
	for (k = 0; k < sys->N_Inter_Types; k++) {
	    k_flag = strcmp(B_DIHEDRAL, sys->Inter_Types[k].inter_type);
	    for (kk = 0; kk < sys->Inter_Types[k].N_pts; kk++) {
		row_k = sys->Inter_Types[k].i_0 + kk;
		if (k_flag == 0) {
		    C[i * N + row_k] = 1.0;
		} else {
		    C[i * N + row_k] = 0.0;
		}
	    }
	}
    }

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	binp[i] = b[i];
    }

    /*  Solve the system using LU Decomposition  */
    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
    time(&start);
    dgglse_(&N_row, &N_col, &Ndih, MT, &LDA, C, &LDC, binp, d, phi, work,
	    &lwork, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgglse_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);
    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	binp[i] = b[i];
    }

    /*  Solve the system using LU Decomposition  */
    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
    time(&start);
    dgglse_(&N_row, &N_col, &Ndih, MT, &LDA, C, &LDC, binp, d, phi, work,
	    &lwork, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgglse_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

    free(MT);
    free(C);
    //free( d ); /* killed by dgglse */
    free(binp);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_constrain_tabdih_2_Bayesian_params(): same as above, but use the
regularization parameters calculated with Bayesian_saveparams() 
*****************************************************************************************/
int solv_lin_eqns_constrain_tabdih_2_Bayesian_params(FILE * fp_log,
						     tW_system * sys,
						     double *alpha,
						     double beta)
{
    int i, j;
    int Ndih = 0;
    for (i = 0; i < sys->N_Inter_Types; i++) {
	if (strcmp(B_DIHEDRAL, sys->Inter_Types[i].inter_type) == 0) {
	    Ndih++;
	}
    }
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_row = N;
    int N_col = N;
    //int N_rhs = 1;             // no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    //int LDB = N;
    int LDC = Ndih;
    //int pivot[N];              // pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT;			// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double row_max[N_row], norm_col[N_col];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    /* variables for constraining dihedrals */
    int k, kk;
    int k_flag, row_k;
    double *C;
    double *d;
    double *binp;
    int lwork = 10 * N;
    double work[lwork];

    fprintf(fp_log, "In solv_lin_eqns.\n");

    MT = (double *) ecalloc(N_row * N_col, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }
    C = (double *) ecalloc(Ndih * N, sizeof(double));
    if (C == NULL) {
	printf("STOP. C ecalloc failed. N: %d\n", N);
	exit(0);
    }
    d = (double *) ecalloc(Ndih, sizeof(double));
    if (d == NULL) {
	printf("STOP. d ecalloc failed. N: %d\n", N);
	exit(0);
    }
    binp = (double *) ecalloc(N, sizeof(double));
    if (binp == NULL) {
	printf("STOP. binp ecalloc failed. N: %d\n", N);
	exit(0);
    }

    /* prepare the matrix using the regularization parameters */
    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    if (i == j) {
		MT[i + N * j] = M[index] + alpha[i] / beta;
		printf("alpha[%d] = %lf \n", i, alpha[i]);
		printf("beta = %lf \n", beta);
	    } else {
		MT[i + N * j] = M[index];
	    }
	}
    }

    for (i = 0; i < Ndih; i++) {
	d[i] = 0.00;
	for (k = 0; k < sys->N_Inter_Types; k++) {
	    k_flag = strcmp(B_DIHEDRAL, sys->Inter_Types[k].inter_type);
	    for (kk = 0; kk < sys->Inter_Types[k].N_pts; kk++) {
		row_k = sys->Inter_Types[k].i_0 + kk;
		if (k_flag == 0) {
		    C[i * N + row_k] = 1.0;
		} else {
		    C[i * N + row_k] = 0.0;
		}
	    }
	}
    }

    /*  Precondition the matrix */
    //Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, FALSE, TRUE, sys);

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	binp[i] = b[i];
    }

    /*  Solve the system using LU Decomposition  */
    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
    time(&start);
    dgglse_(&N_row, &N_col, &Ndih, MT, &LDA, C, &LDC, binp, d, phi, work,
	    &lwork, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgglse_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
    //Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );
    Precondition(N, MT, b, phi, norm_col, row_max, TRUE, TRUE, sys);

    free(MT);
    free(C);
    //free( d ); /* killed by dgglse */
    free(binp);

    return info;
}

/*****************************************************************************************
solv_lin_eqns_rescalephi(): 
JFR - 12.20.13: This routine is a way to right Precondition, but it turns out not to work
well.  I have taken it out of the get_phi_X forks for now, but I also need to remove it
from the input options.
*****************************************************************************************/
int solv_lin_eqns_rescalephi(FILE * fp_log, tW_system * sys)
{
    int i, j;
    int index;
    /* system variables */
    double *b = sys->b;
    double *phi = sys->phi;
    double *M = sys->M;
    /* variables for standard matrix inversion */
    int N = sys->N_coeff;	// no. of lin. eqns.   A is N_row x N_col matrix w/ N=N_row
    int N_rhs = 1;		// no. of {x_i,b_i} pairs  s.t. A x_i = b_i        
    int LDA = N;
    int LDB = N;
    int pivot[N];		// pivot indices for permutation matrix
    int info;			// reports on success of the calculation
    double *MT;			// copy of the matrix, since it will be overwritten.  memory-wise this shouldn't exceed what the calculation already used earlier
    /* variables for preconditioning */
    double /*row_max[N], */ norm_col[N];
    /* variables for the timing things */
    time_t start, end;
    double diff;
    /* variables for force rescaling preconditioning */
    double phi_old[N];
    double tau_phi = sys->RESCALE_var.tau_phi;
    int Nphi = sys->RESCALE_var.Nmax;
    int ctr;
    int flag_conv;
    double relerr;

    fprintf(fp_log, "In solv_lin_eqns.\n");

    MT = (double *) ecalloc(N * N, sizeof(double));
    if (MT == NULL) {
	printf("STOP. MT ecalloc failed. N: %d\n", N);
	exit(0);
    }

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    index = index_Lpacked(i, j, N);
	    MT[i + N * j] = M[index];
	}
    }

    /* Do the normal calculation */

    /*  Precondition the matrix */
//  Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );

    for (i = 0; i < N; i++) {
	/* Copy b -> x which will be overwritten with soln. */
	phi[i] = b[i];
    }

    /*  Solve the system using LU Decomposition  */
    fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
    time(&start);
    dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
    time(&end);
    diff = difftime(end, start);
    if (info != 0) {
	fprintf(fp_log,
		" info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		info);
    }
    fprintf(fp_log,
	    " dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
	    diff, diff / (60.0), diff / (3600.00));

    /* Undo the Preconditioning */
//  Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );

    ctr = 0;
    do {

	/* reset the metric tensor */
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		index = index_Lpacked(i, j, N);
		MT[i + N * j] = M[index];
	    }
	}

	/*  Precondition the matrix normally */
//    Precondition( N, MT, b, phi, norm_col, row_max, FALSE /* => don't restore */, TRUE /* => operate on b and phi */, sys );

	/* calculate the scaling factor such that each force coefficient above 1 will now be 1 */
	for (i = 0; i < N; i++) {
	    if (fabs(phi[i]) > 1.0) {
		norm_col[i] = fabs(phi[i]);
	    } else {
		norm_col[i] = 1.0;
	    }
	}

	/* now do the force rescaling preconditioning */
	for (j = 0; j < N; j++) {
	    for (i = 0; i < N; i++) {
		MT[i + j * N] = MT[i + j * N] * norm_col[j];
	    }
	}

	for (i = 0; i < N; i++) {
	    /* Save the current force */
	    phi_old[i] = phi[i];
	    /* Copy b -> x which will be overwritten with soln. */
	    phi[i] = b[i];
	}

	/*  Solve the system using LU Decomposition  */
	fprintf(fp_log, " Inverting the matrix using LU decomposition \n");
	time(&start);
	dgesv_(&N, &N_rhs, MT, &LDA, pivot, phi, &LDB, &info);
	time(&end);
	diff = difftime(end, start);
	if (info != 0) {
	    fprintf(fp_log,
		    " info = %d after dgesv_().  The inverse condition number may be less than machine precision \n",
		    info);
	}
	fprintf(fp_log,
		" dgesv_() took %lf seconds or %lf minutes or %lf hours \n",
		diff, diff / (60.0), diff / (3600.00));

	/* Undo the normal Preconditioning */
//    Precondition( N, MT, b, phi, norm_col, row_max, TRUE /* => restore */, TRUE /* => operate on b and phi */, sys );

	/* undo the force rescaling preconditioning */
	for (j = 0; j < N; j++) {
	    phi[j] = phi[j] * norm_col[j];
	}

	flag_conv = TRUE;
	for (i = 0; i < N; i++) {
	    relerr = fabs((phi[i] - phi_old[i]) / phi_old[i]);
	    if (relerr > tau_phi) {
		flag_conv = FALSE;
	    }
	}

	ctr++;

    } while ((flag_conv == FALSE) && (ctr < Nphi));

    /* print out the final force scaling terms */
    for (i = 0; i < N; i++) {
	fprintf(fp_log, " scale_f[%d] = %lf \n", i, norm_col[i]);
    }
    if (ctr < Nphi) {
	fprintf(fp_log,
		" Optimization of force scaling parameters converged after %d steps",
		ctr);
    } else {
	fprintf(fp_log,
		" Optimization of force scaling parameters did not converge after Nphi = %d steps",
		Nphi);
    }

    free(MT);

    return info;
}
