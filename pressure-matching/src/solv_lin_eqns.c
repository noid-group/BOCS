/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file solv_lin_eqns.c
@author N. Dunn
@date Apr 24, 2014
@brief Translates pressure matching calculation into LAPACK library

*/

#include <math.h>

#include "safe_mem.h"
#include "solv_lin_eqns.h"


/*! The default solution method for pressure matching.	It uses
 the LAPACK function dgesv to find the basis coefficients. Returns
 the INFO variable from dgesv indicating the success of the calculation:
 *	INFO		(output) INTEGER
 *	= 0:	successful exit
 *	< 0:	if INFO = -i, the i-th argument had an illegal value
 *	> 0:	if INFO = i, U(i,i) is exactly zero.	The factorization
 *	has been completed, but the factor U is exactly
 *	singular, so the solution could not be computed.
*/
int solv_pres_lin_eqns( FILE *fp_log, pres_match_sys *pres_sys )
{
/**
@param fp_log The logfile in which to write documentation
@param pres_sys The main structure for the calculation

@return The INFO integer described above
*/

	int i, j;
	int N = (pres_sys->n_basis); // no. of lin. eqns. A is N_row x N_col matrix w/ N=N_row 
	int N_row = N; 	// matrix is square
	int N_col = N; 	// matrix is square
	int N_rhs = 1;	// no. of {x_i,b_i} pairs  s.t. A x_i = b_i	
	int LDA = pres_sys->n_basis;	// leading dimension of A as stored in AT.  NB A is square
	int pivot[N];			 // pivot indices for permutation matrix
	int info;
	double *MT;
	double *c, *row_max, **Q; //NJD 8/29/12 - changed these arrays to be dynamically allocated
	double * psi = pres_sys->psi;

	c = ecalloc(N_row, sizeof(double));
	row_max = ecalloc(N_col, sizeof(double));

	Q = ecalloc(N, sizeof(double));

	for (i=0; i<N; i++)
	{
		Q[i] = ecalloc(N, sizeof(double));

		for (j=0; j<N; j++)
		{
			Q[i][j] = pres_sys->Q[i][j];
		}
	}


	fprintf( fp_log, "In solv_pres_lin_eqns.\n" );

	MT = ( double * ) ecalloc( N_row*N_col, sizeof( double ) ); 

	/* Do simplest preconditioning by dividing each row by largest element. */
	for ( i=0; i<N_row; i++ )
	{ 
		row_max[i] = 0.; 

		c[i] = pres_sys->b[i]; 

		/* Find max. element. */
		for ( j=0; j<N_col; j++ )
		{ if ( fabs(Q[i][j]) > row_max[i] ) { row_max[i] = fabs(Q[i][j]); } }

		/* Divide through by max. element. */
		c[i] = c[i] / row_max[i];	
		for ( j=0; j<N_col; j++ ) { Q[i][j] = Q[i][j] / row_max[i]; } 
	} 

	for ( j=0; j<N_col; j++ )
	{
		/* Copy c -> x which will be overwritten with soln. */
		psi[j] = c[j];

		for ( i=0; i<N_row; i++ ) { MT[i+N_row*j] = Q[i][j]; }
	}

	fprintf( fp_log, " 1: N	  %d \n", N ); 
	fprintf( fp_log, " 2: N_rhs: %d \n", N_rhs ); 
	fprintf( fp_log, " 4: LDA  : %d \n", LDA ); 
	fprintf( fp_log, " 7: LDB  : %d \n", N ); 
	fprintf( fp_log, "\n" );

	dgesv_( &N, &N_rhs, MT, &LDA, pivot, psi, &N, &info );


	free( MT );
	free( c );
	free( row_max );

	for (i=0; i<N; i++)
	{
		free( Q[i] );
	}
	free( Q );

	return info; 
}


/*! Solves for the pressure coefficients iteratively using the LAPACK function
	dgesvx, which outputs the forward and backward error associated with the 
	solution.	Returns the INFO variable from dgesvx indicating the success 
	of the calculation:
 *	INFO		(output) INTEGER
 *	= 0:	successful exit
 *	< 0:	if INFO = -i, the i-th argument had an illegal value
 *	> 0:	if INFO = i, U(i,i) is exactly zero.	The factorization
 *		has been completed, but the factor U is exactly
 *		singular, so the solution could not be computed.
 */
int solv_pres_lin_eqns_error( FILE *fp_log, pres_match_sys *pres_sys )
{
/**
@param fp_log The logfile in which to write documentation
@param pres_sys The main structure for the calculation

@return The INFO integer described above
*/
	int i, j;
	int N = (pres_sys->n_basis);	 // no. of lin. eqns.	 A is N_row x N_col matrix w/ N=N_row 
	int N_row = N; 	// matrix is square
	int N_col = N; 	// matrix is square
	int N_rhs = 1;	 			 // no. of {x_i,b_i} pairs	s.t. A x_i = b_i	
	int LDA = pres_sys->n_basis;	// leading dimension of A as stored in AT.	NB A is square
	int pivot[N];			 // pivot indices for permutation matrix
	int info;
	double *MT;
	double *c, *row_max, **Q; //NJD 8/29/12 - changed these arrays to be dynamically allocated
	double * psi = pres_sys->psi;

	/* New variables for dgesvx error estimate */
	double *AF;
	int LDAF = pres_sys->n_basis;
	double *R, *C, *FERR, *BERR, *WORK;
	int *IWORK;
	double RCOND;
	char std = 'N';


	R = ecalloc(N, sizeof(double));
	C = ecalloc(N, sizeof(double));
	AF = ecalloc(LDAF*N, sizeof(double));
	FERR = ecalloc(N_rhs, sizeof(double));
	BERR = ecalloc(N_rhs, sizeof(double));
	WORK = ecalloc(4*N, sizeof(double));
	IWORK = ecalloc(N, sizeof(int));


	c = ecalloc(N_row, sizeof(double));

	row_max = ecalloc(N_col, sizeof(double));

	Q = ecalloc(N, sizeof(double));

	for (i=0; i<N; i++)
	{
		Q[i] = ecalloc(N, sizeof(double));

		for (j=0; j<N; j++)
		{
			Q[i][j] = pres_sys->Q[i][j];
		}
	}


	fprintf( fp_log, "In solv_pres_lin_eqns.\n" );

	MT = ( double * ) ecalloc( N_row*N_col, sizeof( double ) ); 

	/* Do simplest preconditioning by dividing each row by largest element. */
	for ( i=0; i<N_row; i++ )
	{ 
		row_max[i] = 0.; 

		c[i] = pres_sys->b[i]; 

		/* Find max. element. */
		for ( j=0; j<N_col; j++ )
		{ if ( fabs(Q[i][j]) > row_max[i] ) { row_max[i] = fabs(Q[i][j]); } }

		/* Divide through by max. element. */
		c[i] = c[i] / row_max[i];	
		for ( j=0; j<N_col; j++ ) { Q[i][j] = Q[i][j] / row_max[i]; } 
	} 

	for ( j=0; j<N_col; j++ )
	{
		/* Copy c -> x which will be overwritten with soln. */
		psi[j] = c[j];

		for ( i=0; i<N_row; i++ )
		{ MT[i+N_row*j] = Q[i][j]; }
	}

	fprintf( fp_log, " 1: N	    %d \n", N ); 
	fprintf( fp_log, " 2: N_rhs: %d \n", N_rhs ); 
	fprintf( fp_log, " 4: LDA  : %d \n", LDA ); 
	fprintf( fp_log, " 7: LDB  : %d \n", N ); 
	fprintf( fp_log, "\n" );

	//dgesv_( &N, &N_rhs, MT, &LDA, pivot, psi, &N, &info );
	dgesvx_(&(std), &(std), &N, &N_rhs, MT, &LDA, AF, &LDAF, pivot, &(std), R, C, c, &N, psi, &N, &RCOND, FERR, BERR, WORK, IWORK, &info );



	fprintf(stdout, "Condition number: %f\nForward Error: %f\nBackward Error: %f\n", 1/RCOND, FERR[0], BERR[0]);

	free( MT );
	free( c );
	free( row_max );

	for (i=0; i<N; i++)
	{
		free( Q[i] );
	}
	free( Q );

	free(R);
	free(C);
	free(AF);
	free(FERR);
	free(BERR);
	free(WORK);
	free(IWORK);

	return info; 
}

/*!	A slight modification of the default solution method for pressure 
 matching.	It uses the LAPACK function dgesv to find the basis 
 coefficients. Prior to this, it weights each bin in the metric tensor
 by the sampling count of that bin. Returns the INFO variable from dgesv
 indicating the success of the calculation:
 *	INFO		(output) INTEGER
 *	= 0:	successful exit
 *	< 0:	if INFO = -i, the i-th argument had an illegal value
 *	> 0:	if INFO = i, U(i,i) is exactly zero.	The factorization
 *		has been completed, but the factor U is exactly
 *		singular, so the solution could not be computed.
 */
int solv_pres_lin_eqns_weighted( FILE *fp_log, pres_match_sys *pres_sys )
{
/**
@param fp_log The logfile in which to write documentation
@param pres_sys The main structure for the calculation

@return The INFO integer described above
*/

	int i, j;
	int N = (pres_sys->n_basis);	 // no. of lin. eqns.	 A is N_row x N_col matrix w/ N=N_row 
	int N_row = N; 	// matrix is square
	int N_col = N; 	// matrix is square
	int N_rhs = 1;	 			 // no. of {x_i,b_i} pairs	s.t. A x_i = b_i	
	int LDA = pres_sys->n_basis;	// leading dimension of A as stored in AT.	NB A is square
	int pivot[N];			 // pivot indices for permutation matrix
	int info;
	double *MT;
	double *MW;
	double *work;
	double *c, *row_max, **Q; //NJD 8/29/12 - changed these arrays to be dynamically allocated
	double * psi = pres_sys->psi;
	double * y;


	c = ecalloc(N_row, sizeof(double));
	row_max = ecalloc(N_col, sizeof(double));

	y = ecalloc (N_row, sizeof(double));
	MW = ecalloc(N_row*N_col, sizeof(double));
	work = ecalloc(N_row*N_col, sizeof(double));

	// Weight each element of the metric tensor by the sampling
	// frequency
	for (i=0; i<N_col; i++)
	{
	for (j=0; j<N_row; j++)
	{
		if (i == j)
		{
			MW[i*N_col + j] = pres_sys->g_cnt[j];
		} else
		{
			MW[i*N_col + j] = 0;
		}
	}
	}

	Q = ecalloc(N, sizeof(double));

	for (i=0; i<N; i++)
	{
		Q[i] = ecalloc(N, sizeof(double));

		for (j=0; j<N; j++)
		{
			Q[i][j] = pres_sys->Q[i][j];
		}
	}


	fprintf( fp_log, "In solv_pres_lin_eqns.\n" );

	MT = ( double * ) ecalloc( N_row*N_col, sizeof( double ) ); 

	/* Do simplest preconditioning by dividing each row by largest element. */
	for ( i=0; i<N_row; i++ )
	{ 
		row_max[i] = 0.; 

		c[i] = pres_sys->b[i]; 

		/* Find max. element. */
		for ( j=0; j<N_col; j++ )
		{ if ( fabs(Q[i][j]) > row_max[i] ) { row_max[i] = fabs(Q[i][j]); } }

		/* Divide through by max. element. */
		c[i] = c[i] / row_max[i];	
		for ( j=0; j<N_col; j++ ) { Q[i][j] = Q[i][j] / row_max[i]; } 
	} 

	for ( j=0; j<N_col; j++ )
	{
		/* Copy c -> x which will be overwritten with soln. */
		psi[j] = c[j];

		for ( i=0; i<N_row; i++ )
		{ MT[i+N_row*j] = Q[i][j]; }
	}

	fprintf( fp_log, " 1: N	    %d \n", N ); 
	fprintf( fp_log, " 2: N_rhs: %d \n", N_rhs ); 
	fprintf( fp_log, " 4: LDA  : %d \n", LDA ); 
	fprintf( fp_log, " 7: LDB  : %d \n", N ); 
	fprintf( fp_log, "\n" );


	dgesv_( &N, &N_rhs, MT, &LDA, pivot, psi, &N, &info );

	for (i=0; i<N; i++)
	{
		fprintf(stdout, "Y: %f, X: %f\n", y[i], psi[i]);
	}


	free( MT );
	free( c );
	free( row_max );

	for (i=0; i<N; i++)
	{
		free( Q[i] );
	}
	free( Q );

	return info; 
}



/*! Returns an array containing the intercept and slope of the line of best fit.
	The solution is determined using the LAPACK function dgelss */
double * get_best_fit(double* dataX, double* dataY, int len )
{
/**
@param dataX Array containing the independent variable in the dataset
@param dataY Array containing the dependent variable in the dataset
@param len The number of points in the dataset

@return A 2-element array of doubles. Element [0] is b and element [1] is m in the equation y=mx+b
*/
	int i;
	int M = len;
	int LDA = M;
	int LDB = M;
	int N = 2; // Linear fit
	int nrhs = 1;
	double *A;
	double *B;
	double *S;
	double rcond = -1.0;
	int rank;
	int lwork = 10 * M;
	int info;
	double *work;


	A = (double*) ecalloc(2*M, sizeof(double));
	B = (double*) ecalloc(M, sizeof(double));
	S = (double*) ecalloc(M, sizeof(double));
	work = (double*) ecalloc(lwork, sizeof(double));


	for (i=0; i<M; i++)
	{
		A[i+0*M] = 1;
		A[i+1*M] = dataX[i];
		B[i] = dataY[i];
	}

	dgelss_( &M, &N, &nrhs, A, &LDA, B, &LDB, S, &rcond, &rank, work, &lwork, &info );

	free(A);
	free(S);
	free(work);

	return B;
}

