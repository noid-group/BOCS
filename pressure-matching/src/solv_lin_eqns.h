/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file solv_lin_eqns.h
@author N. Dunn
@date Apr 30, 2014
@brief Functions for solving the set of pressure matching equations

This file contains headers for the functions that couple to LAPACK
and solve the set of linear equations to determine the pressure
matching coefficients.
 */


#include "pressure_matching_types.h"

#ifndef SOLV_LIN_EQN
#define SOLV_LIN_EQN

/* Coupling functions  */

/*! The default solution method for pressure matching.  It uses
 the LAPACK function dgesv to find the basis coefficients. Returns
 the INFO variable from dgesv indicating the success of the calculation:
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
 *                has been completed, but the factor U is exactly
 *                singular, so the solution could not be computed.
*/
int solv_pres_lin_eqns( FILE *fp_log, pres_match_sys *pres_sys );

/*! Solves for the pressure coefficients iteratively using the LAPACK function
  dgesvx, which outputs the forward and backward error associated with the 
  solution.  Returns the INFO variable from dgesvx indicating the success 
  of the calculation:
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
 *                has been completed, but the factor U is exactly
 *                singular, so the solution could not be computed.
 */
int solv_pres_lin_eqns_error( FILE *fp_log, pres_match_sys *pres_sys );

/*!   A slight modification of the default solution method for pressure 
 matching.  It uses the LAPACK function dgesv to find the basis 
 coefficients. Prior to this, it weights each bin in the metric tensor
 by the sampling count of that bin. Returns the INFO variable from dgesv
 indicating the success of the calculation:
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
 *                has been completed, but the factor U is exactly
 *                singular, so the solution could not be computed.
 */
int solv_pres_lin_eqns_weighted( FILE *fp_log, pres_match_sys *pres_sys );

/*! Returns an array containing the intercept and slope of the line of best fit.
  The solution is determined using the LAPACK function dgelss */
double * get_best_fit(double* dataX, double* dataY, int len );


/* Extern LAPACK functions */

/*!  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGESV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   

    The LU decomposition with partial pivoting and row interchanges is   
    used to factor A as   
       A = P * L * U,   
    where P is a permutation matrix, L is unit lower triangular, and U is   
    upper triangular.  The factored form of A is then used to solve the   
    system of equations A * X = B.

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N coefficient matrix A.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (output) INTEGER array, dimension (N)   
            The pivot indices that define the permutation matrix P;   
            row i of the matrix was interchanged with row IPIV(i).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS matrix of right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization   
                  has been completed, but the factor U is exactly   
                  singular, so the solution could not be computed.   

    =====================================================================  
*/
extern void dgesv_( int *N, int *N_rhs, double AT[], int *LDA,
		int pivot[], double x[], int * N_row, int *info );

/*!  -- LAPACK driver routine (version 3.1) --   
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
       November 2006   


    Purpose   
    =======   

    DGELSS computes the minimum norm solution to a real linear least   
    squares problem:   

    Minimize 2-norm(| b - A*x |).   

    using the singular value decomposition (SVD) of A. A is an M-by-N   
    matrix which may be rank-deficient.   

    Several right hand side vectors b and solution vectors x can be   
    handled in a single call; they are stored as the columns of the   
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix   
    X.   

    The effective rank of A is determined by treating as zero those   
    singular values which are less than RCOND times the largest singular   
    value.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A. N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X. NRHS >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the first min(m,n) rows of A are overwritten with   
            its right singular vectors, stored rowwise.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the M-by-NRHS right hand side matrix B.   
            On exit, B is overwritten by the N-by-NRHS solution   
            matrix X.  If m >= n and RANK = n, the residual   
            sum-of-squares for the solution in the i-th column is given   
            by the sum of squares of elements n+1:m in that column.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,max(M,N)).   

    S       (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The singular values of A in decreasing order.   
            The condition number of A in the 2-norm = S(1)/S(min(m,n)).   

    RCOND   (input) DOUBLE PRECISION   
            RCOND is used to determine the effective rank of A.   
            Singular values S(i) <= RCOND*S(1) are treated as zero.   
            If RCOND < 0, machine precision is used instead.   

    RANK    (output) INTEGER   
            The effective rank of A, i.e., the number of singular values   
            which are greater than RCOND*S(1).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= 1, and also:   
            LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )   
            For good performance, LWORK should generally be larger.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  the algorithm for computing the SVD failed to converge;   
                  if INFO = i, i off-diagonal elements of an intermediate   
                  bidiagonal form did not converge to zero.   

    ===================================================================== 
*/
extern void dgelss_( int *M, int *N, int *NRHS, double A[], int *LDA, 
		double B[], int *LDB, double S[], double *RCOND, 
		int *RANK, double WORK[], int *LWORK, int *INFO ); 

/*!  -- LAPACK driver routine (version 3.1) --   
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
       November 2006   


    Purpose   
    =======   

    DGESVX uses the LU factorization to compute the solution to a real   
    system of linear equations   
       A * X = B,   
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'E', real scaling factors are computed to equilibrate   
       the system:   
          TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B   
          TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B   
          TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B   
       Whether or not the system will be equilibrated depends on the   
       scaling of the matrix A, but if equilibration is used, A is   
       overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')   
       or diag(C)*B (if TRANS = 'T' or 'C').   

    2. If FACT = 'N' or 'E', the LU decomposition is used to factor the   
       matrix A (after equilibration if FACT = 'E') as   
          A = P * L * U,   
       where P is a permutation matrix, L is a unit lower triangular   
       matrix, and U is upper triangular.   

    3. If some U(i,i)=0, so that U is exactly singular, then the routine   
       returns with INFO = i. Otherwise, the factored form of A is used   
       to estimate the condition number of the matrix A.  If the   
       reciprocal of the condition number is less than machine precision,   
       INFO = N+1 is returned as a warning, but the routine still goes on   
       to solve for X and compute error bounds as described below.   

    4. The system of equations is solved for X using the factored form   
       of A.   

    5. Iterative refinement is applied to improve the computed solution   
       matrix and calculate error bounds and backward error estimates   
       for it.   

    6. If equilibration was used, the matrix X is premultiplied by   
       diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so   
       that it solves the original system before equilibration.   

    Arguments   
    =========   

    FACT    (input) CHARACTER*1   
            Specifies whether or not the factored form of the matrix A is   
            supplied on entry, and if not, whether the matrix A should be   
            equilibrated before it is factored.   
            = 'F':  On entry, AF and IPIV contain the factored form of A.   
                    If EQUED is not 'N', the matrix A has been   
                    equilibrated with scaling factors given by R and C.   
                    A, AF, and IPIV are not modified.   
            = 'N':  The matrix A will be copied to AF and factored.   
            = 'E':  The matrix A will be equilibrated if necessary, then   
                    copied to AF and factored.   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Transpose)   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is   
            not 'N', then A must have been equilibrated by the scaling   
            factors in R and/or C.  A is not modified if FACT = 'F' or   
            'N', or if FACT = 'E' and EQUED = 'N' on exit.   

            On exit, if EQUED .ne. 'N', A is scaled as follows:   
            EQUED = 'R':  A := diag(R) * A   
            EQUED = 'C':  A := A * diag(C)   
            EQUED = 'B':  A := diag(R) * A * diag(C).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    AF      (input or output) DOUBLE PRECISION array, dimension (LDAF,N)   
            If FACT = 'F', then AF is an input argument and on entry   
            contains the factors L and U from the factorization   
            A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then   
            AF is the factored form of the equilibrated matrix A.   

            If FACT = 'N', then AF is an output argument and on exit   
            returns the factors L and U from the factorization A = P*L*U   
            of the original matrix A.   

            If FACT = 'E', then AF is an output argument and on exit   
            returns the factors L and U from the factorization A = P*L*U   
            of the equilibrated matrix A (see the description of A for   
            the form of the equilibrated matrix).   

    LDAF    (input) INTEGER   
            The leading dimension of the array AF.  LDAF >= max(1,N).   

    IPIV    (input or output) INTEGER array, dimension (N)   
            If FACT = 'F', then IPIV is an input argument and on entry   
            contains the pivot indices from the factorization A = P*L*U   
            as computed by DGETRF; row i of the matrix was interchanged   
            with row IPIV(i).   

            If FACT = 'N', then IPIV is an output argument and on exit   
            contains the pivot indices from the factorization A = P*L*U   
            of the original matrix A.   

            If FACT = 'E', then IPIV is an output argument and on exit   
            contains the pivot indices from the factorization A = P*L*U   
            of the equilibrated matrix A.   

    EQUED   (input or output) CHARACTER*1   
            Specifies the form of equilibration that was done.   
            = 'N':  No equilibration (always true if FACT = 'N').   
            = 'R':  Row equilibration, i.e., A has been premultiplied by   
                    diag(R).   
            = 'C':  Column equilibration, i.e., A has been postmultiplied   
                    by diag(C).   
            = 'B':  Both row and column equilibration, i.e., A has been   
                    replaced by diag(R) * A * diag(C).   
            EQUED is an input argument if FACT = 'F'; otherwise, it is an   
            output argument.   

    R       (input or output) DOUBLE PRECISION array, dimension (N)   
            The row scale factors for A.  If EQUED = 'R' or 'B', A is   
            multiplied on the left by diag(R); if EQUED = 'N' or 'C', R   
            is not accessed.  R is an input argument if FACT = 'F';   
            otherwise, R is an output argument.  If FACT = 'F' and   
            EQUED = 'R' or 'B', each element of R must be positive.   

    C       (input or output) DOUBLE PRECISION array, dimension (N)   
            The column scale factors for A.  If EQUED = 'C' or 'B', A is   
            multiplied on the right by diag(C); if EQUED = 'N' or 'R', C   
            is not accessed.  C is an input argument if FACT = 'F';   
            otherwise, C is an output argument.  If FACT = 'F' and   
            EQUED = 'C' or 'B', each element of C must be positive.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit,   
            if EQUED = 'N', B is not modified;   
            if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by   
            diag(R)*B;   
            if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is   
            overwritten by diag(C)*B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)   
            If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X   
            to the original system of equations.  Note that A and B are   
            modified on exit if EQUED .ne. 'N', and the solution to the   
            equilibrated system is inv(diag(C))*X if TRANS = 'N' and   
            EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'   
            and EQUED = 'R' or 'B'.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    RCOND   (output) DOUBLE PRECISION   
            The estimate of the reciprocal condition number of the matrix   
            A after equilibration (if done).  If RCOND is less than the   
            machine precision (in particular, if RCOND = 0), the matrix   
            is singular to working precision.  This condition is   
            indicated by a return code of INFO > 0.   

    FERR    (output) DOUBLE PRECISION array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j)   
            is an estimated upper bound for the magnitude of the largest   
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (4*N)   
            On exit, WORK(1) contains the reciprocal pivot growth   
            factor norm(A)/norm(U). The "max absolute element" norm is   
            used. If WORK(1) is much less than 1, then the stability   
            of the LU factorization of the (equilibrated) matrix A   
            could be poor. This also means that the solution X, condition   
            estimator RCOND, and forward error bound FERR could be   
            unreliable. If factorization fails with 0<INFO<=N, then   
            WORK(1) contains the reciprocal pivot growth factor for the   
            leading INFO columns of A.   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, and i is   
                  <= N:  U(i,i) is exactly zero.  The factorization has   
                         been completed, but the factor U is exactly   
                         singular, so the solution and error bounds   
                         could not be computed. RCOND = 0 is returned.   
                  = N+1: U is nonsingular, but RCOND is less than machine   
                         precision, meaning that the matrix is singular   
                         to working precision.  Nevertheless, the   
                         solution and error bounds are computed because   
                         there are a number of situations where the   
                         computed solution can be more accurate than the   
                         value of RCOND would suggest.   

    =====================================================================
*/
extern void dgesvx_( char* fact, char* trans, int* n, int* nrhs,
		double* a, int* lda, double* af, int* ldaf,
		int ipiv[], char* equed, double* r, double* c,
		double* b, int* ldb, double x[], int* ldx,
		double* rcond, double* ferr, double* berr, double* work,
		int* iwork, int *info );



#endif
