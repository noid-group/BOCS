/**
@file wnoid_math.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Functions for basic math operations not basic enough to be in the standard c math
library. 
*/

//c library includes
#include <stdio.h>
#include <math.h>

//local includes
#include "wnoid_math.h"


/*****************************************************************************************
copy_vector(): Copies source into dest. Assumes dvec is an array of DIM elements.
*****************************************************************************************/
int copy_vector( dvec source, dvec dest )
{
  int i;

  for ( i=0; i<DIM; i++ ) { dest[i] = source[i]; }

  return DIM;
}

/*****************************************************************************************
copy_vector(): Copies source into dest. Assumes matrix is DIMxDIM elements.
*****************************************************************************************/
int copy_matrix(matrix source, matrix dest)
{
  int i;

  for (i=0; i<DIM; i++)
  {
    copy_vector(source[i], dest[i]);
  }

  return DIM;
}

/*****************************************************************************************
dot_prod(): Returns the dot product of two dvec's. 
*****************************************************************************************/
double dot_prod( dvec v1, dvec v2 )
{
  int j;
  double sum = 0.0;

  for ( j=0; j<DIM; j++ ) { sum += v1[j] * v2[j]; }

  return sum;
}


/*****************************************************************************************
calc_norm(): Returns the norm of an dvec.
*****************************************************************************************/
double calc_norm( dvec v )
{
  return ( sqrt( dot_prod( v, v ) ) );
}


/*****************************************************************************************
vect_diff(): Calculates diff = v1 - v2.
*****************************************************************************************/
int vect_diff( dvec v1, dvec v2, dvec diff )
{
  int i;

  for ( i=0; i<DIM; i++ ) { diff[i] = v1[i] - v2[i]; }

  return DIM;
}


/*****************************************************************************************
vect_sum(): Calculates sum = v1 + v2.
*****************************************************************************************/
int vect_sum( dvec v1, dvec v2, dvec sum )
{
  int i;

  for ( i=0; i<DIM; i++ ) { sum[i] = v1[i] + v2[i]; }

  return DIM;
}

/*****************************************************************************************
vect_inc(): Calculates v1 = v1 + v2 in place (by incrementing v1)
*****************************************************************************************/
void vect_inc( dvec v1, dvec v2 )
{
  int i;

  for ( i=0; i<DIM; i++ ) { v1[i] = v1[i] + v2[i]; }
}


/*****************************************************************************************
scal_times_vect(): Stores results from scalar-vector multiplication into the variable
result.
*****************************************************************************************/
int scal_times_vect( double sc, dvec vect, dvec result )
{
  int i;

  for ( i=0; i<DIM; i++ ) { result[i] = sc * vect[i]; }

  return DIM;
}


/*****************************************************************************************
cross_product(): Calculates v_cross = v1 x v2.
*****************************************************************************************/
int cross_product( dvec v1, dvec v2, dvec v_cross )
{
  v_cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v_cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v_cross[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return 0;
}


/*****************************************************************************************
calc_centered_diff(): Calculate the centered difference. 
*****************************************************************************************/
int calc_centered_diff( int Nmax, double dx, double *X, double *dX_dx )
{
  int i_x;
  double div = 0.5 / dx;

  for ( i_x=1; i_x<Nmax-1; i_x++ )
  { dX_dx[i_x] = ( X[i_x+1] - X[i_x-1] ) * div; }

  /* Treat end points special. */
  dX_dx[0]      = ( X[1]      - X[0] )      / dx;
  dX_dx[Nmax-1] = ( X[Nmax-1] - X[Nmax-2] ) / dx;

  return 0;
}


/*****************************************************************************************
calc_running_avg(): Calculates a running average.
*****************************************************************************************/
int calc_running_avg( int Nmax, int Nlength, double *X, double *smX )
{
  double norm;
  int  i, j, Nhalf;

  for ( i=0; i<Nmax; i++ )
  { smX[i] = 0.0; }

  Nhalf = ( int ) ( floor ( 0.5 * Nlength ) );

  norm = 1.0 / ( 2.0 * ( ( double ) Nhalf ) + 1.0 );

  for ( i=Nhalf; i<(Nmax-Nhalf); i++ )
  {
    for ( j = -1*Nhalf; j<(Nhalf+1); j++ )
    { smX[i] += norm * X[i+j]; }
  }

  for ( i=0;            i<Nhalf; i++ )  { smX[i] = X[i]; }
  for ( i=(Nmax-Nhalf); i<Nmax;  i++ )  { smX[i] = X[i]; }

  return 0;
}


/*****************************************************************************************
get_difference_unit_vector(): From x_i and x_j, calculates the differenc vector, the 
unit difference vector, and the norm of the difference vector. 
*****************************************************************************************/
int get_difference_unit_vector( dvec x_i, dvec x_j, dvec x_ij, dvec u_ij, double *r_ij )
{
  vect_diff( x_i, x_j, x_ij );                        /* x_ij = x_i - x_j            */
  *r_ij =  calc_norm( x_ij );                         /* r_ij = | x_ij |             */
  scal_times_vect( 1./(*r_ij), x_ij, u_ij );          /* Unit vector along j->i dir. */

  return 0;
}


/*****************************************************************************************
find_max_int(): Returns the maximum of the set {i,j}.
*****************************************************************************************/

int find_max_int( int i, int j )
{
  if ( i > j ) { return i; }

  return j;
}


/*****************************************************************************************
compare_ints():
*****************************************************************************************/
int compare_ints( const void *p, const void *q )
{
  return *(int *)p - *(int *)q;
}


/*****************************************************************************************
calc_running_avg_periodic(): Calculate running averages for periodic functions.
*****************************************************************************************/
int calc_running_avg_periodic( int Nmax, int Nlength, double *X, double *smX )
{
  double norm;
  int  i, j, Nhalf;
  double X_triple[Nmax*3-2];
  int skip_bin = 0;

  for ( i=0; i<Nmax; i++ )
  { smX[i] = 0.0; }

  for ( i=0; i<(3*Nmax); i++ )
  {
    if ( (i==Nmax) || (i==2*Nmax) )
    {
      skip_bin++;
      continue;
    }
    X_triple[(i-skip_bin)] = X[( i - (i/Nmax)*Nmax )];
  }

  Nhalf = ( int ) ( floor ( 0.5 * Nlength ) );

  norm = 1.0 / ( 2.0 * ( ( double ) Nhalf ) + 1.0 );

  for ( i=0; i<Nmax; i++ )
  {
    for ( j = -1*Nhalf; j<(Nhalf+1); j++ )
    { smX[i] += norm * X_triple[(Nmax-1)+i+j]; }
  }

  return 0;
}

/*************************************************************************************************
index_Lpacked(): Calculate the index of element aij in an NxN matrix in Lower-packed storage form.
**************************************************************************************************/
int index_Lpacked( int i, int j, int N )
{
  int index;

  if ( i <= j ) { index = j + i*N - (i*i+i)/2; }
  else { index = i + j*N - (j*j+j)/2; } 

  return index;
}


