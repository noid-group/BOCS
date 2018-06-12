/**
@file table_functions.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Functions for generating table files for gromacs 4 use
*/

//c library includes
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//local includes
#include "cgff_types.h"
#include "wnoid_math.h"
#include "table_types.h"
#include "safe_mem.h"
#include "table_functions.h"



int get_inter_type( char *argv[], int argc )
{
  int type_index;

  if      ( strcmp( argv[3], NB_NAME    ) == 0 ) { type_index = NB_INDEX;    }
  else if ( strcmp( argv[3], BOND_NAME  ) == 0 ) { type_index = BOND_INDEX;  }
  else if ( strcmp( argv[3], ANGLE_NAME ) == 0 ) { type_index = ANGLE_INDEX; }
  else if ( strcmp( argv[3], DIHED_NAME ) == 0 ) { type_index = DIHED_INDEX; }
  else if ( strcmp( argv[3], INTRA_NAME ) == 0 ) { type_index = INTRA_INDEX; }
  else    { print_type_error_message( argv[3] ); exit( EXIT_FAILURE ); }

  if (type_index == BOND_INDEX && argc != 7)
  {
    fprintf(stderr, "ERROR: Wrong # of arguments for bonded force. Terminating.\n"); 
    fprintf(stderr, "Expected: EXE INFILE OUTFILE bond BASIS RMIN RMAX\n\n");
    exit(1);
  } else if (type_index == INTRA_INDEX && argc != 6 )
  {
    fprintf(stderr, "ERROR: Wrong # of arguments for intra force. Terminating.\n");
    fprintf(stderr, "Expected: EXE INPUT OUTPUT intra BASIS RMAX\n\n"); 
    exit(1);
  } else if (type_index == NB_INDEX && argc != 6)
  {
    fprintf(stderr, "ERROR: Wrong # of arguments for nb force. Terminating.\n");                
    fprintf(stderr, "Expected: EXE INPUT OUTPUT nb BASIS RMAX\n\n");
    exit(1);
  } else if (type_index == DIHED_INDEX && argc != 5)
  {
    fprintf(stderr, "ERROR: Wrong # arguments for angle or dihedral force. Terminating.\n"); 
    fprintf(stderr, "Expected: EXE INPUT OUTPUT dihed BASIS\n\n");
    exit(1);
  } else if (type_index == ANGLE_INDEX && argc != 5)
  {
    fprintf(stderr, "ERROR: Wrong # of arguments for angle force. Terminating.\n");
    fprintf(stderr, "Expected: EXE INPUT OUTPUT angle BASIS\n\n");
    exit(1);
  }

  return type_index;
}

int get_basis_type( char *argv[] )
{
  int type_basis;

  if      ( strcmp( argv[4], DELTA_BASIS    ) == 0 ) { type_basis = DELTA_INDEX;  }
  else if ( strcmp( argv[4], LINEAR_BASIS  ) == 0 ) { type_basis = LINEAR_INDEX;  }
  else if ( strcmp( argv[4], BSPLINE_BASIS ) == 0 ) { type_basis = BSPLINE_INDEX; }
  else    { print_type_index_error( type_basis ); exit( EXIT_FAILURE ); }

  return type_basis;
}


void get_table_parameters( Input *input, char *argv[] )
{
  double r_min, r_max;

  switch ( input->type_index )
  {
    case NB_INDEX:       r_min = 0.0;             r_max = atof( argv[5] );   break;
    case BOND_INDEX:     r_min = atof( argv[5] ); r_max = atof( argv[6] );   break;
    case ANGLE_INDEX:    r_min = 0.0;             r_max = 180.0;             break;
    case DIHED_INDEX:    r_min = -180.0;          r_max = 180.0;             break;
    case INTRA_INDEX:    r_min = 0.0;             r_max = atof( argv[5] );   break;
    default:  print_type_index_error( input->type_index ); exit( EXIT_FAILURE );
  }

  input->r_min = r_min;
  input->r_max = r_max;

}


void read_input_forces( int n_pts, Arrays *arrays, FILE *fp_in )
{
  int i;
  int test_sscanf;
  double dr_temp;
  char line[250];

  /* Start from the beginning of the file. */
  rewind( fp_in );

  /* Read in f(r). */
  for ( i=0; i<n_pts; i++ )
  {
    fgets( line, 250, fp_in );
    test_sscanf = sscanf( line, "%lf %lf", &(arrays->r[i]), &(arrays->f[i]) );
    if ( test_sscanf != 2 )  { print_read_input_error1( i ); exit( EXIT_FAILURE ); }
  }

  /* Make sure grid spacing is even. */
  for ( i=0; i<(n_pts-1); i++ )
  {
    dr_temp = arrays->r[i+1] - arrays->r[i];
    if ( fabs( (arrays->r[1] - arrays->r[0]) - dr_temp ) > FLOAT_EPS )
    { print_read_input_error2( i ); exit( EXIT_FAILURE ); }
  }

}

int get_potential( Input *input, Arrays *arrays )
{
  int total_n_pts;

  switch ( input->type_index )
  {
    case NB_INDEX:    total_n_pts = get_nb_potential( input, arrays ); break;
    case BOND_INDEX:  total_n_pts = get_bond_potential( input, arrays ); break;
    case ANGLE_INDEX: total_n_pts = get_angle_potential( input, arrays ); break;
    case DIHED_INDEX: total_n_pts = get_dihedral_potential( input, arrays ); break;
    case INTRA_INDEX: total_n_pts = get_bond_potential( input, arrays ); break;
    default:  print_type_index_error( input->type_index ); exit( EXIT_FAILURE );
  }

  return total_n_pts;
}

int get_nb_potential( Input *input, Arrays *arrays )
{
  int flag;
  int total_n_pts;
  double slope;
  double dr = arrays->r[1] - arrays->r[0];
  double *r_temp, *f_temp, *neg_f;

  /* Allocate memory for temporary arrays. */
  r_temp = ecalloc( input->n_pts, sizeof( double ) );
  f_temp = ecalloc( input->n_pts, sizeof( double ) );

  /* Copy r and f into temporary arrays. */
  copy_array( input->n_pts, arrays->r, r_temp );
  copy_array( input->n_pts, arrays->f, f_temp );

  /* Add the number of points needed to n_pts when interpolating into hard core */
  /* and mark place in resized array with flag.                                 */
  flag = get_total_n_pts( input->n_pts, arrays->r, &total_n_pts, &(input->r_min) );

  /* Allocate memory for full size arrays. */
  arrays->r     = (double  *) erealloc( arrays->r, total_n_pts*sizeof(double) );
  arrays->f     = (double  *) erealloc( arrays->f, total_n_pts*sizeof(double) );
  arrays->u     = (double  *) erealloc( arrays->u, total_n_pts*sizeof(double) );
  neg_f = ecalloc( total_n_pts, sizeof( double ) );

  /* Fill in r. */
  fill_r_nb( total_n_pts, flag, arrays->r, r_temp, dr, input->r_min );

  /* Calculate slope used for linear interpolation. */
  //slope = ( f_temp[1] - f_temp[0] ) / dr; 
  slope = input->slopel;

  /* Fill in f. */
  fill_f_nb( total_n_pts, flag, arrays->f, f_temp, slope, dr );

  /* smooth f */
  smooth_forces( total_n_pts, input, arrays );

  /* Get the negative of f. */
  scal_vect_mult( total_n_pts, -1.0, arrays->f, neg_f );

  /* Calculate the potential. */
  integrate_trapezoid( arrays->r, arrays->u, neg_f, total_n_pts );

  /* Shift u so that it is zero at r[total_n_pts-1]. */
  shift_potential( total_n_pts, arrays->u );

  /* For consistency, caculate the force from the potential. */
  calc_centered_diff( total_n_pts, dr, arrays->u, arrays->f );
  scal_vect_mult( total_n_pts, -1.0, arrays->f, arrays->f );

  return total_n_pts;

}


void shift_potential( int total_n_pts, double *u )
{
  int i;
  double shift;

  shift = u[total_n_pts-1];

  for ( i=0; i<total_n_pts; i++ )
  { u[i] -= shift; }

}


int integrate_trapezoid (double x[], double *int_y, double y[], int no_of_points)
{
  int i;
  double dx;

  dx = x[2] - x[1];
  int_y[0] = 0;
  for (i=1; i<=no_of_points-1; i++)
  { int_y[i] = int_y[i-1] +  0.5 * dx * (y[i] + y[i-1]); }

  return 0;
}


void scal_vect_mult( int total_n_pts, double scalar, double original[], double *result )
{
  int i;
  
  for ( i=0; i<total_n_pts; i++ )
  { result[i] = scalar * original[i]; }

}


void fill_f_nb( int total_n_pts, int flag, double *f, double *f_temp, double slope, double dr )
{
  int i, j;

  /* Extrapolate back. */
  f[flag] = f_temp[0];
  for( i=(flag-1); i>-1; i-- )
  { f[i] = f[i+1] - dr * slope; }

  /* Copy the rest. */
  for ( i=flag, j=0; i<total_n_pts; i++, j++ )
  { f[i] = f_temp[j]; }

}


void fill_r_nb( int total_n_pts, int flag, double *r, double *r_temp, double dr, double r_min )
{
  int i, j;
  double r_test = r_min;

  for ( i=0, j=0; i<total_n_pts; i++ )
  {
    if ( i < flag )
    { r[i] = r_test; }
    else
    { 
      if ( fabs( r_test - r_temp[j] ) > FLOAT_EPS )
      { print_error_fill_r( r_test, r_temp[j] ); exit( EXIT_FAILURE ); }

      r[i] = r_temp[j]; 
      j++;
    }

    r_test += dr;
  }
}


int get_bond_potential( Input *input, Arrays *arrays )
{
  int total_n_pts;
  double dr = arrays->r[1] - arrays->r[0];
  double *r_temp, *f_temp, *neg_f;

  /* Allocate memory for temporary arrays. */
  r_temp = ecalloc( input->n_pts, sizeof( double ) );
  f_temp = ecalloc( input->n_pts, sizeof( double ) );

  /* Copy r and f into temporary arrays. */
  copy_array( input->n_pts, arrays->r, r_temp );
  copy_array( input->n_pts, arrays->f, f_temp );

  /* Find the total number of points. */ 
  total_n_pts = get_total_n_pts_bond( input->n_pts, &(input->r_min), &(input->r_max), dr, arrays->r );

  /* Allocate memory for full size arrays. */
  arrays->r     = (double  *) erealloc( arrays->r, total_n_pts*sizeof(double) );
  arrays->f     = (double  *) erealloc( arrays->f, total_n_pts*sizeof(double) );
  arrays->u     = (double  *) erealloc( arrays->u, total_n_pts*sizeof(double) );
  neg_f         = ecalloc( total_n_pts, sizeof( double ) );

  /* Fill in r. */
  fill_r_bond( total_n_pts, arrays->r, r_temp, input );

  /* Fill in f. */
  fill_f_bond( total_n_pts, arrays->f, f_temp, r_temp, dr, input );

  /* smooth f */
  smooth_forces( total_n_pts, input, arrays );

  /* Get the negative of f. */
  scal_vect_mult( total_n_pts, -1.0, arrays->f, neg_f );

  /* Calculate the potential. */
  integrate_trapezoid( arrays->r, arrays->u, neg_f, total_n_pts );

  /* Shift u so that it is zero at r[total_n_pts-1]. */
  shift_potential_bond( total_n_pts, arrays->u );

  /* For consistency, caculate the force from the potential. */
  calc_centered_diff( total_n_pts, dr, arrays->u, arrays->f );
  scal_vect_mult( total_n_pts, -1.0, arrays->f, arrays->f );

  return total_n_pts;
}


int get_angle_potential( Input *input, Arrays *arrays )
{
  int total_n_pts;
  double dr = arrays->r[1] - arrays->r[0];
  double *r_temp, *f_temp, *neg_f;

  /* Allocate memory for temporary arrays. */
  r_temp = ecalloc( input->n_pts, sizeof( double ) );
  f_temp = ecalloc( input->n_pts, sizeof( double ) );

  /* Copy r and f into temporary arrays. */
  copy_array( input->n_pts, arrays->r, r_temp );
  copy_array( input->n_pts, arrays->f, f_temp );

  /* Find the total number of points. */
  total_n_pts = get_total_n_pts_bond( input->n_pts, &(input->r_min), &(input->r_max), dr, arrays->r );

  /* Make sure that the spacing fits the requirement that r_min and r_max is 0.0 and 180.0. */
  check_angle_min_max( input->r_min, input->r_max );

  /* Allocate memory for full size arrays. */
  arrays->r     = (double  *) erealloc( arrays->r, total_n_pts*sizeof(double) );
  arrays->f     = (double  *) erealloc( arrays->f, total_n_pts*sizeof(double) );
  arrays->u     = (double  *) erealloc( arrays->u, total_n_pts*sizeof(double) );
  neg_f         = ecalloc( total_n_pts, sizeof( double ) );

  /* Fill in r. */
  fill_r_bond( total_n_pts, arrays->r, r_temp, input );

  /* Fill in f. */
  fill_f_bond( total_n_pts, arrays->f, f_temp, r_temp, dr, input );

  /* smooth f */
  smooth_forces( total_n_pts, input, arrays );

  /* Get the negative of f. */
  scal_vect_mult( total_n_pts, -1.0, arrays->f, neg_f );

  /* Calculate the potential. */
  integrate_trapezoid( arrays->r, arrays->u, neg_f, total_n_pts );

  /* Shift u so that it is zero at r[total_n_pts-1]. */
  shift_potential_bond( total_n_pts, arrays->u );

  /* For consistency, caculate the force from the potential. */
  calc_centered_diff( total_n_pts, dr, arrays->u, arrays->f );
  scal_vect_mult( total_n_pts, -1.0, arrays->f, arrays->f );

  return total_n_pts;
}


int get_dihedral_potential( Input *input, Arrays *arrays )
{
  int total_n_pts = input->n_pts;
  double dr = arrays->r[1] - arrays->r[0];
  double *neg_f;

  /* Check the endpoints. */
  check_dihedral_endpoints( input->n_pts, arrays->r );

  /* smooth f */
  smooth_forces_periodic( total_n_pts, input, arrays );

  /* Allocate memory. */
  neg_f = ecalloc( total_n_pts, sizeof( double ) );

  /* Get the negative of f. */
  scal_vect_mult( total_n_pts, -1.0, arrays->f, neg_f );

  /* Calculate the potential. */
  integrate_trapezoid( arrays->r, arrays->u, neg_f, total_n_pts );

  /* Shift u so that it is zero at r[total_n_pts-1]. */
  shift_potential_bond( total_n_pts, arrays->u );

  /* For consistency, caculate the force from the potential. */
  calc_centered_diff( total_n_pts, dr, arrays->u, arrays->f );
  scal_vect_mult( total_n_pts, -1.0, arrays->f, arrays->f );

  return total_n_pts;
}


void check_dihedral_endpoints( int n_pts, double *r )
{
  if ( fabs( r[0] + 180.0 ) > FLOAT_EPS ) 
  { print_check_dihedral_endpoints_error( r[0], -180.0, "r_min" ); exit( EXIT_FAILURE ); }
  if( fabs( r[n_pts-1] - 180.0 ) > FLOAT_EPS )
  { print_check_dihedral_endpoints_error( r[n_pts-1], 180.0, "r_max" ); exit( EXIT_FAILURE ); }
}


void check_angle_min_max( double r_min, double r_max )
{
  if ( fabs( r_min - 0.0 ) > FLOAT_EPS )
  { print_check_angle_min_max_error( r_min, 0.0, "r_min" ); exit( EXIT_FAILURE ); }

  if (fabs( r_max - 180.0 ) > FLOAT_EPS )
  { print_check_angle_min_max_error( r_max, 180.0, "r_max" ); exit( EXIT_FAILURE ); }
}


void shift_potential_bond( int total_n_pts, double *u )
{
  int i;
  double u_min = 1.0e20;
  double shift;

  /* Find the minimum. */
  for ( i=0; i<total_n_pts; i++ )
  { if ( u[i] < u_min ) { u_min = u[i]; } }

  shift = fabs( u_min );

  /* Shift the potential. */
  for ( i=0; i<total_n_pts; i++ )
  { u[i] += shift; }

}


void fill_f_bond( int total_n_pts, double *f, double *f_temp, double *r_temp, double dr, Input *input )
{
  int i, j;
  int flag_min, flag_max;
  double r_test;

  /* Get flags. */
  r_test = input->r_min;
  for ( i=0; i<total_n_pts; i++ )
  {
    if ( fabs( r_test - r_temp[0] )              < FLOAT_EPS ) { flag_min = i; }
    if ( fabs( r_test - r_temp[input->n_pts-1] ) < FLOAT_EPS ) { flag_max = i; }
    r_test += dr;
  }

  /* Copy the given force. */
  for ( i=flag_min, j=0; i<=flag_max; i++, j++ )
  { f[i] = f_temp[j]; }
    
  for ( i=flag_min-1; i>-1; i-- )
  { f[i] = f[i+1] - input->slopel * dr; }

  for ( i=flag_max+1; i<total_n_pts; i++ )
  { f[i] = f[i-1] + input->sloper * dr; }
  
}


int get_total_n_pts_bond( int n_pts, double *r_min, double *r_max, double dr, double *r )
{
  int i;
  int total_n_pts = n_pts;
  double r_test;

  /* For r_min. */
  i = 0;
  r_test = r[0];
  if ( fabs( (r_test - (*r_min) ) > FLOAT_EPS ) )
  {
    do
    {
      r_test -= dr;
      i++;
    } while ( (r_test - FLOAT_EPS) > *r_min );
    total_n_pts += i;
    *r_min = r_test; 
  }
  
  /* For r_max. */
  i = 0;
  r_test = r[n_pts-1];
  if ( fabs( *r_max - r_test ) > FLOAT_EPS )
  {
    do
    {
      r_test += dr;
      i++;
    } while ( (r_test + FLOAT_EPS) < *r_max );
    total_n_pts += i;
    *r_max = r_test;
  }

  return total_n_pts;
}


void fill_r_bond( int total_n_pts, double *r, double *r_temp, Input *input )
{
  int i = 0;
  int j = 0;
  double dr = r[1] - r[0];
  double r_test = input->r_min;

  for ( i=0; i<total_n_pts; i++ )
  {
    r[i] = r_test;
    r_test += dr;
  }
  r_test -= dr;

  if ( fabs( r_test - input->r_max ) > FLOAT_EPS ) 
  { print_fill_r_bond_error( r_test, input->r_max, total_n_pts ); exit( EXIT_FAILURE ); }

}


int get_total_n_pts( int n_pts, double *r, int *total_n_pts, double *r_min )
{
  int i = 0;
  double r_min_temp = r[0];
  double dr = r[1] - r[0];

  do
  {
    /* Quit once r_min_temp is between 0.0 and dr. */
    if ( (r_min_temp + FLOAT_EPS) < dr ) { break; }
    r_min_temp -= dr;
    i++; /* Keeps track of the number of elements to add to the arrays. */
  } while (1);

  *r_min = fabs( r_min_temp );
 
  *total_n_pts = i + n_pts;

  return i;
}


void copy_array( int n_pts, double *array, double *copy )
{
  int i;
  
  for ( i=0; i<n_pts; i++ )
  { copy[i] = array[i]; }
}

void print_table( Input input, Arrays arrays, FILE *fp_out )
{

  switch ( input.type_index )
  {
    case NB_INDEX:    print_nb_table(   input, arrays, fp_out ); break;
    case BOND_INDEX:  print_bond_table( input, arrays, fp_out ); break;
    case ANGLE_INDEX: print_bond_table( input, arrays, fp_out ); break;
    case DIHED_INDEX: print_bond_table( input, arrays, fp_out ); break;
    case INTRA_INDEX: print_nb_table(   input, arrays, fp_out ); break;
    default:  print_type_index_error( input.type_index ); exit( EXIT_FAILURE );
  }

}

void print_nb_table( Input input, Arrays arrays, FILE *fp_out )
{
  int i;
  double dr  = arrays.r[1] - arrays.r[0];
  double r_t;


  for ( i=0; i<input.n_pts; i++ )
  {
    fprintf( fp_out, "%12.6f  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n",
                      arrays.r[i], 0.0, 0.0, 0.0, 0.0, arrays.u[i], arrays.f[i] );
  }

  if ( arrays.r[input.n_pts-1] < input.r_max )
  {
    r_t = arrays.r[input.n_pts-1] + dr;
    do
    {
      fprintf( fp_out, "%12.6f  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n",
                      r_t, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      r_t += dr;
    } while ( r_t < ( input.r_max + FLOAT_EPS ) );
  }

}

void print_bond_table( Input input, Arrays arrays, FILE *fp_out )
{
  int i;

  for ( i=0; i<input.n_pts; i++ )
  {
    fprintf( fp_out, "%12.6f  %12.6e  %12.6e\n", 
             arrays.r[i], arrays.u[i], arrays.f[i] );
  }
}


/***************************************************************************************************/
/*Error Messages************************************************************************************/
/***************************************************************************************************/
void print_type_error_message( char *inter_name )
{
  printf( "\nERROR: Interaction type \"%s\" not supported.\n", inter_name );
  printf( "-----Only supports %s, %s, %s, %s, and %s.\n", NB_NAME, BOND_NAME, ANGLE_NAME, DIHED_NAME, INTRA_NAME );
}

void print_type_index_error( int type_index )
{
  printf( "\nERROR: Basis type index \"%d\" not supported.\n", type_index );  
  fprintf(stderr, "Supported basis types: %s, %s, and %s\n\n",DELTA_BASIS, LINEAR_BASIS, BSPLINE_BASIS);
}


void print_read_input_error1( int i )
{
  printf( "\nERROR: The %dth line does not contain both r and f(r). Check input.\n", i+1 );
}


void print_read_input_error2( int i )
{
  printf( "\nERROR: Grid points for r not evenly spaced in input file.\n" ); 
  printf( "-----Check spacing between elements %d and %d.\n", i+1, i+2 );
}

void print_error_fill_r( double r_test, double r_temp )
{
  printf( "\nERROR: Did not calculate the interpolation of r for nb interaction correctly.\n" );
  printf( "-----r_test: %f  r_temp: %f\n", r_test, r_temp );
}


void print_fill_r_bond_error( double r_test, double r_max, int total_n_pts )
{
  printf( "\nERROR: In fill_r_bond(), the grid spacing did not match quite right.\n" );
  printf( "-----r_test: %f  r_max: %f\n", r_test, r_max );
  printf( "-----total_n_pts: %d\n", total_n_pts );
}


void print_check_angle_min_max_error( double r_test, double r_exact, char *r_test_name )
{
  printf( "\nERROR: The spacing for the tabulated angle potential is not acceptable.\n" );
  printf( "-----%s must be %f but is %f. The current spacing does not allow this.\n", 
                 r_test_name, r_exact, r_test );
}


void print_check_dihedral_endpoints_error( double r_test, double r_exact, char *r_test_name )
{
  printf( "\nERROR: For dihedral force input, %s must be %f.\n", r_test_name, r_exact );
  printf( "-----Your input is %f.\n", r_test );
}

/************* JFR: Functions added 02.27.13 ******************************************************/
int interpolate_forces( Input *input, Arrays *arrays, Arrays *tmparrays )
{
  int i;
  int n_pts = input->interp_fact*(input->n_pts-1) + 1;
  double dr = (arrays->r[1] - arrays->r[0])/((double)input->interp_fact);
  double L_0, L_1;
  double x;
  int ctr, coeff;
  double x0, x1;

  ctr = 0;
  coeff = 0;
  for ( i=0; i<n_pts; i++ )
  {
    tmparrays->r[i] = arrays->r[0] + i*dr;
    x = tmparrays->r[i];
    x0 = arrays->r[coeff];
    x1 = arrays->r[coeff+1];
    L_0 = ( x - x1 ) / ( x0 - x1 );
    L_1 = ( x - x0 ) / ( x1 - x0 );
    tmparrays->f[i] = ( L_0 * arrays->f[coeff] ) + ( L_1 * arrays->f[coeff+1] );
    ctr++;
    if ( (ctr % input->interp_fact) == 0 ) { coeff++; }
  }

  return 0;
}

int smooth_forces( int n_pts, Input *input, Arrays *arrays )
{
  double norm;
  int i, j, Nhalf;
  double tmp[n_pts];
  int sm = input->sm;
  int index;
  double fsm;

  /* set up tmp forces */
  for ( i=0; i<n_pts; i++ ) { tmp[i] = arrays->f[i]; }

  /* do the smoothing. The last sm-1 pts on each side will remain unchanged (this is exact if you have interpolated out of the grid */
  for ( i=(sm-2); i<(n_pts-(sm-2)); i++ )
  {
    fsm = 0.00;
    for ( j = 0; j<sm; j++ ) 
    { 
      index = i + j - ((sm-1)/2);
      fsm += (1.0/((double)sm)) * tmp[index]; 
    }
    arrays->f[i] = fsm;
  }

  return 0;
}

int smooth_forces_periodic( int n_pts, Input *input, Arrays *arrays )
{
  double norm;
  int i, j, Nhalf;
  double tmp[n_pts];
  int sm = input->sm;
  int index;
  double fsm;

  /* set up tmp forces */
  for ( i=0; i<n_pts; i++ ) { tmp[i] = arrays->f[i]; }

  /* do the smoothing */
  for ( i=0; i<n_pts; i++ )
  {
    fsm = 0.00;
    for ( j = 0; j<sm; j++ ) 
    { 
      index = i + j - ((sm-1)/2);
      if ( index < 0 ) { index += (n_pts-1); }
      else if (index >= n_pts ) { index -= (n_pts-1); }
      fsm += (1.0/((double)sm)) * tmp[index]; 
    }
    arrays->f[i] = fsm;
  }

  return 0;
}

int edit_forces( Input *input, Arrays *arrays, Arrays *tmparrays )
{
  int total_n_pts;

  switch ( input->type_index )
  {
    case NB_INDEX:    total_n_pts = edit_nb_force( input, arrays, tmparrays ); break;
    case BOND_INDEX:  total_n_pts = edit_bond_force( input, arrays, tmparrays ); break;
    case ANGLE_INDEX: total_n_pts = edit_angle_force( input, arrays, tmparrays ); break;
    case DIHED_INDEX: /*total_n_pts = edit_dihedral_force( input, arrays, tmparrays );*/ break;
    case INTRA_INDEX: total_n_pts = edit_intra_nb_force( input, arrays, tmparrays ); break;
    default:  print_type_index_error( input->type_index ); exit( EXIT_FAILURE );
  }

  /* interpolate if delta or linear basis */
  interpolate_forces( input, arrays, tmparrays );

  return total_n_pts;
}


int edit_nb_force( Input *input, Arrays *arrays, Arrays *tmparrays )
{
  int total_n_pts;
  double dr = arrays->r[1] - arrays->r[0];
  double *r_temp, *f_temp, *neg_f;
  int n_endpts = ceil(0.02*input->n_pts); /* define the endpts as the last 2% of the data on each end */
  double avg_slope;
  int i_slope;
  int n_switch = floor( (0.1*arrays->r[input->n_pts-1])/dr );
  double dswitch = arrays->r[input->n_pts-1] - arrays->r[input->n_pts-1-n_switch];
  double lam_switch;
  int i;

  /* calc the avg slope for the n_endpts */
  //avg_slope = 0.00;
  //for ( i=0; i<n_endpts; i++ )
  //{
  //  avg_slope += (1.0/((double)n_endpts*dr)) * ( arrays->f[i+1] - arrays->f[i] );
  //}
  //if ( avg_slope < NB_SLOPE ) { input->slopel = avg_slope; }
  //else { input->slopel = NB_SLOPE; }

  /* find the first negative slope */
  avg_slope = 0.00;
  for ( i=0; i<input->n_pts; i++ )
  {
    avg_slope = (1.0/dr) * ( arrays->f[i+1] - arrays->f[i] );
    if ( avg_slope < 0.00 ) { i_slope = i; break; }
  } 
  if ( avg_slope < NB_SLOPE ) { input->slopel = avg_slope; }
  else { input->slopel = NB_SLOPE; }

  /* apply a switching function to the short range end of the force */
  if (i_slope != 0)
  {
    for ( i=0; i<i_slope+1; i++ )
    {
      lam_switch = (cos(3.1459*(arrays->r[i] - arrays->r[0])/(2.0*(arrays->r[i_slope] - arrays->r[0]))));
      lam_switch *= lam_switch;
      arrays->f[i] = /*lam_switch*arrays->f[i] + (1.0-lam_switch)*/( arrays->f[i_slope] + (arrays->r[i]-arrays->r[i_slope])*(input->slopel) );
    }
  }

  /* apply a switching function to the long range end of the force */
  for ( i=input->n_pts-1-n_switch; i<input->n_pts; i++ )
  {
    lam_switch = (cos(3.1459*(arrays->r[i] - arrays->r[input->n_pts-1-n_switch])/(2.0*dswitch)));
    lam_switch *= lam_switch;
    arrays->f[i] *= lam_switch;
  }

  return 0;
}

int edit_bond_force( Input *input, Arrays *arrays, Arrays *tmparrays )
{
  int total_n_pts;
  double dr = arrays->r[1] - arrays->r[0];
  double *r_temp, *f_temp, *neg_f;
  int n_endpts = ceil(0.02*input->n_pts); /* define the endpts as the last 2% of the data on each end */
  double avg_slope;
  int i;

  /* calc the avg slope for the n_endpts */
  avg_slope = 0.00;
  for ( i=0; i<n_endpts; i++ )
  {
    avg_slope += (1.0/((double)n_endpts*dr)) * ( arrays->f[i+1] - arrays->f[i] );
  }
  if ( avg_slope < BOND_SLOPE ) { input->slopel = avg_slope; }
  else { input->slopel = BOND_SLOPE; }

  avg_slope = 0.00;
  for ( i=0; i<n_endpts; i++ )
  {
    avg_slope += (1.0/((double)n_endpts*(-1.0*dr))) * ( arrays->f[input->n_pts-1-(i+1)] - arrays->f[input->n_pts-1-i] );
  }
  if ( avg_slope < BOND_SLOPE ) { input->sloper = avg_slope; }
  else { input->sloper = BOND_SLOPE; }

  return 0;
}

int edit_angle_force( Input *input, Arrays *arrays, Arrays *tmparrays )
{
  int total_n_pts;
  double dr = arrays->r[1] - arrays->r[0];
  double *r_temp, *f_temp, *neg_f;
  int n_endpts = ceil(0.02*input->n_pts); /* define the endpts as the last 2% of the data on each end */
  double avg_slope;
  int i;

  /* calc the avg slope for the n_endpts */
  avg_slope = 0.00;
  for ( i=0; i<n_endpts; i++ )
  {
    avg_slope += (1.0/((double)n_endpts*dr)) * ( arrays->f[i+1] - arrays->f[i] );
  }
  if ( avg_slope < ANGLE_SLOPE ) { input->slopel = avg_slope; }
  else { input->slopel = ANGLE_SLOPE; }

  avg_slope = 0.00;
  for ( i=0; i<n_endpts; i++ )
  {
    avg_slope += (1.0/((double)n_endpts*(-1.0*dr))) * ( arrays->f[input->n_pts-1-(i+1)] - arrays->f[input->n_pts-1-i] );
  }
  if ( (avg_slope < ANGLE_SLOPE) || (arrays->r[input->n_pts-1] > 175) ) { input->sloper = avg_slope; }
  else { input->sloper = ANGLE_SLOPE; }

  return 0;
}

int edit_intra_nb_force( Input *input, Arrays *arrays, Arrays *tmparrays )
{
  int total_n_pts;
  double dr = arrays->r[1] - arrays->r[0];
  double *r_temp, *f_temp, *neg_f;
  int n_endpts = ceil(0.02*input->n_pts); /* define the endpts as the last 2% of the data on each end */
  double avg_slope;
  int i_slope;
  int n_switch = floor( (0.1*arrays->r[input->n_pts-1])/dr );
  double dswitch = arrays->r[input->n_pts-1] - arrays->r[input->n_pts-1-n_switch];
  double lam_switch;
  int i;
  double avg_f;

  /* calc the avg slope for the n_endpts and also the avg value of f */
  //avg_slope = 0.00;
  //avg_f;
  //for ( i=0; i<n_endpts; i++ )
  //{
  //  avg_slope += (1.0/((double)n_endpts*dr)) * ( arrays->f[i+1] - arrays->f[i] );
  //  avg_f += (1.0/((double)n_endpts)) * ( arrays->f[input->n_pts-1-i] );
  //}
  //if ( avg_slope < INTRA_SLOPE ) { input->slopel = avg_slope; }
  //else { input->slopel = INTRA_SLOPE; }

  /* find the first negative slope */
  avg_slope = 0.00;
  for ( i=0; i<input->n_pts; i++ )
  {
    avg_slope = (1.0/dr) * ( arrays->f[i+1] - arrays->f[i] );
    if ( avg_slope < 0.00 ) { i_slope = i; break; }
  }
  if ( avg_slope < NB_SLOPE ) { input->slopel = avg_slope; }
  else { input->slopel = NB_SLOPE; }

  /* apply a switching function to the short range end of the force */
  if (i_slope != 0)
  {
    for ( i=0; i<i_slope+1; i++ )
    {
      lam_switch = (cos(3.1459*(arrays->r[i] - arrays->r[0])/(2.0*(arrays->r[i_slope] - arrays->r[0]))));
      lam_switch *= lam_switch;
      arrays->f[i] = /*lam_switch*arrays->f[i] + (1.0-lam_switch)*/( arrays->f[i_slope] + (arrays->r[i]-arrays->r[i_slope])*(input->slopel) );
    }
  }
//intra-cut version
//  if ( fabs(avg_f) < SWITCH_TOL ) /* Treat the interaction like a nb interaction */
//  { 
    input->sloper = 0.00;

    /* apply a switching function to the long range end of the force */
    for ( i=input->n_pts-1-n_switch; i<input->n_pts; i++ )
    {
      lam_switch = (cos(3.1459*(arrays->r[i] - arrays->r[input->n_pts-1-n_switch])/(2.0*dswitch)));
      lam_switch *= lam_switch;
      arrays->f[i] *= lam_switch;
    }
//  }
//  else /* treat the interaction like a bond */
//  {

    /* calc the avg slope for the right side n_endpts */
//    avg_slope = 0.00;
//    for ( i=0; i<n_endpts; i++ )
//    {
//      avg_slope += (1.0/((double)n_endpts*(-1.0*dr))) * ( arrays->f[input->n_pts-1-(i+1)] - arrays->f[input->n_pts-1-i] );
//    }
//    if ( avg_slope < INTRA_SLOPE ) { input->sloper = avg_slope; }
//    else { input->sloper = INTRA_SLOPE; }
//  }

  return 0;
}

