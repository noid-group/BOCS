/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file tables.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Driver for the tables executable
*/

//c library includes
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//local includes
#include "table_types.h"
#include "safe_mem.h"
#include "read_parameters.h"
#include "table_functions.h"
#include "io_read.h"

int main ( int argc, char *argv[] )
{
  int i;
  int type_index;
  FILE *fp_in; 
  FILE *fp_out; 
  Input input;
  Arrays arrays, smarrays, tmparrays;
  int smn_pts;

  /* Check if input file exists, open input and output files*/
  if ( file_exists(argv[1]) )
  {
  	fp_in  = fopen( argv[1], "r" ); 
  } else
  {
	fprintf(stderr, "ERROR: Input file '%s' not found. Terminating.\n", argv[1]);
	exit(1);
  }
  fp_out = fopen( argv[2], "w" ); 

  /* Index to identify interaction type throughout the code. */
  input.type_index = get_inter_type( argv, argc );

  /* Index to identify tabulated basis set used for the FF calculation. */
  input.type_basis = get_basis_type( argv );

  /* Some interactions need certain parameters. */
  get_table_parameters( &input, argv );

  /* Get the number of grid points for the force. */
  input.n_pts = get_no_points( fp_in );

  /* Allocate memory for the arrays. */
  arrays.r = ecalloc( input.n_pts, sizeof( double ) );
  arrays.f = ecalloc( input.n_pts, sizeof( double ) );

  /* Read in r and f from fp_f_in. */
  read_input_forces( input.n_pts, &arrays, fp_in );

  /* JFR - 02.27.13: interpolate and smooth */
  if ( input.type_basis == DELTA_INDEX ) { input.interp_fact = DELTA_INTERP_FACT; input.sm = DELTA_SM;}
  else if ( input.type_basis == LINEAR_INDEX ) { input.interp_fact = LINEAR_INTERP_FACT; input.sm = LINEAR_SM;}
  else if ( input.type_basis == BSPLINE_INDEX ) { input.interp_fact = BSPLINE_INTERP_FACT; input.sm = BSPLINE_SM;}
  else { printf("BASIS TYPE NOT SUPPORTED. \n"); exit(0); }
  input.smn_pts = input.interp_fact*(input.n_pts - 1) + 1;

  /* Allocate memory for tmp arrays. */
  tmparrays.r = ecalloc( MAX_PTS, sizeof( double ) );
  tmparrays.f = ecalloc( MAX_PTS, sizeof( double ) );

  edit_forces( &input, &arrays, &tmparrays );

  /* free the original arrays */
  free( arrays.r );
  free( arrays.f );

  /* Allocate memory for the smoothed arrays. */
  smarrays.r = ecalloc( input.smn_pts, sizeof( double ) );
  smarrays.f = ecalloc( input.smn_pts, sizeof( double ) );
  smarrays.u = ecalloc( input.smn_pts, sizeof( double ) );

  /* place smoothed data in arrays */
  for ( i=0; i<input.smn_pts; i++ )
  {
    smarrays.r[i] = tmparrays.r[i];
    smarrays.f[i] = tmparrays.f[i];
  }
  /* free tmp memory */
  free( tmparrays.r );
  free( tmparrays.f );

  input.n_pts = input.smn_pts; /* we are dealing with the edited arrays now */

  /* Get the potential as the negative integral from force. */
  input.n_pts = get_potential( &input, &smarrays );

  /* Print GROMACS 4 table. */
  print_table( input, smarrays, fp_out );

  /* free memory */
  free( smarrays.r );
  free( smarrays.f );
  free( smarrays.u );

  return 0;
}

