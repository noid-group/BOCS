/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file pressure_matching_fn.h
@author N. Dunn
@date July 2, 2013
@brief Definitions for functions needed in decoupled pressure matching

*/


#ifndef PMATCH_FN
#define PMATCH_FN

#include "pressure_matching_types.h"



extern const char MD_UNITS[];
extern const char BAR[];
extern const char NM3[];
extern const char LITER[];
extern const char KJ_PER_MOL[];

extern const char DAS_ANDERSEN[];
extern const char DELTA[];
extern const char LINEAR[];

/*! Pulls settings from the config file */
void get_pmatch_settings(pres_match_sys *pres_sys, FILE *settings);

/*! Opens the files needed for pressure matching */
void open_pmatch_files(pres_match_sys *pres_sys, FILE *settings);

/*! Checks the provided files to make sure they have the same # of frames */
void check_file_lengths(pres_match_sys *pres_sys, char filetype[]);

/*! Allocates memory for pressure system variables */
void setup_pres_tables(pres_match_sys *pres_sys);

/*! Reads input files into pres_sys according to the type of calculation */
void read_input_to_sys(pres_match_sys *pres_sys);

/*! Reads input files into pres_sys according to the type of calculation */
void read_input_to_sys(pres_match_sys *pres_sys);

/*! Reads in a reference Fv for use with iterative procedures */
void read_ref_Fv(pres_match_sys *pres_sys);

/*! Calculates the At vector (Paa-Pcg) */
void fill_At_vector(pres_match_sys *pres_sys);

/*! Returns the value of the Das_Andersen pressure correction,
   given a pres_sys that has already been solved. */
double find_DA_pressure(pres_match_sys *pres_sys, int index, double v_CG);

/*! Returns the value of the Das_Andersen potential correction given a
   pres_sys that has already been solved */
double find_DA_potential(pres_match_sys *sys, int index, double v_CG);

/*! Returns the total Das_Andersen pressure correction given a set of psi. */
double find_DA_correction(double *psi, int n_basis, double vavg, int N_sites, double v_CG);

/*! Fills the b vector according to the basis type specified   */
void calc_pres_grids(double *vol, double *pres, int n_points, double *b, int n_basis, char *basis_type, double *g_cnt, \
double *g, double **Q, int N_sites, double Vavg, double vmin, double vmax, double dv );

/*! Trim b and Q based on sampling frequency.  Elements with sampling
below the threshold level of sampling are removed from both b and Q.
Returns the number of bins that were undersampled.  The arrays g_cnt
and b will be reallocated to reflect the sufficiently sampled arrays,
as will the Q matrix.  n_basis will be altered to reflect the new size
of these arrays and matrix.  The zero_list array will remain the length
of the original basis vectors. */
int trim_grids(double **g_cnt, bool **zero_list, double **b, double ***Q, double frac_cutoff, int *n_basis, char *basis_type);

/*!  Returns the arrays in the calculation to the size they were prior
 to trimming to remove undersampled bins. This is required, e.g., when
 the AA and CG fits are done separately as is done in the iterative
 method. Only the size is restored - data that was in the trimmed grids
 cannot be recovered at this point.*/
void restore_grids(double **g_cnt, double **b, double ***Q, int *n_zeroes, int *n_basis);

/*! Finds the value of reference psi for the given volume.  For tabulated volumes,
   will throw an error and halt the program if the volume is outside of the provided
   range. Interpolation between points is linear.   */ 
double get_ref_psi_value(pres_match_sys *pres_sys, double vol);

/*! Releases the memory allocated for the pressure matching system */
void free_pres_match_mem(pres_match_sys *pres_sys);


#endif
