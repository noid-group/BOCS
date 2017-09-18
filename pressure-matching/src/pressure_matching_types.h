/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file pressure_matching_types.h
@author N. Dunn
@date July 2, 2013
@brief Type definitions for pressure matching

This file contains type definitions for structures needed for the pressure matching calculation.  Since this version of pressure matching is decoupled from force matching, these structures are kept separate from those in wnoid_types.h.  A coupled implementation of pressure matching would likely implement these structures in that file.


 */

#ifndef P_MATCH_TYPES
#define P_MATCH_TYPES

#include <stdio.h>
#include <stdbool.h>

typedef struct
{
/**
@struct pres_match_files
@brief Files for use in the pressure matching calculation

*/
	FILE *AA_pres; ///< file containing AA pressures
	FILE *AA_vol;   ///< file containing volumes extracted from AA run
	FILE *CG_pres; ///< file containing CG pressures
	FILE *CG_vol;   ///< file containing volumes from the CG run

	FILE *ref_Fv;  ///< file containing a reference correction for use in iterative processes

	FILE *log;      ///< logfile for any interesting output during the program 
	FILE *psi_out;  ///< file that will contain the results of the calculation
	FILE *Q_out;   ///< file that will contain the correlation matrix Q (in an iterative calc, this is the AA matrix)
	FILE *g_cnt_out; ///< file that will contain the g_cnt array (the aggregate sampling in an iterative calculation)

	/* Iterative-specific files */
	FILE *Q2_out; ///< file that will contain the CG correlation matrix Q in an iterative calc
	FILE *aa_g_cnt_out; ///< file that will contain the g_cnt array for the AA trajectory in an iterative calculation
	FILE *cg_g_cnt_out; ///< file that will contain the g_cnt array for the CG trajectory in an iterative calculation

} pres_match_files;



typedef struct
{
/**
@struct pres_match_sys
@brief Contains the information related to pressure matching

*/

	/* Basic properties of the pressure matching setup */
	int N_aa_frames;  ///< the number of frames considered for pressure matching (AA set in iterative scheme)
	int N_cg_frames; ///< the number of frames in the CG set in the iterative scheme
	int N_sites;   ///< the number of CG molecules in the CG trajectory
	int N_atoms;   ///< the number of molecules in the AA trajectory
	int n_basis; ///< number of basis functions to use for pressure matching
	int n_zeroes; ///< the number of unsampled bins in a delta or linear basis
	int use_kinetic_correction; ///< 0 to ignore the kinetic correction to CG pressure, 1 to include it, 2 for explicit CG KE file
	int iter_flag; ///< 0 for standard matching, 1 for iterative matching
	int ref_Fv_flag; ///< 1 for using a reference correction, 0 otherwise
	int error_est; ///< 1 for estimating error in solution by dgesvx, 0 otherwise


	/* Settings for delta basis and binwise setup */
	double vmax; ///< The maximum volume to consider with the delta basis
	double vmin; ///< The min volume to consider with the delta basis
	double dv; ///< The grid spacing of the the delta basis
	double n_bins; ///< The number of bins to coarsen the input data into
	double frac_cutoff; ///< The minimum fraction of uniform bin population to keep

	double avg_AA_vol; ///< average volume of the atomistic trajectory
	double kB; ///< the value of Boltzmann's constant corresponding to the selected units

	/* Parameters for reference Fv */
	int ref_n_basis; ///< the number of basis functions used in the reference Fv
	double ref_vavg; ///< the reference Vavg for computing the Das Andersen ref Fv
	double ref_N_sites; ///< the reference N_sites for computing the Das Andersen ref Fv
	double * ref_psi; ///< reference pressure matching coefficients
	double ref_vmin; ///< min volume of the reference grid
	double ref_vmax; ///< max volume of the reference grid
	double ref_dv; ///< grid spacing of the reference grid

	/*  Settings stored as strings */
	char *pressure_units; ///< the units used for pressure
	char *volume_units; ///< the units used for volume
	char *energy_units; ///< the units used in the input kinetic energies
	char *basis_type; ///< the type of basis function used to represent F(V)

	/* Information gathered from AA trajectory */
	double *AA_times;   ///< a vector of the time indices for the trajectory
	double *AA_pres;  ///< a vector of atomistic pressure for each frame
	double *AA_volumes; ///< a vector of volumes for each frame, from the AA trajectory


	/* Information gathered from mapped trajectory */
	double *CG_times; ///< framewise vector of times from the CG simulation
	double *CG_volumes; ///< framewise vector of volumes from the CG simulation
	double *CG_pres; ///< framewise vector of CG pressures


	/* Calculated quantities */
	double * At;       ///< vector of difference in CG and AA pressures
	double * b;  ///< c vector for pressure matching
	double * g;	  ///< g basis vector
	double * g_cnt;   ///< counts the samples per bin
	bool * zero_list; ///< entry is true for insufficient sampling, false otherwise
	double ** Q;      ///< Q matrix for pressure matching
	double * psi;     ///< pressure matching coefficients

	/* Files for the system */
	pres_match_files *files; ///< the files needed for input/output for the calculation

} pres_match_sys;




#endif
