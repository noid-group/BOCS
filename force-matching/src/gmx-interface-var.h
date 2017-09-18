/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file gmx-interface-var.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Definitions and structures needed for the gromacs interface
*/

#ifndef INTR_VAR
#define INTR_VAR

#ifndef bool
#define bool int
#endif

#ifndef DIM
#define DIM 3
typedef double dvec[DIM];
#endif

#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif


typedef struct {
    bool b_Gromacs;		// struct. stored in gmx files ?
    bool b_Forces;		// forces present ?
    bool b_PBC;			// pbc ?
    double box[DIM][DIM];	// box for pbc
    bool b_Forces_1;		// forces for any struct?
    bool b_Forces_N;		// forces for all struct?
} tW_gmx_info;


#endif
