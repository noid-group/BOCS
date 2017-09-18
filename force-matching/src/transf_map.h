/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file transf_map.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef TRANSF_MAP
#define TRANSF_MAP

#include "cgff_types.h"
#include "gmx-interface.h"




int setup_CG_trr(int N_sites, tW_gmx_topology *top, tW_gmx_trxframe *fr_aa, tW_gmx_trxframe *fr_cg);
/* copies information from fr_aa into fr_cg; allocates memory for x,v,f arrays; when present
 * copies atom information from a .tpr file */

int map_config(int N_sites, tW_site_map CG_map[], tW_gmx_trxframe *fr_aa, tW_gmx_trxframe * fr_cg);

void print_gro_frame(FILE *out_gro, int N_sites, tW_site_map *CG_map, tW_gmx_trxframe *fr); 

#endif
