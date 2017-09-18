/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file read_map.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef READ_MAP
#define READ_MAP

#include "cgff_types.h"

void print_site_mapping(tW_site_map map);
/*  prints out information stored in map */

void print_mapping_N(int N_sites, tW_site_map map[]);
/* prints out all mapping information. */

tW_site_map *get_CG_map(FILE * map_top, int *N_sites);
/* Reads mapping file pointed to by fp and generates and stores map into
   an array - CG_map - with N_sites elements. N_sites is assigned within
   the function. */

#endif
