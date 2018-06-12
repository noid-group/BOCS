/**
@file read_map.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
*/

#ifndef READ_MAP
#define READ_MAP

#include "cgff_types.h"

void print_site_mapping(tW_site_map map);
/*  prints out information stored in map */

void print_mapping_N(int N_sites, tW_site_map map[]);
/* prints out all mapping information. */

tW_site_map *get_CG_map(FILE * map_top, int *N_sites);

tW_site_map *NEW_get_CG_map(FILE * map_top, int *N_sites);
/* this was my rewrite of it but we patched the original instead */


#endif
