/**
@file transf_map.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Functions related to mapping to a cg representation
*/

//c library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//local includes
#include "transf_map.h"
#include "safe_mem.h"


int setup_CG_trr(int N_sites, tW_gmx_topology * top, tW_gmx_trxframe *fr_aa, tW_gmx_trxframe *fr_cg)
{
    copy_trxframe_info(fr_aa, fr_cg);

    fr_cg->set_natoms(fr_cg, N_sites);

    fr_cg->set_atom_labels(fr_cg, top);

     return 0;
}

int map_config(int N_sites, tW_site_map *CG_map, tW_gmx_trxframe *fr_aa, tW_gmx_trxframe *fr_cg)
{

    copy_trxframe_info(fr_aa, fr_cg);

    map_trxframe(fr_aa, fr_cg, CG_map);

    return 0;
}


/* To alleviate the issue of having to generate a CG gro file before mapping,
 * we will print the first frame of the trajectory */
void print_gro_frame(FILE *out_gro, int N_sites, tW_site_map *CG_map, tW_gmx_trxframe *fr) 
{
    int mol, site;
    tW_word current_molname, next_molname;
    dvec x, v;
    matrix box;

    fprintf(out_gro, "first cg frame\n");
    fprintf(out_gro, "%5d\n", N_sites);

    mol = 0;
    for (site=0; site<N_sites; site++)
    {
	fr->get_pos_of(fr, site, x);
	if (fr->contents->bV) {fr->get_vel_of(fr, site, v); } // if MRD 12/14/2017

	//res# resnm atmnm atom# posx posy posz velx vely velz
        //"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
	fprintf(out_gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", 
	  CG_map[site].res_no+1, CG_map[site].mol_name, CG_map[site].cg_name, site+1, x[0], x[1], x[2], v[0], v[1], v[2] );
    }

    fr->get_box(fr, box);

    fprintf(out_gro, "%8.5f %8.5f %8.5f\n", box[0][0], box[1][1], box[2][2]);
}
