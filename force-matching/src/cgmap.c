/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file cgmap.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Driver for the cgmap executable
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id$";

//c library includes
#include <math.h>

//local includes
#include "wnoid_math.h"
#include "read_map.h"
#include "transf_map.h"
#include "safe_mem.h"
#include "io_read.h"
#include "gmx-interface.h"

#define DEBUG_setup FALSE

int main(int argc, char *argv[])
{
    int N_sites;

    FILE *fp_map, *single_gro;

    tW_site_map *CG_map;

    tW_gmx_topology *CG_top;
    tW_gmx_trxframe *fr_in, *fr_out;
    tW_gmx_cgmap_input *input;
    tW_gmx_output *output;

    CG_top = init_tW_gmx_topology();

    fr_in = init_tW_gmx_trxframe();
    fr_out = init_tW_gmx_trxframe();

    input = init_tW_gmx_cgmap_input(argc, argv);

    output = init_tW_gmx_output(stderr, argv[0]);

    output->copyright(output);

    if (file_exists( input->get_filename(input, "-p") )) {
	fp_map = fopen(input->get_filename(input, "-p"), "r");	// for reading map file
    } else {
	printf("ERROR: map.top file '%s' not found.\n", input->get_filename(input, "-p"));
	exit(1);
    }

    if ( !input->arg_is_set(input, "-c") && !input->arg_is_set(input, "-o") )
    {
	fprintf(stderr, "You haven't provided a -c or -o argument, so there's nothing to do.\n");
	exit(1);
    }


    if ( input->arg_is_set(input, "-s") )
    {
        if ( file_exists( input->get_filename(input, "-s") ) )
	{
	    read_tpr(CG_top, input->get_filename(input, "-s") );

	} else {
	    printf("ERROR: tpr file '%s' not found.\n", input->get_filename(input, "-s") );
	    exit(1);

	} 
    }

    /* The first time we read data is a little special */
    read_first_trxframe(fr_in, input->get_filename(input,"-f"));

    /* read CG mapping */
    CG_map = get_CG_map(fp_map, &N_sites);

    /* setup CG fr_out */
    setup_CG_trr(N_sites, CG_top, fr_in, fr_out);


    /* map first atomistic fr_in into cg fr_out */
    map_config(N_sites, CG_map, fr_in, fr_out);	// HEREIAM

    /* If the user asks for it, produce the first frame as a gro file */
    if ( input->arg_is_set(input, "-c") )
    {
        single_gro = fopen( input->get_filename(input, "-c"),  "w");
	print_gro_frame(single_gro, N_sites, CG_map, fr_out); 
    }


    /* If the user asks for it, map the trajectory to the specified file */
    if ( input->arg_is_set(input, "-o") )
    {
        open_new_trxfile(fr_out, input->get_filename(input, "-o") );

        do {
    	    /* map atomistic fr_in into cg fr_out */
    	    map_config(N_sites, CG_map, fr_in, fr_out);	// HEREIAM 
 
            write_trxframe_to_file(fr_out);  	  

        }
          while (read_next_trxframe(fr_in));

    }

    output->thanks(output);

    return 0;
}
