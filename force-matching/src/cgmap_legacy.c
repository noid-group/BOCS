/**
@file cgmap.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
@brief Driver for the cgmap executable
**/

/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id$";

//c library includes
#include <math.h>
#include <string.h>

//local includes
#include "wnoid_math.h"
#include "read_map.h"
#include "transf_map.h"
#include "safe_mem.h"
#include "io_read.h"
#include "gromacs_topology.h"

#define DEBUG_setup FALSE

int main(int argc, char *argv[])
{
    copyright();
    int N_sites;

    FILE *fp_map, *single_gro;

    tW_site_map *CG_map;

    tW_gmx_topology *CG_top;
    tW_gmx_trxframe *fr_in, *fr_out;
    tW_gmx_cgmap_input *input;
    CG_top = init_tW_gmx_topology();

    fr_in = init_tW_gmx_trxframe();
    fr_out = init_tW_gmx_trxframe();

    input = init_tW_gmx_cgmap_input(argc, argv);


    if (file_exists( input->get_filename(input, "-p") )) {
        fp_map = fopen(input->get_filename(input, "-p"), "r");  // for reading map file
    } else {
        printf("ERROR: map.top file '%s' not found.\n", input->get_filename(input, "-p"));
        exit(1);
    }

// eFileType must be set so the correct functions are called from read_(first/next)_frame or write_frame
    fr_in->eFileType = input->file_types[eAA].value;
    if ( input->arg_is_set(input,"-o") ) { fr_out->eFileType = input->file_types[eCG].value; }

    if (file_exists( input->get_filename(input, "-p") )) {
	fp_map = fopen(input->get_filename(input, "-p"), "r");	// for reading map file
    } else {
	printf("ERROR: map.top file '%s' not found.\n", input->get_filename(input, "-p"));
	exit(1);
    }

    if ( input->arg_is_set(input, "-s") )
    {
        if ( file_exists( input->get_filename(input, "-s") ) )
	{
// eFileType must be set so the correct function is called from read_topology
	    CG_top->eFileType = input->file_types[eCGTPR].value;
	    CG_top->b_tpr = read_topology(CG_top,input->get_filename(input,"-s"));
	} else {
	    printf("ERROR: topology file '%s' not found.\n", input->get_filename(input, "-s") );
	    exit(1);
	} 
    }

    /* The first time we read data is a little special */
    read_first_frame(fr_in,input->get_filename(input,"-f"));
    read_next_frame(fr_in,TRUE);

    /* read CG mapping */
    CG_map = get_CG_map(fp_map, &N_sites);
    /* setup CG fr_out */
    setup_CG_trr(N_sites, CG_top, fr_in, fr_out);


    /* map first atomistic fr_in into cg fr_out */
    map_config(N_sites, CG_map, fr_in, fr_out);

    /* If the user asks for it, produce the first frame as a gro file */
    if ( input->arg_is_set(input, "-c") )
    {
        single_gro = fopen( input->get_filename(input, "-c"),  "w");
	print_gro_frame(single_gro, N_sites, CG_map, fr_out); 
    }

    /* If the user asks for it, map the trajectory to the specified file */

    if (! input->arg_is_set(input,"-o"))
    {
	fprintf(stderr,"no output file set with -o. nothing to do!\n");
	return 0;
    }

    bool bB, bE, bDT;

    bB = input->arg_is_set(input,"-b");
    bE = input->arg_is_set(input,"-e");
    bDT = input->arg_is_set(input,"-dt");
    bool bThisFrame, bEND = FALSE;

    open_write_trajectory(fr_out, input->get_filename(input,"-o"));

    do 
    {
    /* map atomistic fr_in into cg fr_out */
	bThisFrame = TRUE;
	if (bB) { if (fr_in->contents->time < input->times[eBTIME].value) { bThisFrame = FALSE; }}
	if (bDT) 
	{ 
	    float test = fr_in->contents->time / input->times[eDELTAT].value; 
	    if ((int) test != test) { bThisFrame = FALSE;}
	}
	if (bE) { if (fr_in->contents->time > input->times[eETIME].value) { bThisFrame = FALSE; bEND = TRUE; }}
	if (bThisFrame) 
	{
	    map_config(N_sites, CG_map, fr_in, fr_out);	    	    
	    write_frame(fr_out);
	}   
    } while (read_next_frame(fr_in,TRUE) && (! bEND));
    fclose(fr_out->fp);

    return 0;
}
