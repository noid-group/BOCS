/**
@file cgmap.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
@brief Driver for the cgmap executable

MRD 09/17/2018: 
I modified the user input to match the analysis programs I wrote.
There's less command line options required now.
We infer what trajectory file types are desired for both the
input and output files based on the file extensions (and manually add a .trr
extension if none of the standard extensions are found).
Note, this new version can't do -b -e -dt like the old version.
If you only want a subset of the trajectory mapped, then trim down the atomistic
.trr file itself. Or, map the entire atomistic trajectory, and then trim the cg one.
**/

#include <math.h>
#include <string.h>

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
    int i;
    FILE *fp_map, *single_gro;

    tW_site_map *CG_map;

    tW_gmx_topology *CG_top;
    tW_gmx_trxframe *fr_in, *fr_out;
    CG_top = init_tW_gmx_topology();
    tW_word top_fnm, fri_fnm, fro_fnm;


    fr_in = init_tW_gmx_trxframe();
    fr_out = init_tW_gmx_trxframe();


    int n_params = 5;
    int n_arg_found;
    t_pargs *params = (t_pargs *) ecalloc(n_params,sizeof(t_pargs));

    init_arg(&(params[0]),"-p",etSTRING,"mapping topology filename (req.)",TRUE);
    init_arg(&(params[1]),"-f",etSTRING,"input AA trajectory filename (req.)",TRUE);
    init_arg(&(params[2]),"-o",etSTRING,"output CG trajectory filename (opt.)",FALSE);
    init_arg(&(params[3]),"-c",etSTRING,"output CG coordinate (.gro) file (opt.)",FALSE);
    init_arg(&(params[4]),"-s",etSTRING,"input CG bocs topology file (opt.)",FALSE);

    n_arg_found = get_command_line_args(argc, argv, n_params, params);

    if (print_arg_table(n_params,params,(n_arg_found > 0 ? FALSE : TRUE)))
    {
      return 0;
    }

    /* Check for mandatory args */
    check_mand_args(params,n_params);


    fp_map = open_file(params[0].value, 'r');

    build_filename(fri_fnm, params[1].value, eTRAJ, ".trr");
    if (params[2].bSet) { build_filename(fro_fnm, params[2].value, eTRAJ, ".trr"); }

    /* The first time we read data is a little special */
    read_first_frame(fr_in, (const char *) fri_fnm);
    read_next_frame(fr_in,TRUE);

    /* read CG mapping */
    CG_map = get_CG_map(fp_map, &N_sites);
    /* setup CG fr_out */
    setup_CG_trr(N_sites, CG_top, fr_in, fr_out);


    /* map first atomistic fr_in into cg fr_out */
    map_config(N_sites, CG_map, fr_in, fr_out);

    /* If the user asks for it, produce the first frame as a gro file */
    if ( params[3].bSet )
    {
        single_gro = open_file(params[3].value,'w');
        print_gro_frame(single_gro, N_sites, CG_map, fr_out);
    }

    if (! params[2].bSet)
    {
        return 0;
    }

    open_write_trajectory(fr_out, (char *)fro_fnm);
    rewind(fr_in->fp);

    if (params[4].bSet) { CG_top->b_tpr = read_topology(CG_top,params[4].value); }
    if ((fr_out->eFileType == eGRO) || (fr_out->eFileType == eLMPTRJ) || (fr_out->eFileType == eLMPDATA))
    {
      if (! params[4].bSet)
      {
        fprintf(stderr,"ERROR: unable to map entire trajectory to .gro .lmp or .data file formats without a bocs topology file (provided with -s)\n");
        return 1;
      }
    }


    while (read_next_frame(fr_in,TRUE))
    {
        map_config(N_sites, CG_map, fr_in, fr_out);
        write_frame(fr_out, CG_top);
    } ;
    fclose(fr_out->fp);
    fclose(fr_in->fp);
    return 0;
}

