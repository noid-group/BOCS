#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "omp.h"

#include "safe_mem.h"
#include "io_read.h"
#include "gromacs_topology.h"
#include "LDD.h"

int main(int argc, char * argv[])
{
  int i, j;
  tW_gmx_trxframe *fr_in, *fr_out;
  tW_gmx_topology *top;
  top = init_tW_gmx_topology();
  fr_in = init_tW_gmx_trxframe();
  fr_out = init_tW_gmx_trxframe();
  bool bSelf = FALSE;
  tW_word top_fnm, fr_in_fnm, fr_out_fnm, lmp_info_fnm;  
  double **hiLDs, **loLDs;

  int n_params = 5;
  int n_arg_found;
  t_pargs *params = (t_pargs *) ecalloc(n_params,sizeof(t_pargs));
  
  init_arg(&(params[0]),"-s",etSTRING,"input topology filename",TRUE);
  init_arg(&(params[1]),"-f",etSTRING,"input trajectory filename",TRUE);
  init_arg(&(params[2]),"-o",etSTRING,"output trajectory filename",TRUE);
  init_arg(&(params[3]),"-i",etSTRING,"local density info filename",TRUE);
  init_arg(&(params[4]),"-self",etBOOL,"include self-density",FALSE);
  n_arg_found = get_command_line_args(argc, argv, n_params, params);

  if (print_arg_table(n_params,params,(n_arg_found > 0 ? FALSE : TRUE)))
  {
    return 0;
  }

  /* Check for mandatory args */
  check_mand_args(params,n_params);

  /* Construct filenames */
  build_filename(top_fnm, params[0].value, eTOP, ".btp");
  build_filename(fr_in_fnm, params[1].value, eTRAJ, ".trr");
  build_filename(fr_out_fnm, params[2].value, eTRAJ, ".btj");
  build_filename(lmp_info_fnm, params[3].value, 9, ".txt");
  if (params[4].bSet) { bSelf = TRUE; }

  if (! read_topology(top,(const char *) top_fnm))
  {
    fprintf(stderr,"ERROR: unable to read topology file: %s \n",top_fnm);
    return 1;
  } 

  int n_atoms = read_first_frame(fr_in, (const char *) fr_in_fnm);
  open_write_trajectory(fr_out, (char *) fr_out_fnm);
  
  copy_trxframe_info(fr_in, fr_out); 
  set_natoms(fr_out,fr_in->contents->natoms);
  init_LD_info(fr_out,top,lmp_info_fnm);
 
  hiLDs = (double **) ecalloc(top->n_atomtypes, sizeof(double *));
  loLDs = (double **) ecalloc(top->n_atomtypes, sizeof(double *));

  for (i = 0; i < top->n_atomtypes; ++i)
  {
    hiLDs[i] = (double *) ecalloc(top->n_atomtypes, sizeof(double));
    loLDs[i] = (double *) ecalloc(top->n_atomtypes, sizeof(double));
    for (j = 0; j < top->n_atomtypes; ++j)
    {
      hiLDs[i][j] = 0.0;
      loLDs[i][j] = 9e20;
    }
  }  
 
  int n_threads;
  #pragma omp parallel
  {
    n_threads = omp_get_num_threads();
  } 


  fprintf(stderr,"bX: %d %d  bV: %d %d  bF: %d %d \n",fr_in->contents->bX,fr_out->contents->bX,
                                                      fr_in->contents->bV,fr_out->contents->bV,
                                                      fr_in->contents->bF,fr_out->contents->bF);
  while (read_next_frame(fr_in,FALSE))
  {
    ++(fr_out->counter);
    copy_matrix(fr_in->contents->box, fr_out->contents->box);
    copy_xvf(fr_in,fr_out,TRUE,TRUE,TRUE);
    calculate_local_densities(top,fr_out,bSelf,hiLDs,loLDs);
    write_frame(fr_out, top);
  }
  
  fclose(fr_in->fp);
  fclose(fr_out->fp);

  FILE *fp = open_file("hilo_LDs.dat",'w');
  fprintf(fp,"!%5s  %5s  %10s  %10s \n","itype","jtype","lo_dens","hi_dens");
  for (i = 0; i < fr_out->n_atomtypes; ++i)
  {
    for (j = 0; j < fr_out->n_atomtypes; ++j)
    {
      fprintf(fp," %5d  %5d  %10.7f  %10.7f \n",i,j,loLDs[i][j], hiLDs[i][j]);
    }
  }
  fclose(fp);

  return 0;
}
