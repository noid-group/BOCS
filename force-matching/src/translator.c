/**
 * @file translator.c
 * @authors Michael DeLyser
 * @brief This is a driver for the translator exec. Generates .btp files 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cgff_types.h"
#include "safe_mem.h"
#include "io_read.h"
#include "gromacs_topology.h"


int main(int argc, char * argv[])
{

   // First some preliminary set up [This is standard to any BOCS executable driver]
  tW_word fnmi, fnmo, fnms; // Strings for filenames we'll use
  tW_gmx_trxframe *fri = init_tW_gmx_trxframe(); // Functions for allocating pointers to files 
  tW_gmx_trxframe *fro = init_tW_gmx_trxframe();
  tW_gmx_topology *top = init_tW_gmx_topology();
  
  int n_params = 3; // num params driver accomodates from user
  int n_arg_found; // num params we actually find
  int i = 0;
  t_pargs *params = (t_pargs *) ecalloc(n_params, sizeof(t_pargs)); // allocates memory for param structs
  
  init_arg(&(params[0]),"-f",etSTRING,"input filename (req.)",TRUE); // Here we specify info about the param [is it req?]
  init_arg(&(params[1]),"-o",etSTRING,"output filename (req.)",TRUE);
  init_arg(&(params[2]),"-s",etSTRING,"ref. topology filename",FALSE);

  n_arg_found = get_command_line_args(argc, argv, n_params, params); // Extracts no. args received from user to fill param, returns -9999 when -h is passed
  
  print_arg_table(n_params,params,(n_arg_found > 0 ? FALSE : TRUE));
  
  if (n_arg_found <= 0) { return 0; } // break the program if user asked help

  check_mand_args(params,n_params); // break the program if mand args not passed
  
 // Now this part of the code is specific to the driver function [translator.c]
  if (strstr(params[1].value,".btp") != NULL)
  {
    fprintf(stderr,"Found .btp file extension in output file. Assuming you're translating a dumped tpr file.\n");
    
    top->eFileType = eDUMP; // This sets the appropriate mode in read_topology
    top->b_tpr = read_topology(top,params[0].value);
    if (top->b_tpr) { write_bocs_top(params[1].value,top); }
    else
    {
      fprintf(stderr,"ERROR: there was a problem trying to read the dumped topology file %s \n",params[0].value);
      return 1;
    } 
    return 0;
  }

  build_filename(fnmi, params[0].value, eTRAJ, ".trr"); 
  build_filename(fnmo, params[1].value, eTRAJ, ".trr");
  
  open_write_trajectory(fro,fnmo);
  read_first_frame(fri,fnmi);
  if ((fro->eFileType == eGRO) || (fro->eFileType == eLMPDATA) || (fro->eFileType == eLMPTRJ))
  {
    if (! params[2].bSet)
    {
      fprintf(stderr,"ERROR: this output format requires inclusion of a corresponding topology file!\n");
      fprintf(stderr,"Please specify a .btp file with -s\n");
      return 1;
    }
    build_filename(fnms, params[2].value, eTOP, ".btp");
    top->b_tpr = read_topology(top,fnms);
    set_atom_labels(fri,top);
    if (top->contents->atoms.nr != fri->contents->natoms)
    {
      fprintf(stderr,"ERROR: topology file %s contains information on %d atoms\n",fnms,top->contents->atoms.nr);
      fprintf(stderr,"However, trajectory file %s contains information on %d atoms\n",fnmi,fri->contents->natoms);
      return 1;
    }
  }

  fro->contents = fri->contents;
  if (fri->bLDs)
  {
    fro->bLDs = TRUE;
    fro->n_atomtypes = fri->n_atomtypes;
    fprintf(stderr,"n_atomtypes: %d \n",fri->n_atomtypes);
    fro->local_densities = (double **) ecalloc(fro->contents->natoms, sizeof(double *));
    for (i = 0; i < fro->contents->natoms; ++i)
    {
      fro->local_densities[i] = &(fri->linear_lds[i*(fri->n_atomtypes)]);
    }
    if (fri->bLD_Grads)
    {
      fro->bLD_Grads = TRUE;
      fro->local_density_gradients = (dvec **) ecalloc(fro->contents->natoms, sizeof(dvec *));
      for (i = 0; i < fro->contents->natoms; ++i)
      {
        fro->local_density_gradients[i] = &(fri->linear_ldgs[i*(fri->n_atomtypes)]);
      }
    }
  }



  if (fro->eFileType == eLMPDATA)
  {
    read_next_frame(fri,TRUE);
    copy_trxframe_info(fri,fro);
    write_lammps_data(fro,top);
    return 0;
  }
  while (read_next_frame(fri,TRUE))
  {
    fro->counter = fri->counter;
    do_PBC(fri);
    write_frame(fro,top);
  }
  fprintf(stderr,"\n");
  fclose(fri->fp);
  fclose(fro->fp);

  return 0;
}
