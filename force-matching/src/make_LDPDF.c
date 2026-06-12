/*
 *  Calculate the LDPDF 
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "safe_mem.h"
#include "io_read.h"
#include "gromacs_topology.h"
#include "LDD.h"

int main(int argc, char * argv[])
{
  /* Generic indices */
  int i, j, k;
  double dr, minr, maxr, zlo, zhi;

  /* I assume you'll be analyzing at least one trajectory. If you're looking at multiple, add them */
  tW_gmx_trxframe *fr;
  tW_gmx_topology *top;

  /* Initialize the data structures */
  top = init_tW_gmx_topology();
  fr = init_tW_gmx_trxframe();

  /* Make words for their filenames */
  tW_word top_fnm, fr_fnm;

  /* Change n_params to be however many parameters you're going to need input on the command line */
  int n_params = 9;
  int n_arg_found;
  t_pargs *params = (t_pargs *) ecalloc(n_params,sizeof(t_pargs));

  
  init_arg(&(params[0]),"-s",etSTRING,"input topology filename (req.)",TRUE);
  init_arg(&(params[1]),"-f",etSTRING,"input trajectory filename (req.)",TRUE);
  init_arg_def(&(params[2]),"-lo",etREAL,"low LD value","0.0");
  init_arg_def(&(params[3]),"-hi",etREAL,"hi LD value","100.0");
  init_arg_def(&(params[4]),"-dr",etREAL,"drho","0.1");
  init_arg_def(&(params[5]),"-zlo",etREAL,"only consider particles with z > zlo","-9999999.999999");
  init_arg_def(&(params[6]),"-zhi",etREAL,"only consider particles with z < zhi","9999999.999999");
  init_arg(&(params[7]),"-var",etBOOL,"calculate variances",FALSE);
  init_arg(&(params[8]), "-SG", etBOOL, "calculate SG distro instead of LD distro", FALSE);

  n_arg_found = get_command_line_args(argc, argv, n_params, params);

  if (print_arg_table(n_params,params,(n_arg_found > 0 ? FALSE : TRUE)))
  {
    return 0;
  }

  /* Check for mandatory args */
  check_mand_args(params,n_params);

  /* Construct filenames */
  build_filename(top_fnm, params[0].value, eTOP, ".btp");
  build_filename(fr_fnm, params[1].value, eTRAJ, ".btj");

  /* Try to read the topology file */
  if (! read_topology(top,(const char *) top_fnm))
  {
    fprintf(stderr,"ERROR: unable to read topology file: %s \n",top_fnm);
    return 1;
  } 

  read_first_frame(fr,(const char *)fr_fnm);
  
  minr = atof(params[2].value);
  maxr = atof(params[3].value);
  dr = atof(params[4].value);
  zlo = atof(params[5].value);
  zhi = atof(params[6].value);
  bool bCalcVar = params[7].bSet;
  bool bSGcalc = params[8].bSet;

  double qty_to_bin = 0; // This will hold either LDs or SGs 
  if (!fr->bLD_Grads && bSGcalc) 
  {
  fprintf(stderr, "ERROR: LDPDF can only calculate gradients when info is stored in the trajectory (lammps ldd format or btj files only)\n");
  exit(1);
  }
  /* Set Up Tallying for PDF in proper LD or SG Mode */
  tally_func_t tally; 
  if (bSGcalc) {tally = sg_tally_impl;}
  else {tally = ld_tally_impl;}


  double **averages = (double **) ecalloc(top->n_atomtypes, sizeof(double *));
  double **variances = (double **) ecalloc(top->n_atomtypes, sizeof(double *));
  double **nPoints = (double **) ecalloc(top->n_atomtypes, sizeof(double *));

  double ***LDPDFs = (double ***) ecalloc(top->n_atomtypes, sizeof(double **));
  int n_bins = (int)((maxr - minr) / dr);
  int itype, bin_idx;
  double **norms = (double **) ecalloc(top->n_atomtypes, sizeof(double *));
  for (i = 0; i < top->n_atomtypes; ++i)
  {
    averages[i] = (double *) ecalloc(top->n_atomtypes, sizeof(double));
    variances[i] = (double *) ecalloc(top->n_atomtypes, sizeof(double));
    nPoints[i] = (double *) ecalloc(top->n_atomtypes, sizeof(double));

    norms[i] = (double *) ecalloc(top->n_atomtypes, sizeof(double));
    LDPDFs[i] = (double **) ecalloc(top->n_atomtypes, sizeof(double *));
    for (j = 0; j < top->n_atomtypes; ++j)
    {
      averages[i][j] = 0.0;
      variances[i][j] = 0.0;
      nPoints[i][j] = 0.0;
      norms[i][j] = 0.0;
      LDPDFs[i][j] = (double *) ecalloc(n_bins, sizeof(double));
    }
  }

  /* loop over frames */
  while (read_next_frame(fr,TRUE))
  {
    /* Write some function to perform the analysis you need done at each frame and call it here */
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      if ((zlo <= fr->contents->x[i][2]) && (zhi >= fr->contents->x[i][2]))
      {
        itype = get_type(*(top->contents->atoms.atomtype[i]),top);
        for (j = 0; j < top->n_atomtypes; ++j)
        {
	  qty_to_bin = tally(fr, i, j);
          bin_idx = (int) ((qty_to_bin - minr) / dr);
          if (bin_idx < n_bins && bin_idx >= 0) 
          { 
            ++LDPDFs[itype][j][bin_idx]; 
            norms[itype][j] += 1.0;
            averages[itype][j] += qty_to_bin;
            nPoints[itype][j] += 1.0;
          }
	  if (bin_idx > n_bins && !bSGcalc)
          { 
            fprintf(stderr,"WARNING: LD %g is higher than maxld: %g\n",qty_to_bin,maxr); 
            fprintf(stderr,"particle idx: %d   frame: %d \n",i,fr->counter);
          }
	  else if ((bin_idx > n_bins && bSGcalc))
	  {
	    fprintf(stderr,"ERROR: SG %g is higher than maxld: %g\n",qty_to_bin,maxr);
            fprintf(stderr,"particle idx: %d   frame: %d \n",i,fr->counter);
	    exit(1);
	  }
        }
      }
    }
  }
  for (i = 0; i < top->n_atomtypes; ++i)
  {
    for (j = 0; j < top->n_atomtypes; ++j)
    {
      if (nPoints[i][j] > 0.0)
      {
        averages[i][j] /= nPoints[i][j];
      }
    }
  }

  if (bCalcVar)
  {
    fprintf(stderr, "Now calculating distribution variances\n");
    fr->counter = -1; // see init_gmx_trx_frame -1 is start
    rewind(fr->fp);
    while (read_next_frame(fr,TRUE))
    {
      for (i = 0; i < fr->contents->natoms; ++i)
      {
        if ((zlo <= fr->contents->x[i][2]) && (zhi >= fr->contents->x[i][2]))
        {
          itype = get_type(*(top->contents->atoms.atomtype[i]),top);
          for (j = 0; j < top->n_atomtypes; ++j)
          {
            qty_to_bin = tally(fr, i, j);
            bin_idx = (int) ((qty_to_bin - minr) / dr);
            if (bin_idx < n_bins && bin_idx >= 0)
            {
              variances[itype][j] += pow((qty_to_bin - averages[itype][j]),2);
            }
            else if (bin_idx >= 0)
            {
              fprintf(stderr,"WARNING: LD %g is higher than maxld: %g\n",qty_to_bin,maxr);
              fprintf(stderr,"particle idx: %d   frame: %d \n",i,fr->counter);
           }
          }
        }
      }
    }
    for (i = 0; i < top->n_atomtypes; ++i)
    {
      for (j = 0; j < top->n_atomtypes; ++j)
      {
        if (nPoints[i][j] > 1.0)
        {
          variances[i][j] /= (nPoints[i][j]-1.0);
        }
      }
    }
  }
  fclose(fr->fp);

  /* Save the analysis here */
  int *ntypes = (int *) ecalloc(top->n_atomtypes, sizeof(int));
  for (i = 0; i < fr->contents->natoms; ++i)
  {
    itype = get_type(*(top->contents->atoms.atomtype[i]),top);
    ++ntypes[itype];
  }
  
  for (i = 0; i < top->n_atomtypes; ++i)
  {
    for (j = 0; j < top->n_atomtypes; ++j)
    {
      fprintf(stderr,"\nAverage");
      if (bSGcalc) {fprintf(stderr, " SG");}
      else {fprintf(stderr, " LD");}
      fprintf(stderr," %d %d: %g \n",i,j,averages[i][j]);
      if (bCalcVar) { fprintf(stderr,"Variance: %g  StdDev: %g \n",variances[i][j],sqrt(variances[i][j])); }

      char *fnm = (char *) ecalloc(50,sizeof(char));
      if (bSGcalc) {sprintf(fnm,"SGPDF_%d_%d.dat",i,j);}
      else {sprintf(fnm,"LDPDF_%d_%d.dat",i,j);}
      FILE *fp = open_file(fnm,'w');
      
      fprintf(fp,"# Made with command: ");
      for (k = 0; k < argc; ++k)
      {
        fprintf(fp,"%s ",argv[k]);
      }
      fprintf(fp,"\n");
      fprintf(fp,"# Average: %g ",averages[i][j]);
      if (bCalcVar) { fprintf(fp," Variance: %g  StdDev: %g ",variances[i][j],sqrt(variances[i][j])); }
      fprintf(fp,"\n");
      for (k = 0; k < n_bins; ++k)
      {
//        fprintf(fp,"%f %.20f \n",minr + dr * k, (double) ((double)LDPDFs[i][j][k] / (double)(ntypes[i] * fr->counter * dr)));
        fprintf(fp,"%f %.20f \n",minr + dr * k, (double) ((double)LDPDFs[i][j][k] / (double)(norms[i][j] * dr)));
      }
      fclose(fp);
      efree(fnm);
    }
  }

  return 0;
}
