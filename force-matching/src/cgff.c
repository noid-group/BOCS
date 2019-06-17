/**
@file cgff.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
@brief MPI driver for the cgff calculation 
*/

//c library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//openmpi library includes
#include "mpi.h"

//local includes
#include "safe_mem.h"
//#include "gmx-interface.h"
#include "gromacs_topology.h"
#include "cgff_fn.h"
#include "calc_grids.h"
#include "read_parameters.h"
#include "io_output.h"
#include "calc_ref_potential.h"

#define DEBUG_tpr              FALSE
#define DEBUG_setup_CG_struct  FALSE

#define DEBUG            FALSE
#define DEBUG_NFRAMES    100



int main (int argc, char *argv[])
{
  int i, j;
  int N_coeff;
  int N_pack;
  int N_sites;
  int n_frames;
  int n_frames_local;
  int np;
  int local_rank;
  double w_local;
  tW_word tag;                        /* Used to label output files */
  tW_files files;                /* Collection of files                */
  tW_system sys;                /* Contains stuff for total system    */
 // tW_system *sys_top;                /* Copy for single topology           *///JFR - 04.13.11: changed this to a pointer so I can choose not to make a sys cpy
 // tW_system *sys_top_global;        //JFR - 04.13.11: changed this to a pointer so I can choose not to make a sys cpy
  tW_system sys_global;
  tW_gmx_info info;                /* Flags for GROMACS.                 */
  tW_CG_site *CG_struct;        /* Stuff for each CG site             */
  tW_ref_potential ref_potential;

  tW_word tpr_filename;                // JFR - added 04.06.12: for over-writing default TPR filenames
  strcpy (tpr_filename, "");

  /* GROMACS specific variables */
  bool bGromacs = FALSE;        /* Logical variable re using GROMACS      */
  bool bF = FALSE;                /* Logical variable re presence of forces */

  tW_gmx_trxframe *fr, *fr_ref;
  //fr = init_tW_gmx_trxframe();
  //fr_ref = init_tW_gmx_trxframe();

  tW_gmx_topology *top;
  top = init_tW_gmx_topology();

  tW_gmx_info info_ref;

  /* End GROMACS specific variables */

  /* Initialize mpi processes. */
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &local_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np);

if (local_rank == 0) { copyright(); }

  /* Open files. */
  files.fp_par = fopen("par.txt", "r");
  files.fp_log = fopen("log.txt", "w");

  /* Read parameters. *//* JFR - 01.27.14: read params first */
  read_par (&files, &sys, &ref_potential);

  /* JFR - added 04.11.12:  Check to make sure that the number of structure files match the number of tpr files if explicitly specified */
  if ((sys.TPR_var.flag_TPR == TRUE) && (files.N_struct != sys.TPR_var.N_TPR))
  {
    fprintf(stderr, "ERROR: The number of .btp files specified in par.txt does not match the number of structure files.\n");
    exit (EXIT_FAILURE);
  }

  /* Allocate memory for trr output info structure *///4.5.3
  //snew (oenv, 1);
  //snew (oenv_ref, 1);

  /* set program name, command line, and default values for output options */// 4.5.3
  //output_env_init_default (oenv);
  //output_env_init_default (oenv_ref);

  /* JFR - added 04.19.12: calculate memory estimate and exit if in TEST_INP mode */
  if (sys.CalcMODE_var.CalcMODE == ITEST_INP)
  {
    estimate_memory_usage (files, &sys);
    exit (0);
  }

  /* Allocate memory for arrays. */
  N_coeff = setup_tables (&sys);
  N_pack = (N_coeff * N_coeff + N_coeff) / 2;        //JFR - 04.16.12: Number of elements in the matrix array in packed storage form

  /* Build interaction lists for each site type */
  get_all_iList(&sys);

  /* Write input summary to log file. */
  summarize_input (files, sys, ref_potential);

  /* JFR - 08.02.12: Turn on reference flag in the case that you are passing in a reference trr file */
  if (sys.REF_var.flag_reftrr == TRUE)
  {
    sys.flag_ref_potential = TRUE;
  }

  /* Make the global copy. */
  if (local_rank == 0)
  {
    initialize_sys(&sys_global);
    setup_sys_copy (sys, &sys_global, TRUE);
  }

  if (strcmp (files.mode, GMX_MODE) == 0)
  {
    /* Initialize flags for GROMACS. */
    init_gmx_info (&info);
    if (sys.REF_var.flag_reftrr == TRUE)
    {
      init_gmx_info (&info_ref);
    }

      /* Loop through CG structures. */
    for (i = 0; i < files.N_struct; i++)
    {
      /* JFR - 07.16.12: For more efficient I/O (e.g., for large trajectories), let each processor deal with a seperate file */
      if ((sys.REF_var.flag_splitfiles == TRUE) && (i % np != local_rank)) { continue; }

      /* Get the ith topology information. */
      /* JFR - added 04.06.12:  Pass variables to over-write default TPR filenames */
      if (sys.TPR_var.flag_TPR == TRUE) { strcpy (tpr_filename, sys.TPR_var.TPR_files[i]); }
      else // MRD 11.06.2017
      { 
        strcpy(tpr_filename,files.structures[i]); // Put trr filename in topology filename
        char *file_ext = strstr(tpr_filename,".trr"); // Find the file extension
        if (file_ext == NULL) // Check to make sure file extension found
        {
          fprintf(stderr,"ERROR: Did not find file extension \".trr\" in trajectory filename: %s\n",files.structures[i]);
          return 1;
        }
        strcpy(file_ext,".btp"); // Replace .trr extension with .btp
      }        
      bGromacs = read_topology(top,tpr_filename); // MRD 11.06.2017

      // Allocate memory for the number of sites in this trajectory
      CG_struct = (tW_CG_site *) ecalloc (top->get_natoms(top), sizeof (tW_CG_site));

      /* Is tpr read correctly? */
      if (DEBUG_tpr)
      {
        FILE *fp_tpr;
        fp_tpr = fopen("Print_tpr_file.txt","w");
        print_tpr_file (fp_tpr, top);
        fclose(fp_tpr);
      }        

      /* GROMACS structures? */
      if (bGromacs)
      {
        /* Setup the CG structure. */
        N_sites = setup_CG_struct (&sys, top, CG_struct, sys.Bonded_Inter_Types);

        if (DEBUG_setup_CG_struct)
        {
//        print_CG_struct (files.fp_log, N_sites, CG_struct, sys); 
          print_Bond_Types (files.fp_log, sys); //NJD 5-24-16 - changed from sys_top to sys
          return 0;
        }

        if (sys.CalcMODE_var.CalcMODE != ISECOND_HALF)
        {                /* JFR - added 04.06.12: If in SECOND_HALF MODE, skip the trr loop */
          /* Open trr file for current topology. */
          //open_trr_file (files.structures[i], oenv, &status, &info, &fr);        // 4.5.3
          fr = init_tW_gmx_trxframe();
          fr_ref = init_tW_gmx_trxframe();
                      
          // read_first_trxframe(fr, files.structures[i]); MRD 11.09.17
          if (read_first_frame(fr, files.structures[i]) == -1)
          {
            fprintf(stderr,"There was an error reading file: %s\n",files.structures[i]);
            return 1;
          }
          else { read_next_frame(fr,FALSE); }
          // Originally, read_first_trxframe read the header for the first frame, rewound
          // the file, and then got the box matrix and the xvf vectors for the first frame
          // In the new read_first_frame function, only the header of the first frame
          // gets read. To actually get the box matrix and the xvf vectors, you have to call
          // read_next_frame.

          if (sys.REF_var.flag_reftrr == TRUE)
          {
            //open_trr_file (sys.REF_var.reftrr_fnm[i], oenv_ref, &status_ref, &info_ref, &fr_ref);
            //read_first_trxframe(fr_ref, sys.REF_var.reftrr_fnm[i]); MRD
            //if (read_first_trxframe(fr_ref,sys.REF_var.reftrr_fnm[i]) == -1)
            if (read_first_frame(fr_ref,sys.REF_var.reftrr_fnm[i]) == -1)
            {
              fprintf(stderr,"There was an error reading file: %s\n",sys.REF_var.reftrr_fnm[i]);
              return 1;
            }
            else { read_next_frame(fr_ref,FALSE); }
          }

          // MRD
          fr->set_atom_labels(fr,top);
          // end MRD

          /* Copy trr information to CG_struct. */
          bF = copy_trr_2_CGstruct (fr, CG_struct);
          update_info_trr(&info, fr);
                      
          if (sys.REF_var.flag_reftrr == TRUE)
          {
            copy_trr_2_CGstruct_ref (fr_ref, CG_struct);
          }

          /* Allocate memory if skipping triple loop */
          if (sys.SKIP_TRIPLE_LOOP)
          {
            sys.linear_half_matrix = (dvec *) ecalloc( N_sites * sys.N_coeff, sizeof(dvec) );
            sys.half_matrix = (dvec **) ecalloc(N_sites, sizeof(dvec *));
            sys.bm_half_mat = (bitMask *) ecalloc(GET_N_SPOTS(N_sites * sys.N_coeff), sizeof(bitMask));
            for (j = 0; j < N_sites; ++j)
            {
              sys.half_matrix[i] = &(sys.linear_half_matrix[j * sys.N_coeff]);
            }
          }

          /* Set frame counter to zero for each topology. */
          n_frames = 0;
          n_frames_local = 0;

          /* NJD - set normalization weight to zero for each topology */
          if (sys.FRAMEWEIGHT_var.flag_FRAMEWEIGHT == TRUE)
          {
            sys.wt_norm = 0;
          }

          // set up iLists for all site types here

          /* Loop over all frames */
          do
          {

            if (DEBUG && (n_frames >= DEBUG_NFRAMES))
            {
              fprintf (files.fp_log, "Exiting loop at n_frames: %d.\n", n_frames);
              break;
            }

            n_frames++;
            sys_global.REG_var.Nframes++;

            /* JFR - 07.16.12: If the splitfiles option is off, split the frames between processes, opening one file at a time */
            if ((sys.REF_var.flag_splitfiles == TRUE)|| (n_frames % np == local_rank))
            {
              n_frames_local++;

              if ((sys.flag_ref_potential == TRUE) && (sys.REF_var.flag_reftrr == FALSE))
              {
                get_ref_forces (files.fp_log, N_sites, CG_struct, info, top, ref_potential);
              }
              if (sys.SKIP_TRIPLE_LOOP) { calc_grids3(files.fp_log, info, N_sites, CG_struct, &sys); }
              else { calc_grids2 (files.fp_log, info, N_sites, CG_struct, &sys); }

            }

            if (sys.REF_var.flag_reftrr == TRUE)
            {
              read_trr_2_CGstruct_ref (&info_ref, fr_ref, CG_struct);
            }
          } while (read_trr_2_CGstruct (&info, fr, CG_struct));        // 4.5.3
          fprintf (files.fp_log, "Read %d frames for %s.\n\n", n_frames, files.structures[i]);

          /* NJD - Here, we copy back the G_wt etc., to G etc.,
             while normalizing by the weight for this processor/
             trajcectory. */
          if (sys.FRAMEWEIGHT_var.flag_FRAMEWEIGHT == TRUE)
          {
            sys.wt_norm  = sys.wt_norm / n_frames_local;

            printf ("Avg wt_factor is %lg, n_frames is %d\n", sys.wt_norm, n_frames_local);

            for (j=0; j<sys.N_coeff; j++)
            {
              sys.b[j] = sys.b_wt[j] / sys.wt_norm;
              sys.b_wt[j] = 0.0;
            }
              
            for (j=0; j<sys.N_pack; j++)
            {
              sys.M[j] = sys.M_wt[j] / sys.wt_norm;
              sys.M2[j] = sys.M2_wt[j] / sys.wt_norm;
              sys.M_wt[j] = 0.0;
              sys.M2_wt[j] = 0.0;
            }
          }

          /* Normalize arrays for ith topology. */
          normalize_arrays_top (top->get_natoms(top), n_frames_local, &sys);

          /* Evaluate w_local. */
          w_local = ((double) n_frames_local) / ((double) n_frames);

          /* Weight the sums from each local proc */
          weight_local_top (&sys, w_local, N_coeff);

          // If we have more than one topology, we must collect all information for this topology from all
          // processors, weight it, then add the weighted info to sys_global
          if ( (files.N_struct > 1) && (sys.REF_var.flag_splitfiles == FALSE) ) 
          {
//                        fprintf(stderr, "About to reduce arrays for first topology\n");

            // Do an in-place MPI reduce onto the arrays of the mother processor
            if (local_rank == 0)
            {
              MPI_Reduce (MPI_IN_PLACE, sys.b, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (MPI_IN_PLACE, sys.b_ref, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (MPI_IN_PLACE, sys.g, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (MPI_IN_PLACE, sys.L, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (MPI_IN_PLACE, sys.g_cnt, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (MPI_IN_PLACE, sys.M, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (MPI_IN_PLACE, sys.M2, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (sys.M_cnt != NULL)
              {
                MPI_Reduce (MPI_IN_PLACE, sys.M_cnt, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              }
              /* JFR - 06.16.12: Chi2 */
              MPI_Reduce (MPI_IN_PLACE, &sys.Chi2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (MPI_IN_PLACE, sys.d2b, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (MPI_IN_PLACE, sys.d2M, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            } else {
              MPI_Reduce (sys.b, sys.b, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (sys.b_ref, sys.b_ref, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (sys.g, sys.g, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (sys.L, sys.L, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (sys.g_cnt, sys.g_cnt, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (sys.M, sys.M, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (sys.M2, sys.M2, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (sys.M_cnt != NULL)
              {
                MPI_Reduce (sys.M_cnt, sys.M_cnt, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              }
              /* JFR - 06.16.12: Chi2 */
              MPI_Reduce (&sys.Chi2, &sys.Chi2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (sys.d2b, sys.d2b, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce (sys.d2M, sys.d2M, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
//                        fprintf(stderr, "Just reduced arrays for first topology\n");

            if  (local_rank == 0)
            {
              // Add weighted contents of sys to sys_global
              update_total_arrays(sys, files.p_struct[i], &sys_global);
            }
          } else {        /* We only have one topology (or one topology for this processor), so we can MPI_Reduce right to sys_global*/
            /* Weight the arrays by the topology weighting */
            weight_local_top (&sys, files.p_struct[i], N_coeff); 

            /* Combine weighted arrays for current top straight into sys_global on root */
            MPI_Reduce (sys.b, sys_global.b, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (sys.b_ref, sys_global.b_ref, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (sys.g, sys_global.g, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (sys.L, sys_global.L, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (sys.g_cnt, sys_global.g_cnt, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (sys.M, sys_global.M, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (sys.M2, sys_global.M2, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (sys.M_cnt != NULL)
            {
              MPI_Reduce (sys.M_cnt, sys_global.M_cnt, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            MPI_Reduce (&sys.Chi2, &sys_global.Chi2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (sys.d2b, sys_global.d2b, N_coeff, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce (sys.d2M, sys_global.d2M, N_pack, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
          }

          if (files.N_struct > 1)
          {
            // Print results for the current topology if the user requested we do so
            if (sys.MEM_var.flag_mult_top == TRUE)
            {
              if (local_rank == 0)
              {
                get_top_tag (files.structures[i], tag);
                if (sys.SKIP_TRIPLE_LOOP) { process_G_matrix(&sys); }
                get_results (files.fp_log, &sys, info.b_Forces_N, tag); 
                print_output (info.b_Forces_N, sys, tag);

                if (sys.REF_var.flag_calcbref == TRUE)
                {
                  //get_top_tag (files.structures[i], tag);
                  print_bref (sys.N_coeff, sys.b_ref, tag);
                }

                // Reset sys properties for the next topology
                //free_sys_copy(&sys); 
                //initialize_sys(&sys);
                setup_sys_copy(sys_global, &sys, TRUE); // reset all properties to initial values
              }
            }
          }
          efree(CG_struct);
          efree(fr);
          efree(fr_ref);    

          /* Free memory if skipping triple loop */
          if (sys.SKIP_TRIPLE_LOOP)
          {
            efree(sys.half_matrix);
            efree(sys.linear_half_matrix);
            efree(sys.bm_half_mat);
          }

        }                /* End if CalcMODE != SECOND_HALF */
        else { printf("CalcMODE = %s, skipping loop over .trr files \n", SECOND_HALF);}
      } /* End if ( bGromacs ). */
      else
      {
        printf ("ERROR: Trouble reading GROMACS topology.\n");
        exit (EXIT_FAILURE);
      }

      /* Reset flags for next topology. */
      reset_gmx_info (&info);
      if (sys.REF_var.flag_reftrr == TRUE)
      {
        reset_gmx_info (&info_ref);        /*close_trj(status_ref); */
      }
    }                        /* End the loop over topologies. */
  }                                /* End GROMACS LOOP. */
  else
  {
    printf ("PDB mode only available in serial \n");
    exit (0);
  }


  /* Get final results */
  if (local_rank == 0)
  {
    if (sys.REF_var.flag_calcbref == FALSE)
    {                        /* JFR - 07.16.12: do the normal stuff */
      /* Get final results */
      if (sys.CalcMODE_var.CalcMODE == ISECOND_HALF)
      {                        /* JFR - added 04.06.12: If in SECOND_HALF MODE, read in save state */
        read_save_state (files.fp_log, &sys_global);        /* JFR - added 04.06.12: If in SECOND_HALF MODE, read in the save state */
      }
      else
      {
        if (sys_global.SKIP_TRIPLE_LOOP) { process_G_matrix(&sys_global); }
      }
      if (sys_global.MT_var.flag_print == TRUE)
      {
        print_M_matrix2 (&sys_global);
      }                        // JFR - added 04.06.12: check if the user wants to print the matrix
      get_results (files.fp_log, &sys_global, info.b_Forces_N, "total");
      print_output (info.b_Forces_N, sys_global, "total");
    }
    else
    { /* JFR - 07.13.12: The triple loop was skipped, so just print out a file with bref */
      print_bref (sys_global.N_coeff, sys_global.b_ref, "total");
    }
  }

  printf ("\nGoodbye from proc. %d.\n", local_rank);

  MPI_Finalize ();

  return 0;
}
