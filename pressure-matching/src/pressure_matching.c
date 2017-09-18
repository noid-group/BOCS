/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file pressure_matching.c
@author N. Dunn
@date Apr 24, 2011
@brief Driver for pressure matching algorithm

*/

#include <string.h>

#include "safe_mem.h"
#include "distro.h"
#include "io.h"
#include "pressure_matching_fn.h"
#include "solv_lin_eqns.h"

/*! Driver function for calculating the required pressure correction for simulating a
CG system that will reproduce a corresponding AA volume distribution. */
int main( int argc,char *argv[] )
{
	int i;
	FILE *settings;
	pres_match_sys pres_sys;
	int nargs = 1; // the only argument should be the path to the settings file

	if (argc != nargs + 1)
	{
		fprintf(stderr, "Error, %d arguments expected and %d given. Aborting.\n", nargs, argc-1);
		exit(1);
	}

	settings = safe_open_file( argv[1], "r");

	/*===================== Perform setup ======================================*/

	get_pmatch_settings(&pres_sys, settings);
	fprintf(stdout, "checkpoint 0: settings read\n");

	open_pmatch_files(&pres_sys, settings);
	fprintf(stdout, "checkpoint 1: files opened\n");

	setup_pres_tables(&pres_sys);
	fprintf(stdout, "checkpoint 2: structure memory allocated\n");

	read_input_to_sys(&pres_sys);
	fprintf(stdout, "checkpoint 3: data read in\n");

	/*====================== Perform calculations ==============================*/

	pres_sys.avg_AA_vol = get_avg( pres_sys.AA_volumes, pres_sys.N_aa_frames );
	fprintf(stdout, "avg AA volume is %f\n", pres_sys.avg_AA_vol );

	if (pres_sys.iter_flag == 0)
	{
		fill_At_vector(&pres_sys);

		calc_pres_grids(pres_sys.AA_volumes, pres_sys.At, pres_sys.N_aa_frames, pres_sys.b,\
			pres_sys.n_basis, pres_sys.basis_type, pres_sys.g_cnt, pres_sys.g, pres_sys.Q,\
			pres_sys.N_sites, pres_sys.avg_AA_vol, pres_sys.vmin, pres_sys.vmax, pres_sys.dv );


		//if tabulated basis, trim zeroes from calculation based on cutoff
		if (strcmp(DELTA, pres_sys.basis_type) == 0 )
		{
			pres_sys.n_zeroes =  trim_grids(&(pres_sys.g_cnt), &(pres_sys.zero_list),\
				&(pres_sys.b), &(pres_sys.Q), pres_sys.frac_cutoff, &(pres_sys.n_basis), pres_sys.basis_type);
		} else
		{
			pres_sys.n_zeroes = 0;
		}

		//solve linear equations based on choice of solution method
		if ( pres_sys.error_est == 0 )
		{
			solv_pres_lin_eqns( pres_sys.files->log, &pres_sys );
		} else
		{
			solv_pres_lin_eqns_error( pres_sys.files->log, &pres_sys );
		}

	} else if (pres_sys.iter_flag == 1)
	{
		//fit both aa and cg P vs V with the basis function separately, then store the difference
		//for printing
		double *psi_aa, *psi_cg; //the trimmed solutions for the aa and cg fits
		bool *zeroes_aa, *zeroes_cg;  //the zero lists for the aa and cg trajectories
		double *fit_aa, *fit_cg;  //each is a length 2 array, with fit[0] = b, fit[1] = m in  the equation p = v * m + b
		double *cnt_aa, *cnt_cg; //the trimmed g_cnt arrays for the aa and cg data
		int aa_index, cg_index; //tracks the index within a trimmed array
		double vol; //the volume corresponding to a specific grid point

		zeroes_aa = ecalloc(pres_sys.n_basis, sizeof(bool));
		zeroes_cg = ecalloc(pres_sys.n_basis, sizeof(bool));

		/***************************************************/
		// Fit the AA P-V data to the basis function

		calc_pres_grids(pres_sys.AA_volumes, pres_sys.AA_pres, pres_sys.N_aa_frames, pres_sys.b,\
			pres_sys.n_basis, pres_sys.basis_type, pres_sys.g_cnt, pres_sys.g, pres_sys.Q,\
			pres_sys.N_sites, pres_sys.avg_AA_vol, pres_sys.vmin, pres_sys.vmax, pres_sys.dv );

		//if tabulated basis, trim zeroes from calculation based on cutoff
		if (strcmp(DELTA, pres_sys.basis_type) == 0  )
		{
			pres_sys.n_zeroes = trim_grids(&(pres_sys.g_cnt), &(pres_sys.zero_list),\
				&(pres_sys.b), &(pres_sys.Q), pres_sys.frac_cutoff, &(pres_sys.n_basis), pres_sys.basis_type);
		} else
		{
			pres_sys.n_zeroes = 0;
		}

		//solve linear equations based on choice of solution method
		if ( pres_sys.error_est == 0 )
		{
			solv_pres_lin_eqns( pres_sys.files->log, &pres_sys );
		} else
		{
			solv_pres_lin_eqns_error( pres_sys.files->log, &pres_sys );
		}

		psi_aa = ecalloc(pres_sys.n_basis, sizeof(double));
		cnt_aa = ecalloc(pres_sys.n_basis, sizeof(double));

		// Copy the AA results into temporary storage
		aa_index = 0;
		for (i=0; i<(pres_sys.n_basis + pres_sys.n_zeroes); i++)
		{
			if (!pres_sys.zero_list[i])
			{
				psi_aa[aa_index] = pres_sys.psi[aa_index];
				cnt_aa[aa_index] = pres_sys.g_cnt[aa_index];
				zeroes_aa[i] = false;
				aa_index++;
			} else
			{
				zeroes_aa[i] = true;
			}
		}

		//reset array and matrix sizes for CG calculation
		restore_grids(&(pres_sys.g_cnt), &(pres_sys.b), &(pres_sys.Q), &(pres_sys.n_zeroes), &(pres_sys.n_basis));


		/***************************************************/
		// Fit the CG P-V data to the basis function	

		calc_pres_grids(pres_sys.CG_volumes, pres_sys.CG_pres, pres_sys.N_cg_frames, pres_sys.b,\
			pres_sys.n_basis, pres_sys.basis_type, pres_sys.g_cnt, pres_sys.g, pres_sys.Q,\
			pres_sys.N_sites, pres_sys.avg_AA_vol, pres_sys.vmin, pres_sys.vmax, pres_sys.dv );

		//if tabulated basis, trim zeroes from calculation based on cutoff
		if (strcmp(DELTA, pres_sys.basis_type) == 0 )
		{
			pres_sys.n_zeroes = trim_grids(&(pres_sys.g_cnt), &(pres_sys.zero_list),\
				&(pres_sys.b), &(pres_sys.Q), pres_sys.frac_cutoff, &(pres_sys.n_basis), pres_sys.basis_type);
		} else
		{
			pres_sys.n_zeroes = 0;
		}

		//solve linear equations based on choice of solution method
		if ( pres_sys.error_est == 0 )
		{
			solv_pres_lin_eqns( pres_sys.files->log, &pres_sys );
		} else
		{
			solv_pres_lin_eqns_error( pres_sys.files->log, &pres_sys );
		}


		psi_cg = ecalloc(pres_sys.n_basis, sizeof(double));
		cnt_cg = ecalloc(pres_sys.n_basis, sizeof(double));

		// Copy the CG results into temporary storage
		cg_index = 0;
		for (i=0; i<(pres_sys.n_basis + pres_sys.n_zeroes); i++)
		{
			if (!pres_sys.zero_list[i])
			{
				psi_cg[cg_index] = pres_sys.psi[cg_index];
				cnt_cg[cg_index] = pres_sys.g_cnt[cg_index];
				zeroes_cg[i] = false;
				cg_index++;
			} else
			{
				zeroes_cg[i] = true;
			}
		}

		//reset array and matrix sizes for containing the results of the calculation
		restore_grids(&(pres_sys.g_cnt), &(pres_sys.b), &(pres_sys.Q), &(pres_sys.n_zeroes), &(pres_sys.n_basis));


		/***************************************************/
		// Merge the AA and CG data based on their overlap

		fit_aa = get_best_fit(pres_sys.AA_volumes, pres_sys.AA_pres, pres_sys.N_aa_frames);
		fit_cg = get_best_fit(pres_sys.CG_volumes, pres_sys.CG_pres, pres_sys.N_cg_frames);

		aa_index = 0;
		cg_index = 0;
		for (i=0; i<pres_sys.n_basis; i++)
		{
			// For bins with direct overlap, Fv is exactly the difference of the values in
			// the bins, and the sampling frequency is the average of the AA and CG.
			// Note: for an analytic basis, only this first case applies, since all 'bins'
			// are sampled by every frame
			if (!zeroes_cg[i] && !zeroes_aa[i])
			{
				pres_sys.psi[i] = psi_aa[aa_index] - psi_cg[cg_index];
				pres_sys.g_cnt[i] = (cnt_aa[aa_index] + cnt_cg[cg_index])/2;
				aa_index++;
				cg_index++;

			// For bins without direct overlap, where either the AA or CG system is sampled
			// but the other is not, use the linear fit to the parent dataset of the missing
			// value to approximate it.  The sampling frequency is half that of the sampled
			// bin.
			} else if  (zeroes_cg[i] && !zeroes_aa[i])
			{
				vol = get_bin_val(i, pres_sys.dv, pres_sys.vmin);
				pres_sys.psi[i] = psi_aa[aa_index] - (vol * fit_cg[1] + fit_cg[0]);
				pres_sys.g_cnt[i] = cnt_aa[aa_index]/2;
				aa_index++;
			} else if (!zeroes_cg[i] && zeroes_aa[i])
			{
				vol = get_bin_val(i, pres_sys.dv, pres_sys.vmin);
				pres_sys.psi[i] =  (vol * fit_aa[1] + fit_aa[0]) - psi_cg[cg_index];
				pres_sys.g_cnt[i] = cnt_cg[cg_index]/2;
				cg_index++;
			} else
			{
				pres_sys.psi[i] = (vol * fit_aa[1] + fit_aa[0]) - (vol * fit_cg[1] + fit_cg[0]);
				pres_sys.g_cnt[i] = 1;
			}
		}

		// Print AA and CG g_cnt copies
		print_array_d(cnt_aa, aa_index, pres_sys.files->aa_g_cnt_out);
		print_array_d(cnt_cg, cg_index, pres_sys.files->cg_g_cnt_out);


		// If tabulated basis, trim zeroes from calculation based on cutoff
		// Note: this trimming procedure replaces the b array with the psi array,
		if (strcmp(DELTA, pres_sys.basis_type) == 0 )
		{
			pres_sys.n_zeroes = trim_grids(&(pres_sys.g_cnt), &(pres_sys.zero_list),\
				&(pres_sys.psi), &(pres_sys.Q), 0.0, &(pres_sys.n_basis), "delta");
		} else
		{
			pres_sys.n_zeroes = 0;
		}


		free(fit_aa);
		free(fit_cg);
		free(psi_aa);
		free(psi_cg);
		free(zeroes_aa);
		free(zeroes_cg);
		free(cnt_aa);
		free(cnt_cg);
	}

	fprintf(stdout, "psi has been solved for\n");

	/*====================== Write output ======================================*/

	print_psi_table(&pres_sys);

	print_Q(&pres_sys);

	print_array_d(pres_sys.g_cnt, pres_sys.n_basis, pres_sys.files->g_cnt_out);


	fprintf(stdout, "psi table has been printed\n");


	/*====================== Close out files and free mem ======================*/

	free_pres_match_mem(&pres_sys);

	fprintf(stdout, "pres_match_sys memory has been freed\n");

	fclose(settings);



	return 0;
}
