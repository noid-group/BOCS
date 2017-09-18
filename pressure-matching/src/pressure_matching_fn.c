/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file pressure_matching_fn.c
@author N. Dunn
@date Apr 24, 2014
@brief Functions needed in decoupled pressure matching

*/

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "io.h"
#include "safe_mem.h"
#include "solv_lin_eqns.h"
#include "basis.h"
#include "distro.h"

#include "pressure_matching_fn.h"

const char MD_UNITS[] = "md"; ///< String used for comparison against input settings.  This one indicates pressure units of kJ/mol*nm^3
const char BAR[] = "bar"; ///< String used for comparison against input settings. This one indicates pressure units of bar.
const char NM3[] = "nm3"; ///< String used for comparison against input settings. This one indicates volume units of nm^3.
const char LITER[] = "liter"; ///< String used for comparison against input settings. This one indicates volume units of liters.
const char KJ_PER_MOL[] = "kj_per_mol"; ///< String used for comparison against input settings. This one indicates energy units of kJ/mol.

const char DAS_ANDERSEN[] = "das_andersen"; ///< String used for comparison against input settings. This one indicates an analytic Das Andersen basis type as in MSCG V.
const char DELTA[] = "delta"; ///< String used for comparison against input settings. This one indicates a delta basis representation.



/*! Pulls settings from the config file */
void get_pmatch_settings(pres_match_sys *pres_sys, FILE *settings)
{
/**
@arg @c pres_sys The main structure for the calculation - it will contain all settings
@arg @c settings The configuration file for the calculation
*/

	//get #cg_sites, #atoms and #pres_coeff from the command line args
	char pressure_units[128], volume_units[128], energy_units[128], basis_type[128];

	pres_sys->n_zeroes = 0;

	get_string_setting("Pressure_units:", settings, pressure_units);
	get_string_setting("Volume_units:", settings, volume_units);
	get_string_setting("Energy_units:", settings, energy_units);
	get_string_setting("Basis_type:", settings, basis_type);

	pres_sys->iter_flag = get_setting_i("Iterative:", settings);
	pres_sys->ref_Fv_flag = get_setting_i("Use_ref_Fv:", settings);
	pres_sys->error_est = get_setting_i("Error_estimate:", settings);

	//initialize all numeric variables to error flag variables
	pres_sys->vmax = -1;
	pres_sys->vmin = -1;
	pres_sys->dv = -1;
	pres_sys->frac_cutoff = -1;
	pres_sys->n_bins = -1;
	pres_sys->n_basis = -1;
	pres_sys->avg_AA_vol = -1;



	if ( strcmp(DAS_ANDERSEN, basis_type) == 0 )
	{
		pres_sys->basis_type = ecalloc( sizeof(DAS_ANDERSEN), sizeof(char));
		strcpy(pres_sys->basis_type, DAS_ANDERSEN);

		pres_sys->n_basis = get_setting_i("N_pres_coeff:", settings);


	} else if ( strcmp(DELTA, basis_type) == 0 )
	{
		pres_sys->basis_type = ecalloc( sizeof(DELTA), sizeof(char));
		strcpy(pres_sys->basis_type, DELTA);

		pres_sys->vmax = get_setting_d("vmax:", settings);
		pres_sys->vmin = get_setting_d("vmin:", settings);
		pres_sys->dv = get_setting_d("dv:", settings);
		pres_sys->frac_cutoff = get_setting_d("frac_cutoff:", settings);

		pres_sys->n_basis = (int) (pres_sys->vmax - pres_sys->vmin) / pres_sys->dv;
		pres_sys->n_bins = pres_sys->n_basis;

		fprintf(stdout, "Number of grids: %d\n", pres_sys->n_basis);

	}  else
	{
		fprintf(stderr, "Basis type \'%s\' not recognized.  Current options are %s, %s\n", basis_type, DAS_ANDERSEN, DELTA);
		exit(1);
	}

	pres_sys->N_atoms = get_setting_i("N_atoms:", settings);
	pres_sys->N_sites = get_setting_i("N_sites:", settings);


	//decide pressure units (kB here assumes cubic nm volumes)
	if ( strcmp(MD_UNITS, pressure_units) == 0 ) //should be its own function for checking the pressure units
	{
		pres_sys->kB = 8.314510E-03;
		pres_sys->pressure_units = ecalloc( sizeof(MD_UNITS), sizeof(char));
		strcpy(pres_sys->pressure_units, MD_UNITS);

	} else if ( strcmp(BAR, pressure_units) == 0 )
	{
		pres_sys->kB = 1.3806488E-01; //double check this value
		pres_sys->pressure_units = ecalloc( sizeof(BAR), sizeof(char));
		strcpy(pres_sys->pressure_units, BAR);
	} else
	{
		fprintf(stderr, "Pressure units \'%s\' are not recognized. Options are \'md\' and \'bar\'\n",\
		pressure_units );
		exit(1);
	}

	//decide volume units
	if ( strcmp(NM3, volume_units) == 0 )
	{
		pres_sys->volume_units = ecalloc( sizeof(NM3), sizeof(char));
		strcpy(pres_sys->volume_units, NM3);

	} else if ( strcmp(LITER, volume_units) == 0 )
	{
		pres_sys->kB = (pres_sys->kB / 1.0E+24);
		pres_sys->volume_units = ecalloc( sizeof(LITER), sizeof(char));
		strcpy(pres_sys->volume_units, LITER);

	} else
	{
		fprintf(stderr, "Volume units \'%s\' are not recognized. Options are \'nm3\' and \'liter\'\n",\
		volume_units );
		exit(1);
	}

	//decide energy units
	if ( strcmp(KJ_PER_MOL, energy_units) == 0 )
	{
		pres_sys->energy_units = ecalloc( sizeof(KJ_PER_MOL), sizeof(char));
		strcpy(pres_sys->energy_units, KJ_PER_MOL);
	} else
	{
		fprintf(stderr, "Energy units \'%s\' are not recognized. Only \'kj_per_mol\' is supported.\n",\
		volume_units );
		exit(1);
	}

	fprintf(stdout, "There are %d basis functions\n", pres_sys->n_basis);

}


/*! Opens the files needed for pressure matching */
void open_pmatch_files(pres_match_sys *pres_sys, FILE *settings)
{
/**
@arg @c pres_sys The main structure for the calculation - it will contain all settings
@arg @c settings The configuration file for the calculation

*/
	char filename[128];

	pres_sys->files = (pres_match_files *) emalloc( sizeof(pres_match_files) );

	get_string_setting("AA_pressures:", settings, filename);
	pres_sys->files->AA_pres = safe_open_file(filename, "r");

	get_string_setting("CG_pressures:", settings, filename);
	pres_sys->files->CG_pres = safe_open_file(filename, "r");

	get_string_setting("AA_volumes:", settings, filename);
	pres_sys->files->AA_vol = safe_open_file(filename, "r");

	if (pres_sys->iter_flag == 1)
	{
		/* If this isn't an iterative calculation, then it's a mapped trajectory
		   and the CG volumes are 1:1 with the AA volumes */
		get_string_setting("CG_volumes:", settings, filename);
		pres_sys->files->CG_vol = safe_open_file(filename, "r");

		get_string_setting("CG_Q_output:", settings, filename);
		pres_sys->files->Q2_out = safe_open_file(filename,"w");

		get_string_setting("AA_g_cnt_output:", settings, filename);
		pres_sys->files->aa_g_cnt_out = safe_open_file(filename,"w");

		get_string_setting("CG_g_cnt_output:", settings, filename);
		pres_sys->files->cg_g_cnt_out = safe_open_file(filename,"w");
	}

	//Check the input files to make sure the frames will line up
	check_file_lengths(pres_sys, "xvg");

	if ( pres_sys->ref_Fv_flag )
	{
		get_string_setting("ref_Fv_file:", settings, filename);
		pres_sys->files->ref_Fv = safe_open_file(filename,"r");
	}

	get_string_setting("log_file:", settings, filename);
	pres_sys->files->log = safe_open_file(filename,"w");

	get_string_setting("psi_output:", settings, filename);
	pres_sys->files->psi_out = safe_open_file(filename, "w");

	get_string_setting("Q_output:", settings, filename);
	pres_sys->files->Q_out = safe_open_file(filename,"w");

	get_string_setting("g_cnt_output:", settings, filename);
	pres_sys->files->g_cnt_out = safe_open_file(filename,"w");
}

/*! Checks the provided files to make sure they have the same # of frames */
void check_file_lengths(pres_match_sys *pres_sys, char filetype[])
{
/**
@arg @c pres_sys The main structure for the calculation
@arg @c filetype The three-letter file extension (currently only xvg is supported)

*/
	FILE **files;
	int N_files = 3;

	//Iterative calculations have an additional file: cg_volumes	
	if (pres_sys->iter_flag)
	{
		N_files++;
	}

	files = ecalloc(N_files, sizeof(FILE*));
	files[0] = pres_sys->files->AA_pres;
	files[1] = pres_sys->files->AA_vol;
	files[2] = pres_sys->files->CG_pres;

	if (pres_sys->iter_flag)
	{
		files[3] = pres_sys->files->CG_vol;
	}

	if (strcmp("xvg", filetype) == 0)
	{
		// Iterative calculations don't require AA files to have the same length as CG files
		// so long as both AA and both CG files match
		if (pres_sys->iter_flag)
		{
			pres_sys->N_aa_frames = check_xvg_lengths(files, 2);
			pres_sys->N_cg_frames = check_xvg_lengths(&(files[2]), 2);
		} else
		{
			pres_sys->N_aa_frames = check_xvg_lengths(files, 3);
			pres_sys->N_cg_frames = pres_sys->N_aa_frames;
		}
	} else
	{
		fprintf(stderr, "ERROR: Filetype %s not recognized by check_file_lengths()\n", filetype);
	}


	free(files);

	if ( (pres_sys->N_aa_frames == -1) || (pres_sys->N_cg_frames == -1) )
	{

		fprintf(stderr, "ERROR: There is a mismatch in the length of the input files\n");
		exit(1);
	}

}



/*! Allocates memory for pressure system variables */
void setup_pres_tables(pres_match_sys *pres_sys)
{
/**
@arg @c pres_sys The main structure for the calculation
*/
	int i;
	int n_frames;
	int cg_frames;
	int N_pres_coeff;

	//One-off declarations
	n_frames = pres_sys->N_aa_frames;
	cg_frames = pres_sys->N_cg_frames;
	N_pres_coeff = pres_sys->n_basis;

	/* Each of these arrays needs an entry per frame */
	pres_sys->AA_times = ecalloc(n_frames, sizeof(double));

	pres_sys->At = ecalloc(n_frames, sizeof(double));

	pres_sys->AA_volumes = ecalloc(n_frames, sizeof(double));

	pres_sys->AA_pres = ecalloc(n_frames, sizeof(double));

	pres_sys->CG_pres = ecalloc(cg_frames, sizeof(double));

	pres_sys->CG_volumes = ecalloc(cg_frames, sizeof(double));

	pres_sys->CG_times = ecalloc(cg_frames, sizeof(double));

	/* These arrays are limited to orders of N_pres_coeff */
	pres_sys->b = ecalloc(N_pres_coeff, sizeof(double));

	/* These arrays are limited to orders of N_pres_coeff */
	pres_sys->zero_list = ecalloc(N_pres_coeff, sizeof(double));

	pres_sys->g = ecalloc(N_pres_coeff, sizeof(double));

	pres_sys->g_cnt = ecalloc(N_pres_coeff, sizeof(double));

	pres_sys->psi = ecalloc(N_pres_coeff, sizeof(double));

	pres_sys->Q = ecalloc(N_pres_coeff, sizeof(double));

	for (i=0; i<N_pres_coeff; i++)
	{
		pres_sys->Q[i] = ecalloc(N_pres_coeff, sizeof(double));
	}

}


/*! Reads input files into pres_sys according to the type of calculation */
void read_input_to_sys(pres_match_sys *pres_sys)
{
/**
@arg @c pres_sys The main structure for the calculation
*/
	if (pres_sys->iter_flag)
	{
		read_xvg_to_array(pres_sys->files->AA_vol, pres_sys->AA_volumes, 1);
		read_xvg_to_array(pres_sys->files->AA_pres, pres_sys->AA_times, 0);
		read_xvg_to_array(pres_sys->files->AA_pres, pres_sys->AA_pres, 1)
;
		read_xvg_to_array(pres_sys->files->CG_pres, pres_sys->CG_pres, 1);
		read_xvg_to_array(pres_sys->files->CG_vol, pres_sys->CG_volumes, 1);
		read_xvg_to_array(pres_sys->files->CG_pres, pres_sys->CG_times, 0);

	} else
	{

		read_xvg_to_array(pres_sys->files->AA_vol, pres_sys->AA_volumes, 1);
		read_xvg_to_array(pres_sys->files->AA_pres, pres_sys->AA_times, 0);
		read_xvg_to_array(pres_sys->files->AA_pres, pres_sys->AA_pres, 1);

		read_xvg_to_array(pres_sys->files->CG_pres, pres_sys->CG_pres, 1);

	}


	if (pres_sys->ref_Fv_flag == 1)
	{
		read_ref_Fv(pres_sys);
	}


}

/* Reads in a reference Fv for use with iterative procedures */
void read_ref_Fv(pres_match_sys *pres_sys)
{
/**
@arg @c pres_sys The main structure for the calculation
*/

	int i;
	pres_match_sys ref_sys;

	pres_sys->ref_vavg = -1;
	pres_sys->ref_n_basis = -1;
	pres_sys->ref_vmin = -1;
	pres_sys->ref_vmax = -1;
	pres_sys->ref_dv = -1;

	if ( strcmp(DAS_ANDERSEN, pres_sys->basis_type) == 0 )
	{
		read_psi_table(&ref_sys, pres_sys->files->ref_Fv);

		pres_sys->ref_n_basis = ref_sys.n_basis;
		pres_sys->ref_vavg = ref_sys.avg_AA_vol;
		pres_sys->ref_N_sites = ref_sys.N_sites;

		pres_sys->ref_psi = ecalloc(pres_sys->ref_n_basis, sizeof(double));

		for (i=0; i<pres_sys->ref_n_basis; i++)
		{
			pres_sys->ref_psi[i] = ref_sys.psi[i];
		}

		free(ref_sys.psi);

	} else if ( strcmp(DELTA, pres_sys->basis_type) == 0 )
	{
		read_tabulated_basis(&ref_sys, pres_sys->files->ref_Fv);

		pres_sys->ref_n_basis = ref_sys.n_basis;

		pres_sys->ref_dv = ref_sys.dv;
		pres_sys->ref_vmin = ref_sys.vmin;
		pres_sys->ref_vmax = ref_sys.vmax;

		pres_sys->ref_psi = ecalloc(pres_sys->ref_n_basis, sizeof(double));

		for (i=0; i<pres_sys->ref_n_basis; i++)
		{
			pres_sys->ref_psi[i] = ref_sys.psi[i];
		}

		free(ref_sys.psi);
	}  else
	{
		fprintf(stderr, "Basis type \'%s\' not recognized.  Current options are %s, %s\n", pres_sys->basis_type, DAS_ANDERSEN, DELTA);
		exit(1);
	}
}





/*! Calculates the At vector (Paa-Pcg) for a 1:1 mapped
calculation. */
void fill_At_vector(pres_match_sys *pres_sys)
{
/**
@arg @c pres_sys The main structure for the calculation
*/
	int i;

	for (i=0; i<pres_sys->N_aa_frames; i++)
	{
		pres_sys->At[i] = pres_sys->AA_pres[i] - pres_sys->CG_pres[i];
	}

}


/*! Returns the value of a term in the Das_Andersen pressure correction,
   given a pres_sys that has already been solved. */
double find_DA_pressure(pres_match_sys *pres_sys, int index, double v_CG)
{
/**
@arg @c pres_sys The main data structure of the calculation, post-solution.
@arg @c index The basis function index D
@arg @c v_CG The CG volume to find the correction for
*/
	return pres_sys->psi[index-1] * get_pres_basis_force(index, pres_sys->N_sites, pres_sys->avg_AA_vol, v_CG);
}

/*! Returns the value of a term in the Das_Andersen potential correction given a
   pres_sys that has already been solved */
double find_DA_potential(pres_match_sys *pres_sys, int index, double v_CG)
{
/**
@arg @c pres_sys The main data structure of the calculation, post-solution.
@arg @c index The basis function index D
@arg @c v_CG The CG volume to find the correction for
*/
	return pres_sys->psi[index-1] * get_pres_basis_potential(index, pres_sys->N_sites, pres_sys->avg_AA_vol, v_CG);

}

/*! Returns the total Das_Andersen pressure correction given a set of psi. */
double find_DA_correction(double *psi, int n_basis, double vavg, int N_sites, double v_CG)
{
/**
@arg @c psi The set of pressure correction coefficients
@arg @c n_basis The number of terms in the Das_Andersen correction
@arg @c vavg The average volume of the corresponding AA simulation
@arg @c N_sites The number of CG sites in the system
@arg @c v_CG The volume for which a correction is needed
*/
	int i;
	double corr;

	corr = 0;

	for (i=1; i<=n_basis; i++)
	{
		corr += psi[i-1] * get_pres_basis_force(i, N_sites, vavg, v_CG);
	}


	return corr;
}

/*! Fills the b vector according to the basis type specified   */
void calc_pres_grids(double *vol, double *pres, int n_points, double *b, int n_basis, char *basis_type, double *g_cnt, \
double *g, double **Q, int N_sites, double Vavg, double vmin, double vmax, double dv )
{
/**
@arg @c vol The ordered array of sampled volumes
@arg @c pres The ordered array of sampled pressures
@arg @c n_points The number of entries in both vol and pres
@arg @c b Pointer to the vector being filled (must already be allocated)
@arg @c n_basis The number of entries in b
@arg @c basis_type The type of basis function being employed to represent the pressure
@arg @c g_cnt A vector that keeps track of the hits per bin in a binned basis
@arg @c Q The correlation matrix
@arg @c N_sites The number of CG sites in the system
@arg @c Vavg The average volume of the atomistic simulation
@arg @c vmin The minimum volume of a binned basis
@arg @c vmax The maximum volume of a binned basis
@arg @c dv The bin spacing of a binned basis
*/

	int i, j, k;
	double basis;

	for (i=0; i < n_points; i++) //loop over frames
	{

		// Zeroes g between frames
		for (j=0; j<n_basis; j++)
		{
			g[j] = 0;
		}

		// Populates, b, g, g_cnt for this frame
		for (j=0; j < n_basis; j++) //loop over basis functions
		{
			if ( strcmp(DAS_ANDERSEN, basis_type) == 0 )
			{
				basis = get_pres_basis_force(j+1, N_sites, Vavg, vol[i]);

				b[j] += (basis * pres[i]);
				g[j] += basis;

			} else if (strcmp(DELTA, basis_type) == 0)
			{
				basis = get_delta_basis(j, vmin, vmax, dv, vol[i]);

				if (basis == 1) 
				{
					b[j] += (basis * pres[i]);
					g_cnt[j] += 1; 
					g[j] += basis;
				}

			} 
		}

		// Populates Q for this frame
		for (j=0; j<n_basis; j++)
		{
			for (k=0; k<n_basis; k++)
			{
				Q[j][k] += (g[j] * g[k]);
			}
		}
	}
}


/*! Trim b and Q based on sampling frequency.  Elements with sampling
below the threshold level of sampling are removed from both b and Q.
Returns the number of bins that were undersampled.  The arrays g_cnt
and b will be reallocated to reflect the sufficiently sampled arrays,
as will the Q matrix.  n_basis will be altered to reflect the new size
of these arrays and matrix.  The zero_list array will remain the length
of the original basis vectors. */
int trim_grids(double **g_cnt, bool **zero_list, double **b, double ***Q, double frac_cutoff, int *n_basis, char *basis_type)
{
/**
@arg @c g_cnt The array of hit counts per bin.  Will be reallocated in the function.
@arg @c zero_list The array keeping track of undersampled bins
@arg @c b The array of average pressures. Will be reallocated by the function
@arg @c Q The correlation matrix.  Will be reallocated by the function.
@arg @c frac_cutoff The fraction of a uniformly sampled bin that is required as a minimum \
	amount of sampling
@arg @c n_basis The number of gridpoints.  Will be updated to reflect the number of \
	well-sampled grid points.
@arg @c basis_type String that defines the type of basis function being used
*/
	int i, j;
	int index1, index2;
	int n_zeroes = 0;
	int new_n_basis;
	double cutoff;
	double total_hits = 0;
	double *temp_b;
	double *temp_g_cnt;
	double **temp_Q;
	int first_bin = 0;
	int last_bin = (*n_basis)-1;



 
	for (i=first_bin; i<=last_bin; i++)
	{
		total_hits += (*g_cnt)[i];
	}

	// cutoff is a some frac_cutoff times the number of hits you would
	// expect if the data were uniformly distributed
	cutoff = ((total_hits) / ((double) (last_bin - first_bin + 1) ) ) * frac_cutoff;


	// The bins in the delta basis are evaluated based only on their own population
	for (i=first_bin; i<=last_bin; i++)
	{
		if ( (*g_cnt)[i] == 0 || (*g_cnt)[i] < cutoff )
		{
			(*zero_list)[i] = true;
			n_zeroes++;
		} else
		{
			(*zero_list)[i] = false;
		}
	}

	// Allocate memory for temporary storage
	new_n_basis = (*n_basis) - n_zeroes;

	temp_b = ecalloc(new_n_basis, sizeof(double));
	temp_g_cnt = ecalloc(new_n_basis, sizeof(double));
	temp_Q = ecalloc(new_n_basis, sizeof(double *));

	for (i=0; i<new_n_basis; i++)
	{
		temp_Q[i] = ecalloc(new_n_basis, sizeof(double));
	}

	// Copy the data from sampled bins into temporary storage
	index1 = 0;
	for (i=0; i<(*n_basis); i++)
	{
		if ( !(*zero_list)[i] )
		{
			temp_b[index1] = (*b)[i];
			temp_g_cnt[index1] = (*g_cnt)[i];

			index2 = 0;
			for (j=0; j<(*n_basis); j++)
			{
				if ( !(*zero_list)[j] )
				{
					temp_Q[index1][index2] = (*Q)[i][j];
					index2++;
				}
			}
			index1++;
		}
	}

	// Get rid of the old data, their pointers are then directed to the 
	// new arrays containing data from only the sampled bins
	free((*b));
	free((*g_cnt));

	for (i=0; i<(*n_basis); i++)
	{
		free((*Q)[i]);
	}
	free((*Q));

	(*b) = temp_b;
	(*g_cnt) = temp_g_cnt;
	(*Q) = temp_Q;
	(*n_basis) = new_n_basis;

	return n_zeroes;
}


/*!  Returns the arrays in the calculation to the size they were prior
 to trimming to remove undersampled bins. This is required, e.g., when
 the AA and CG fits are done separately as is done in the iterative
 method. Only the size is restored - data that was in the trimmed grids
 cannot be recovered at this point.
*/
void restore_grids(double **g_cnt, double **b, double ***Q, int *n_zeroes, int *n_basis)
{
/**
@arg @c g_cnt The array of hit counts per bin.  Will be reallocated in the function.
@arg @c b The array of average pressures. Will be reallocated by the function
@arg @c Q The correlation matrix.  Will be reallocated by the function.
@arg @c n_zeroes The number of undersampled bins previously removed from the \
	calculation.  Will be set to zero in this function.
@arg @c n_basis The number of gridpoints.  Will be updated to reflect the number of \
	well-sampled grid points.

*/
	int i;
	int total_bins;

	total_bins = (*n_zeroes) + (*n_basis);

	free(*g_cnt);
	free(*b);
	free(*Q);

	(*g_cnt) = ecalloc(total_bins, sizeof(double));
	(*b) = ecalloc(total_bins, sizeof(double));
	(*Q) = ecalloc(total_bins, sizeof(double *));

	for (i=0; i<total_bins; i++)
	{
		(*Q)[i] = ecalloc(total_bins, sizeof(double));
	}

	(*n_zeroes) = 0;
	(*n_basis) = total_bins;
}


/*! Finds the value of reference psi for the given volume.  For tabulated volumes,
   will throw an error and halt the program if the volume is outside of the provided
   range. Interpolation between points is linear.   */ 
double get_ref_psi_value(pres_match_sys *pres_sys, double vol)
{
/**
@arg @c The main data structure of the calculation
@arg @c vol The volume to which the returned psi value corresponds
*/
	int grid;
	double ref_psi;
	double dx;
	double point, next_point;
	double slope;

	if ( strcmp(DAS_ANDERSEN, pres_sys->basis_type) == 0 )
	{
		find_DA_correction(pres_sys->ref_psi, pres_sys->ref_n_basis, pres_sys->ref_vavg, pres_sys->ref_N_sites, vol);
	} else if (  (strcmp(DELTA, pres_sys->basis_type) == 0)  )
	{
		grid = get_bin_index(vol, pres_sys->ref_dv, pres_sys->ref_vmin);

		if ( (grid+1) >= pres_sys->ref_n_basis)
		{
			fprintf(stderr, "ERROR: grid point outside of tabulated reference psi range\n");
			exit(1);
		}

		dx = vol - (pres_sys->ref_vmin + grid * pres_sys->ref_dv);

		point = pres_sys->ref_psi[grid];
		next_point = pres_sys->ref_psi[grid+1];

		slope = (next_point - point) / pres_sys->ref_dv;

		ref_psi = (slope * dx) + point;

	} else
	{
		fprintf(stderr, "ERROR: Basis type not recognized in 'get_ref_psi_value'\n");
		exit(1);
	}

	return ref_psi;
}


/*! Releases the memory allocated for the pressure matching system */
void free_pres_match_mem(pres_match_sys *pres_sys)
{
/**
@arg @c pres_sys The main structure for the calculation
*/
	int i;

	/* Close files */
	fclose(pres_sys->files->AA_pres);
	fclose(pres_sys->files->CG_pres);
	fclose(pres_sys->files->AA_vol);

	if(pres_sys->iter_flag) { fclose(pres_sys->files->CG_vol); }
	if(pres_sys->ref_Fv_flag) { fclose(pres_sys->files->ref_Fv); }

	fclose(pres_sys->files->log);
	fclose(pres_sys->files->psi_out);
	fclose(pres_sys->files->Q_out);
	fclose(pres_sys->files->g_cnt_out);

	free(pres_sys->files);

	/* Free Arrays */
	free(pres_sys->pressure_units);
	free(pres_sys->volume_units);
	free(pres_sys->energy_units);

	free(pres_sys->basis_type);

	free(pres_sys->AA_times);
	free(pres_sys->AA_pres);
	free(pres_sys->AA_volumes);

	free(pres_sys->CG_pres);
	free(pres_sys->CG_volumes);
	free(pres_sys->CG_times);

	free(pres_sys->At);
	free(pres_sys->b);
	free(pres_sys->zero_list);
	free(pres_sys->g);
	free(pres_sys->g_cnt);

	free(pres_sys->psi);
	if (pres_sys->ref_Fv_flag) { free(pres_sys->ref_psi); }


	for (i=0; i<pres_sys->n_basis; i++)
	{
		free(pres_sys->Q[i]);
	}
	free(pres_sys->Q);

}


