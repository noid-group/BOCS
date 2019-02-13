/**
@file calc_grids.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
@brief Main driver of the cgff calculation 
*/

//c library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//gromacs includes
//#include "physics.h" // for definition of BOLTZ - used only here and only once

//local includes
#include "safe_mem.h"
#include "cgff_types.h"
#include "wnoid_math.h"
#include "calc_grids.h"
#include "solv_lin_eqns.h"
#include "io_read.h"
#include "io_output.h"


/* DEBUG calc_grids. */
#define DEBUG_eval_bond_basis FALSE
#define DEBUG_calc_grids      FALSE

/*****************************************************************************************
setup_box_dimensions(): Stores the box vectors for a cubic box from the diagonal
elements of matrix. Returns False if the box is not cubic.
*****************************************************************************************/
bool setup_box_dimensions(dvec box, double matrix[DIM][DIM])
{
    int i, j;

    for (i = 0; i < DIM; i++) {
	for (j = 0; j < DIM; j++) {
	    if (i == j) {
		continue;
	    }
	    if (matrix[i][j] > FLOAT_EPS) {
		printf("\nERROR: Only rectangular boxes supported.\n");
		printf("i: %d   j: %d  matrix: %f  \t FLOAT_EPS: %f \n",
		       i, j, matrix[i][j], FLOAT_EPS);
		return FALSE;
	    }
	}
    }

    for (i = 0; i < DIM; i++) {
	box[i] = matrix[i][i];
    }

    return TRUE;
}

/*****************************************************************************************
get_all_iList():
*****************************************************************************************/
void get_all_iList(tW_system *sys)
{
	int i, j;
	int k;
	int N_i;
	int b_Match1, b_Match2;

	sys->Inter_Map = calloc(sys->N_Site_Types, sizeof(tW_word *));
	sys->Inter_Map_Len = calloc(sys->N_Site_Types, sizeof(int));
	sys->Inter_iMap = calloc(sys->N_Site_Types, sizeof(int *));

	for (i=0; i<sys->N_Site_Types; i++)
	{
		fprintf(stderr, "Site type %s: ", sys->Site_Types[i]);

		N_i = 0;
		for (j=0; j<sys->N_Inter2_Types; j++)
		{
			b_Match1 = strcmp(sys->Site_Types[i], sys->Inter2_Type_List[j].name1);
			b_Match2 = strcmp(sys->Site_Types[i], sys->Inter2_Type_List[j].name2);

			if (b_Match1 == 0 || b_Match2 == 0)
			{
				N_i++;
			}
		}

		sys->Inter_Map_Len[i] = N_i;

		sys->Inter_Map[i] = calloc(N_i, sizeof(tW_word));
		sys->Inter_iMap[i] = calloc(N_i, sizeof(int));

		k = 0;
		for (j=0; j<sys->N_Inter2_Types; j++)
		{
			fprintf(stderr, "Checking %s-%s\n", sys->Inter2_Type_List[j].name1, sys->Inter2_Type_List[j].name2);
			b_Match1 = strcmp(sys->Site_Types[i], sys->Inter2_Type_List[j].name1);
			b_Match2 = strcmp(sys->Site_Types[i], sys->Inter2_Type_List[j].name2);

			if (b_Match1 == 0) 
			{
				strcpy(sys->Inter_Map[i][k], sys->Inter2_Type_List[j].name2);
				sys->Inter_iMap[i][k] = j;

				fprintf(stderr, "%s %d.%d ", sys->Inter_Map[i][k], j, k);

				k++;
			} else if (b_Match2 == 0 )
			{
				strcpy(sys->Inter_Map[i][k], sys->Inter2_Type_List[j].name1);
				sys->Inter_iMap[i][k] = j;

				fprintf(stderr, "%s %d.%d ", sys->Inter_Map[i][k], j, k);

				k++;
			}
		}
		fprintf(stderr, "\n");
	}
}



/*****************************************************************************************
get_iList():
*****************************************************************************************/
int get_iList(tW_word name, int N_Inter_Types,
	      tW_type_inter2 InterList[], int *i1_list, tW_word * n1_list)
{
    int i, j;
    int i1_inter = 0;
    int b_Match1, b_Match2;

    /* Loop over all pair interaction types. */
    for (i = 0; i < N_Inter_Types; i++) {
	/* Does the site "name" appear in the ith interaction? */
	b_Match1 = strcmp(name, InterList[i].name1);
	b_Match2 = strcmp(name, InterList[i].name2);

	/* If so, copy interaction index and second site name in lists. */
	if (b_Match1 == 0) {
	    i1_list[i1_inter] = i;
	    strcpy(n1_list[i1_inter], InterList[i].name2);
	    i1_inter++;
	} else if (b_Match2 == 0) {
	    i1_list[i1_inter] = i;
	    strcpy(n1_list[i1_inter], InterList[i].name1);
	    i1_inter++;
	}
    }

    return i1_inter;
}


/*****************************************************************************************
det_min_image():
*****************************************************************************************/
int det_min_image(dvec box, dvec Origin, dvec x)
/*
 * routine determines image of prtcl at R closest to the
 * particle at the Origin.  routine returns the no. of transformations
 * NB the function ASSUMES rectangular box, s.t. box_j = L_j
 */
{
    int j;
    double r_j, d_j, L_j, sign;
    int result = 0;

    for (j = 0; j < DIM; j++) {
	r_j = x[j] - Origin[j];	/* Signed distance.         */
	d_j = fabs((double) r_j);	/* Unsgnd distance.         */
	L_j = box[j];		/* ASSUMES rectangular box. */

	if (d_j > 0.5 * L_j) {
	    sign = d_j / r_j;	/* Determines direction of shift. */
	    x[j] = x[j] - sign * L_j;	/* Calculates new coord.          */
	    result++;
	}
    }

    return result;
}


/*****************************************************************************************
normalize_arrays_top(): Normalizes the arrays.
*****************************************************************************************/
int normalize_arrays_top(int N_sites, int N_frames, tW_system * sys_top)
{
    int i;
    int N_coeff = sys_top->N_coeff;
    int N_pack = sys_top->N_pack;
    double norm = (double) (3 * N_sites * N_frames);

    for (i = 0; i < N_coeff; i++) {
	sys_top->b[i] /= norm;
	sys_top->b_ref[i] /= norm;
	sys_top->g[i] /= norm;
	sys_top->L[i] /= norm;

	sys_top->g_cnt[i] /= N_frames;

	sys_top->d2b[i] /= norm;

    }

    // JFR - added 04.11.12: put the matrix in packed form
    for (i = 0; i < N_pack; i++) {
	sys_top->M[i] /= norm;
	sys_top->M2[i] /= norm;
	sys_top->d2M[i] /= norm;

    }

    /* JFR - 06.27.12: Chi2 */
    sys_top->Chi2 /= norm;

    return 0;
}


/*****************************************************************************************
update_N_instances_nb(): Currently, this function is not used.
*****************************************************************************************/
int update_N_instances_nb(tW_system sys_top_global, tW_system * sys_global)
{
    int i;

    for (i = 0; i < sys_top_global.N_Inter2_Types; i++) {
	sys_global->Inter2_Type_List[i].N_inter +=
	    sys_top_global.Inter2_Type_List[i].N_inter;
    }

    return 0;
}


/*****************************************************************************************
update_total_arrays(): Update the master grids with the data from the current topology.
*****************************************************************************************/
int update_total_arrays(tW_system sys_top, double pr_top, tW_system * sys)
{
    int i;
    int N_coeff = sys_top.N_coeff;
    int N_pack = sys_top.N_pack;

    for (i = 0; i < N_coeff; i++) {
	sys->b[i] += pr_top * sys_top.b[i];
	sys->b_ref[i] += pr_top * sys_top.b_ref[i];
	sys->g[i] += pr_top * sys_top.g[i];
	sys->g_cnt[i] += pr_top * sys_top.g_cnt[i];
	sys->L[i] += pr_top * sys_top.L[i];
	sys->d2b[i] += pr_top * sys_top.d2b[i];
	///*if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) {*/ sys->d2b[i] += pr_top * sys_top.d2b[i]; /*}*/
    }

    // JFR - added 04.11.12: put the matrix in packed form
    for (i = 0; i < N_pack; i++) {
	sys->M[i] += pr_top * sys_top.M[i];
	sys->M2[i] += pr_top * sys_top.M2[i];
	if (sys->M_cnt != NULL) {
	    sys->M_cnt[i] += pr_top * sys_top.M_cnt[i];
	}			//JFR - added 04.13.12: check if in LOWMEM mode
	sys->d2M[i] += pr_top * sys_top.d2M[i];
	///*if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 ) {*/ sys->d2M[i] += pr_top * sys_top.d2M[i]; /*}*/
    }

    /* JFR - 06.27.12: Chi2 */
    sys->Chi2 += pr_top * sys_top.Chi2;

    return 0;
}


/*****************************************************************************************
get_results(): (1) Calculate and store b_forces and b_struct. (2) Remove rows of 
zero from the matrix M. (3) Calculate the coeff. phi.
*****************************************************************************************/
int get_results(FILE * fp_log, tW_system * sys, int forces, tW_word tag)
{
    FILE *fp;
    int dummy, test_sscanf, i, j, k, l;
    tW_line inp_line;
    double bref;
/* MRD 02.13.2019 Only allocate memory for this later if needed */
    double *M2; 
	//double M2[sys->N_pack];
    int index;

    if (sys->CalcMODE_var.CalcMODE != ISECOND_HALF) {
	get_b_struct_and_b_forces(fp_log, sys, tag);
	calc_d2b(sys);
	calc_d2M(sys);
	///*if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) {*/ calc_d2b( sys ); /*}*/
	///*if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 ) {*/ calc_d2M( sys ); /*}*/
	get_rescale(sys);
	print_save_state(fp_log, sys);	/* JFR - added 04.06.12: If not SECOND_HALF MODE, print the save state */
    } else if (sys->REF_var.flag_readbref == TRUE) {	/* JFR - 02.26.13: read in and subtract off precomputed bref in SECONDHALF */
	if (!file_exists("bref.dat")) {
	    printf
		("ERROR: Unable to open the file 'bref.dat' for reading.\n");
	    exit(EXIT_FAILURE);
	}
	fp = fopen("bref.dat", "r");
	for (i = 0; i < sys->N_coeff; i++) {
	    dummy = get_next_line(fp, inp_line);
	    test_sscanf = sscanf(inp_line, "%lf", &bref);
	    sys->b_struct[i] -= bref;
	    sys->b_forces[i] -= bref;
	}
	fclose(fp);
    }

    /* JFR - 12.03.13 */
    if (sys->ITER_var.flag_AAM2 == FALSE) {
	/* put the normal matrix back together */
	for (i = 0; i < sys->N_pack; i++) {
	    sys->M[i] += sys->M2[i];
	}
    } else {
	fp = fopen(sys->ITER_var.AAM2_fnm, "r");
	M2 = (double *) ecalloc(sys->N_pack, sizeof(double));
	for (i = 0; i < sys->N_pack; i++) {
	    dummy = get_next_line(fp, inp_line);
	    test_sscanf = sscanf(inp_line, "%lf", &M2[i]);
	}
	fclose(fp);

	for (j = 0; j < sys->N_Inter_Types; j++) {	/* loop over interaction types */
	    for (l = 0; l < sys->Inter_Types[j].N_coeff; l++) {	/* loop over the coefficients in order */
		i = sys->Inter_Types[j].i_0 + l;
		for (k = 0; k < sys->N_coeff; k++) {	/* second loop over coefficients (rows) */
		    index = index_Lpacked(k, i, sys->N_coeff);
		    if ((strcmp(sys->Inter_Types[j].inter_type, NB_PAIR) ==
			 0)
			||
			(strcmp
			 (sys->Inter_Types[j].inter_type,
			  B_NB_PAIR_BOND) == 0)
			||
			(strcmp(sys->Inter_Types[j].inter_type, B_DIHEDRAL)
			 == 0)) {
			sys->M[i] += sys->M2[i];	/* These are CG 2-body correlations */
		    } else {
			sys->M2[index] = M2[index];	/* These are AA 2-body correlations */
		    }
		}
	    }
	}

	/* put the modified matrix together */
	for (i = 0; i < sys->N_pack; i++) {
	    sys->M[i] += sys->M2[i];
	}
	efree(M2);
    }


    trim_Mb(fp_log, sys);

    if (sys->CalcMODE_var.CalcMODE != IFIRST_HALF) {
	get_phi_struct_and_phi_forces(fp_log, sys, forces);

	if (sys->MFD_var.flag_MFD == TRUE)	// JFR - 04.06.12: added option to shut this off
	{
	    get_b_soln(fp_log, sys);	//JFR-04.04.11 <- this is when I originally implemented this
	}

	if (sys->ITER_var.flag_bsolnerr == TRUE)	// JFR - 12.13.13
	{
	    get_b_soln_err(fp_log, sys);	// older version without diff from AA
	    //get_b_soln_errAA( fp_log, sys ); // the diff from AA is not what I really wanted but may be useful later
	}
    }

    return 0;
}


/*****************************************************************************************
get_b_struct_and_b_forces(): (1) Allocates memory for sys->b_forces and sys->b_struct. 
(2) Copies sys->b to sys->b_forces. (3) Calculates b_struct and stores it in 
sys->b_struct. (4) If present, subtracts sys->b_ref from sys->b_struct. (5) Clears 
sys->b.
*****************************************************************************************/
int get_b_struct_and_b_forces(FILE * fp_log, tW_system * sys, tW_word tag)
{
    int i;
    int N_coeff = sys->N_coeff;
    double temperature = sys->Temperature;

    /* JFR - 07.16.12: variables for reading in bref */
    FILE *fp;
    int dummy, test_sscanf;
    tW_line inp_line;
    double bref;

    fprintf(fp_log, "\nIn get_b_struct_and_b_forces.\n");

    sys->b_forces = (double *) ecalloc(N_coeff, sizeof(double));
    sys->b_struct = (double *) ecalloc(N_coeff, sizeof(double));

    for (i = 0; i < N_coeff; i++) {
	sys->b_forces[i] = sys->b[i];
    }

    for (i = 0; i < N_coeff; i++) {
	sys->b[i] = 0.0;
    }

    calc_b_struct(sys, temperature, tag);
    for (i = 0; i < N_coeff; i++) {
	sys->b_struct[i] = sys->b[i];
    }

    if (sys->REF_var.flag_readbref == TRUE) {	/* JFR - 07.16.12: read in and subtract off precomputed bref */
	fp = fopen("bref.dat", "r");
	for (i = 0; i < N_coeff; i++) {
	    dummy = get_next_line(fp, inp_line);
	    test_sscanf = sscanf(inp_line, "%lf", &bref);
	    sys->b_struct[i] -= bref;
	    sys->b_forces[i] -= bref;
	}
	fclose(fp);
    } else if (sys->flag_ref_potential) {	/* JFR - 07.16.12: else do the normal stuff */
	for (i = 0; i < N_coeff; i++) {
	    sys->b_struct[i] -= sys->b_ref[i];
	    // JFR - added 06.27.12: use b ref
	    sys->b_forces[i] -= sys->b_ref[i];
	}
    }

    for (i = 0; i < N_coeff; i++) {
	sys->b[i] = 0.0;
    }

    free(sys->b_ref);		//JFR - added 04.13.12: this is the last time we need b_ref

    return 0;
}


/*****************************************************************************************
calc_dg_dz(): Calculates the derivative of g and stores it into sm_dg_dz. Prints g, g_cnt
dg_dz, and sm_dg_dz. Special treatment for dihedral interaction to make sure the -180 and
180 grid points are the same and to account for the periodic nature of this coordinate.
Currently, the last grid point for dihedrals are cleared since it duplicates the first
grid point. Smoothing of dg_dz is only performed if set in par.txt. 
*****************************************************************************************/
void calc_dg_dz(tW_Inter_Types * inter, double *dg_dz, double *sm_dg_dz,
		tW_word tag)
{
    int i;
    double r = inter->R_min;
    double dr = inter->dr;
    FILE *fp;
    tW_word fname;


    if (strcmp(inter->inter_type, B_DIHEDRAL) == 0) {
	inter->ptr_g[inter->N_pts - 1] = inter->ptr_g[0];
    }

    calc_centered_diff(inter->N_pts, inter->dr, inter->ptr_g, dg_dz);

    if (strcmp(inter->inter_type, B_DIHEDRAL) == 0) {
	dg_dz[0] =
	    (1.0 / (2.0 * inter->dr)) * (inter->ptr_g[1] -
					 inter->ptr_g[inter->N_pts - 2]);
	dg_dz[inter->N_pts - 1] = dg_dz[0];
    }

    if (inter->n_smooth > 0) {
	if (strcmp(inter->inter_type, B_DIHEDRAL) == 0) {
	    calc_running_avg_periodic(inter->N_pts, inter->n_smooth, dg_dz,
				      sm_dg_dz);
	} else {
	    calc_running_avg(inter->N_pts, inter->n_smooth, dg_dz,
			     sm_dg_dz);
	}
    } else {
	for (i = 0; i < inter->N_pts; i++) {
	    sm_dg_dz[i] = dg_dz[i];
	}
    }

    if (strcmp(inter->inter_type, B_DIHEDRAL) == 0) {
	inter->ptr_g[inter->N_pts - 1] = 0.0;
	dg_dz[inter->N_pts - 1] = 0.0;
	sm_dg_dz[inter->N_pts - 1] = 0.0;
    }

    sprintf(fname, "dg_dz.%s.%s.%s.dat", tag, inter->inter_type,
	    inter->inter_name);
    fp = fopen(fname, "w");
    for (i = 0; i < inter->N_pts; i++) {
	fprintf(fp, "%f  %f  %f  %f  %f\n", r, inter->ptr_g[i],
		inter->ptr_g_cnt[i], dg_dz[i], sm_dg_dz[i]);
	r += dr;
    }
    fclose(fp);

}


/*****************************************************************************************
calc_b_struct(): Calculate b_struct from g and L.
*****************************************************************************************/
int calc_b_struct(tW_system * sys, double temperature, tW_word tag)
{
    int i, j, k;
    double kT = BOLTZ * temperature;
    double *dg_dz;
    double *sm_dg_dz;
    tW_Inter_Types *inter;

     /*JFR*/ double sign;

    for (i = 0; i < sys->N_Inter_Types; i++) {
	inter = &(sys->Inter_Types[i]);

	if (inter->i_basis == DELTA_BASIS_INDEX) {
	    dg_dz = (double *) ecalloc(inter->N_pts, sizeof(double));
	    sm_dg_dz = (double *) ecalloc(inter->N_pts, sizeof(double));

	    calc_dg_dz(inter, dg_dz, sm_dg_dz, tag);

	    for (j = 0; j < inter->N_pts; j++) {
		inter->ptr_b[j] = kT * (sm_dg_dz[j] - inter->ptr_L[j]);
	    }

	    free(dg_dz);
	    free(sm_dg_dz);
	}
	/* START JFR */
	else if (inter->i_basis == LINEAR_BASIS_INDEX) {

	    dg_dz = (double *) ecalloc(inter->N_pts, sizeof(double));
	    sm_dg_dz = (double *) ecalloc(inter->N_pts, sizeof(double));

	    if (LINEAR_STRUCT == 0) {	/* calculate b same as delta basis */
		/* Calculate dg_dz. */
		calc_dg_dz(inter, dg_dz, sm_dg_dz, tag);
		sign = 1.0;
	    } else {		/* calc b by taking the derivative of the linear basis analytically */

		sign = -1.0;
		if (inter->n_smooth > 0) {
		    if (strcmp(inter->inter_type, B_DIHEDRAL) == 0) {
			calc_running_avg_periodic(inter->N_pts,
						  inter->n_smooth,
						  inter->ptr_g, sm_dg_dz);
		    } else {
			calc_running_avg(inter->N_pts, inter->n_smooth,
					 inter->ptr_g, sm_dg_dz);
		    }
		} else {
		    for (k = 0; k < inter->N_pts; k++) {
			sm_dg_dz[k] = inter->ptr_g[k];
		    }
		}
	    }

	    /* Calculation b_struct. */
	    for (j = 0; j < inter->N_pts; j++) {
		inter->ptr_b[j] =
		    kT * (sign * sm_dg_dz[j] - inter->ptr_L[j]);
	    }

	    free(dg_dz);
	    free(sm_dg_dz);

	}
	/* End linear basis function. */
	/* END JFR */
	else if (inter->i_basis == BSPLINE_BASIS_INDEX) {	/* JFR - 07.22.12 */

	    dg_dz = (double *) ecalloc(inter->N_pts, sizeof(double));
	    sm_dg_dz = (double *) ecalloc(inter->N_pts, sizeof(double));

	    if (BSPLINE_STRUCT == 0) {	/* calculate b same as delta basis */
		/* Calculate dg_dz. */
		calc_dg_dz(inter, dg_dz, sm_dg_dz, tag);
		sign = 1.0;
	    } else {		/* calc b by taking the derivative of the linear basis analytically */

		sign = -1.0;
		if (inter->n_smooth > 0) {
		    if (strcmp(inter->inter_type, B_DIHEDRAL) == 0) {
			calc_running_avg_periodic(inter->N_pts,
						  inter->n_smooth,
						  inter->ptr_g, sm_dg_dz);
		    } else {
			calc_running_avg(inter->N_pts, inter->n_smooth,
					 inter->ptr_g, sm_dg_dz);
		    }
		} else {
		    for (k = 0; k < inter->N_pts; k++) {
			sm_dg_dz[k] = inter->ptr_g[k];
		    }
		}
	    }

	    /* Calculation b_struct. */
	    for (j = 0; j < inter->N_pts; j++) {
		inter->ptr_b[j] =
		    kT * (sign * sm_dg_dz[j] - inter->ptr_L[j]);
	    }

	    free(dg_dz);
	    free(sm_dg_dz);

	} /* End Bspline basis function. */
	else if (inter->i_basis == HARMONIC_BASIS_INDEX) {
	    get_b_harmonic_bond(inter, kT);
	}

	else if (inter->i_basis == RYCKAERT_BELLEMANS_BASIS_INDEX) {
	    get_b_rb_dihedral(inter, kT);
	}

	else if (inter->i_basis == TOY_DIHED_INDEX) {
	    get_b_TOY_dihedral(inter, kT);
	}

	else if (inter->i_basis == POWER_INDEX) {
	    get_b_power(inter, kT);
	}

	else {
	    printf
		("\nERROR: Bond basis function not supported in calc_b_struct_Bond_Inter().\n");
	    exit(EXIT_FAILURE);
	}

    }

    return 0;
}


/*****************************************************************************************
trim_Mb(): Eliminates rows from sys->M that contain only zeros (i.e. the basis functions
were not sampled in the configurations. (1) Compile the Zero_list by finding which basis
functions were not sampled during the calculation. (2) 
*****************************************************************************************/
int trim_Mb(FILE * fp, tW_system * sys)
{
    int N_coeff = sys->N_coeff;
    int N_coeff0 = sys->N_coeff;
    int N_zero;
    int Zero_list[N_coeff0];

    fprintf(fp, "\n\nIn trim_Mb.\n");

    initialize_Zero_list(N_coeff0, Zero_list);

    /* Determine Zero_list, N_zero, and N_coeff. */
    N_zero = get_Zero_list(Zero_list, &N_coeff, sys->g_cnt, /*JFR*/ sys);

    remove_zero_rows_mem(fp, N_coeff0, N_zero, Zero_list, sys);

    resize_arrays_trim_Mb(sys, N_coeff);

    set_pointers_trim_Mb(sys);

    fprintf(fp, "\nExiting trim_Mb.\n");

    return N_coeff;
}


/*****************************************************************************************
initialize_Zero_list(): Initialize elements of Zero_list to -1.
*****************************************************************************************/
void initialize_Zero_list(int N_coeff0, int *Zero_list)
{
    int i;

    for (i = 0; i < N_coeff0; i++) {
	Zero_list[i] = -1;
    }
}


/*****************************************************************************************
get_Zero_list(): Determine Zero_list, N_zero, and N_coeff. (1) Zero_list[] contains the
indices to elements of g_cnt that are zero. (2) N_zero is the number of indices stored 
in Zero_list. (3) N_coeff is the number of coeff. to be determined after removing those
corresponding to the Zero_list.
*****************************************************************************************/
int get_Zero_list(int *Zero_list, int *N_coeff, double *g_cnt,
		  /*JFR*/ tW_system * sys)
{
    int i;
    int N_coeff0 = *N_coeff;
    int N_coeff_trim = *N_coeff;
    int N_zero = 0;

    double g_cnt_total;
    double Nhits_uniform;

    int end_flag = FALSE;
    int j, l;

    /* JFR - 02.06.13:  I have rewritten this entire function */
    for (j = 0; j < sys->N_Inter_Types; j++) {	/* loop over interaction types */
	g_cnt_total = 0.00;
	for (l = 0; l < sys->Inter_Types[j].N_coeff; l++) {	/* calculate the total number of hits for each interaction */
	    i = sys->Inter_Types[j].i_0 + l;
	    g_cnt_total += g_cnt[i];
	}
	/* calculate the number of hits per bin if the hits were uniformly distributed */
	Nhits_uniform =
	    ((sys->Inter_Types[j].dr * g_cnt_total) /
	     (sys->Inter_Types[j].R_max - sys->Inter_Types[j].R_min));

	for (l = 0; l < sys->Inter_Types[j].N_coeff; l++) {	/* loop over the coefficients in order */
	    i = sys->Inter_Types[j].i_0 + l;
	    /* trim bins that are sampled l.t. a specified fraction (sys->TRIM_var.FE) of the uniform value */
	    if (g_cnt[i] < (sys->TRIM_var.FE * Nhits_uniform)) {
		end_flag = TRUE;
	    }
	    /* do the same for the sampling of the input atomistic simulation (for iter-FMing) */
	    else if (sys->rescale[i] < (sys->TRIM_var.FE * Nhits_uniform)) {
		end_flag = TRUE;
	    }
	    /* for the linear basis, trim the first and last bins no matter what for proper supports */
	    else if (sys->Inter_Types[j].i_basis == LINEAR_BASIS_INDEX) {
		if (strcmp(sys->Inter_Types[j].inter_type, B_DIHEDRAL) ==
		    0) {
		    continue;
		}
		if ((l == 0) || (l == (sys->Inter_Types[j].N_coeff - 1))) {
		    end_flag = TRUE;
		}
	    }
	    /* for the Bspline basis, trim the first and last k/2+1 bins no matter what for proper supports */
	    else if (sys->Inter_Types[j].i_basis == BSPLINE_BASIS_INDEX) {
		if (l == (sys->Inter_Types[j].N_coeff - 1)) {
		    if (strcmp(sys->Inter_Types[j].inter_type, B_DIHEDRAL)
			== 0) {
			if (fabs
			    ((sys->Inter_Types[j].R_min +
			      l * sys->Inter_Types[j].dr) - M_PI) <
			    FLOAT_EPS) {
			    end_flag = TRUE;
			}
		    }
		}
		//if( (l<=sys->Inter_Types[j].kspline/2) || (l >= sys->Inter_Types[j].N_coeff-1-sys->Inter_Types[j].kspline/2 ) ) { end_flag = TRUE; }
	    }

	    /* update the trimming variables */
	    if (end_flag == TRUE) {
		Zero_list[N_zero] = i;	/* Record index of zero element of g_cnt.                       */
		N_zero++;	/* Keep track of the number of zero elements.                   */
		N_coeff_trim--;	/* For each zero element, there's one less coeff. to determine. */
	    }
	    end_flag = FALSE;
	}
    }

    *N_coeff = N_coeff_trim;

    return N_zero;
}


///*****************************************************************************************
//setup_arrays_trim_Mb(): Allocate memory for temporary storage arrays. 
//*****************************************************************************************/
//void setup_arrays_trim_Mb(tW_system * sys, int N_coeff)
//{
//    int N_pack = (N_coeff * N_coeff + N_coeff) / 2;
//
//    sys->x = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->b = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->b_forces = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->b_struct = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->g = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->g_cnt = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->L = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->phi = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->phi_forces = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->phi_struct = (double *) ecalloc(N_coeff, sizeof(double));
//    // JFR - added 04.11.12: put the matrix in packed form
//    sys->M = (double *) ecalloc(N_pack, sizeof(double));
//    sys->M2 = (double *) ecalloc(N_pack, sizeof(double));
//    if (sys->M_cnt != NULL) {
//	sys->M_cnt = (double *) ecalloc(N_pack, sizeof(double));
//    }
//    sys->d2b = (double *) ecalloc(N_coeff, sizeof(double));
//    sys->d2M = (double *) ecalloc(N_pack, sizeof(double));
//
//    sys->rescale = (double *) ecalloc(N_coeff, sizeof(double));
//}


///*****************************************************************************************
//free_arrays_trim_Mb(): Free the memory for arrays in sys.
//*****************************************************************************************/
//void free_arrays_trim_Mb(tW_system * sys, int N_coeff)
//{
//
//    free(sys->b);
//    free(sys->b_forces);
//    free(sys->b_struct);
//    free(sys->g);
//    free(sys->g_cnt);
//    free(sys->L);
//    // JFR - added 04.11.12: put the matrix in packed form
//    free(sys->M);
//    free(sys->M2);
//    if (sys->M_cnt != NULL) {
//	free(sys->M_cnt);
//    }
//    free(sys->d2b);
//    free(sys->d2M);
//    ///*if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) {*/ free( sys->d2b  ); /*}*/
//    ///*if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 ) {*/ free( sys->d2M  ); /*}*/
//    free(sys->rescale);
//
//    free(sys->x);
//    free(sys->phi);
//    free(sys->phi_forces);
//    free(sys->phi_struct);
//
//}


///*****************************************************************************************
//remove_zero_rows(): Removes rows of zero from matrix M and the corresponding elements in
//the other arrays by copying only the elements that correspond to the non-zero g_cnt into
//temporary storage.
//*****************************************************************************************/
//int remove_zero_rows(FILE * fp, int N_coeff0, int N_zero, int *Zero_list,
//		     tW_system * sys, tW_system * sys_temp)
//{
//    int i_inter;
//    int i0, j0, i0_inter;
//    int i1, j1, i1_inter;
//    int i_zero, j_zero;
//    double r, dr;
//    tW_Inter_Types *inter_ptr;
//
//    int index0, index1;
//    int N0 = sys->N_coeff;
//    int N1 = sys->N_coeff - N_zero;
//
//    i0 = 0;
//    i1 = 0;
//    i_zero = 0;
//    for (i_inter = 0; i_inter < sys->N_Inter_Types; i_inter++) {
//	inter_ptr = &(sys->Inter_Types[i_inter]);
//
//	inter_ptr->i_0 = i1;
//
//	dr = inter_ptr->dr;
//	r = inter_ptr->R_min - dr;
//
//	i1_inter = 0;
//	for (i0_inter = 0; i0_inter < inter_ptr->N_coeff; i0_inter++) {
//	    r += dr;
//	    fprintf(fp,
//		    "i0: %6d    i0_inter: %6d    sys->g_cnt[%6d]: %15.6lf    ",
//		    i0, i0_inter, i0, sys->g_cnt[i0]);
//
//	    if (i0 != Zero_list[i_zero]) {
//		fprintf(fp, "i1: %6d    i1_inter: %6d\n", i1, i1_inter);
//
//		sys_temp->x[i1] = r;
//		sys_temp->b[i1] = sys->b[i0];
//		sys_temp->b_forces[i1] = sys->b_forces[i0];
//		sys_temp->b_struct[i1] = sys->b_struct[i0];
//		sys_temp->g[i1] = sys->g[i0];
//		sys_temp->g_cnt[i1] = sys->g_cnt[i0];
//		sys_temp->L[i1] = sys->L[i0];
//
//
//		sys_temp->d2b[i1] = sys->d2b[i0];
//		///*if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) {*/ sys_temp->d2b[i1] = sys->d2b[i0]; /*}*/
//		sys_temp->rescale[i1] = sys->rescale[i0];
//
//		j1 = 0;
//		j_zero = 0;
//		for (j0 = 0; j0 < N_coeff0; j0++) {
//		    if (j0 == Zero_list[j_zero]) {
//			if (j_zero == N_zero) {
//			    exit(EXIT_FAILURE);
//			}
//			j_zero++;
//			continue;
//		    }
//		    // JFR - added 04.11.12: put the matrix in packed form
//		    index0 = index_Lpacked(i0, j0, N0);
//		    index1 = index_Lpacked(i1, j1, N1);
//		    //sys_temp->M[i1][j1] = sys->M[i0][j0];
//		    //sys_temp->M_cnt[i1][j1] = sys->M_cnt[i0][j0];
//		    sys_temp->M[index1] = sys->M[index0];
//		    sys_temp->M2[index1] = sys->M2[index0];
//		    if (sys_temp->M_cnt != NULL) {
//			sys_temp->M_cnt[index1] = sys->M_cnt[index0];
//		    }
//		    sys_temp->d2M[index1] = sys->d2M[index0];
//		    ///*if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 ) {*/ sys_temp->d2M[index1] = sys->d2M[index0]; /*}*/
//
//		    j1++;
//		}
//
//		i1++;
//		i1_inter++;
//	    } else {
//		fprintf(fp, "looping.\n");
//		i_zero++;
//	    }
//
//	    i0 += 1;
//	}
//
//	inter_ptr->N_coeff = i1_inter;
//
//	if ((inter_ptr->i_basis == DELTA_BASIS_INDEX)
//	    || (inter_ptr->i_basis == LINEAR_BASIS_INDEX)
//	    || (inter_ptr->i_basis == BSPLINE_BASIS_INDEX)) {
//	    inter_ptr->N_pts = i1_inter;
//	}
//	/* JFR */
//    }
//
//    return 0;
//}


/********************************************************************************************
remove_zero_rows_mem(): Removes rows of zero from matrix M and the corresponding elements in
the other arrays by moving the elements that correspond to the non-zero g_cnt into the lowest
available spots in the array.  This function is tricky and hard to read, so I've tried to
put extensive comments below. - JFR 04.13.12
********************************************************************************************/
int remove_zero_rows_mem(FILE * fp, int N_coeff0, int N_zero,
			 int *Zero_list, tW_system * sys)
{
    int i_inter;  // the index within the interaction list
    int i0, j0;       // the row and column indices in the untrimmed arrays
    int i1, j1;       // the row and column indices in the trimmed arrays
    int i1_inter; // index within the grid of a particular trimmed interaction
    int i0_inter; // index within the grid of a particular untrimmed interaction

    int i_zero, j_zero; // track the grid position in the zero list for rows and columns, respectively
    int j_zero_start = 0; // the starting position to look at the zero list for columns

    double r, dr; // the grid value and grid spacing for this interaction

    tW_Inter_Types *inter_ptr; // pointer to the current interaction

    int index0, index1; // indices in the packed matrices for the untrimmed and trimmed systems
    int N0 = N_coeff0; // the number of coefficients in the original system
    int N1 = N_coeff0 - N_zero; // the number of coefficients in the trimmed system

    /* This doesn't have memory yet since I didn't make a copy of the system structure */
    sys->x = (double *) ecalloc(N1, sizeof(double));	/* We only need this array to be of the new length */


    i0 = 0;
    i1 = 0;
    i_zero = 0;
    for (i_inter = 0; i_inter < sys->N_Inter_Types; i_inter++) {	/* Loop over the interactions since we need to get r for each gridpt */
	inter_ptr = &(sys->Inter_Types[i_inter]);

	inter_ptr->i_0 = i1;

	dr = inter_ptr->dr;
	r = inter_ptr->R_min - dr;

	i1_inter = 0;
	for (i0_inter = 0; i0_inter < inter_ptr->N_coeff; i0_inter++) {	
        /* This loop + the one above is equivalent to looping over all the coefficients */
        /* i0 keeps track of the coefficient */
	    r += dr;

	    fprintf(fp, "i0: %6d    i0_inter: %6d    sys->g_cnt[%6d]: %15.6lf    \n", i0, i0_inter, i0, sys->g_cnt[i0]);



	    if (i_zero == N_zero + 1) {
		printf
		    ("ERROR: i_zero=%d is larger than N_zero=%d trim_Mb( remove_zeros() ) \n", i_zero, N_zero);
		exit(0);
	    }
	    /* double check the implementation */
	    if (i0 != Zero_list[i_zero]) {	
            /* i0 == Zero_list[i_zero] => i0 is a zero and you don't want to copy anything */
            /* Zero_list has N_coeff0 elements, all elements past N_zero have values -1 */
		fprintf(fp, "i1: %6d    i1_inter: %6d\n", i1, i1_inter);

		/* i0 >= i1, so we don't need a copy here */
		if (i0 < i1) {
		    printf("ERROR: i0=%d<i1=%d in trim_Mb( remove_zeros() ) \n", i0, i1);
		    exit(0);
		}		/* double check the implementation */
		sys->x[i1] = r;
		sys->b[i1] = sys->b[i0];
		sys->b_forces[i1] = sys->b_forces[i0];
		sys->b_struct[i1] = sys->b_struct[i0];
		sys->g[i1] = sys->g[i0];
		sys->g_cnt[i1] = sys->g_cnt[i0];
		sys->L[i1] = sys->L[i0];
	        sys->d2b[i1] = sys->d2b[i0];
		sys->rescale[i1] = sys->rescale[i0];

		j1 = i1;
		j_zero = j_zero_start;
		for (j0 = i0; j0 < N_coeff0; j0++) {	/* Now, for each i0, loop over all the coefficients, j0 */

		    if (j_zero == N_zero + 1) {
			printf("ERROR: j_zero=%d is larger than N_zero=%d trim_Mb( remove_zeros() ) \n", j_zero, N_zero);
			exit(0);
		    }
		    /* double check the implementation */
		    if (j0 == Zero_list[j_zero]) {	
                    /* j0 == j_cmp => j0 is a zero and you don't want to copy anything */
                    /* Zero_list has N_coeff0 elements, all elements past N_zero have values -1 */
			j_zero++;
			continue;
		    }
		    // JFR - added 04.11.12: put the matrix in packed form
		    /*  This should search the matrix array in order such that j0 >= j1, so we don't need a copy here either */
		    if (j0 < j1) {
			printf("ERROR: j0=%d<j1=%d in trim_Mb( remove_zeros() ) \n", j0, j1);
			exit(0);
		    }		/* double check the implementation */
		    index0 = index_Lpacked(i0, j0, N0);
		    index1 = index_Lpacked(i1, j1, N1);

		    sys->M[index1] = sys->M[index0];
		    sys->M2[index1] = sys->M2[index0];
		    if (sys->M_cnt != NULL) {
			sys->M_cnt[index1] = sys->M_cnt[index0];
		    }
		    sys->d2M[index1] = sys->d2M[index0];

		    j1++;
		}

		i1++;
		i1_inter++;
	    } else {
		fprintf(fp, "looping.\n");
		i_zero++;
		j_zero_start++;	/* If you pass a zero in the i variable, you should begin the j search in the next index of Zero_list */
	    }

	    i0 += 1;

	}

	inter_ptr->N_coeff = i1_inter;

	if ((inter_ptr->i_basis == DELTA_BASIS_INDEX)
	    || (inter_ptr->i_basis == LINEAR_BASIS_INDEX)
	    || (inter_ptr->i_basis == BSPLINE_BASIS_INDEX)) {
	    inter_ptr->N_pts = i1_inter;
	}
	/* JFR */
    }

   // sys->N_coeff = N1;
   // sys->N_pack = (N1 * N1 + N1) / 2;

   // /*  NJD - Moved solution array allocation here, to where we first know their dimension */
   // sys->phi = (double *) ecalloc(sys->N_coeff, sizeof(double));
   // sys->phi_forces = (double *) ecalloc(sys->N_coeff, sizeof(double));
   // sys->phi_struct = (double *) ecalloc(sys->N_coeff, sizeof(double));

    return 0;
}


/*****************************************************************************************
resize_arrays_trim_Mb(): erealloc arrays to proper size after trimming.
*****************************************************************************************/
int resize_arrays_trim_Mb(tW_system * sys, int N_coeff)
{

    sys->N_coeff = N_coeff;
    sys->N_pack = (N_coeff * N_coeff + N_coeff) / 2;

    //sys->x = (double *) erealloc ( sys->x, N_coeff*sizeof(double) ); /* no longer need to resize this */
    //if ( sys->x == NULL ) { printf( "ERROR: erealloc failed for sys->x in resize_arrays_trim_Mb \n" ); exit(0); }

    sys->b = (double *) erealloc(sys->b, N_coeff * sizeof(double));
    if (sys->b == NULL) {
        printf
            ("ERROR: erealloc failed for sys->b in resize_arrays_trim_Mb \n");
        exit(0);
    }

    sys->b_forces =
        (double *) erealloc(sys->b_forces, N_coeff * sizeof(double));
    if (sys->b_forces == NULL) {
        printf
            ("ERROR: erealloc failed for sys->b_forces in resize_arrays_trim_Mb \n");
        exit(0);
    }

    sys->b_struct =
        (double *) erealloc(sys->b_struct, N_coeff * sizeof(double));
    if (sys->b_struct == NULL) {
        printf
            ("ERROR: erealloc failed for sys->b_struct in resize_arrays_trim_Mb \n");
        exit(0);
    }

    sys->g = (double *) erealloc(sys->g, N_coeff * sizeof(double));
    if (sys->g == NULL) {
        printf
            ("ERROR: erealloc failed for sys->g in resize_arrays_trim_Mb \n");
        exit(0);
    }

    sys->g_cnt = (double *) erealloc(sys->g_cnt, N_coeff * sizeof(double));
    if (sys->g_cnt == NULL) {
        printf
            ("ERROR: erealloc failed for sys->g_cnt in resize_arrays_trim_Mb \n");
        exit(0);
    }

    sys->L = (double *) erealloc(sys->L, N_coeff * sizeof(double));
    if (sys->L == NULL) {
        printf
            ("ERROR: erealloc failed for sys->L in resize_arrays_trim_Mb \n");
        exit(0);
    }

    sys->M = (double *) erealloc(sys->M, sys->N_pack * sizeof(double));
    if (sys->M == NULL) {
        printf
            ("ERROR: erealloc failed for sys->M in resize_arrays_trim_Mb \n");
        exit(0);
    }

    sys->M2 = (double *) erealloc(sys->M2, sys->N_pack * sizeof(double));
    if (sys->M2 == NULL) {
        printf
            ("ERROR: erealloc failed for sys->M2 in resize_arrays_trim_Mb \n");
        exit(0);
    }

    sys->M_cnt =
        (double *) erealloc(sys->M_cnt, sys->N_pack * sizeof(double));
    if (sys->M_cnt == NULL) {
        printf
            ("ERROR: erealloc failed for sys->M_cnt in resize_arrays_trim_Mb \n");
        exit(0);
    }


    /*  NJD - Moved solution array allocation here, to where we first know their dimension */
    sys->phi = (double *) ecalloc(sys->N_coeff, sizeof(double));
    sys->phi_forces = (double *) ecalloc(sys->N_coeff, sizeof(double));
    sys->phi_struct = (double *) ecalloc(sys->N_coeff, sizeof(double));

    return 0;
}

///*****************************************************************************************
//copy_arrays_trim_Mb(): Copy arrays from sys_temp to sys.
//*****************************************************************************************/
//int copy_arrays_trim_Mb(tW_system * sys, tW_system sys_temp, int N_coeff)
//{
//    int i;
//    int N_pack = (N_coeff * N_coeff + N_coeff) / 2;
//
//    for (i = 0; i < N_coeff; i++) {
//	sys->x[i] = sys_temp.x[i];
//	sys->b[i] = sys_temp.b[i];
//	sys->b_forces[i] = sys_temp.b_forces[i];
//	sys->b_struct[i] = sys_temp.b_struct[i];
//	sys->g[i] = sys_temp.g[i];
//	sys->g_cnt[i] = sys_temp.g_cnt[i];
//	sys->L[i] = sys_temp.L[i];
//	sys->d2b[i] = sys_temp.d2b[i];
//	///*if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) {*/ sys->d2b[i] = sys_temp.d2b[i]; /*}*/
//	sys->rescale[i] = sys_temp.rescale[i];
//    }
//
//    // JFR - added 04.11.12: put the matrix in packed form
//    for (i = 0; i < N_pack; i++) {
//	sys->M[i] = sys_temp.M[i];
//	sys->M2[i] = sys_temp.M2[i];
//	if (sys->M_cnt != NULL) {
//	    sys->M_cnt[i] = sys_temp.M_cnt[i];
//	}
//	sys->d2M[i] = sys_temp.d2M[i];
//	///*if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 ) {*/ sys->d2M[i] = sys_temp.d2M[i]; /*}*/
//    }
//
//    sys->N_coeff = N_coeff;
//    sys->N_pack = N_pack;
//
//    return 0;
//}


/*****************************************************************************************
set_pointers_trim_Mb(): Set pointers to master array.
*****************************************************************************************/
int set_pointers_trim_Mb(tW_system * sys)
{
    int i;
    int i_0;
    tW_Inter_Types *inter;

    /* Loop over all interaction types. */
    for (i = 0; i < sys->N_Inter_Types; i++) {
	inter = &(sys->Inter_Types[i]);
	i_0 = inter->i_0;

	/* Set pointers to master grid. */
	inter->ptr_x = &(sys->x[i_0]);
	inter->ptr_b = &(sys->b[i_0]);
	inter->ptr_b_forces = &(sys->b_forces[i_0]);
	inter->ptr_b_struct = &(sys->b_struct[i_0]);
	inter->ptr_g = &(sys->g[i_0]);
	inter->ptr_g_cnt = &(sys->g_cnt[i_0]);
	inter->ptr_L = &(sys->L[i_0]);
	inter->ptr_phi = &(sys->phi[i_0]);
	inter->ptr_phi_forces = &(sys->phi_forces[i_0]);
	inter->ptr_phi_struct = &(sys->phi_struct[i_0]);

    }




    return 0;
}


/*****************************************************************************************
get_phi_struct_and_phi_forces(): Calculates the coefficients using b_forces if forces
were available and b_struct in all cases.
*****************************************************************************************/
int get_phi_struct_and_phi_forces(FILE * fp_log, tW_system * sys,
				  int forces)
{
    fprintf(fp_log, "\nIn get_phi_struct_and_phi_forces.\n");

    if (forces) {
	get_phi_forces(fp_log, sys);
    }

    get_phi_struct(fp_log, sys);

    /* JFR - 12.20.13: evaluate options for calculating the eigenspectrum of MT */

    if (sys->SVD_var.flag_printSV == TRUE) {	/* Explicitly calculate and print the singular values */
	/* Explicitly calculate and print the singular values */
	/* eigenvectors can be printed with eigen(), below */
	SVD(fp_log, sys);
    }

    if (sys->Eigen_var.flag_Eigen == TRUE) {
	eigen(fp_log, sys, FALSE);	/* calculate the eigenvalues and eigenvectors of the full metric tensor */

	if (sys->Eigen_var.flag_Gbar == TRUE) {
	    eigen(fp_log, sys, TRUE);	/* calculate the eigenvalues and eigenvectors of the ``3-body'' components only */
	}
    }

    return 0;
}


/*****************************************************************************************
get_phi_forces(): Calculates the coefficients when forces are available.
*****************************************************************************************/
int get_phi_forces(FILE * fp_log, tW_system * sys)
{
    int i;
    int N_coeff = sys->N_coeff;

    fprintf(fp_log, "\nIn get_phi_forces.\n");

    for (i = 0; i < N_coeff; i++) {
	sys->b[i] = sys->b_forces[i];
    }

    direct_solve_lin_eqns(fp_log, sys, "forces");

    for (i = 0; i < N_coeff; i++) {
	sys->phi_forces[i] = sys->phi[i];
    }

    for (i = 0; i < N_coeff; i++) {
	sys->b[i] = 0.0;
	sys->phi[i] = 0.0;
    }

    return 0;
}

/*****************************************************************************************
get_phi_struct(): Calculates the coefficients using b_struct.
*****************************************************************************************/
int get_phi_struct(FILE * fp_log, tW_system * sys)
{
    int i;
    int N_coeff = sys->N_coeff;

    fprintf(fp_log, "\nIn get_phi_struct.\n");

    for (i = 0; i < N_coeff; i++) {
	sys->b[i] = sys->b_struct[i];
    }

    direct_solve_lin_eqns(fp_log, sys, "struct");

    for (i = 0; i < N_coeff; i++) {
	sys->phi_struct[i] = sys->phi[i];
    }

    for (i = 0; i < N_coeff; i++) {
	sys->b[i] = 0.0;
	sys->phi[i] = 0.0;
    }

    return 0;
}

/*****************************************************************************************
direct_solve_lin_eqns(): JFR - 12.20.13: takes care of deciding which matrix inversion
routine to use for both the force and struct calc
*****************************************************************************************/
int direct_solve_lin_eqns(FILE * fp_log, tW_system * sys, tW_word info)
{
    int N_coeff = sys->N_coeff;

    /* regularization parameters for bayes method with equality constraints */
    double alpha[N_coeff];
    double beta = 0.00;

    /* JFR - 12.03.13: new forks */
    if ((sys->REG_var.flag_REG == FALSE) && (sys->CONSTRAIN_var.flag_CONSTRAIN == FALSE)) 	/* No regularization or constraints */
    {
	/* do the normal calculations */
	if ((sys->MEM_var.flag_LOWMEM == TRUE) && (strcmp(sys->MEM_var.info, info) == 0)) 
	{
	    solv_lin_eqns_symm(fp_log, sys);

	} else if (sys->MEM_var.flag_LOWMEM == FALSE) 
	{
	    solv_lin_eqns(fp_log, sys);
	}
    } else {
	if ((sys->REG_var.flag_REG == TRUE)) {	/* Regularize */
	    if (strcmp(sys->REG_var.type, "BAYES") == 0) 
	    {	/* Regularize with Bayesian inference method */
		solv_lin_eqns_Bayesian_saveparams(fp_log, sys, alpha, &beta);

	    } else if (strcmp(sys->REG_var.type, "UNCERT") == 0) 
	    {	/* Regularize with simple scaling of MT uncertainties */
		beta = sys->REG_var.tau_beta;

	    } else 
	    {
		printf("ERROR: Regularization type incorrect");
		exit(0);
	    }
	} else 
	{
	    beta = 0.00;
	}			/* No regularization */

	solv_lin_eqns_constrain_tabdih_regularize(fp_log, sys, alpha, beta);	/* Solve with regularization and/or constraints */
    }


    /* Solve by Perturbation Theory (Iterative Inversion) */
    if (sys->PT_var.flag_PT == TRUE) 
    {
	solv_lin_eqns_PT(fp_log, sys);	/* JFR- solves the system of eqns up to N_PT orders by matrix perturbation theory */
    }

    return 0;
}

/*****************************************************************************************
eval_bond_basis_vectors(): Evaluates the basis vectors for bonded interactions for site
i_site. Indices for each basis vector are stored in D_b[] and each vector is stored in 
calcG_b[]. The functions g, g_cnt, and L are also evaluated for each interaction.
NOTE: For calc_grids2(), this function gets called twice, but g, g_cnt, and L should only
be evaluated once. Therefore, I add i_flag to indicate if these functions should be
evaluated.
*****************************************************************************************/
int eval_bond_basis_vectors(FILE * fp, tW_CG_site * CG_struct,
			    tW_Bonded_Inter * Bonded_Inter_Types,
			    int i_site, int *D_b, dvec * calG_b,
			    bool i_flag,
			    int *D_b_inter_index
			    /* JFR - 02.26.13: for sep M2 */ )
{
    int i;
    int bond_coeff_ctr = 0;
    int tmp_coeff_ctr, j;	/* JFR - 02.26.13: for sep M2 */
    tW_CG_site *i_site_ptr = &(CG_struct[i_site]);

    /* JFR - 02.26.13: I changed a few things concerning bond_coeff_ctr below for sep M2 */

    if (DEBUG_eval_bond_basis) {
	fprintf(fp, "\nIn eval_bond_basis_vectors() for Site: %d.\n",
		i_site);
    }

    /* Loop over all bond interactions found in the GROMACS topology. */
   for (i = 0; i < i_site_ptr->nr_bonds; i++) {
//fprintf(stderr,"i_site_ptr->bond_type[i]: %d  Bonded_Inter_Types[i_site_ptr->bond_type[i]].name: %s\n",i_site_ptr->bond_type[i],Bonded_Inter_Types[i_site_ptr->bond_type[i]].name); // MRD
	if (strcmp
	    (B_BOND_STRETCH,
	     Bonded_Inter_Types[i_site_ptr->bond_type[i]].name) == 0) {
	    if (DEBUG_eval_bond_basis) {
		print_debug_eval_bond_basis(i_site_ptr, Bonded_Inter_Types,
					    i, fp);
	    }
	    tmp_coeff_ctr =
		get_BondStretch_basis_vectors(i_site, i, bond_coeff_ctr,
					      i_site_ptr,
					      Bonded_Inter_Types, D_b,
					      calG_b, CG_struct, i_flag);

	} else
	    if (strcmp
		(B_ANGLE,
		 Bonded_Inter_Types[i_site_ptr->bond_type[i]].name) == 0) {
	    if (DEBUG_eval_bond_basis) {
		print_debug_eval_bond_basis(i_site_ptr, Bonded_Inter_Types,
					    i, fp);
	    }
	    tmp_coeff_ctr =
		get_Angle_basis_vectors(i_site, i, bond_coeff_ctr,
					i_site_ptr, Bonded_Inter_Types,
					D_b, calG_b, CG_struct, i_flag);
	} else
	    if (strcmp
		(B_DIHEDRAL,
		 Bonded_Inter_Types[i_site_ptr->bond_type[i]].name) == 0) {
	    if (DEBUG_eval_bond_basis) {
		print_debug_eval_bond_basis(i_site_ptr, Bonded_Inter_Types,
					    i, fp);
	    }
	    tmp_coeff_ctr =
		get_Dihedral_basis_vectors(i_site, i, bond_coeff_ctr,
					   i_site_ptr, Bonded_Inter_Types,
					   D_b, calG_b, CG_struct, i_flag);
	} else
	    if (strcmp
		(B_NB_PAIR_BOND,
		 Bonded_Inter_Types[i_site_ptr->bond_type[i]].name) == 0) {
	    if (DEBUG_eval_bond_basis) {
		print_debug_eval_bond_basis(i_site_ptr, Bonded_Inter_Types,
					    i, fp);
	    }
	    tmp_coeff_ctr =
		get_IntraMolec_NB_Pair_basis_vectors(i_site, i,
						     bond_coeff_ctr,
						     i_site_ptr,
						     Bonded_Inter_Types,
						     D_b, calG_b,
						     CG_struct, i_flag);
	} else {
	    printf
		("\nERROR: Problem reading inter. name in eval_bond_basis_vectors() for Site: %d.\n",
		 i_site);
	    exit(EXIT_FAILURE);
	}

	for (j = 0; j < tmp_coeff_ctr; j++) {
	    D_b_inter_index[bond_coeff_ctr + j] = i;	/* JFR - 02.26.13: for sep M2 */
	}
	bond_coeff_ctr += tmp_coeff_ctr;

    }				/* End loop over total_no_bonds for i_site. */

    return bond_coeff_ctr;
}


/*****************************************************************************************
get_IntraMolec_NB_Pair_basis_vectors():
*****************************************************************************************/
int get_IntraMolec_NB_Pair_basis_vectors(int i_site, int i_inter,
					 int bond_coeff_ctr,
					 tW_CG_site * i_site_ptr,
					 tW_Bonded_Inter *
					 Bonded_Inter_Types, int *D_b,
					 dvec * calG_b,
					 tW_CG_site * CG_struct,
					 bool i_flag)
{
    int i;
    int n_coeff;
    int i_bond_type;
    int i_basis_type;
    int *sites;
    double norm_ij;
    double power;
    double exponent;
    double r_exp;
    int *powers;
    dvec u_ij;

    i_bond_type = i_site_ptr->bond_type[i_inter];
    i_basis_type = Bonded_Inter_Types[i_bond_type].i_basis;
    sites = i_site_ptr->bond_site[i_inter];

    norm_ij = get_IntraMolec_NB_Pair_info(i_site, sites, CG_struct, u_ij);

    if (Bonded_Inter_Types[i_bond_type].i_basis == DELTA_BASIS_INDEX) {
	n_coeff =
	    get_delta_IntraMolec_NB_Pair_basis_vector(bond_coeff_ctr,
						      i_bond_type, D_b,
						      calG_b,
						      Bonded_Inter_Types,
						      norm_ij, u_ij,
						      i_flag);
    }
    /* START JFR */
    else if (Bonded_Inter_Types[i_bond_type].i_basis == LINEAR_BASIS_INDEX) {
	n_coeff =
	    get_linear_IntraMolec_NB_Pair_basis_vector(bond_coeff_ctr,
						       i_bond_type, D_b,
						       calG_b,
						       Bonded_Inter_Types,
						       norm_ij, u_ij,
						       i_flag);
    }
    /* END JFR */
    else if (Bonded_Inter_Types[i_bond_type].i_basis ==
	     BSPLINE_BASIS_INDEX) {
	n_coeff =
	    get_Bspline_IntraMolec_NB_Pair_basis_vector(bond_coeff_ctr,
							i_bond_type, D_b,
							calG_b,
							Bonded_Inter_Types,
							norm_ij, u_ij,
							i_flag);
    } else if (Bonded_Inter_Types[i_bond_type].i_basis == POWER_INDEX) {
	powers = Bonded_Inter_Types[i_bond_type].powers;
	for (i = 0; i < Bonded_Inter_Types[i_bond_type].N_powers; i++) {
	    exponent = powers[i] + 1.0;
	    power = pow(norm_ij, exponent);
	    scal_times_vect(powers[i] / power, u_ij,
			    calG_b[bond_coeff_ctr + i]);
	    D_b[bond_coeff_ctr + i] =
		Bonded_Inter_Types[i_bond_type].i_0 + i;

	    exponent = powers[i] + 2.0;
	    r_exp = pow(norm_ij, exponent);

	    if (i_flag == TRUE) {
		Bonded_Inter_Types[i_bond_type].ptr_g[i] +=
		    (-1.0) * (powers[i]) * (powers[i] - 1.0) / r_exp;
		Bonded_Inter_Types[i_bond_type].ptr_g_cnt[i] += 1.0;
//        Bonded_Inter_Types[i_bond_type].ptr_g[i]     += (-1.0) * (powers[i]) * (powers[i]+1.0) / r_exp;
//        Bonded_Inter_Types[i_bond_type].ptr_L[i]     += 2.0 * (powers[i]) / r_exp;
	    }
	}
	n_coeff = Bonded_Inter_Types[i_bond_type].N_powers;
    } else {
	printf
	    ("\nERROR: Unsupported basis in get_IntraMolec_NB_Pair_basis_vectors().\n");
	exit(EXIT_FAILURE);
    }


    return n_coeff;
}


/*****************************************************************************************
get_delta_IntraMolec_NB_Pair_basis_vector():
*****************************************************************************************/
int get_delta_IntraMolec_NB_Pair_basis_vector(int bond_coeff_ctr,
					      int i_bond_type, int *D_b,
					      dvec * calG_b,
					      tW_Bonded_Inter *
					      Bonded_Inter_Types,
					      double norm_ij, dvec u_ij,
					      bool i_flag)
{
    double dr = Bonded_Inter_Types[i_bond_type].dr;
    double R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    double R_max = Bonded_Inter_Types[i_bond_type].R_max;
    int i_0 = Bonded_Inter_Types[i_bond_type].i_0;
    dvec zero = { 0.0, 0.0, 0.0 };

    if ((norm_ij >= R_0) && (norm_ij < (R_max + 0.5 * dr))) {
	D_b[bond_coeff_ctr] =
	    get_grid_index_for_delta_basis(norm_ij, i_0, dr, R_0);
	copy_vector(u_ij, calG_b[bond_coeff_ctr]);

	if (i_flag == TRUE) {
	    Bonded_Inter_Types[i_bond_type].ptr_g[D_b[bond_coeff_ctr] -
						  i_0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[D_b[bond_coeff_ctr] -
						      i_0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr] -
						  i_0] += 2.0 / norm_ij;
	}
    } else {
	return 0;
    }

    return 1;			/* Number of coeff. to determine. */
}

/* START JFR */
/*****************************************************************************************
get_linear_IntraMolec_NB_Pair_basis_vector():
*****************************************************************************************/
int get_linear_IntraMolec_NB_Pair_basis_vector(int bond_coeff_ctr,
					       int i_bond_type, int *D_b,
					       dvec * calG_b,
					       tW_Bonded_Inter *
					       Bonded_Inter_Types,
					       double norm_ij, dvec u_ij,
					       bool i_flag)
{
    int n_coeff;

    n_coeff =
	get_linear_bond_basis_vector(bond_coeff_ctr, i_bond_type, D_b,
				     calG_b, Bonded_Inter_Types, norm_ij,
				     u_ij, i_flag);


    return n_coeff;		/* Number of coeff. to determine. */
}

/* END JFR */

/*****************************************************************************************
get_Bspline_IntraMolec_NB_Pair_basis_vector(): JFR - 07.23.12: Bspline
*****************************************************************************************/
int get_Bspline_IntraMolec_NB_Pair_basis_vector(int bond_coeff_ctr,
						int i_bond_type, int *D_b,
						dvec * calG_b,
						tW_Bonded_Inter *
						Bonded_Inter_Types,
						double norm_ij, dvec u_ij,
						bool i_flag)
{
    int n_coeff;

    n_coeff =
	get_Bspline_bond_basis_vector(bond_coeff_ctr, i_bond_type, D_b,
				      calG_b, Bonded_Inter_Types, norm_ij,
				      u_ij, i_flag);


    return n_coeff;		/* Number of coeff. to determine. */
}

/*****************************************************************************************
get_IntrarMolec_NB_Pair_info():
*****************************************************************************************/
double get_IntraMolec_NB_Pair_info(int i_site, int *sites,
				   tW_CG_site * CG_struct, dvec u_ij)
{
    int j_site;
    double norm_ij;
    dvec r_ij;

    if (i_site == sites[0]) {
	j_site = sites[1];
    } else {
	j_site = sites[0];
	if (i_site != sites[1]) {
	    printf("\nError in get_IntraMolec_NB_Pair_info().\n");
	    exit(EXIT_FAILURE);
	}
    }

    get_difference_unit_vector(CG_struct[i_site].r, CG_struct[j_site].r,
			       r_ij, u_ij, &norm_ij);

    return norm_ij;
}


/*****************************************************************************************
get_Dihedral_basis_vectors(): 
*****************************************************************************************/
int get_Dihedral_basis_vectors(int i_site, int i_inter, int bond_coeff_ctr,
			       tW_CG_site * i_site_ptr,
			       tW_Bonded_Inter * Bonded_Inter_Types,
			       int *D_b, dvec * calG_b,
			       tW_CG_site * CG_struct, bool i_flag)
{
    int n_coeff;
    int i_bond_type;
    int i_basis_type;
    int *sites;
    tW_dihedral dihedral;

    i_bond_type = i_site_ptr->bond_type[i_inter];
    i_basis_type = Bonded_Inter_Types[i_bond_type].i_basis;
    dihedral.dihedral_quartet = i_site_ptr->bond_site[i_inter];

    dihedral.atom_position =
	get_atom_position_dihedral(dihedral.dihedral_quartet, i_site);

    eval_grad_B(&dihedral, CG_struct);


    if (Bonded_Inter_Types[i_bond_type].i_basis ==
	RYCKAERT_BELLEMANS_BASIS_INDEX) {
	n_coeff =
	    get_RB_dihedral_basis_vector(bond_coeff_ctr, i_bond_type, D_b,
					 calG_b, Bonded_Inter_Types,
					 dihedral, i_flag);
    } else if (Bonded_Inter_Types[i_bond_type].i_basis ==
	       DELTA_BASIS_INDEX) {
	n_coeff =
	    get_delta_dihedral_basis_vector(bond_coeff_ctr, i_bond_type,
					    D_b, calG_b,
					    Bonded_Inter_Types, dihedral,
					    i_flag);
    }
    /* START JFR */
    else if (Bonded_Inter_Types[i_bond_type].i_basis == LINEAR_BASIS_INDEX) {
	n_coeff =
	    get_linear_dihedral_basis_vector(bond_coeff_ctr, i_bond_type,
					     D_b, calG_b,
					     Bonded_Inter_Types, dihedral,
					     i_flag);
    }
    /* END JFR */
    else if (Bonded_Inter_Types[i_bond_type].i_basis == BSPLINE_BASIS_INDEX) {	/* JFR - 07.23.12: Bspline */
	n_coeff =
	    get_Bspline_dihedral_basis_vector(bond_coeff_ctr, i_bond_type,
					      D_b, calG_b,
					      Bonded_Inter_Types, dihedral,
					      i_flag);
    } else if (Bonded_Inter_Types[i_bond_type].i_basis == TOY_DIHED_INDEX) {
	n_coeff =
	    get_TOY_dihedral_basis_vector(bond_coeff_ctr, i_bond_type, D_b,
					  calG_b, Bonded_Inter_Types,
					  dihedral, i_flag);
    } else {
	printf
	    ("\nERROR: Unsupported basis in get_Dihedral_basis_vectors().\n");
	exit(EXIT_FAILURE);
    }

    return n_coeff;
}


/*****************************************************************************************
eval_grad_B(): Code adopted from pages 278-283 of Rapaport's "The Art of MD Sim." This
gives the same result as Noid's notes, but this implementation is more concise.
*****************************************************************************************/
void eval_grad_B(tW_dihedral * dihedral, tW_CG_site * CG_struct)
{
    int *quartet = dihedral->dihedral_quartet;
    double C_jj, C_jk, C_jl, C_kk, C_kl, C_ll;
    double t1, t2, t3, t4, t5, t6;
    double sqrt_q, p;
    double K, L;
    double ipr;
    double sign;
    dvec r_i, r_j, r_k, r_l;
    dvec b_j, b_k, b_l;
    dvec f_i, f_j, f_k, f_l;
    dvec b_j_temp, b_k_temp, b_l_temp;
    dvec f_i_temp, f_l_temp;
    dvec n_jkl;


    /* Copy coordinates into local variables. */
    copy_vector(CG_struct[quartet[0]].r, r_i);
    copy_vector(CG_struct[quartet[1]].r, r_j);
    copy_vector(CG_struct[quartet[2]].r, r_k);
    copy_vector(CG_struct[quartet[3]].r, r_l);

    /* Calculate bond vectors. */
    vect_diff(r_j, r_i, b_j);
    vect_diff(r_k, r_j, b_k);
    vect_diff(r_l, r_k, b_l);

    /* Calculate the appropriate dot products of bond vectors. */
    C_jj = dot_prod(b_j, b_j);
    C_jk = dot_prod(b_j, b_k);
    C_jl = dot_prod(b_j, b_l);
    C_kk = dot_prod(b_k, b_k);
    C_kl = dot_prod(b_k, b_l);
    C_ll = dot_prod(b_l, b_l);

    /* Calculate the scalars. */
    t1 = (C_jl * C_kk) - (C_jk * C_kl);
    t2 = (C_jj * C_kl) - (C_jk * C_jl);
    t3 = (C_jk * C_jk) - (C_jj * C_kk);
    t4 = (C_kk * C_ll) - (C_kl * C_kl);
    t5 = (C_jl * C_kl) - (C_jk * C_ll);
    t6 = -t1;
    sqrt_q = sqrt(-1.0 * t3 * t4);

    /* Calculate f_i. */
    scal_times_vect(t1, b_j, b_j_temp);
    scal_times_vect(t2, b_k, b_k_temp);
    scal_times_vect(t3, b_l, b_l_temp);
    vect_sum(b_j_temp, b_k_temp, f_i);
    vect_sum(b_l_temp, f_i, f_i);
    K = C_kk / (-1.0 * sqrt_q * t3);
    scal_times_vect(K, f_i, f_i);

    /* Calculate f_l. */
    scal_times_vect(t4, b_j, b_j_temp);
    scal_times_vect(t5, b_k, b_k_temp);
    scal_times_vect(t6, b_l, b_l_temp);
    vect_sum(b_j_temp, b_k_temp, f_l);
    vect_sum(b_l_temp, f_l, f_l);
    K = C_kk / (sqrt_q * t4);
    scal_times_vect(K, f_l, f_l);

    /* Calculate f_j. */
    K = -1.0 * (1.0 + C_jk / C_kk);
    scal_times_vect(K, f_i, f_i_temp);
    L = C_kl / C_kk;
    scal_times_vect(L, f_l, f_l_temp);
    vect_sum(f_i_temp, f_l_temp, f_j);

    /* Calculate f_k. */
    K = C_jk / C_kk;
    scal_times_vect(K, f_i, f_i_temp);
    L = -1.0 * (1.0 + C_kl / C_kk);
    scal_times_vect(L, f_l, f_l_temp);
    vect_sum(f_i_temp, f_l_temp, f_k);

    /* Save grad_B. */
    copy_vector(f_i, dihedral->grad_B[0]);
    copy_vector(f_j, dihedral->grad_B[1]);
    copy_vector(f_k, dihedral->grad_B[2]);
    copy_vector(f_l, dihedral->grad_B[3]);

    /* Calculate the dihedral angle. */
    p = C_jl * C_kk - C_jk * C_kl;
    dihedral->cos_angle = -1.0 * p / sqrt_q;
    if (fabs(dihedral->cos_angle) > 1.0) {
	if (dihedral->cos_angle > 1.0) {
	    dihedral->cos_angle = 1.0;
	} else {
	    dihedral->cos_angle = -1.0;
	}
    }
    dihedral->angle = acos(dihedral->cos_angle);
    cross_product(b_k, b_l, n_jkl);
    ipr = dot_prod(b_j, n_jkl);
    sign = (ipr < 0.0) ? -1.0 : 1.0;
    dihedral->angle = (dihedral->angle) * sign;
    dihedral->sin_angle = sin(dihedral->angle);

}



/*****************************************************************************************
get_TOY_dihedral_basis_vector():
*****************************************************************************************/
int get_TOY_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				  int *D_b, dvec * calG_b,
				  tW_Bonded_Inter * Bonded_Inter_Types,
				  tW_dihedral dihedral, bool i_flag)
{
    double psi = dihedral.angle;
    double g_temp;
    dvec grad_dihedral;
    dvec grad_B;


    copy_vector(dihedral.grad_B[dihedral.atom_position - 1], grad_B);

    D_b[bond_coeff_ctr] = Bonded_Inter_Types[i_bond_type].i_0;
    D_b[bond_coeff_ctr + 1] = Bonded_Inter_Types[i_bond_type].i_0 + 1;
    D_b[bond_coeff_ctr + 2] = Bonded_Inter_Types[i_bond_type].i_0 + 2;

    if (fabs(dihedral.sin_angle) > DIHED_PREC) {
	scal_times_vect((1.0), grad_B, calG_b[bond_coeff_ctr]);
	scal_times_vect((3.0 * sin(3 * psi) / dihedral.sin_angle), grad_B,
			calG_b[bond_coeff_ctr + 1]);
	scal_times_vect((1.0 * sin(psi + 0.785398163) /
			 dihedral.sin_angle), grad_B,
			calG_b[bond_coeff_ctr + 2]);

	scal_times_vect(1.0 / (dihedral.sin_angle), grad_B, grad_dihedral);
	g_temp = get_g_dihedral_angle(grad_dihedral);

	if (i_flag == TRUE) {
	    Bonded_Inter_Types[i_bond_type].ptr_g[0] +=
		dihedral.cos_angle * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g[1] +=
		9.0 * cos(3.0 * psi) * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g[2] +=
		cos(psi + 0.785398163) * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[0] += 1;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[1] += 1;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[2] += 1;
	    Bonded_Inter_Types[i_bond_type].ptr_L[0] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[1] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[2] += 0.0;
	}

    } else {
	if (i_flag == TRUE) {
	    Bonded_Inter_Types[i_bond_type].ptr_g[0] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g[1] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g[2] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[0] += 1;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[1] += 1;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[2] += 1;
	    Bonded_Inter_Types[i_bond_type].ptr_L[0] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[1] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[2] += 0.0;
	}
    }

    return 3;			/* Number of coeff. to determine. */
}


/*****************************************************************************************
get_RB_dihedral_basis_vector():
*****************************************************************************************/
int get_RB_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				 int *D_b, dvec * calG_b,
				 tW_Bonded_Inter * Bonded_Inter_Types,
				 tW_dihedral dihedral, bool i_flag)
{
    double phi;
    double cos1, cos2, cos3, cos4, cos5, sin1, sin2;
    double g_temp;
    dvec grad_B;
    dvec grad_dihedral;

    copy_vector(dihedral.grad_B[dihedral.atom_position - 1], grad_B);

    D_b[bond_coeff_ctr] = Bonded_Inter_Types[i_bond_type].i_0;
    D_b[bond_coeff_ctr + 1] = Bonded_Inter_Types[i_bond_type].i_0 + 1;
    D_b[bond_coeff_ctr + 2] = Bonded_Inter_Types[i_bond_type].i_0 + 2;
    D_b[bond_coeff_ctr + 3] = Bonded_Inter_Types[i_bond_type].i_0 + 3;
    D_b[bond_coeff_ctr + 4] = Bonded_Inter_Types[i_bond_type].i_0 + 4;

    phi = dihedral.angle - M_PI;
    cos1 = cos(phi);
    cos2 = cos1 * cos1;
    cos3 = cos2 * cos1;
    cos4 = cos3 * cos1;
    cos5 = cos4 * cos1;
    sin1 = sin(phi);
    sin2 = sin1 * sin1;

    scal_times_vect((-1.0), grad_B, calG_b[bond_coeff_ctr]);
    scal_times_vect((-2.0 * cos1), grad_B, calG_b[bond_coeff_ctr + 1]);
    scal_times_vect((-3.0 * cos2), grad_B, calG_b[bond_coeff_ctr + 2]);
    scal_times_vect((-4.0 * cos3), grad_B, calG_b[bond_coeff_ctr + 3]);
    scal_times_vect((-5.0 * cos4), grad_B, calG_b[bond_coeff_ctr + 4]);
    scal_times_vect(1.0 / (dihedral.sin_angle), grad_B, grad_dihedral);

    /* For now, avoid the 1 / sin( phi ) term where it blows up! */
    if (fabs(dihedral.sin_angle) > 0.0) {
	/* The Laplacian term is zero. */
	if (i_flag == TRUE) {
	    g_temp = get_g_dihedral_angle(grad_dihedral);
	    Bonded_Inter_Types[i_bond_type].ptr_g[0] += cos1 * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g[1] +=
		2 * (cos2 - sin2) * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g[2] +=
		(3 * cos3 - 6 * cos1 * sin2) * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g[3] +=
		(4 * cos4 - 12 * cos2 * sin2) * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g[4] +=
		(5 * cos5 - 20 * cos3 * sin2) * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[1] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[2] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[3] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[4] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[0] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[1] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[2] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[3] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[4] += 0.0;
	}
    } else {
	if (i_flag == TRUE) {
	    Bonded_Inter_Types[i_bond_type].ptr_g[0] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g[1] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g[2] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g[3] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g[4] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[1] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[2] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[3] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[4] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[0] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[1] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[2] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[3] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[4] += 0.0;
	}
    }

    return 5;			/* Number of coeff. to determine. */
}


/*****************************************************************************************
get_delta_dihedral_basis_vector(): Range is not checked since the R_min and R_max are
assumed to be -180 and 180. For now, everyting is set to zero if sin(phi) gets
too close to 0. When the 180 grid point is updated, then the index is switched to the 
-180 grid point since these are the same point. 
*****************************************************************************************/
int get_delta_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				    int *D_b, dvec * calG_b,
				    tW_Bonded_Inter * Bonded_Inter_Types,
				    tW_dihedral dihedral, bool i_flag)
{
    int index;
    int N_pts = Bonded_Inter_Types[i_bond_type].N_pts;
    double dr = Bonded_Inter_Types[i_bond_type].dr;
    double R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    double R_max = Bonded_Inter_Types[i_bond_type].R_max;
    int i_0 = Bonded_Inter_Types[i_bond_type].i_0;
    dvec grad_B;


    copy_vector(dihedral.grad_B[dihedral.atom_position - 1], grad_B);

    if ((dihedral.angle >= R_0) && (dihedral.angle < (R_max + 0.5 * dr))) {
	if (fabs(dihedral.sin_angle) > DIHED_PREC) {
	    index =
		get_grid_index_for_delta_basis(dihedral.angle, i_0, dr,
					       R_0);

	    if (index == (i_0 + N_pts - 1)) {
		index = i_0;
	    }

	    D_b[bond_coeff_ctr] = index;
	    scal_times_vect(1.0 / (dihedral.sin_angle), grad_B,
			    calG_b[bond_coeff_ctr]);

	    if (i_flag == TRUE) {
		Bonded_Inter_Types[i_bond_type].ptr_g[D_b[bond_coeff_ctr] -
						      i_0] +=
		    get_g_dihedral_angle(calG_b[bond_coeff_ctr]);
		Bonded_Inter_Types[i_bond_type].
		    ptr_g_cnt[D_b[bond_coeff_ctr] - i_0] += 1.0;
		Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr] -
						      i_0] += 0.0;
	    }
	} else {
	    D_b[bond_coeff_ctr] =
		get_grid_index_for_delta_basis(dihedral.angle, i_0, dr,
					       R_0);
	    scal_times_vect(0.0, grad_B, calG_b[bond_coeff_ctr]);
	}
    } else {
	return 0;
    }

    return 1;			/* Number of coeff. to determine. */
}

/* START JFR */
/*****************************************************************************************
get_linear_dihedral_basis_vector(): Range is not checked since the R_min and R_max are
assumed to be -180 and 180. For now, everyting is set to zero if sin(phi) gets
too close to 0. When the 180 grid point is updated, then the index is switched to the 
-180 grid point since these are the same point. 
*****************************************************************************************/
int get_linear_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				     int *D_b, dvec * calG_b,
				     tW_Bonded_Inter * Bonded_Inter_Types,
				     tW_dihedral dihedral, bool i_flag)
{
    int N_pts = Bonded_Inter_Types[i_bond_type].N_pts;
    double dr = Bonded_Inter_Types[i_bond_type].dr;
    double R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    double R_max = Bonded_Inter_Types[i_bond_type].R_max;
    int i_0 = Bonded_Inter_Types[i_bond_type].i_0;
    dvec grad_B;
    double A, B, Ap, Bp;

     /*JFR*/ dvec gtemp;


    copy_vector(dihedral.grad_B[dihedral.atom_position - 1], grad_B);

    if (fabs(dihedral.sin_angle) > DIHED_PREC) {
	A = calc_linear_spline_A(bond_coeff_ctr, D_b, dr, R_0, i_0, N_pts,
				 dihedral.angle, &Ap, &Bp, TRUE);
	B = 1 - A;
	/* A and B are weighting factors which depend on where the angle 
	   is relative to the two nearest gridpoints */

	scal_times_vect(A / (dihedral.sin_angle), grad_B,
			calG_b[bond_coeff_ctr]);
	scal_times_vect(B / (dihedral.sin_angle), grad_B,
			calG_b[bond_coeff_ctr + 1]);

	if (i_flag == TRUE) {
	    if (LINEAR_STRUCT == 0) {
		scal_times_vect(1.0 / (dihedral.sin_angle), grad_B, gtemp);
		Bonded_Inter_Types[i_bond_type].ptr_g[D_b[bond_coeff_ctr] -
						      i_0] +=
		    A * get_g_dihedral_angle(gtemp);
		scal_times_vect(1.0 / (dihedral.sin_angle), grad_B, gtemp);
		Bonded_Inter_Types[i_bond_type].
		    ptr_g[D_b[bond_coeff_ctr + 1] - i_0] +=
		    B * get_g_dihedral_angle(gtemp);
	    } else {
		scal_times_vect(1.0 / (dihedral.sin_angle), grad_B, gtemp);
		Bonded_Inter_Types[i_bond_type].ptr_g[D_b[bond_coeff_ctr] -
						      i_0] +=
		    Ap * get_g_dihedral_angle(gtemp);
		scal_times_vect(1.0 / (dihedral.sin_angle), grad_B, gtemp);
		Bonded_Inter_Types[i_bond_type].
		    ptr_g[D_b[bond_coeff_ctr + 1] - i_0] +=
		    Bp * get_g_dihedral_angle(gtemp);
	    }
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[D_b[bond_coeff_ctr] -
						      i_0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].
		ptr_g_cnt[D_b[bond_coeff_ctr + 1] - i_0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr] -
						  i_0] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr + 1] -
						  i_0] += 0.0;

	}

	return 2;		/* Number of coeff. to determine. */
    }

    return 0;			/* Number of coeff. to determine. */

}

/* END JFR */

/*****************************************************************************************
get_Bspline_dihedral_basis_vector(): Range is not checked since the R_min and R_max are
assumed to be -180 and 180. For now, everyting is set to zero if sin(phi) gets
too close to 0. When the 180 grid point is updated, then the index is switched to the 
-180 grid point since these are the same point. 
*****************************************************************************************/
int get_Bspline_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				      int *D_b, dvec * calG_b,
				      tW_Bonded_Inter * Bonded_Inter_Types,
				      tW_dihedral dihedral, bool i_flag)
{
    int i, j;
    int coeff_ctr = 0;

    int N_pts = Bonded_Inter_Types[i_bond_type].N_pts;
    double dr = Bonded_Inter_Types[i_bond_type].dr;
    double R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    double R_max = Bonded_Inter_Types[i_bond_type].R_max;
    int i_0 = Bonded_Inter_Types[i_bond_type].i_0;
    int k = Bonded_Inter_Types[i_bond_type].kspline;
    dvec grad_B;
    double B, Bp;
    double Nik[MAX_BSPLINE_COEFF];
    double Npik[MAX_BSPLINE_COEFF];

     /*JFR*/ dvec gtemp;

    copy_vector(dihedral.grad_B[dihedral.atom_position - 1], grad_B);

    if (fabs(dihedral.sin_angle) > DIHED_PREC) {

	i = get_grid_index_for_linear_basis(dihedral.angle, i_0, dr,
					    R_0) - i_0;

	B = calc_Bspline(dr, R_0, i_0, N_pts, dihedral.angle, i, k, Nik,
			 Npik);
	for (j = 0; j < k; j++) {
	    /* Account for periodicity */
	    if ((i - j + k / 2) < 0) {
		D_b[bond_coeff_ctr + coeff_ctr] =
		    (N_pts - 1) + i - j + k / 2 + i_0;
	    } else if ((i - j + k / 2) >= (N_pts - 1)) {
		D_b[bond_coeff_ctr + coeff_ctr] =
		    i - j + k / 2 - (N_pts - 1) + i_0;
	    } else {
		D_b[bond_coeff_ctr + coeff_ctr] = i - j + k / 2 + i_0;
	    }			/* k/2 shifts the supports up to the proper place, but this only works for even order Bsplines! */

	    B = Nik[k - 1 - j];
	    Bp = calc_Bspline_deriv(dr, R_0, i_0, N_pts, dihedral.angle, i,
				    j, k, Nik, Npik);

	    scal_times_vect(B / (dihedral.sin_angle), grad_B,
			    calG_b[bond_coeff_ctr + coeff_ctr]);

	    if (i_flag == TRUE) {
		scal_times_vect(1.0 / (dihedral.sin_angle), grad_B, gtemp);
		if (BSPLINE_STRUCT == 0) {
		    Bonded_Inter_Types[i_bond_type].
			ptr_g[D_b[bond_coeff_ctr + coeff_ctr] - i_0] +=
			B * get_g_dihedral_angle(gtemp);
		} else {
		    Bonded_Inter_Types[i_bond_type].
			ptr_g[D_b[bond_coeff_ctr + coeff_ctr] - i_0] +=
			Bp * get_g_dihedral_angle(gtemp);
		}
		Bonded_Inter_Types[i_bond_type].
		    ptr_g_cnt[D_b[bond_coeff_ctr + coeff_ctr] - i_0] +=
		    1.0;
		Bonded_Inter_Types[i_bond_type].
		    ptr_L[D_b[bond_coeff_ctr + coeff_ctr] - i_0] += 0.0;

	    }
	    coeff_ctr++;
	}
	return coeff_ctr;	/* Number of coeff. to determine. */
    }

    return 0;			/* Number of coeff. to determine. */

}

/*****************************************************************************************
get_Angle_basis_vectors(): 
*****************************************************************************************/
int get_Angle_basis_vectors(int i_site, int i_inter, int bond_coeff_ctr,
			    tW_CG_site * i_site_ptr,
			    tW_Bonded_Inter * Bonded_Inter_Types, int *D_b,
			    dvec * calG_b, tW_CG_site * CG_struct,
			    bool i_flag)
{
    int n_coeff;
    int i_bond_type;
    int i_basis_type;
    int *sites;

    i_bond_type = i_site_ptr->bond_type[i_inter];
    i_basis_type = Bonded_Inter_Types[i_bond_type].i_basis;
    sites = i_site_ptr->bond_site[i_inter];

    if (Bonded_Inter_Types[i_bond_type].i_basis == HARMONIC_BASIS_INDEX) {
	n_coeff =
	    get_harmonic_angle_basis_vector(i_site, bond_coeff_ctr,
					    i_bond_type, sites, D_b,
					    calG_b, Bonded_Inter_Types,
					    CG_struct, i_flag);
    } else if (Bonded_Inter_Types[i_bond_type].i_basis ==
	       DELTA_BASIS_INDEX) {
	n_coeff =
	    get_delta_angle_basis_vector(i_site, bond_coeff_ctr,
					 i_bond_type, sites, D_b, calG_b,
					 Bonded_Inter_Types, CG_struct,
					 i_flag);
    }
    /* START JFR */
    else if (Bonded_Inter_Types[i_bond_type].i_basis == LINEAR_BASIS_INDEX) {
	n_coeff =
	    get_linear_angle_basis_vector(i_site, bond_coeff_ctr,
					  i_bond_type, sites, D_b, calG_b,
					  Bonded_Inter_Types, CG_struct,
					  i_flag);
    }
    /* END JFR */
    else if (Bonded_Inter_Types[i_bond_type].i_basis ==
	     BSPLINE_BASIS_INDEX) {
	n_coeff =
	    get_Bspline_angle_basis_vector(i_site, bond_coeff_ctr,
					   i_bond_type, sites, D_b, calG_b,
					   Bonded_Inter_Types, CG_struct,
					   i_flag);
    } else {
	printf
	    ("\nERROR: Unsupported basis in get_Angle_basis_vectors().\n");
	exit(EXIT_FAILURE);
    }

    return n_coeff;
}


/*****************************************************************************************
get_Angle_info():
*****************************************************************************************/
int get_Angle_info(int i_site, int *sites, tW_CG_site * CG_struct,
		   dvec L_2, dvec L_3, double *cos_theta,
		   double *sin_theta, double *theta, double *norm_12,
		   double *norm_13)
{
    int atom_position;
    dvec r_12, r_13;
    dvec u_12, u_13;

    atom_position = get_atom_position_angle(sites, i_site);
    get_difference_unit_vector(CG_struct[sites[1]].r,
			       CG_struct[sites[0]].r, r_12, u_12, norm_12);
    get_difference_unit_vector(CG_struct[sites[1]].r,
			       CG_struct[sites[2]].r, r_13, u_13, norm_13);
    *cos_theta = dot_prod(u_12, u_13);
    *theta = acos(*cos_theta);
    *sin_theta = sin(*theta);
    get_L2_L3_angles(u_12, u_13, *cos_theta, *norm_12, *norm_13, L_2, L_3);

    return atom_position;
}


/*****************************************************************************************
get_harmonic_angle_basis_vector():
*****************************************************************************************/
int get_harmonic_angle_basis_vector(int i_site, int bond_coeff_ctr,
				    int i_bond_type, int *sites, int *D_b,
				    dvec * calG_b,
				    tW_Bonded_Inter * Bonded_Inter_Types,
				    tW_CG_site * CG_struct, bool i_flag)
{
    int atom_position;
    double norm_12, norm_13;
    double cos_theta, sin_theta, theta;
    double g_temp, L_temp;
    dvec L_2, L_3;
    dvec sum;

    atom_position = get_Angle_info(i_site, sites, CG_struct, L_2, L_3,
				   &cos_theta, &sin_theta, &theta,
				   &norm_12, &norm_13);

    D_b[bond_coeff_ctr] = Bonded_Inter_Types[i_bond_type].i_0;
    D_b[bond_coeff_ctr + 1] = Bonded_Inter_Types[i_bond_type].i_0 + 1;

    /* avoid dividing by zero when sin_theta = 0: JFR -> I added this 08.24.10 */
    if (fabs(sin_theta) < FLOAT_EPS) {
	if (sin_theta > 0.0) {
	    sin_theta += FLOAT_EPS;
	}
	if (sin_theta < 0.0) {
	    sin_theta -= FLOAT_EPS;
	}
//DEBUG JFR
//printf( " sin_theta = %lf \n", sin_theta );
    }

    if (atom_position == 0) {	/* end atom */
	scal_times_vect(-1.0 * theta / (sin_theta), L_2,
			calG_b[bond_coeff_ctr]);
	scal_times_vect(1.0 / (sin_theta), L_2,
			calG_b[bond_coeff_ctr + 1]);
    } else if (atom_position == 1) {	/* central atom */
	vect_sum(L_2, L_3, sum);
	scal_times_vect(theta / (sin_theta), sum, calG_b[bond_coeff_ctr]);
	scal_times_vect(-1.0 / (sin_theta), sum,
			calG_b[bond_coeff_ctr + 1]);
    } else if (atom_position == 2) {	/* end atom */
	scal_times_vect(-1.0 * theta / (sin_theta), L_3,
			calG_b[bond_coeff_ctr]);
	scal_times_vect(1.0 / (sin_theta), L_3,
			calG_b[bond_coeff_ctr + 1]);
    }

    if (atom_position == 0) {

	if (i_flag == TRUE) {
	    g_temp = get_g_theta_angle(cos_theta, norm_12, norm_13);
	    L_temp =
		get_laplacian_theta_angle(cos_theta, sin_theta, norm_12,
					  norm_13);
//DEBUG JFR
//if( L_temp > 1e5 )
//{
//  printf( " sin_theta = %lf \n", sin_theta );
//}
	    Bonded_Inter_Types[i_bond_type].ptr_g[0] += -1.0 * g_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_g[1] += 0.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[1] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[0] += -theta * L_temp;
	    Bonded_Inter_Types[i_bond_type].ptr_L[1] += L_temp;
	}
    }

    return 2;			/* Number of coeff. to determine. */
}


/*****************************************************************************************
get_delta_angle_basis_vector():
*****************************************************************************************/
int get_delta_angle_basis_vector(int i_site, int bond_coeff_ctr,
				 int i_bond_type, int *sites, int *D_b,
				 dvec * calG_b,
				 tW_Bonded_Inter * Bonded_Inter_Types,
				 tW_CG_site * CG_struct, bool i_flag)
{
    int atom_position;
    int i_0;
    double norm_12, norm_13;
    double cos_theta, sin_theta, theta;
    double g_temp, L_temp;
    double dr, R_0, R_max;
    dvec L_2, L_3;
    dvec sum;
    dvec zero = { 0.0, 0.0, 0.0 };

    atom_position = get_Angle_info(i_site, sites, CG_struct, L_2, L_3,
				   &cos_theta, &sin_theta, &theta,
				   &norm_12, &norm_13);

    dr = Bonded_Inter_Types[i_bond_type].dr;
    R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    R_max = Bonded_Inter_Types[i_bond_type].R_max;
    i_0 = Bonded_Inter_Types[i_bond_type].i_0;

    /* avoid dividing by zero when sin_theta = 0: JFR -> I added this 08.24.10 */
    if (fabs(sin_theta) < FLOAT_EPS) {
	if (sin_theta > 0.0) {
	    sin_theta += FLOAT_EPS;
	}
	if (sin_theta < 0.0) {
	    sin_theta -= FLOAT_EPS;
	}
    }

    /* Test if angle is within the range */
    if ((theta >= R_0) && (theta < (R_max + 0.5 * dr))) {
	if (atom_position == 0) {	/* end atom (-L_2) */
	    scal_times_vect(1. / (sin_theta), L_2, calG_b[bond_coeff_ctr]);
	} else if (atom_position == 1) {	/* central atom (L_2 + L_3) */
	    vect_sum(L_2, L_3, sum);
	    scal_times_vect(-1. / (sin_theta), sum,
			    calG_b[bond_coeff_ctr]);
	} else if (atom_position == 2) {	/* end atom (-L_3) */
	    scal_times_vect(1. / (sin_theta), L_3, calG_b[bond_coeff_ctr]);
	}

	D_b[bond_coeff_ctr] =
	    get_grid_index_for_delta_basis(theta, i_0, dr, R_0);

	/* Only need to evaluate this once for each angle in this form */
	if (atom_position == 0) {
	    if (i_flag == TRUE) {
		Bonded_Inter_Types[i_bond_type].ptr_g[D_b[bond_coeff_ctr] -
						      i_0]
		    += get_g_theta_angle(cos_theta, norm_12, norm_13);
		Bonded_Inter_Types[i_bond_type].
		    ptr_g_cnt[D_b[bond_coeff_ctr] - i_0]
		    += 1.0;
		Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr] -
						      i_0]
		    += get_laplacian_theta_angle(cos_theta, sin_theta,
						 norm_12, norm_13);
	    }
	}
    }
    /* Not sure if this is needed. */
    else {
	return 0;
    }

    return 1;			/* Number of coeff. to determine. */
}

/* START JFR */
/*****************************************************************************************
get_linear_angle_basis_vector():
*****************************************************************************************/
int get_linear_angle_basis_vector(int i_site, int bond_coeff_ctr,
				  int i_bond_type, int *sites, int *D_b,
				  dvec * calG_b,
				  tW_Bonded_Inter * Bonded_Inter_Types,
				  tW_CG_site * CG_struct, bool i_flag)
{
    int atom_position;
    int i_0, N_pts;
    double norm_12, norm_13;
    double cos_theta, sin_theta, theta;
    double g_temp, L_temp;
    double dr, R_0, R_max, A, B, Ap, Bp;
    dvec L_2, L_3;
    dvec sum;

    atom_position = get_Angle_info(i_site, sites, CG_struct, L_2, L_3,
				   &cos_theta, &sin_theta, &theta,
				   &norm_12, &norm_13);

    N_pts = Bonded_Inter_Types[i_bond_type].N_pts;
    dr = Bonded_Inter_Types[i_bond_type].dr;
    R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    R_max = Bonded_Inter_Types[i_bond_type].R_max;
    i_0 = Bonded_Inter_Types[i_bond_type].i_0;

    /* avoid dividing by zero when sin_theta = 0: JFR -> I added this 08.24.10 */
    if (fabs(sin_theta) < FLOAT_EPS) {
	if (sin_theta > 0.0) {
	    sin_theta += FLOAT_EPS;
	}
	if (sin_theta < 0.0) {
	    sin_theta -= FLOAT_EPS;
	}
    }

    /* Test if angle is within the range */
    if ((theta >= R_0) && (theta < R_max)) {
	A = calc_linear_spline_A(bond_coeff_ctr, D_b, dr, R_0, i_0, N_pts,
				 theta, &Ap, &Bp, FALSE);
	B = 1 - A;
	/* A and B are weighting factors which depend on where the angle is relative to the two nearest gridpoints */

	if (atom_position == 0) {	/* end atom (-L_2) */
	    scal_times_vect(A / (sin(theta)), L_2, calG_b[bond_coeff_ctr]);
	    scal_times_vect(B / (sin(theta)), L_2,
			    calG_b[bond_coeff_ctr + 1]);
	} else if (atom_position == 1) {	/* central atom (L_2 + L_3) */
	    vect_sum(L_2, L_3, sum);
	    scal_times_vect(-A / (sin(theta)), sum,
			    calG_b[bond_coeff_ctr]);
	    scal_times_vect(-B / (sin(theta)), sum,
			    calG_b[bond_coeff_ctr + 1]);
	} else if (atom_position == 2) {	/* end atom (-L_3) */
	    scal_times_vect(A / (sin(theta)), L_3, calG_b[bond_coeff_ctr]);
	    scal_times_vect(B / (sin(theta)), L_3,
			    calG_b[bond_coeff_ctr + 1]);
	} else {
	    printf("ERROR: atom position not valid!");
	    printf("See function: get_linear_angle_basis_vector");
	    exit(0);
	}

	/* Only need to evaluate this once for each angle in this form */
	if (atom_position == 0) {
	    if (i_flag == TRUE) {
		if (LINEAR_STRUCT == 0) {
		    Bonded_Inter_Types[i_bond_type].
			ptr_g[D_b[bond_coeff_ctr] - i_0]
			+= A * get_g_theta_angle(cos_theta, norm_12,
						 norm_13);
		    Bonded_Inter_Types[i_bond_type].
			ptr_g[D_b[bond_coeff_ctr + 1] - i_0]
			+= B * get_g_theta_angle(cos_theta, norm_12,
						 norm_13);
		} else {
		    Bonded_Inter_Types[i_bond_type].
			ptr_g[D_b[bond_coeff_ctr] - i_0]
			+= Ap * get_g_theta_angle(cos_theta, norm_12,
						  norm_13);
		    Bonded_Inter_Types[i_bond_type].
			ptr_g[D_b[bond_coeff_ctr + 1] - i_0]
			+= Bp * get_g_theta_angle(cos_theta, norm_12,
						  norm_13);
		}
		Bonded_Inter_Types[i_bond_type].
		    ptr_g_cnt[D_b[bond_coeff_ctr] - i_0]
		    += 1.0;
		Bonded_Inter_Types[i_bond_type].
		    ptr_g_cnt[D_b[bond_coeff_ctr + 1] - i_0]
		    += 1.0;
		Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr] -
						      i_0]
		    += A * get_laplacian_theta_angle(cos_theta, sin_theta,
						     norm_12, norm_13);
		Bonded_Inter_Types[i_bond_type].
		    ptr_L[D_b[bond_coeff_ctr + 1] - i_0]
		    += B * get_laplacian_theta_angle(cos_theta, sin_theta,
						     norm_12, norm_13);

	    }
	}
	return 2;		/* Number of coeff. to determine. */
    }

    return 0;			/* Number of coeff. to determine. */
}

/* END JFR */

/*****************************************************************************************
get_Bspline_angle_basis_vector(): JFR - 07.23.12
*****************************************************************************************/
int get_Bspline_angle_basis_vector(int i_site, int bond_coeff_ctr,
				   int i_bond_type, int *sites, int *D_b,
				   dvec * calG_b,
				   tW_Bonded_Inter * Bonded_Inter_Types,
				   tW_CG_site * CG_struct, bool i_flag)
{
    int atom_position;
    int i_0, N_pts;
    double norm_12, norm_13;
    double cos_theta, sin_theta, theta;
    double g_temp, L_temp;
    double dr, R_0, R_max, B, Bp;
    dvec L_2, L_3;
    dvec sum;

    int i, j, k;
    int coeff_ctr = 0;

    double Nik[MAX_BSPLINE_COEFF];
    double Npik[MAX_BSPLINE_COEFF];

    atom_position = get_Angle_info(i_site, sites, CG_struct, L_2, L_3,
				   &cos_theta, &sin_theta, &theta,
				   &norm_12, &norm_13);

    N_pts = Bonded_Inter_Types[i_bond_type].N_pts;
    dr = Bonded_Inter_Types[i_bond_type].dr;
    R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    R_max = Bonded_Inter_Types[i_bond_type].R_max;
    i_0 = Bonded_Inter_Types[i_bond_type].i_0;
    k = Bonded_Inter_Types[i_bond_type].kspline;

    /* avoid dividing by zero when sin_theta = 0: JFR -> I added this 08.24.10 */
    if (fabs(sin_theta) < FLOAT_EPS) {
	if (sin_theta > 0.0) {
	    sin_theta += FLOAT_EPS;
	}
	if (sin_theta < 0.0) {
	    sin_theta -= FLOAT_EPS;
	}
    }

    /* Test if angle is within the range */
    if ((theta >= R_0) && (theta < R_max)) {

	i = get_grid_index_for_linear_basis(theta, i_0, dr, R_0) - i_0;

	B = calc_Bspline(dr, R_0, i_0, N_pts, theta, i, k, Nik, Npik);
	for (j = 0; j < k; j++) {
	    D_b[bond_coeff_ctr + coeff_ctr] = i - j + k / 2 + i_0;	/* k/2 shifts the supports up to the proper place, but this only works for even order Bsplines! */

	    if (((i - j + k / 2) < 0) || ((i - j + k / 2) >= N_pts - 1)) {
		continue;
	    }
	    /* You are beyond the grid, there are no supports here */
	    B = Nik[k - 1 - j];
	    Bp = calc_Bspline_deriv(dr, R_0, i_0, N_pts, theta, i, j, k,
				    Nik, Npik);

	    /* A and B are weighting factors which depend on where the angle is relative to the two nearest gridpoints */
	    if (atom_position == 0) {	/* end atom (-L_2) */
		scal_times_vect(B / (sin(theta)), L_2,
				calG_b[bond_coeff_ctr + coeff_ctr]);
	    } else if (atom_position == 1) {	/* central atom (L_2 + L_3) */
		vect_sum(L_2, L_3, sum);
		scal_times_vect(-B / (sin(theta)), sum,
				calG_b[bond_coeff_ctr + coeff_ctr]);
	    } else if (atom_position == 2) {	/* end atom (-L_3) */
		scal_times_vect(B / (sin(theta)), L_3,
				calG_b[bond_coeff_ctr + coeff_ctr]);
	    } else {
		printf("ERROR: atom position not valid!");
		printf("See function: get_linear_angle_basis_vector");
		exit(0);
	    }

	    /* Only need to evaluate this once for each angle in this form */
	    if (atom_position == 0) {
		if (i_flag == TRUE) {
		    if (BSPLINE_STRUCT == 0) {
			Bonded_Inter_Types[i_bond_type].
			    ptr_g[D_b[bond_coeff_ctr + coeff_ctr] - i_0]
			    += B * get_g_theta_angle(cos_theta, norm_12,
						     norm_13);
		    } else {
			Bonded_Inter_Types[i_bond_type].
			    ptr_g[D_b[bond_coeff_ctr + coeff_ctr] - i_0]
			    += Bp * get_g_theta_angle(cos_theta, norm_12,
						      norm_13);
		    }
		    Bonded_Inter_Types[i_bond_type].
			ptr_g_cnt[D_b[bond_coeff_ctr + coeff_ctr] - i_0]
			+= 1.0;
		    Bonded_Inter_Types[i_bond_type].
			ptr_L[D_b[bond_coeff_ctr + coeff_ctr] - i_0]
			+= B * get_laplacian_theta_angle(cos_theta,
							 sin_theta,
							 norm_12, norm_13);
		}
	    }
	    coeff_ctr++;
	}
	return coeff_ctr;	/* Number of coeff. to determine. */
    }

    return 0;			/* Number of coeff. to determine. */
}

/*****************************************************************************************
get_laplacian_theta_angle(): 
*****************************************************************************************/
double get_laplacian_theta_angle(double cos_theta, double sin_theta,
				 double r_12, double r_13)
{
    double R = (r_12 / r_13) + (r_13 / r_12);
    double C = 2.0 / (r_12 * r_13 * sin_theta);
    return (C * (R * cos_theta - 1.0));
}


/*****************************************************************************************
get_g_theta_angle():
*****************************************************************************************/
double get_g_theta_angle(double cos_theta, double r_12, double r_13)
{
    double R = (r_12 / r_13) + (r_13 / r_12);

    return ((2.0 / (r_12 * r_13)) * (R - cos_theta));
}


/*****************************************************************************************
get_L2_L3_angles():
*****************************************************************************************/
int get_L2_L3_angles(dvec u_12, dvec u_13, double cos_theta,
		     double norm_12, double norm_13, dvec L_2, dvec L_3)
{
    dvec temp_12, temp_13, temp_A, temp_B;

    scal_times_vect(cos_theta, u_12, temp_12);
    vect_diff(u_13, temp_12, temp_A);
    scal_times_vect(1 / norm_12, temp_A, L_2);

    scal_times_vect(cos_theta, u_13, temp_13);
    vect_diff(u_12, temp_13, temp_B);
    scal_times_vect(1 / norm_13, temp_B, L_3);

    return 0;
}


/*****************************************************************************************
get_BondStretch_basis_vectors():
*****************************************************************************************/
int get_BondStretch_basis_vectors(int i_site, int i_inter,
				  int bond_coeff_ctr,
				  tW_CG_site * i_site_ptr,
				  tW_Bonded_Inter * Bonded_Inter_Types,
				  int *D_b, dvec * calG_b,
				  tW_CG_site * CG_struct, bool i_flag)
{
    int n_coeff;
    int i_bond_type;
    int i_basis_type;
    int *sites;
    double norm_ij;
    dvec u_ij;

    i_bond_type = i_site_ptr->bond_type[i_inter];
    i_basis_type = Bonded_Inter_Types[i_bond_type].i_basis;
    sites = i_site_ptr->bond_site[i_inter];

    norm_ij = get_BondStretch_info(i_site, sites, CG_struct, u_ij);

    if (Bonded_Inter_Types[i_bond_type].i_basis == HARMONIC_BASIS_INDEX) {
	n_coeff =
	    get_harmonic_bond_basis_vector(bond_coeff_ctr, i_bond_type,
					   D_b, calG_b, Bonded_Inter_Types,
					   norm_ij, u_ij, i_flag);
    } else if (Bonded_Inter_Types[i_bond_type].i_basis ==
	       DELTA_BASIS_INDEX) {
	n_coeff =
	    get_delta_bond_basis_vector(bond_coeff_ctr, i_bond_type, D_b,
					calG_b, Bonded_Inter_Types,
					norm_ij, u_ij, i_flag);
    }
    /* START JFR */
    else if (Bonded_Inter_Types[i_bond_type].i_basis == LINEAR_BASIS_INDEX) {
	n_coeff =
	    get_linear_bond_basis_vector(bond_coeff_ctr, i_bond_type, D_b,
					 calG_b, Bonded_Inter_Types,
					 norm_ij, u_ij, i_flag);
    }
    /* END JFR */
    else if (Bonded_Inter_Types[i_bond_type].i_basis == BSPLINE_BASIS_INDEX) {	/* JFR - 07.23.12: Bspline */
	n_coeff =
	    get_Bspline_bond_basis_vector(bond_coeff_ctr, i_bond_type, D_b,
					  calG_b, Bonded_Inter_Types,
					  norm_ij, u_ij, i_flag);
    } else {
	printf
	    ("\nERROR: Unsupported basis in get_BondStretch_basis_vectors().\n");
	exit(EXIT_FAILURE);
    }

    return n_coeff;
}


/*****************************************************************************************
get_BondStretch_info(): Calculates and stores data needed to calculated the BondStretch
basis vectors.
*****************************************************************************************/
double get_BondStretch_info(int i_site, int *sites, tW_CG_site * CG_struct,
			    dvec u_ij)
{
    int j_site;
    double norm_ij;
    dvec r_ij;

    if (i_site == sites[0]) {
	j_site = sites[1];
    } else {
	j_site = sites[0];
	if (i_site != sites[1]) {
	    printf("\nError in get_BondStretch_info().\n");
	    exit(EXIT_FAILURE);
	}
    }

    get_difference_unit_vector(CG_struct[i_site].r, CG_struct[j_site].r,
			       r_ij, u_ij, &norm_ij);

    return norm_ij;
}


/*****************************************************************************************
get_delta_bond_basis_vector():
*****************************************************************************************/
int get_delta_bond_basis_vector(int bond_coeff_ctr, int i_bond_type,
				int *D_b, dvec * calG_b,
				tW_Bonded_Inter * Bonded_Inter_Types,
				double norm_ij, dvec u_ij, bool i_flag)
{
    int i_0;
    double dr, R_0, R_max;
    dvec zero = { 0.0, 0.0, 0.0 };

    dr = Bonded_Inter_Types[i_bond_type].dr;
    R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    R_max = Bonded_Inter_Types[i_bond_type].R_max;
    i_0 = Bonded_Inter_Types[i_bond_type].i_0;

    if ((norm_ij >= R_0) && (norm_ij < (R_max + 0.5 * dr))) {
	D_b[bond_coeff_ctr] =
	    get_grid_index_for_delta_basis(norm_ij, i_0, dr, R_0);
	copy_vector(u_ij, calG_b[bond_coeff_ctr]);

	if (i_flag == TRUE) {
	    Bonded_Inter_Types[i_bond_type].ptr_g[D_b[bond_coeff_ctr] -
						  i_0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[D_b[bond_coeff_ctr] -
						      i_0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr] -
						  i_0] += 2.0 / norm_ij;
	}
    } else {
	return 0;
    }

    return 1;			/* Number of coeff. to determine. */
}

/* START JFR */
/*****************************************************************************************
get_linear_bond_basis_vector(): 
*****************************************************************************************/
int get_linear_bond_basis_vector(int bond_coeff_ctr, int i_bond_type,
				 int *D_b, dvec * calG_b,
				 tW_Bonded_Inter * Bonded_Inter_Types,
				 double norm_ij, dvec u_ij, bool i_flag)
{
    int i_0, N_pts;
    double dr, R_0, R_max, A, B, Ap, Bp;
    dvec u_temp;

    N_pts = Bonded_Inter_Types[i_bond_type].N_pts;
    dr = Bonded_Inter_Types[i_bond_type].dr;
    R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    R_max = Bonded_Inter_Types[i_bond_type].R_max;
    i_0 = Bonded_Inter_Types[i_bond_type].i_0;

    /* Test if bond is within the range */
    if ((norm_ij >= R_0) && (norm_ij < R_max)) {
	A = calc_linear_spline_A(bond_coeff_ctr, D_b, dr, R_0, i_0, N_pts,
				 norm_ij, &Ap, &Bp, FALSE);
	B = 1.0 - A;
	/* A and B are weighting factors which depend on where the distance 
	   is relative to the two nearest gridpoints */

	scal_times_vect(A, u_ij, u_temp);
	copy_vector(u_temp, calG_b[bond_coeff_ctr]);

	scal_times_vect(B, u_ij, u_temp);
	copy_vector(u_temp, calG_b[bond_coeff_ctr + 1]);

	if (i_flag == TRUE) {
	    if (LINEAR_STRUCT == 0) {
		Bonded_Inter_Types[i_bond_type].ptr_g[D_b[bond_coeff_ctr] -
						      i_0] += A;
		Bonded_Inter_Types[i_bond_type].
		    ptr_g[D_b[bond_coeff_ctr + 1] - i_0] += B;
	    } else {
		Bonded_Inter_Types[i_bond_type].ptr_g[D_b[bond_coeff_ctr] -
						      i_0] += Ap;
		Bonded_Inter_Types[i_bond_type].
		    ptr_g[D_b[bond_coeff_ctr + 1] - i_0] += Bp;
	    }
	    Bonded_Inter_Types[i_bond_type].ptr_g_cnt[D_b[bond_coeff_ctr] -
						      i_0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].
		ptr_g_cnt[D_b[bond_coeff_ctr + 1] - i_0] += 1.0;
	    Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr] -
						  i_0] +=
		(A * 2.0) / norm_ij;
	    Bonded_Inter_Types[i_bond_type].ptr_L[D_b[bond_coeff_ctr + 1] -
						  i_0] +=
		(B * 2.0) / norm_ij;
	}

	return 2;		/* Number of coeff. to determine. */

    }

    return 0;			/* Number of coeff. to determine. */

}

/* END JFR */

/* START JFR */
/*****************************************************************************************
calc_linear_spline_A(): 
*****************************************************************************************/
double calc_linear_spline_A(int bond_coeff_ctr, int *D_b, double dr,
			    double R_0, int i_0, int N_pts, double R,
			    double *Ap, double *Bp, int periodic)
{
    double A, B, r1, r2;

    D_b[bond_coeff_ctr] = get_grid_index_for_linear_basis(R, i_0, dr, R_0);

    D_b[bond_coeff_ctr + 1] = D_b[bond_coeff_ctr] + 1;

    if ((periodic == TRUE)
	&& (D_b[bond_coeff_ctr + 1] == (i_0 + N_pts - 1))) {
	D_b[bond_coeff_ctr + 1] = i_0;
    }

    r1 = dr * (D_b[bond_coeff_ctr] - i_0);
    r2 = dr * (D_b[bond_coeff_ctr] + 1 - i_0);

    A = (r2 - (R - R_0)) / (r2 - r1);
    B = 1.0 - A;
    /* A and B are weighting factors which depend on where the distance 
       is relative to the two nearest gridpoints */

    *Bp = 1.0 / (r2 - r1);
    *Ap = -1.0 * (*Bp);
    /* Ap and Bp are the derivatives of A and B */

    //DEBUG JFR
    //if( A < 0.0 || A > 1.0 || B < -0.0000 || B >1.0 )
    //{
    //printf("ERROR:  A, B not in [0,1]\n");
    //printf(" A=%lf; B=%lf\n", A, B);
    //printf(" D_b[%d]=%d\n",bond_coeff_ctr, D_b[bond_coeff_ctr]);
    //printf(" r1=%lf; r2=%lf\n", r1, r2);
    //printf(" R=%lf\n", R);
    //printf(" dr=%lf\n", dr);
    //exit(0);
    //}

    return A;

}

/* END JFR */

/*****************************************************************************************
get_Bspline_bond_basis_vector(): 
*****************************************************************************************/
int get_Bspline_bond_basis_vector(int bond_coeff_ctr, int i_bond_type,
				  int *D_b, dvec * calG_b,
				  tW_Bonded_Inter * Bonded_Inter_Types,
				  double norm_ij, dvec u_ij, bool i_flag)
{
    int i_0, N_pts;
    double dr, R_0, R_max, B, Bp;
    dvec u_temp;

    int i, j, k;
    int coeff_ctr = 0;

    double Nik[MAX_BSPLINE_COEFF];
    double Npik[MAX_BSPLINE_COEFF];

    N_pts = Bonded_Inter_Types[i_bond_type].N_pts;
    dr = Bonded_Inter_Types[i_bond_type].dr;
    R_0 = Bonded_Inter_Types[i_bond_type].R_0;
    R_max = Bonded_Inter_Types[i_bond_type].R_max;
    i_0 = Bonded_Inter_Types[i_bond_type].i_0;
    k = Bonded_Inter_Types[i_bond_type].kspline;

    /* Test if bond is within the range */
    if ((norm_ij >= R_0) && (norm_ij < R_max)) {

	i = get_grid_index_for_linear_basis(norm_ij, i_0, dr, R_0) - i_0;

	B = calc_Bspline(dr, R_0, i_0, N_pts, norm_ij, i, k, Nik, Npik);
	for (j = 0; j < k; j++) {
	    D_b[bond_coeff_ctr + coeff_ctr] = i - j + k / 2 + i_0;	/* k/2 shifts the supports up to the proper place, but this only works for even order Bsplines! */

	    if (((i - j + k / 2) < 0) || ((i - j + k / 2) >= N_pts - 1)) {
		continue;
	    }
	    /* You are beyond the grid, there are no supports here */
	    B = Nik[k - 1 - j];
	    Bp = calc_Bspline_deriv(dr, R_0, i_0, N_pts, norm_ij, i, j, k,
				    Nik, Npik);

	    scal_times_vect(B, u_ij, u_temp);
	    copy_vector(u_temp, calG_b[bond_coeff_ctr + coeff_ctr]);

	    if (i_flag == TRUE) {
		if (BSPLINE_STRUCT == 0) {
		    Bonded_Inter_Types[i_bond_type].
			ptr_g[D_b[bond_coeff_ctr + coeff_ctr] - i_0] += B;
		} else {
		    Bonded_Inter_Types[i_bond_type].
			ptr_g[D_b[bond_coeff_ctr + coeff_ctr] - i_0] += Bp;
		}
		Bonded_Inter_Types[i_bond_type].
		    ptr_g_cnt[D_b[bond_coeff_ctr + coeff_ctr] - i_0] +=
		    1.0;
		Bonded_Inter_Types[i_bond_type].
		    ptr_L[D_b[bond_coeff_ctr + coeff_ctr] - i_0] +=
		    (B * 2.0) / norm_ij;
	    }
	    coeff_ctr++;
	}
	return coeff_ctr;	/* Number of coeff. to determine. */

    }

    return 0;			/* Number of coeff. to determine. */

}

/************************************************************************************************************
calc_Bspline(): JFR - 07.16.12: This is the De Boor iterative method for calculating normalized Bsplines.
*************************************************************************************************************/
double calc_Bspline(double dr, double R_0, int i_0, int N_pts, double R,
		    int i, int k, double *Nik, double *Npik)
{
    int r, s;
    double DP, DM;
    double M;
    double Nik_tmp[MAX_BSPLINE_COEFF][MAX_BSPLINE_COEFF];

    /* De Boor's algorithm */
    Npik[0] = 0.00;
    Nik_tmp[0][0] = 1.0;
    for (s = 0; s < k - 1; s++) {
	Nik_tmp[0][s + 1] = 0.0;
	for (r = 0; r <= s; r++) {
	    DP = R_0 + (i + 1 + r) * dr - R;
	    DM = R - (R_0 + (i - s + r) * dr);
	    M = Nik_tmp[r][s] / (DP + DM);
	    Nik_tmp[r][s + 1] += DP * M;
	    Nik_tmp[r + 1][s + 1] = DM * M;
	    if (s == k - 3) {
		Npik[r + 1] = Nik_tmp[r][s + 1];
		if (r == s) {
		    Npik[r + 2] = Nik_tmp[r + 1][s + 1];
		}
	    }
	    if (s == k - 2) {
		Nik[r] = Nik_tmp[r][s + 1];
		if (r == s) {
		    Nik[r + 1] = Nik_tmp[r + 1][s + 1];
		}
	    }
	}
    }

    return Nik[k - 1];

}

/************************************************************************************************************
calc_Bspline_deriv(): JFR - 07.16.12: This is the De Boor iterative method for calculating normalized Bsplines.
*************************************************************************************************************/
double calc_Bspline_deriv(double dr, double R_0, int i_0, int N_pts,
			  double R, int i, int j, int k, double *Nik,
			  double *Npik)
{
    double DP;
    double Bp = 0.00;

//     Bp = Npik[k-1-j];
/* Try calculating the derivative here */
//if ( j == k-1 ) 
//{ 
//Bp = (R_0 + (i - j + k/2 + k)*dr - r_ij)*( (1.0-((double)k))*Nik[k-1-j][k-1] ); 
//printf( "Bp = %lf , Nik[k-1-j][k-1] = %lf \n", Bp, Nik[k-1-j][k-1] );
//}
//else 
//{ 
//Bp = (R_0 + (i - j + k/2 + k)*dr - r_ij)*( (1.0-((double)k))*Nik[k-1-j][k-1] + ((double)k)*Nik[k-2-j][k-2] ); 
//printf( "Bp = %lf , Nik[k-1-j][k-1] = %lf, Nik[k-2-j][k-2] = %lf \n", Bp, Nik[k-1-j][k-1], Nik[k-2-j][k-2] );
//}
//if ( j == k-1 ) 
//{
//Bp = ( ((double)k) / (r_ij - (R_0 + (i - j + k/2 + k)*dr)) )*Nik[k-1-j][k-1];
//}
//else 
//{ 
//Bp = ( ((double)k) / (r_ij - (R_0 + (i - j + k/2 + k)*dr)) )*Nik[k-1-j][k-1] +  (1.0/dr)*(r_ij - (R_0 + (i - j + k/2 + k)*dr) + ((double)k)*((double)k)*dr)*Nik[k-2-j][k-2]/( (1.0 -((double)k))*(R_0 + (i - j + k/2 + k)*dr - r_ij) ); 
//}
//if ( j == k-1 )
//{
//Bp = (1.0 / ((R_0 + (i - j + k/2 + k)*dr) - r_ij))*( (1.0-((double)k))*Nik[k-1-j][k-1] );
//}
//else
//{
//Bp = (1.0 / ((R_0 + (i - j + k/2 + k)*dr) - r_ij))*( (1.0-((double)k))*Nik[k-1-j][k-1] + (2.0*((double)k)-1.0)*Nik[k-2-j][k-2] );
//}
    /* MRD 02.05.2019 got rid of the k / 2
    DP = 1.0 / (R_0 + (i - j + k / 2 + k) * dr - R); */
    DP = 1.0 / (R_0 + (i - j + k) * dr - R);
    if (j == k - 1) {
//Bp = DP*(1.0-((double)k))*Nik[k-1-j][k-1];
	Bp = DP * (1.0 - ((double) k)) * Nik[k - 1 - j];
    } else {
//Bp = DP*( (1.0-((double)k))*Nik[k-1-j][k-1] + ((double)k)*Nik[k-2-j][k-2] );
	Bp = DP * ((1.0 - ((double) k)) * Nik[k - 1 - j] + ((double) k) * Npik[k - 1 - j]);
		   /*((double) k) * Npik[k - 2 - j]); MRD 02.05.2019 changed the 2 to a 1*/
    }

    return Bp;

}

/*****************************************************************************************
get_harmonic_bond_basis_vector():
*****************************************************************************************/
int get_harmonic_bond_basis_vector(int bond_coeff_ctr, int i_bond_type,
				   int *D_b, dvec * calG_b,
				   tW_Bonded_Inter * Bonded_Inter_Types,
				   double norm_ij, dvec u_ij, bool i_flag)
{

    D_b[bond_coeff_ctr] = Bonded_Inter_Types[i_bond_type].i_0;
    D_b[bond_coeff_ctr + 1] = Bonded_Inter_Types[i_bond_type].i_0 + 1;

    scal_times_vect(-1.0 * norm_ij, u_ij, calG_b[bond_coeff_ctr]);
    copy_vector(u_ij, calG_b[bond_coeff_ctr + 1]);

    if (i_flag == TRUE) {
	Bonded_Inter_Types[i_bond_type].ptr_g[0] += -1.0;
	Bonded_Inter_Types[i_bond_type].ptr_g[1] += 0.0;
	Bonded_Inter_Types[i_bond_type].ptr_g_cnt[0] += 1.0;
	Bonded_Inter_Types[i_bond_type].ptr_g_cnt[1] += 1.0;
	Bonded_Inter_Types[i_bond_type].ptr_L[0] += -2.0;
	Bonded_Inter_Types[i_bond_type].ptr_L[1] += 2.0 / norm_ij;
    }

    return 2;			/* Number of coeff. to determine. */
}


/*****************************************************************************************
print_debug_eval_bond_basis():
*****************************************************************************************/
void print_debug_eval_bond_basis(tW_CG_site * i_site_ptr,
				 tW_Bonded_Inter * Bonded_Inter_Types,
				 int i, FILE * fp)
{
    int j;
    int N_Int_Sites;

    fprintf(fp, "  BondInter: %d of %d  %s  Sites: ",
	    i + 1, i_site_ptr->nr_bonds,
	    Bonded_Inter_Types[i_site_ptr->bond_type[i]].name);
    fflush(fp);
    N_Int_Sites = Bonded_Inter_Types[i_site_ptr->bond_type[i]].N_Int_Sites;
    for (j = 0; j < N_Int_Sites; j++) {
	fprintf(fp, "%d ", i_site_ptr->bond_site[i][j]);
	fflush(fp);
    }
    fprintf(fp, "\n");
    fflush(fp);
}


/*****************************************************************************************
get_grid_index_for_delta_basis():
*****************************************************************************************/
int get_grid_index_for_delta_basis(double r, int i_0, double dr,
				   double R_0)
{
    int index, coeff;
    double resid;

    coeff = (int) floor((r - R_0) / dr);	/* ID last grid point before r.                            */
    resid = r - (R_0 + coeff * dr);	/* Distance from prev. grid point.                         */
    if (resid > 0.5 * dr) {
	coeff += 1;
    }				/* If r closer to next grid_point update the next grid pt. */
    index = i_0 + coeff;

    return index;
}

/* START JFR */
/*****************************************************************************************
get_grid_index_for_linear_basis():
*****************************************************************************************/
int get_grid_index_for_linear_basis(double r, int i_0, double dr,
				    double R_0)
{
    int index, coeff;

    coeff = (int) floor((r - R_0) / dr);	/* ID last grid point before r. */
    index = i_0 + coeff;	/* index = lower gridpoint */

    return index;
}

/* END JFR */

/*****************************************************************************************
get_atom_position_angle():
*****************************************************************************************/
int get_atom_position_angle(int angle_triple[], int i_site)
{
    int i;
    int atom_position = -2;

    for (i = 0; i < 3; i++) {
	if (i_site == angle_triple[i]) {
	    atom_position = i;
	    break;
	} else {
	    atom_position = -1;
	}
    }
    if (atom_position == -1) {
	printf("ERROR: position of site %d in angle not found.\n", i_site);
	exit(EXIT_FAILURE);
    }

    return atom_position;
}


/*****************************************************************************************
get_atom_position_dihedral():
*****************************************************************************************/
int get_atom_position_dihedral(int dihedral_quartet[], int i_site)
{
    int i;
    int atom_position;

    /* Get site i's position in the dihedral gromacs < 1234 > */
    for (i = 0; i < 4; i++) {
	if (i_site == dihedral_quartet[i]) {
	    atom_position = i + 1;
	    break;
	} else {
	    atom_position = -1;
	}
    }
    if (atom_position == -1) {
	printf("\nERROR: position of site %d in dihedral not found.\n",
	       i_site);
	exit(EXIT_FAILURE);
    }

    return atom_position;
}



/*****************************************************************************************
get_g_dihedral_angle():
*****************************************************************************************/
double get_g_dihedral_angle(dvec basis)
{
    double g_temp;

    g_temp = dot_prod(basis, basis);

    return g_temp;
}


/*****************************************************************************************
eval_bonds_M_b():
*****************************************************************************************/
int eval_bonds_M_b(tW_system * sys, int nr_bond_coeff, int D_b[],
		   dvec calG_b[], bool b_F, dvec f_i, dvec f_i_ref,
		   int *D_b_inter_index /* JFR - 02.26.13: for sep M2 */ )
{
    int i1_bond_coeff, i2_bond_coeff;
    int D_i1_bond, D_i2_bond;
    double inner_prod;
    int N = sys->N_coeff;
    int index;

    for (i1_bond_coeff = 0; i1_bond_coeff < nr_bond_coeff; i1_bond_coeff++) {
	D_i1_bond = D_b[i1_bond_coeff];

	/* Evaluate b */
	if (b_F) {
	    sys->b[D_i1_bond] += dot_prod(f_i, calG_b[i1_bond_coeff]);
	}
	if (sys->flag_ref_potential) {
	    sys->b_ref[D_i1_bond] +=
		dot_prod(f_i_ref, calG_b[i1_bond_coeff]);
	}

	/* If necessary, keep track of <b^2>, for variance calculation */
	//if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 )
	//{
	if (b_F) {
	    sys->d2b[D_i1_bond] +=
		dot_prod(f_i, calG_b[i1_bond_coeff]) * dot_prod(f_i,
								calG_b
								[i1_bond_coeff]);
	}
	//}

	/* Evaluate M */
	// JFR - added 04.11.12: put the matrix in packed form
	//for ( i2_bond_coeff=0; i2_bond_coeff<nr_bond_coeff; i2_bond_coeff++ )
	for (i2_bond_coeff = 0; i2_bond_coeff <= i1_bond_coeff;
	     i2_bond_coeff++) {
	    D_i2_bond = D_b[i2_bond_coeff];
	    inner_prod =
		dot_prod(calG_b[i1_bond_coeff], calG_b[i2_bond_coeff]);
	    // JFR - added 04.11.12: put the matrix in packed form
	    //sys->M[D_i1_bond][D_i2_bond] += inner_prod;
	    //sys->M_cnt[D_i1_bond][D_i2_bond] += 1.0;
	    index = index_Lpacked(D_i1_bond, D_i2_bond, N);
	    if (D_b_inter_index[i1_bond_coeff] ==
		D_b_inter_index[i2_bond_coeff]
		/* JFR - 02.26.13: for sep M2 */ ) {
		sys->M2[index] += inner_prod;
	    } else {
		sys->M[index] += inner_prod;
	    }
	    if (sys->M_cnt != NULL) {
		sys->M_cnt[index] += 1.0;
	    }			//JFR - added 04.13.12: check if in LOWMEM mode

	    /* If necessary, keep track of <M^2>, for variance calculation */
	    //if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 )
	    //{
	    sys->d2M[index] += inner_prod * inner_prod;
	    //}
	}

    }

    return 0;
}


/*****************************************************************************************
eval_M_nbPair_bonds():
*****************************************************************************************/
int eval_M_nbPair_bonds(tW_system * sys, int n_basis_ij, int ij_index[],
			dvec ij_basis[], int nr_i_bond_coeff, int D_b[],
			dvec calG_b[])
{
    int i, j;
    double inner_prod;
    int N = sys->N_coeff;
    int index;

    for (i = 0; i < n_basis_ij; i++) {
	for (j = 0; j < nr_i_bond_coeff; j++) {
	    inner_prod = dot_prod(ij_basis[i], calG_b[j]);
	    // JFR - added 04.11.12: put the matrix in packed form
	    index = index_Lpacked(ij_index[i], D_b[j], N);

	    sys->M[index] += inner_prod;
	    if (sys->M_cnt != NULL) {
		sys->M_cnt[index] += 1.0;
	    }			//JFR - added 04.13.12: check if in LOWMEM mode

	    /* If necessary, keep track of <M^2>, for variance calculation */
	    //if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 )
	    //{
	    sys->d2M[index] += inner_prod * inner_prod;
	    //}
	}
    }

    return 0;
}


/*****************************************************************************************
get_nb_pair_info():
*****************************************************************************************/
void get_nb_pair_info(double *r_ij, dvec x_i, dvec x_j, dvec x_ij,
		      int b_PBC, dvec box)
{
    /* Determine closest j-site if periodic boundary conditions. */
    if (b_PBC) {
	det_min_image(box, x_i, x_j);
    }
    vect_diff(x_i, x_j, x_ij);
    *r_ij = calc_norm(x_ij);
}


/*****************************************************************************************
get_nb_pair_inter_ptr():
*****************************************************************************************/
tW_type_inter2 *get_nb_pair_inter_ptr(tW_CG_site * j_site, tW_system * sys,
				      int N_i1, tW_word *nList1,
				      int *iList1)
{
    int j_iList1;
    int ij_type;

    /* ID interaction between site i and j. */
    /* ASSUMES that no two particles are involved in more than 1 interaction. */
    j_iList1 = match_word(N_i1, j_site->name, nList1);

    /* Do i and j interact? */
    if (j_iList1 < 0) {
	return NULL;
    }

    /* If so, get the index for sys->Inter2_Type_List. */
    ij_type = iList1[j_iList1];

    return &(sys->Inter2_Type_List[ij_type]);
}


/*****************************************************************************************
check_inter_range():
*****************************************************************************************/
int check_inter_range(double r_ij, double R_0, double R_max, double dr)
{
    if (r_ij < R_0) {
	return -1;
    }
    if (r_ij > R_max + 0.5 * dr) {
	return -1;
    }

    return 0;
}

/* START JFR */
/*****************************************************************************************
check_linear_inter_range():  Changed r_ij > R_max + 0.5*dr to r_ij > R_max 
*****************************************************************************************/
int check_linear_inter_range(double r_ij, double R_0, double R_max,
			     double dr)
{
    if (r_ij < R_0) {
	return -1;
    }
    if (r_ij >= R_max) {
	return -1;
    }

    return 0;
}

/* END JFR */

/*****************************************************************************************
check_Bspline_inter_range():  JFR - 07.22.12: same as linear right now
*****************************************************************************************/
int check_Bspline_inter_range(double r_ij, double R_0, double R_max,
			      double dr)
{
    if (r_ij < R_0) {
	return -1;
    }
    if (r_ij >= R_max) {
	return -1;
    }

    return 0;
}

/*****************************************************************************************
update_nb_pair_grids(): Currently, not using L for this basis.
*****************************************************************************************/
void update_nb_pair_grids(int n_basis_ij, int ij_index[], dvec ij_basis[],
			  bool b_F, double r_ij, tW_system * sys,
			  tW_CG_site site_i, tW_type_inter2 * ij_inter)
{
    int i;
    double exponent;
    double r_exp;
    int *powers = ij_inter->powers;
    /* START JFR */
    double A, B, Ap, Bp;
    double dr = ij_inter->dr;
    double R_0 = ij_inter->R_0;
    int i_0 = ij_inter->i_0;
    int N_pts = ij_inter->N_pts;
    /* END JFR */
    int k = ij_inter->kspline;	/* JFR - 07.22.12 */

    for (i = 0; i < n_basis_ij; i++) {
	if (ij_inter->i_basis == DELTA_BASIS_INDEX) {
	    sys->g[ij_index[i]] += 1.0;
	    sys->g_cnt[ij_index[i]] += 1.0;
	    sys->L[ij_index[i]] += 2.0 / r_ij;
	}
	/* START JFR */
	else if (ij_inter->i_basis == LINEAR_BASIS_INDEX) {
	    if (i == 0) {	/* Update everything at once */

		A = calc_linear_spline_A(i, ij_index, dr, R_0, i_0, N_pts,
					 r_ij, &Ap, &Bp, FALSE);
		B = 1 - A;

		if (LINEAR_STRUCT == 0) {
		    sys->g[ij_index[i]] += A;
		    sys->g[ij_index[i + 1]] += B;
		} else {
		    sys->g[ij_index[i]] += Ap;
		    sys->g[ij_index[i + 1]] += Bp;
		}
		sys->g_cnt[ij_index[i]] += 1.0;
		sys->g_cnt[ij_index[i + 1]] += 1.0;
		sys->L[ij_index[i]] += (A * 2.0) / r_ij;
		sys->L[ij_index[i + 1]] += (B * 2.0) / r_ij;
	    }
	}
	/* END JFR */
	else if (ij_inter->i_basis == POWER_INDEX) {
	    exponent = powers[i] + 2.0;
	    r_exp = pow(r_ij, exponent);
	    sys->g_cnt[ij_index[i]] += 1.0;
	    sys->g[ij_index[i]] +=
		(-1.0) * (powers[i]) * (powers[i] - 1.0) / r_exp;
	}

	/* If forces are present, update b. */
	if (b_F) {
	    sys->b[ij_index[i]] += dot_prod(site_i.f, ij_basis[i]);
	}
	if (sys->flag_ref_potential) {
	    sys->b_ref[ij_index[i]] += dot_prod(site_i.ref_f, ij_basis[i]);
	}
    }
}


/*****************************************************************************************
get_b_harmonic_bond():
*****************************************************************************************/
int get_b_harmonic_bond(tW_Inter_Types * inter, double kT)
{

    inter->ptr_b[0] = -kT * (inter->ptr_g[0] + inter->ptr_L[0]);
    inter->ptr_b[1] = -kT * (inter->ptr_g[1] + inter->ptr_L[1]);

    return 0;
}


/*****************************************************************************************
get_b_rb_dihedral():
*****************************************************************************************/
int get_b_rb_dihedral(tW_Inter_Types * inter, double kT)
{

    inter->ptr_b[0] = -kT * (inter->ptr_g[0]);
    inter->ptr_b[1] = -kT * (inter->ptr_g[1]);
    inter->ptr_b[2] = -kT * (inter->ptr_g[2]);
    inter->ptr_b[3] = -kT * (inter->ptr_g[3]);
    inter->ptr_b[4] = -kT * (inter->ptr_g[4]);

    return 0;
}


/*****************************************************************************************
skip_excl():
*****************************************************************************************/
bool skip_excl(int nr_excl, int *excl_list, int j)
{
    int i;

    for (i = 0; i < nr_excl; i++) {
	if (j == excl_list[i]) {
	    return TRUE;
	}
    }

    return FALSE;
}


/*****************************************************************************************
subtract_ref_forces(): Subtracts the force calculated from a reference potential from 
the force read from the trajectory.
*****************************************************************************************/
void subtract_ref_forces(int N_sites, tW_CG_site * CG_struct)
{
    int i, j;

    for (i = 0; i < N_sites; i++) {
	for (j = 0; j < 3; j++) {
	    CG_struct[i].f[j] -= CG_struct[i].ref_f[j];
	}
    }
}


/*****************************************************************************************
get_b_TOY_dihedral():
*****************************************************************************************/
void get_b_TOY_dihedral(tW_Inter_Types * inter, double kT)
{
    inter->ptr_b[0] = -kT * (inter->ptr_g[0]);
    inter->ptr_b[1] = -kT * (inter->ptr_g[1]);
    inter->ptr_b[2] = -kT * (inter->ptr_g[2]);
}


/*****************************************************************************************
eval_delta_basis_vectors():
*****************************************************************************************/
int eval_delta_basis_vectors(tW_type_inter2 * ij_inter, double n_basis_ij,
			     dvec u_ij, dvec * ij_basis, double r_ij,
			     int *ij_index, tW_system * sys,
			     tW_CG_site site_i, bool b_F)
{
    /* JFR - 07.23.12: update the grids */
    sys->g[ij_index[0]] += 1.0;
    sys->g_cnt[ij_index[0]] += 1.0;
    sys->L[ij_index[0]] += 2.0 / r_ij;

    /* If forces are present, update b. */
    if (b_F) {
	sys->b[ij_index[0]] += dot_prod(site_i.f, ij_basis[0]);
    }
    if (sys->flag_ref_potential) {
	sys->b_ref[ij_index[0]] += dot_prod(site_i.ref_f, ij_basis[0]);
    }

    if (b_F) {
	sys->d2b[ij_index[0]] +=
	    dot_prod(site_i.f, ij_basis[0]) * dot_prod(site_i.f,
						       ij_basis[0]);
    }
    /* end update the grids */

    return 0;

}

/* START JFR */
/*****************************************************************************************
eval_linear_basis_vectors():
*****************************************************************************************/
int eval_linear_basis_vectors(tW_type_inter2 * ij_inter, double n_basis_ij,
			      dvec u_ij, dvec * ij_basis, double r_ij,
			      int *ij_index, tW_system * sys,
			      tW_CG_site site_i, bool b_F, int flag_grids)
{
    int i;
    int i_0, N_pts;
    double dr, A, B, R_0, Ap, Bp;
    int coeff_ctr = 0;

    dr = ij_inter->dr;
    R_0 = ij_inter->R_0;
    i_0 = ij_inter->i_0;
    N_pts = ij_inter->N_pts;

    A = calc_linear_spline_A(coeff_ctr, ij_index, dr, R_0, i_0, N_pts,
			     r_ij, &Ap, &Bp, FALSE);
    B = 1 - A;

    scal_times_vect(A, u_ij, ij_basis[0]);
    scal_times_vect(B, u_ij, ij_basis[1]);

    if (flag_grids == TRUE) {	/* JFR - 07.23.12: update the grids */
	if (LINEAR_STRUCT == 0) {
	    sys->g[ij_index[0]] += A;
	    sys->g[ij_index[1]] += B;
	} else {
	    sys->g[ij_index[0]] += Ap;
	    sys->g[ij_index[1]] += Bp;
	}
	sys->g_cnt[ij_index[0]] += 1.0;
	sys->g_cnt[ij_index[1]] += 1.0;
	sys->L[ij_index[0]] += (A * 2.0) / r_ij;
	sys->L[ij_index[1]] += (B * 2.0) / r_ij;
	/* If forces are present, update b. */
	if (b_F) {
	    sys->b[ij_index[0]] += dot_prod(site_i.f, ij_basis[0]);
	    sys->b[ij_index[1]] += dot_prod(site_i.f, ij_basis[1]);
	}
	if (sys->flag_ref_potential) {
	    sys->b_ref[ij_index[0]] += dot_prod(site_i.ref_f, ij_basis[0]);
	    sys->b_ref[ij_index[1]] += dot_prod(site_i.ref_f, ij_basis[1]);
	}
	//if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) 
	//{ 
	if (b_F) {
	    sys->d2b[ij_index[0]] +=
		dot_prod(site_i.f, ij_basis[0]) * dot_prod(site_i.f,
							   ij_basis[0]);
	    sys->d2b[ij_index[1]] +=
		dot_prod(site_i.f, ij_basis[1]) * dot_prod(site_i.f,
							   ij_basis[1]);
	}
	//}
    }
    /* end update the grids */
    return 0;

}

/* END JFR */

/*****************************************************************************************
eval_Bspline_basis_vectors(): JFR - 07.22.12
*****************************************************************************************/
int eval_Bspline_basis_vectors(tW_type_inter2 * ij_inter, dvec u_ij,
			       dvec * ij_basis, double r_ij, int *ij_index,
			       tW_system * sys, tW_CG_site site_i,
			       bool b_F, int flag_grids)
{
    int i, j, l;
    int i_0, N_pts, k;
    double dr, B, R_0, Bp;

    dr = ij_inter->dr;
    R_0 = ij_inter->R_0;
    i_0 = ij_inter->i_0;
    N_pts = ij_inter->N_pts;
    k = ij_inter->kspline;
    int coeff_ctr = 0;

    double Nik[MAX_BSPLINE_COEFF];
    double Npik[MAX_BSPLINE_COEFF];

    i = get_grid_index_for_linear_basis(r_ij, i_0, dr, R_0) - i_0;

    B = calc_Bspline(dr, R_0, i_0, N_pts, r_ij, i, k, Nik, Npik);
    for (j = 0; j < k; j++) {
	ij_index[coeff_ctr] = i - j + k / 2 + i_0;	/* The k/2 shifts the supports up to the proper place, but this only works for even order Bsplines! */

	if (((i - j + k / 2) < 0) || ((i - j + k / 2) >= N_pts - 1)) {
	    continue;
	}			/* You are beyond the grid, there are no supports here */
	B = Nik[k - 1 - j];
	Bp = calc_Bspline_deriv(dr, R_0, i_0, N_pts, r_ij, i, j, k, Nik,
				Npik);

	scal_times_vect(B, u_ij, ij_basis[coeff_ctr]);

	if (flag_grids == TRUE) {	/* JFR - 07.23.12: update the grids */
	    if (BSPLINE_STRUCT == 0) {
		sys->g[ij_index[coeff_ctr]] += B;
	    } else {
		sys->g[ij_index[coeff_ctr]] += Bp;
	    }
	    sys->g_cnt[ij_index[coeff_ctr]] += 1.0;
	    sys->L[ij_index[coeff_ctr]] += (B * 2.0) / r_ij;
	    /* If forces are present, update b. */
	    if (b_F) {
		sys->b[ij_index[coeff_ctr]] +=
		    dot_prod(site_i.f, ij_basis[coeff_ctr]);
	    }
	    if (sys->flag_ref_potential) {
		sys->b_ref[ij_index[coeff_ctr]] +=
		    dot_prod(site_i.ref_f, ij_basis[coeff_ctr]);
	    }
	    //if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) 
	    //{ 
	    if (b_F) {
		sys->d2b[ij_index[coeff_ctr]] +=
		    dot_prod(site_i.f,
			     ij_basis[coeff_ctr]) * dot_prod(site_i.f,
							     ij_basis
							     [coeff_ctr]);
	    }
	    //}  

	}			/* end update the grids */
	coeff_ctr++;
    }

    return coeff_ctr;
}

/*****************************************************************************************
eval_power_basis_vectors():
*****************************************************************************************/
int eval_power_basis_vectors(tW_type_inter2 * ij_inter, dvec u_ij,
			     dvec * power_ij_basis, double r_ij,
			     int *ij_index, tW_system * sys,
			     tW_CG_site site_i, bool b_F, int flag_grids)
{
    int i;
    int N_coeff = ij_inter->N_coeff;
    double power;
    double exponent;

    for (i = 0; i < N_coeff; i++) {
	exponent = ij_inter->powers[i] + 1.0;
	power = pow(r_ij, exponent);
	scal_times_vect(ij_inter->powers[i] / power, u_ij,
			power_ij_basis[i]);
	ij_index[i] = ij_inter->i_0 + i;

	if (flag_grids == TRUE) {	/* JFR - 07.23.12: update the grids */
	    sys->g_cnt[ij_index[i]] += 1.0;
	    sys->g[ij_index[i]] +=
		(-1.0) * (ij_inter->powers[i]) * (ij_inter->powers[i] -
						  1.0) / power;
	    /* If forces are present, update b. */
	    if (b_F) {
		sys->b[ij_index[i]] +=
		    dot_prod(site_i.f, power_ij_basis[i]);
	    }
	    if (sys->flag_ref_potential) {
		sys->b_ref[ij_index[i]] +=
		    dot_prod(site_i.ref_f, power_ij_basis[i]);
	    }
	    //if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) 
	    //{ 
	    if (b_F) {
		sys->d2b[ij_index[i]] +=
		    dot_prod(site_i.f,
			     power_ij_basis[i]) * dot_prod(site_i.f,
							   power_ij_basis
							   [i]);
	    }
	    //}
	}			/* end update the grids */
    }

    return 0;
}


/*****************************************************************************************
eval_M_nb_pair_inter():
*****************************************************************************************/
void eval_M_nb_pair_inter(int n_basis_ij, int ij_index[], dvec ij_basis[],
			  int n_basis_ik, int ik_index[], dvec ik_basis[],
			  tW_system * sys, int flag_M2)
{

    int i, j;
    double inner_prod;
    int N = sys->N_coeff;
    int index;

    for (i = 0; i < n_basis_ij; i++) {
	for (j = 0; j < n_basis_ik; j++) {
	    inner_prod = dot_prod(ij_basis[i], ik_basis[j]);
	    // JFR - added 04.11.12: put the matrix in packed form
	    if (ik_index[j] <= ij_index[i]) {
		index = index_Lpacked(ij_index[i], ik_index[j], N);
		if (flag_M2 == TRUE) {
		    sys->M2[index] += inner_prod;
		} else {
		    sys->M[index] += inner_prod;
		}
		if (sys->M_cnt != NULL) {
		    sys->M_cnt[index] += 1.0;
		}		//JFR - added 04.13.12: check if in LOWMEM mode
		sys->d2M[index] += inner_prod * inner_prod;
		///*if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 ) {*/ sys->d2M[index] += inner_prod*inner_prod; /*}*/
	    }
	}
    }
}


/*****************************************************************************************
get_b_power():
*****************************************************************************************/
void get_b_power(tW_Inter_Types * inter, double kT)
{
    int i;
    int N_coeff = inter->N_coeff;

    for (i = 0; i < N_coeff; i++) {
	inter->ptr_b[i] = -kT * (inter->ptr_g[i] + inter->ptr_L[i]);
    }

}


/*****************************************************************************************
calc_grids2(): The fast loops.
*****************************************************************************************/
int calc_grids2(FILE * fp, tW_gmx_info info, int N_sites, tW_CG_site CG_struct[], tW_system * sys)
{
    int i, j, k;
    int n_basis_ij, n_basis_ji, n_basis_ik, n_basis_jk;
    int N_i1, N_j1;
    int ij_index[MAX_BSPLINE_COEFF], ji_index[MAX_BSPLINE_COEFF], ik_index[MAX_BSPLINE_COEFF], jk_index[MAX_BSPLINE_COEFF];	/* JFR - 07.22.12: POWER -> BSPLINE for more coeff */
    int nr_i_bond_coeff, nr_j_bond_coeff;
 //   int iList1[Max_Num_Prot_Inter1], jList1[Max_Num_Prot_Inter1];
    int *iList1, *jList1;
    int iSite, jSite;
    double r_ij, r_ik, r_jk;
    bool case1_flag, case2_flag, case3_flag;
    bool b_F = info.b_Forces;
    bool b_PBC = info.b_PBC;
    dvec box;
    //tW_word nList1_i[Max_Num_Prot_Inter1], nList1_j[Max_Num_Prot_Inter1];
    tW_word *nList1_i, *nList1_j;
    int D_b_i[MAX_NUM_BOND_COEFF], D_b_j[MAX_NUM_BOND_COEFF],
	D_b_inter_index[MAX_NUM_BOND_COEFF]
	/* JFR - 02.26.13: for sep M2 */ ;
    dvec calG_b_i[MAX_NUM_BOND_COEFF], calG_b_j[MAX_NUM_BOND_COEFF];
    tW_type_inter2 *ij_inter, *ik_inter, *jk_inter;
    dvec x_i, x_j, x_k;
    dvec x_ij, x_ji, x_ik, x_jk;
    dvec u_ij, u_ji, u_ik, u_jk;
    dvec ij_basis[MAX_BSPLINE_COEFF], ji_basis[MAX_BSPLINE_COEFF], ik_basis[MAX_BSPLINE_COEFF], jk_basis[MAX_BSPLINE_COEFF];	/* JFR - 07.22.12: POWER -> BSPLINE for more coeff */




    /* Make sure rectangular box was used for periodic boundary conditions. */
    if (b_PBC) {
	if (!setup_box_dimensions(box, info.box)) {
	    exit(EXIT_FAILURE);
	}
    }

    /* For each site. */
    for (i = 0; i < N_sites; i++) {
	/* JFR - 06.27.12: Chi2 */
	sys->Chi2 += dot_prod(CG_struct[i].f, CG_struct[i].f);

	/* Get indices and basis vectors for bonded interactions for site i. */
	nr_i_bond_coeff =
	    eval_bond_basis_vectors(fp, CG_struct, sys->Bonded_Inter_Types,
				    i, D_b_i, calG_b_i, TRUE,
				    D_b_inter_index
				    /* JFR - 02.26.13: for sep M2 */ );

	/* Eval b for bond basis vectors and also M for correlations b/w bonds for site i. */
	eval_bonds_M_b(sys, nr_i_bond_coeff, D_b_i, calG_b_i, b_F,
		       CG_struct[i].f, CG_struct[i].ref_f,
		       D_b_inter_index /* JFR - 02.26.13: for sep M2 */ );


	/* Get list of NB pair interactions involving site i, list of site_names involved, and no. interactions. */
//	N_i1 = get_iList(CG_struct[i].name, sys->N_Inter2_Types, sys->Inter2_Type_List, iList1, nList1_i);

	iSite = match_word(sys->N_Site_Types, CG_struct[i].name, sys->Site_Types);
	N_i1 = sys->Inter_Map_Len[iSite];
	iList1 = sys->Inter_iMap[iSite];
	nList1_i = sys->Inter_Map[iSite];

	/* Note: j != i. */
	for (j = 0; j < i; j++) {
	    /* If i and j do not have nb interaction, skip the innter loops. */
	    if (skip_excl(CG_struct[i].nr_excl, CG_struct[i].excl_list, j)) {
		continue;
	    }

	    /* Get list of NB pair interactions involving site j, list of site_names involved, and no. interactions. */
//	    N_j1 = get_iList(CG_struct[j].name, sys->N_Inter2_Types, sys->Inter2_Type_List, jList1, nList1_j);

	    jSite = match_word(sys->N_Site_Types, CG_struct[j].name, sys->Site_Types);
	    N_j1 = sys->Inter_Map_Len[jSite];
	    jList1 = sys->Inter_iMap[jSite];
	    nList1_j = sys->Inter_Map[jSite];

	    /* Get indices and basis vectors for bonded interactions for site i */
	    nr_j_bond_coeff =
		eval_bond_basis_vectors(fp, CG_struct,
					sys->Bonded_Inter_Types, j, D_b_j,
					calG_b_j, FALSE,
					D_b_inter_index
					/* JFR - 02.26.13: for sep M2 */ );

	    /* This is the ij_inter and the ji_inter. */
	    ij_inter =
		get_nb_pair_inter_ptr(&(CG_struct[j]), sys, N_i1, nList1_i,
				      iList1);
	    if (ij_inter == NULL) {
		continue;
	    }

	    /* Get pair info. for the ij and ji pair. */
	    copy_vector(CG_struct[i].r, x_i);
	    copy_vector(CG_struct[j].r, x_j);
	    get_nb_pair_info(&r_ij, x_i, x_j, x_ij, b_PBC, box);
	    get_nb_pair_info(&r_ij, x_j, x_i, x_ji, b_PBC, box);

	    /* Evaluate the basis vectors. */
	    if (ij_inter->i_basis == DELTA_BASIS_INDEX) {
		n_basis_ij = 1;
		n_basis_ji = 1;

		if (check_inter_range
		    (r_ij, ij_inter->R_0, ij_inter->R_max,
		     ij_inter->dr) != 0) {
		    continue;
		}

		ij_inter->N_inter += 2;

		ij_index[0] =
		    get_grid_index_for_delta_basis(r_ij, ij_inter->i_0,
						   ij_inter->dr,
						   ij_inter->R_0);
		ji_index[0] = ij_index[0];

		scal_times_vect(1.0 / r_ij, x_ij, ij_basis[0]);
		scal_times_vect(1.0 / r_ij, x_ji, ji_basis[0]);

		eval_delta_basis_vectors(ij_inter, n_basis_ij, u_ij,
					 ij_basis, r_ij, ij_index, sys,
					 CG_struct[i], b_F);
		eval_delta_basis_vectors(ij_inter, n_basis_ji, u_ji,
					 ji_basis, r_ij, ji_index, sys,
					 CG_struct[j], b_F);
	    }
	    /* START JFR */
	    else if (ij_inter->i_basis == LINEAR_BASIS_INDEX) {
		n_basis_ij = 2;
		n_basis_ji = 2;

		if (check_linear_inter_range
		    (r_ij, ij_inter->R_0, ij_inter->R_max,
		     ij_inter->dr) != 0) {
		    continue;
		}

		ij_inter->N_inter += 2;	/* need to check this -> looks good to me JFR 11/29/10 */

		scal_times_vect(1.0 / r_ij, x_ij, u_ij);
		scal_times_vect(1.0 / r_ij, x_ji, u_ji);

		eval_linear_basis_vectors(ij_inter, n_basis_ij, u_ij,
					  ij_basis, r_ij, ij_index, sys,
					  CG_struct[i], b_F, TRUE);
		eval_linear_basis_vectors(ij_inter, n_basis_ji, u_ji,
					  ji_basis, r_ij, ji_index, sys,
					  CG_struct[j], b_F, TRUE);

	    }
	    /* END JFR */
	    else if (ij_inter->i_basis == BSPLINE_BASIS_INDEX) {	/* JFR - 07.22.12 */

	   // The number of basis vectors affected depends on your distance from the endpoints
	   //        n_basis_ij  = 2*ij_inter->kspline;
	   //        n_basis_ji  = 2*ij_inter->kspline;

		if (check_Bspline_inter_range
		    (r_ij, ij_inter->R_0, ij_inter->R_max,
		     ij_inter->dr) != 0) {
		    continue;
		}

		ij_inter->N_inter += 2;

		scal_times_vect(1.0 / r_ij, x_ij, u_ij);
		scal_times_vect(1.0 / r_ij, x_ji, u_ji);

		n_basis_ij =
		    eval_Bspline_basis_vectors(ij_inter, u_ij, ij_basis,
					       r_ij, ij_index, sys,
					       CG_struct[i], b_F, TRUE);
		n_basis_ji =
		    eval_Bspline_basis_vectors(ij_inter, u_ji, ji_basis,
					       r_ij, ji_index, sys,
					       CG_struct[j], b_F, TRUE);

	    } else if (ij_inter->i_basis == POWER_INDEX) {
		if (r_ij <= ij_inter->R_max) {
		    n_basis_ij = ij_inter->N_powers;
		    n_basis_ji = ij_inter->N_powers;

		    scal_times_vect(1.0 / r_ij, x_ij, u_ij);
		    scal_times_vect(1.0 / r_ij, x_ji, u_ji);

		    eval_power_basis_vectors(ij_inter, u_ij, ij_basis,
					     r_ij, ij_index, sys,
					     CG_struct[i], b_F, TRUE);
		    eval_power_basis_vectors(ij_inter, u_ji, ji_basis,
					     r_ij, ji_index, sys,
					     CG_struct[j], b_F, TRUE);
		} else {
		    continue;
		}
	    } else {
		printf
		    ("ERROR: No other basis functions are implemented.\n");
		exit(EXIT_FAILURE);
	    }
	    /* JFR - 07.23.12: update the grids in eval_----_basis_vectors() now */
	    //update_nb_pair_grids( n_basis_ij, ij_index, ij_basis, b_F, r_ij, sys, CG_struct[i], ij_inter );
	    //update_nb_pair_grids( n_basis_ji, ji_index, ji_basis, b_F, r_ij, sys, CG_struct[j], ij_inter );
	    eval_M_nbPair_bonds(sys, n_basis_ij, ij_index, ij_basis,
				nr_i_bond_coeff, D_b_i, calG_b_i);
	    eval_M_nbPair_bonds(sys, n_basis_ji, ji_index, ji_basis,
				nr_j_bond_coeff, D_b_j, calG_b_j);


	    if (sys->REF_var.flag_calcbref == FALSE) {	/* JFR - 07.16.12: if you just want to calculate bref, skip the inner loops */
		/* k loop 1. */
		for (k = 0; k <= j; k++) {
		    case1_flag = FALSE;
		    case2_flag = FALSE;

		    copy_vector(CG_struct[k].r, x_k);
		    ik_inter =
			get_nb_pair_inter_ptr(&(CG_struct[k]), sys, N_i1,
					      nList1_i, iList1);
		    jk_inter =
			get_nb_pair_inter_ptr(&(CG_struct[k]), sys, N_j1,
					      nList1_j, jList1);

		    /* Case 1: ij:ik, since k<=j<i, k!=i in this loop. */
		    if ((ik_inter != NULL)
			&&
			!(skip_excl
			  (CG_struct[i].nr_excl, CG_struct[i].excl_list,
			   k))) {
			get_nb_pair_info(&r_ik, x_i, x_k, x_ik, b_PBC,
					 box);

			if (ik_inter->i_basis == DELTA_BASIS_INDEX) {
			    n_basis_ik = 1;

			    if (check_inter_range
				(r_ik, ik_inter->R_0, ik_inter->R_max,
				 ik_inter->dr) == 0) {
				ik_index[0] =
				    get_grid_index_for_delta_basis(r_ik,
								   ik_inter->
								   i_0,
								   ik_inter->
								   dr,
								   ik_inter->
								   R_0);
				scal_times_vect(1.0 / r_ik, x_ik,
						ik_basis[0]);
				case1_flag = TRUE;
			    }
			}
			/* START JFR */
			else if (ik_inter->i_basis == LINEAR_BASIS_INDEX) {

			    n_basis_ik = 2;

			    if (check_linear_inter_range
				(r_ik, ik_inter->R_0, ik_inter->R_max,
				 ik_inter->dr) == 0) {
				scal_times_vect(1.0 / r_ik, x_ik, u_ik);
				eval_linear_basis_vectors(ik_inter,
							  n_basis_ik, u_ik,
							  ik_basis, r_ik,
							  ik_index, sys,
							  CG_struct[k],
							  b_F, FALSE);
				case1_flag = TRUE;
			    }
			}
			/* END JFR */
			else if (ik_inter->i_basis == BSPLINE_BASIS_INDEX) {	/* JFR - 07.22.12 */

			// The number of basis vectors affected depends on your distance from the endpoints
			//              n_basis_ik  = 2*ik_inter->kspline;

			    if (check_Bspline_inter_range
				(r_ik, ik_inter->R_0, ik_inter->R_max,
				 ik_inter->dr) == 0) {
				scal_times_vect(1.0 / r_ik, x_ik, u_ik);
				n_basis_ik =
				    eval_Bspline_basis_vectors(ik_inter,
							       u_ik,
							       ik_basis,
							       r_ik,
							       ik_index,
							       sys,
							       CG_struct
							       [k], b_F,
							       FALSE);
				case1_flag = TRUE;
			    }
			} else if (ik_inter->i_basis == POWER_INDEX) {
			    if (r_ik <= ik_inter->R_max) {
				n_basis_ik = ik_inter->N_powers;
				scal_times_vect(1.0 / r_ik, x_ik, u_ik);
				eval_power_basis_vectors(ik_inter, u_ik,
							 ik_basis, r_ik,
							 ik_index, sys,
							 CG_struct[i], b_F,
							 FALSE);
				case1_flag = TRUE;
			    }
			} else {
			    printf
				("ERROR: No other basis functions are implemented.\n");
			    exit(EXIT_FAILURE);
			}
			if (case1_flag == TRUE) {
			    if (j == k) {	/* This is a 2-body contribution */
				eval_M_nb_pair_inter(n_basis_ij, ij_index,
						     ij_basis, n_basis_ik,
						     ik_index, ik_basis,
						     sys, TRUE);
				eval_M_nb_pair_inter(n_basis_ik, ik_index,
						     ik_basis, n_basis_ij,
						     ij_index, ij_basis,
						     sys, TRUE);
			    } else {	/* This is a 3-body contribution */

				eval_M_nb_pair_inter(n_basis_ij, ij_index,
						     ij_basis, n_basis_ik,
						     ik_index, ik_basis,
						     sys, FALSE);
				eval_M_nb_pair_inter(n_basis_ik, ik_index,
						     ik_basis, n_basis_ij,
						     ij_index, ij_basis,
						     sys, FALSE);
			    }

			}
		    }

		    /* End case 1 k excl if. */
		    /* Case 2: ji:jk, at this point i!=j, but j may equal k so I need to take care. */
		    if ((jk_inter != NULL)
			&&
			!(skip_excl
			  (CG_struct[j].nr_excl, CG_struct[j].excl_list,
			   k)) && (j != k)) {
			get_nb_pair_info(&r_jk, x_j, x_k, x_jk, b_PBC,
					 box);

			if (jk_inter->i_basis == DELTA_BASIS_INDEX) {
			    n_basis_jk = 1;

			    if (check_inter_range
				(r_jk, jk_inter->R_0, jk_inter->R_max,
				 jk_inter->dr) == 0) {
				jk_index[0] =
				    get_grid_index_for_delta_basis(r_jk,
								   jk_inter->
								   i_0,
								   jk_inter->
								   dr,
								   jk_inter->
								   R_0);
				scal_times_vect(1.0 / r_jk, x_jk,
						jk_basis[0]);
				case2_flag = TRUE;
			    }
			}
			/* START JFR */
			else if (jk_inter->i_basis == LINEAR_BASIS_INDEX) {
			    n_basis_jk = 2;

			    if (check_linear_inter_range
				(r_jk, jk_inter->R_0, jk_inter->R_max,
				 jk_inter->dr) == 0) {
				scal_times_vect(1.0 / r_jk, x_jk, u_jk);
				eval_linear_basis_vectors(jk_inter,
							  n_basis_jk, u_jk,
							  jk_basis, r_jk,
							  jk_index, sys,
							  CG_struct[i],
							  b_F, FALSE);
				case2_flag = TRUE;
			    }
			}
			/* END JFR */
			else if (jk_inter->i_basis == BSPLINE_BASIS_INDEX) {	/* JFR - 07.22.12 */
			// The number of basis vectors affected depends on your distance from the endpoints
			//              n_basis_jk  = 2*jk_inter->kspline;

			    if (check_Bspline_inter_range
				(r_jk, jk_inter->R_0, jk_inter->R_max,
				 jk_inter->dr) == 0) {
				scal_times_vect(1.0 / r_jk, x_jk, u_jk);
				n_basis_jk =
				    eval_Bspline_basis_vectors(jk_inter,
							       u_jk,
							       jk_basis,
							       r_jk,
							       jk_index,
							       sys,
							       CG_struct
							       [i], b_F,
							       FALSE);
				case2_flag = TRUE;
			    }
			} else if (jk_inter->i_basis == POWER_INDEX) {
			    if (r_jk <= jk_inter->R_max) {
				n_basis_jk = jk_inter->N_powers;
				scal_times_vect(1.0 / r_jk, x_jk, u_jk);
				eval_power_basis_vectors(jk_inter, u_jk,
							 jk_basis, r_jk,
							 jk_index, sys,
							 CG_struct[i], b_F,
							 FALSE);
				case2_flag = TRUE;
			    }
			} else {
			    printf
				("ERROR: No other basis functions are implemented.\n");
			    exit(EXIT_FAILURE);
			}
			if (case2_flag == TRUE) {
			    if (j == k) {	/* This is a 2-body contribution */
				eval_M_nb_pair_inter(n_basis_ji, ji_index,
						     ji_basis, n_basis_jk,
						     jk_index, jk_basis,
						     sys, TRUE);
				eval_M_nb_pair_inter(n_basis_jk, jk_index,
						     jk_basis, n_basis_ji,
						     ji_index, ji_basis,
						     sys, TRUE);
			    } else {	/* This is a 3-body contribution */

				eval_M_nb_pair_inter(n_basis_ji, ji_index,
						     ji_basis, n_basis_jk,
						     jk_index, jk_basis,
						     sys, FALSE);
				eval_M_nb_pair_inter(n_basis_jk, jk_index,
						     jk_basis, n_basis_ji,
						     ji_index, ji_basis,
						     sys, FALSE);
			    }

			}
		    }
		    /* End case 2 k excl if. */
		}		/* end k loop 1. */

		/* k loop 2. */
		for (k = j + 1; k < i; k++) {
		    case3_flag = FALSE;

		    copy_vector(CG_struct[k].r, x_k);
		    jk_inter =
			get_nb_pair_inter_ptr(&(CG_struct[k]), sys, N_j1,
					      nList1_j, jList1);

		    /* Case 3: ji:jk, at this point i!=j and k!=j becase j+1<= k <i. */
		    if ((jk_inter != NULL)
			&&
			!(skip_excl
			  (CG_struct[j].nr_excl, CG_struct[j].excl_list,
			   k)) && (j != k)) {
			get_nb_pair_info(&r_jk, x_j, x_k, x_jk, b_PBC,
					 box);

			if (jk_inter->i_basis == DELTA_BASIS_INDEX) {
			    n_basis_jk = 1;

			    if (check_inter_range
				(r_jk, jk_inter->R_0, jk_inter->R_max,
				 jk_inter->dr) == 0) {
				jk_index[0] =
				    get_grid_index_for_delta_basis(r_jk,
								   jk_inter->
								   i_0,
								   jk_inter->
								   dr,
								   jk_inter->
								   R_0);
				scal_times_vect(1.0 / r_jk, x_jk,
						jk_basis[0]);
				case3_flag = TRUE;
			    }
			}
			/* START JFR */
			else if (jk_inter->i_basis == LINEAR_BASIS_INDEX) {

			    n_basis_jk = 2;

			    if (check_linear_inter_range
				(r_jk, jk_inter->R_0, jk_inter->R_max,
				 jk_inter->dr) == 0) {
				scal_times_vect(1.0 / r_jk, x_jk, u_jk);
				eval_linear_basis_vectors(jk_inter,
							  n_basis_jk, u_jk,
							  jk_basis, r_jk,
							  jk_index, sys,
							  CG_struct[i],
							  b_F, FALSE);
				case3_flag = TRUE;
			    }
			} else if (jk_inter->i_basis == BSPLINE_BASIS_INDEX) {	/* JFR - 07.22.12 */
			// The number of basis vectors affected depends on your distance from the endpoints
			//              n_basis_jk  = 2*jk_inter->kspline;

			    if (check_Bspline_inter_range
				(r_jk, jk_inter->R_0, jk_inter->R_max,
				 jk_inter->dr) == 0) {
				scal_times_vect(1.0 / r_jk, x_jk, u_jk);
				n_basis_jk =
				    eval_Bspline_basis_vectors(jk_inter,
							       u_jk,
							       jk_basis,
							       r_jk,
							       jk_index,
							       sys,
							       CG_struct
							       [i], b_F,
							       FALSE);
				case3_flag = TRUE;
			    }
			}
			/* END JFR */
			else if (jk_inter->i_basis == POWER_INDEX) {
			    if (r_jk <= jk_inter->R_max) {
				n_basis_jk = jk_inter->N_powers;
				scal_times_vect(1.0 / r_jk, x_jk, u_jk);
				eval_power_basis_vectors(jk_inter, u_jk,
							 jk_basis, r_jk,
							 jk_index, sys,
							 CG_struct[i], b_F,
							 FALSE);
				case3_flag = TRUE;
			    }
			} else {
			    printf
				("ERROR: No other basis functions are implemented.\n");
			    exit(EXIT_FAILURE);
			}
			if (case3_flag == TRUE) {
			    if (i == k) {	/* This is a 2-body contribution */
				eval_M_nb_pair_inter(n_basis_ji, ji_index,
						     ji_basis, n_basis_jk,
						     jk_index, jk_basis,
						     sys, TRUE);
				eval_M_nb_pair_inter(n_basis_jk, jk_index,
						     jk_basis, n_basis_ji,
						     ji_index, ji_basis,
						     sys, TRUE);
			    } else {	/* This is a 3-body contribution */

				eval_M_nb_pair_inter(n_basis_ji, ji_index,
						     ji_basis, n_basis_jk,
						     jk_index, jk_basis,
						     sys, FALSE);
				eval_M_nb_pair_inter(n_basis_jk, jk_index,
						     jk_basis, n_basis_ji,
						     ji_index, ji_basis,
						     sys, FALSE);
			    }
			}
		    } /* end if k excl. */

		}		/* end k loop 2. */

	    } /* end if sys->REF_var.flag_calcbref == FALSE */

	}			/* end j loop. */

    }				/* end i loop. */


    /* NJD 
     * To reweight in a framewise manner, we should do the following here: 
     * 1) determine the weight of the frame by whatever mechanism
     * 2) to the arrays G_wt, b_wt, etc., we add G*scale, b*scale,
     * 3) clear the current G, b, etc. so they can be freshly calculated and weighted next frame
     * Note that this only applies to the force calculation in its simplest form - no reference
     * forces or structure calculations have yet been implemented in this manner.
     * */
    if (sys->FRAMEWEIGHT_var.flag_FRAMEWEIGHT == TRUE ) {
	//printf("Weighting frame by weight factor \n");
	double wt_factor = 1.0;

	if ( strcmp(sys->FRAMEWEIGHT_var.FRAMEWEIGHT, "NPT") == 0 )
	{
		wt_factor = box[0]*box[0]; 
	}

	sys->wt_norm += wt_factor;

	for (i=0; i<sys->N_coeff; i++)
	{
		sys->b_wt[i] += sys->b[i] * wt_factor;
		sys->b[i] = 0;
	}

	for (i=0; i<sys->N_pack; i++)
	{
		sys->M_wt[i] += sys->M[i] * wt_factor;
		sys->M2_wt[i] += sys->M2[i] * wt_factor;

		sys->M[i] = 0;
		sys->M2[i] = 0;
	}

    }

    return 0;
}


/*****************************************************************************************
weight_local_top(): Weight the sums from each local proc. 
*****************************************************************************************/
int weight_local_top(tW_system * sys_top, double w_local, int N_coeff)
{
    int i;
    int N_pack = (N_coeff * N_coeff + N_coeff) / 2;

    for (i = 0; i < N_coeff; i++) {
	sys_top->b[i] *= w_local;
	sys_top->b_ref[i] *= w_local;
	sys_top->g[i] *= w_local;
	sys_top->L[i] *= w_local;
	sys_top->g_cnt[i] *= w_local;

        sys_top->d2b[i] *= w_local;

        //*if ( strcmp( sys_top->PC_var.LPC, "bvar" ) == 0 ) { */ sys_top->d2b[i] *= w_local; /*} *//* JFR - 01.31.13 */
    }

    // JFR - added 04.11.12: put the matrix in packed form
    for (i = 0; i < N_pack; i++) {
	sys_top->M[i] *= w_local;
	sys_top->M2[i] *= w_local;
	if (sys_top->M_cnt != NULL) {
	    sys_top->M_cnt[i] *= w_local;
	}			//JFR - added 04.13.12: check if in LOWMEM mode
	sys_top->d2M[i] *= w_local;

	//*if ( strcmp( sys_top->PC_var.RPC, "MTvar" ) == 0 ) {*/ sys_top->d2M[i] *= w_local; /*}*/ /* JFR - 01.31.13 */
    }

    /* JFR - 06.27.12: Chi2 */
    sys_top->Chi2 *= w_local;

    return 0;
}

/*****************************************************************************************
get_b_soln(): multiply the obtained forces with the original matrix to obtain the solns, b
*****************************************************************************************/
int get_b_soln(FILE * fp_log, tW_system * sys)
{
    int i, j, k, l;
    double b_soln_struct = 0.00;
    double b_soln_forces = 0.00;
    double *phi_struct = sys->phi_struct;
    double *phi_forces = sys->phi_forces;
    double *M = sys->M;
    double *g = sys->g;
    FILE *fp_soln_struct, *fp_soln_forces, *fp_soln_struct_diag,
	*fp_soln_forces_diag;
    tW_word fname;

    int index;
    int N = sys->N_coeff;

    fprintf(fp_log, "In get_b_soln.\n");

    //JFR - Separate the contributions to b
    for (i = 0; i < sys->N_Inter_Types; i++) {

	sprintf(fname, "%s.%s.dat", "b_soln_struct_diag",
		sys->Inter_Types[i].inter_name);
	fp_soln_struct_diag = fopen(fname, "w");

	sprintf(fname, "%s.%s.dat", "b_soln_forces_diag",
		sys->Inter_Types[i].inter_name);
	fp_soln_forces_diag = fopen(fname, "w");

	for (j = 0; j < sys->N_Inter_Types; j++) {

	    sprintf(fname, "%s.%s.%s.dat", "b_soln_struct",
		    sys->Inter_Types[i].inter_name,
		    sys->Inter_Types[j].inter_name);
	    fp_soln_struct = fopen(fname, "w");

	    sprintf(fname, "%s.%s.%s.dat", "b_soln_forces",
		    sys->Inter_Types[i].inter_name,
		    sys->Inter_Types[j].inter_name);
	    fp_soln_forces = fopen(fname, "w");

	    for (k = 0; k < sys->Inter_Types[i].N_pts; k++) {

		if ((i == j)) {

		    fprintf(fp_soln_struct_diag, "%10.5lf    %15.5lf \n",
			    sys->Inter_Types[i].ptr_x[k],
			    (g[sys->Inter_Types[i].i_0 + k] *
			     phi_struct[sys->Inter_Types[i].i_0 +
					k]) / sys->Inter_Types[i].dr);
		    fprintf(fp_soln_forces_diag, "%10.5lf    %15.5lf \n",
			    sys->Inter_Types[i].ptr_x[k],
			    (g[sys->Inter_Types[i].i_0 + k] *
			     phi_forces[sys->Inter_Types[i].i_0 +
					k]) / sys->Inter_Types[i].dr);

		}

		for (l = 0; l < sys->Inter_Types[j].N_pts; l++) {

		    // JFR - added 04.11.12: put the matrix in packed form
		    index =
			index_Lpacked((sys->Inter_Types[i].i_0 + k),
				      (sys->Inter_Types[j].i_0 + l), N);

		    if ((i == j) && (k == l)) {
			b_soln_struct +=
			    (M[index] -
			     g[sys->Inter_Types[i].i_0 +
			       k]) * phi_struct[sys->Inter_Types[j].i_0 +
						l];
			b_soln_forces +=
			    (M[index] -
			     g[sys->Inter_Types[i].i_0 +
			       k]) * phi_forces[sys->Inter_Types[j].i_0 +
						l];
		    } else {
			b_soln_struct +=
			    M[index] * phi_struct[sys->Inter_Types[j].i_0 +
						  l];
			b_soln_forces +=
			    M[index] * phi_forces[sys->Inter_Types[j].i_0 +
						  l];

		    }

		}

		fprintf(fp_soln_struct, "%10.5lf    %15.5lf \n",
			sys->Inter_Types[i].ptr_x[k],
			b_soln_struct / sys->Inter_Types[i].dr);

		fprintf(fp_soln_forces, "%10.5lf    %15.5lf \n",
			sys->Inter_Types[i].ptr_x[k],
			b_soln_forces / sys->Inter_Types[i].dr);

		b_soln_struct = 0.00;
		b_soln_forces = 0.00;
	    }
	    fclose(fp_soln_struct);
	    fclose(fp_soln_forces);

	}
	fclose(fp_soln_struct_diag);
	fclose(fp_soln_forces_diag);

    }

    return 0;
}

/*****************************************************************************************
get_b_soln_err(): multiply the obtained forces with the original matrix to obtain the solns, b
*****************************************************************************************/
int get_b_soln_err(FILE * fp_log, tW_system * sys)
{
    int i, j, k, l, dummy, test_sscanf;
    double b_soln_i = 0.00;
    double b_soln_imin1 = 0.00;
    double b_soln_err;
    double *phi_struct = sys->phi_struct;
    double *phi_forces = sys->phi_forces;
    double *M = sys->M;
    double *g = sys->g;
    double *b = sys->b_forces;
    FILE *fp_bsoln_err, *fp_bAA, *fp_forces_i, *fp_forces_imin1;
    tW_word fname;
    tW_line inp_line;
    int index;
    int N = sys->N_coeff;

    int f_imin1_pts;
    double f_imin1_Rmin;
    double phi_imin1[sys->N_coeff];
    int N_front, N_mid, N_back, sign_front, sign_back;

    fprintf(fp_log, "In get_b_soln.\n");

    //JFR - Separate the contributions to b
    for (i = 0; i < sys->N_Inter_Types; i++) {

	// JFR - write the current force vector to a file with the # of pts and min gridpt in the header
	sprintf(fname, "%s.%s.dat", "force_i",
		sys->Inter_Types[i].inter_name);
	fp_forces_i = fopen(fname, "w");

	fprintf(fp_forces_i, "%d \n", sys->Inter_Types[i].N_coeff);
	fprintf(fp_forces_i, "%15.5lf \n",
		sys->x[sys->Inter_Types[i].i_0]);
	for (k = 0; k < sys->Inter_Types[i].N_coeff; k++) {
	    fprintf(fp_forces_i, "%15.5lf \n",
		    phi_forces[sys->Inter_Types[i].i_0 + k]);
	}
	fclose(fp_forces_i);

	// JFR - read the force vector from the previous iteration
	sprintf(fname, "%s.%s.dat", "force_imin1",
		sys->Inter_Types[i].inter_name);
	fp_forces_imin1 = fopen(fname, "r");
	if (fp_forces_imin1 == NULL) {
	    printf("ERROR: fp_force_imin1 points to null");
	    exit(0);
	}
	// get Npts
	dummy = get_next_line(fp_forces_imin1, inp_line);
	test_sscanf = sscanf(inp_line, "%d", &f_imin1_pts);
	printf("f_imin1_pts = %d\n", f_imin1_pts);
	// get R_min
	dummy = get_next_line(fp_forces_imin1, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &f_imin1_Rmin);
	printf("f_imin1_Rmin = %lf\n", f_imin1_Rmin);
	if (sys->Inter_Types[i].dr > FLOAT_EPS) {
	    N_front =
		(int) ((fabs
			((sys->x[sys->Inter_Types[i].i_0] -
			  f_imin1_Rmin) / sys->Inter_Types[i].dr)) +
		       FLOAT_EPS);
	    N_back =
		abs(abs(sys->Inter_Types[i].N_coeff - f_imin1_pts) -
		    N_front);
	} else {
	    N_front = 0;
	    N_back = 0;
	}

	if (sys->x[sys->Inter_Types[i].i_0] <= f_imin1_Rmin) {
	    sign_front = 1;
	} else {
	    sign_front = -1;
	}

	if ( abs(sys->Inter_Types[i].N_coeff - f_imin1_pts) >= N_front ) {
	    sign_back = 1;
	} else {
	    sign_back = -1;
	}

	printf("sys->Inter_Types[i].dr = %lf\n", sys->Inter_Types[i].dr);
	printf("N_front = %d\n", N_front);
	printf("N_back = %d\n", N_back);
	printf("sys->N_coeff = %d\n", sys->N_coeff);
	printf("sys->Inter_Types[i].N_pts = %d\n",
	       sys->Inter_Types[i].N_pts);
	printf("sys->Inter_Types[i].N_coeff = %d\n",
	       sys->Inter_Types[i].N_coeff);

	switch (sign_front) {
	case 1:
	    for (k = 0; k < N_front; k++) {
		phi_imin1[sys->Inter_Types[i].i_0 + k] =
		    phi_forces[sys->Inter_Types[i].i_0 + k];
	    }
	    switch (sign_back) {
	    case 1:
		N_mid = sys->Inter_Types[i].N_coeff - N_front - N_back;
		printf("N_mid = %d\n", N_mid);
		for (k = N_front; k < (N_mid + N_front); k++) {
		    dummy = get_next_line(fp_forces_imin1, inp_line);
		    test_sscanf =
			sscanf(inp_line, "%lf",
			       &phi_imin1[sys->Inter_Types[i].i_0 + k]);
		}
		for (k = (N_mid + N_front);
		     k < sys->Inter_Types[i].N_coeff; k++) {
		    phi_imin1[sys->Inter_Types[i].i_0 + k] =
			phi_forces[sys->Inter_Types[i].i_0 + k];
		}
		break;
	    case -1:
		N_mid =
		    (sys->Inter_Types[i].N_coeff + f_imin1_pts - N_front -
		     N_back) / 2;
		printf("N_mid = %d\n", N_mid);
		for (k = N_front; k < (N_mid + N_front); k++) {
		    dummy = get_next_line(fp_forces_imin1, inp_line);
		    test_sscanf =
			sscanf(inp_line, "%lf",
			       &phi_imin1[sys->Inter_Types[i].i_0 + k]);
		}
		break;
	    }
	    break;
	case -1:
	    for (k = 0; k < N_front; k++) {
		dummy = get_next_line(fp_forces_imin1, inp_line);
	    }
	    switch (sign_back) {
	    case -1:
		N_mid = f_imin1_pts - N_front - N_back;
		printf("N_mid = %d\n", N_mid);
		for (k = 0; k < N_mid; k++) {
		    dummy = get_next_line(fp_forces_imin1, inp_line);
		    test_sscanf =
			sscanf(inp_line, "%lf",
			       &phi_imin1[sys->Inter_Types[i].i_0 + k]);
		}
		break;
	    case 1:
		N_mid =
		    (sys->Inter_Types[i].N_coeff + f_imin1_pts - N_front -
		     N_back) / 2;
		printf("N_mid = %d\n", N_mid);
		for (k = 0; k < N_mid; k++) {
		    dummy = get_next_line(fp_forces_imin1, inp_line);
		    test_sscanf =
			sscanf(inp_line, "%lf",
			       &phi_imin1[sys->Inter_Types[i].i_0 + k]);
		}
		for (k = (N_mid + N_front);
		     k < sys->Inter_Types[i].N_coeff; k++) {
		    phi_imin1[sys->Inter_Types[i].i_0 + k] =
			phi_forces[sys->Inter_Types[i].i_0 + k];
		}
		break;
	    }
	    break;
	}
	fclose(fp_forces_imin1);

	//JFR - check the two forces
	sprintf(fname, "%s.%s.dat", "force_imin1_check",
		sys->Inter_Types[i].inter_name);
	fp_forces_i = fopen(fname, "w");

	fprintf(fp_forces_i, "%d \n", sys->Inter_Types[i].N_coeff);
	fprintf(fp_forces_i, "%15.5lf \n",
		sys->x[sys->Inter_Types[i].i_0]);
	for (k = 0; k < sys->Inter_Types[i].N_coeff; k++) {
	    fprintf(fp_forces_i, "%15.5lf \n",
		    phi_imin1[sys->Inter_Types[i].i_0 + k]);
	}
	fclose(fp_forces_i);


    }

    //JFR - /* now calc bsoln */
    for (i = 0; i < sys->N_Inter_Types; i++) {

	sprintf(fname, "%s.%s.dat", "b_soln_err",
		sys->Inter_Types[i].inter_name);
	fp_bsoln_err = fopen(fname, "w");

	sprintf(fname, "%s.%s.dat", "bAA", sys->Inter_Types[i].inter_name);
	fp_bAA = fopen(fname, "w");

	for (k = 0; k < sys->Inter_Types[i].N_coeff; k++) {

	    for (j = 0; j < sys->N_coeff; j++) {
		// JFR - added 04.11.12: put the matrix in packed form
		index = index_Lpacked((sys->Inter_Types[i].i_0 + k), j, N);

		b_soln_i += M[index] * phi_forces[j];
		b_soln_imin1 += M[index] * phi_imin1[j];

	    }

	    b_soln_err = b_soln_i - b_soln_imin1;

	    fprintf(fp_bsoln_err, "%15.5lf \n", b_soln_err);

	    fprintf(fp_bAA, "%15.5lf \n", b[sys->Inter_Types[i].i_0 + k]);

	    b_soln_i = 0.00;
	    b_soln_imin1 = 0.00;
	}
	fclose(fp_bsoln_err);
	fclose(fp_bAA);

    }

    return 0;
}

/*****************************************************************************************
get_b_soln_errAA(): multiply the obtained forces with the original matrix to obtain the solns, b
*****************************************************************************************/
int get_b_soln_errAA(FILE * fp_log, tW_system * sys)
{
    int i, j, k, l, dummy, test_sscanf;
    double b_soln_i = 0.00;
    double b_soln_imin1 = 0.00;
    double b_soln_err;
    double b_soln_errAA;
    double *phi_struct = sys->phi_struct;
    double *phi_forces = sys->phi_forces;
    double *M = sys->M;
    double *g = sys->g;
    double *b = sys->b_forces;
    FILE *fp_bsoln_err, *fp_bsoln_errAA, *fp_bAA, *fp_forces_i,
	*fp_forces_imin1;
    tW_word fname;
    tW_line inp_line;
    int index;
    int N = sys->N_coeff;

    int f_imin1_pts;
    double f_imin1_Rmin;
    double phi_imin1[sys->N_coeff];
    int N_front, N_mid, N_back, sign_front, sign_back;

    fprintf(fp_log, "In get_b_soln.\n");

    //JFR - Separate the contributions to b
    for (i = 0; i < sys->N_Inter_Types; i++) {

	// JFR - write the current force vector to a file with the # of pts and min gridpt in the header
	sprintf(fname, "%s.%s.dat", "force_i",
		sys->Inter_Types[i].inter_name);
	fp_forces_i = fopen(fname, "w");

	fprintf(fp_forces_i, "%d \n", sys->Inter_Types[i].N_coeff);
	fprintf(fp_forces_i, "%15.5lf \n",
		sys->x[sys->Inter_Types[i].i_0]);
	for (k = 0; k < sys->Inter_Types[i].N_coeff; k++) {
	    fprintf(fp_forces_i, "%15.5lf \n",
		    phi_forces[sys->Inter_Types[i].i_0 + k]);
	}
	fclose(fp_forces_i);

	// JFR - read the force vector from the previous iteration
	sprintf(fname, "%s.%s.dat", "force_imin1",
		sys->Inter_Types[i].inter_name);
	fp_forces_imin1 = fopen(fname, "r");
	if (fp_forces_imin1 == NULL) {
	    printf("ERROR: fp_force_imin1 points to null");
	    exit(0);
	}
	// get Npts
	dummy = get_next_line(fp_forces_imin1, inp_line);
	test_sscanf = sscanf(inp_line, "%d", &f_imin1_pts);
	printf("f_imin1_pts = %d\n", f_imin1_pts);
	// get R_min
	dummy = get_next_line(fp_forces_imin1, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &f_imin1_Rmin);
	printf("f_imin1_Rmin = %lf\n", f_imin1_Rmin);
	if (sys->Inter_Types[i].dr > FLOAT_EPS) {
	    N_front =
		(int) ((fabs
			((sys->x[sys->Inter_Types[i].i_0] -
			  f_imin1_Rmin) / sys->Inter_Types[i].dr)) +
		       FLOAT_EPS);
	    N_back =
		abs(abs(sys->Inter_Types[i].N_coeff - f_imin1_pts) -
		    N_front);
	} else {
	    N_front = 0;
	    N_back = 0;
	}

	if (sys->x[sys->Inter_Types[i].i_0] <= f_imin1_Rmin) {
	    sign_front = 1;
	} else {
	    sign_front = -1;
	}

	if ( abs(sys->Inter_Types[i].N_coeff - f_imin1_pts) >= N_front ) {
	    sign_back = 1;
	} else {
	    sign_back = -1;
	}

//printf("sys->Inter_Types[i].dr = %lf\n",sys->Inter_Types[i].dr);
//printf("N_front = %d\n",N_front);
//printf("N_back = %d\n",N_back);
//printf("sys->N_coeff = %d\n",sys->N_coeff);
//printf("sys->Inter_Types[i].N_pts = %d\n",sys->Inter_Types[i].N_pts);
//printf("sys->Inter_Types[i].N_coeff = %d\n",sys->Inter_Types[i].N_coeff);

	switch (sign_front) {
	case 1:
	    for (k = 0; k < N_front; k++) {
		phi_imin1[sys->Inter_Types[i].i_0 + k] =
		    phi_forces[sys->Inter_Types[i].i_0 + k];
	    }
	    switch (sign_back) {
	    case 1:
		N_mid = sys->Inter_Types[i].N_coeff - N_front - N_back;
//printf("N_mid = %d\n",N_mid);
		for (k = N_front; k < (N_mid + N_front); k++) {
		    dummy = get_next_line(fp_forces_imin1, inp_line);
		    test_sscanf =
			sscanf(inp_line, "%lf",
			       &phi_imin1[sys->Inter_Types[i].i_0 + k]);
		}
		for (k = (N_mid + N_front);
		     k < sys->Inter_Types[i].N_coeff; k++) {
		    phi_imin1[sys->Inter_Types[i].i_0 + k] =
			phi_forces[sys->Inter_Types[i].i_0 + k];
		}
		break;
	    case -1:
		N_mid =
		    (sys->Inter_Types[i].N_coeff + f_imin1_pts - N_front -
		     N_back) / 2;
//printf("N_mid = %d\n",N_mid);
		for (k = N_front; k < (N_mid + N_front); k++) {
		    dummy = get_next_line(fp_forces_imin1, inp_line);
		    test_sscanf =
			sscanf(inp_line, "%lf",
			       &phi_imin1[sys->Inter_Types[i].i_0 + k]);
		}
		break;
	    }
	    break;
	case -1:
	    for (k = 0; k < N_front; k++) {
		dummy = get_next_line(fp_forces_imin1, inp_line);
	    }
	    switch (sign_back) {
	    case -1:
		N_mid = f_imin1_pts - N_front - N_back;
//printf("N_mid = %d\n",N_mid);
		for (k = 0; k < N_mid; k++) {
		    dummy = get_next_line(fp_forces_imin1, inp_line);
		    test_sscanf =
			sscanf(inp_line, "%lf",
			       &phi_imin1[sys->Inter_Types[i].i_0 + k]);
		}
		break;
	    case 1:
		N_mid =
		    (sys->Inter_Types[i].N_coeff + f_imin1_pts - N_front -
		     N_back) / 2;
//printf("N_mid = %d\n",N_mid);
		for (k = 0; k < N_mid; k++) {
		    dummy = get_next_line(fp_forces_imin1, inp_line);
		    test_sscanf =
			sscanf(inp_line, "%lf",
			       &phi_imin1[sys->Inter_Types[i].i_0 + k]);
		}
		for (k = (N_mid + N_front);
		     k < sys->Inter_Types[i].N_coeff; k++) {
		    phi_imin1[sys->Inter_Types[i].i_0 + k] =
			phi_forces[sys->Inter_Types[i].i_0 + k];
		}
		break;
	    }
	    break;
	}
	fclose(fp_forces_imin1);

	//JFR - check the two forces
	sprintf(fname, "%s.%s.dat", "force_imin1_check",
		sys->Inter_Types[i].inter_name);
	fp_forces_i = fopen(fname, "w");

	fprintf(fp_forces_i, "%d \n", sys->Inter_Types[i].N_coeff);
	fprintf(fp_forces_i, "%15.5lf \n",
		sys->x[sys->Inter_Types[i].i_0]);
	for (k = 0; k < sys->Inter_Types[i].N_coeff; k++) {
	    fprintf(fp_forces_i, "%15.5lf \n",
		    phi_imin1[sys->Inter_Types[i].i_0 + k]);
	}
	fclose(fp_forces_i);


    }

    //JFR - /* now calc bsoln */
    for (i = 0; i < sys->N_Inter_Types; i++) {

	sprintf(fname, "%s.%s.dat", "b_soln_err",
		sys->Inter_Types[i].inter_name);
	fp_bsoln_err = fopen(fname, "w");

	sprintf(fname, "%s.%s.dat", "b_soln_errAA",
		sys->Inter_Types[i].inter_name);
	fp_bsoln_errAA = fopen(fname, "w");

	sprintf(fname, "%s.%s.dat", "bAA", sys->Inter_Types[i].inter_name);
	fp_bAA = fopen(fname, "w");

	for (k = 0; k < sys->Inter_Types[i].N_coeff; k++) {

	    for (j = 0; j < sys->N_coeff; j++) {

		// JFR - added 04.11.12: put the matrix in packed form
		index = index_Lpacked((sys->Inter_Types[i].i_0 + k), j, N);

		b_soln_i += M[index] * phi_forces[j];
		b_soln_imin1 += M[index] * phi_imin1[j];

	    }

	    b_soln_err = b_soln_i - b_soln_imin1;

	    b_soln_errAA = b_soln_i - b[sys->Inter_Types[i].i_0 + k];

	    fprintf(fp_bsoln_err, "%15.5lf \n", b_soln_err);

	    fprintf(fp_bsoln_errAA, "%15.5lf \n", b_soln_errAA);

	    fprintf(fp_bAA, "%15.5lf \n", b[sys->Inter_Types[i].i_0 + k]);

	    b_soln_i = 0.00;
	    b_soln_imin1 = 0.00;
	}
	fclose(fp_bsoln_err);
	fclose(fp_bsoln_errAA);
	fclose(fp_bAA);

    }

    return 0;
}

/*****************************************************************************************
calc_d2b(): calculate the variance of b
*****************************************************************************************/
int calc_d2b(tW_system * sys)
{
    int i;

    for (i = 0; i < sys->N_coeff; i++) {
	sys->d2b[i] -= sys->b_forces[i] * sys->b_forces[i];
    }

    return 0;
}

/*****************************************************************************************
calc_d2M(): calculate the variance of M
*****************************************************************************************/
int calc_d2M(tW_system * sys)
{
    int i;

    for (i = 0; i < sys->N_pack; i++) {
	sys->d2M[i] -= sys->M[i] * sys->M[i];
    }

    return 0;
}

/*****************************************************************************************
get_rescale(): copy g_cnt to rescale
*****************************************************************************************/
int get_rescale(tW_system * sys)
{
    int i;

    for (i = 0; i < sys->N_coeff; i++) {
	sys->rescale[i] = sys->g_cnt[i];
    }

    return 0;
}
