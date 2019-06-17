/**
@file read_parameters.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Functions related to reading in settings for the cgff calculation
*/

//c library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//local includes
#include "cgff_types.h"
#include "read_parameters.h"
#include "io_read.h"
#include "safe_mem.h"

/*****************************************************************************************
read_par(): Reads the information from the parameter file and files listed in the param-
eter file. Stores this information in the data structures files, sys, and ref_potential.
A non-blank line is read from the parameter file and is passed to get_parameters() to 
read and store information on associated line for each directive ([...]). After all lines 
have been read, check_input() is called to make sure the input is consistent.
*****************************************************************************************/
int read_par(tW_files * files, tW_system * sys,
	     tW_ref_potential * ref_potential)
{
    tW_line inp_line;
    Par_Flags flags;

    if (files->fp_par == NULL) {
	printf("ERROR: 'par.txt' not found.\n");
	exit(EXIT_FAILURE);
    }

    initialize_sys(sys);
    initialize_files(files);
    initialize_par_flags(&flags);

    do {
	if (get_next_line(files->fp_par, inp_line) == -1) {
	    break;
	}

	get_parameters(files->fp_par, inp_line, sys, files, ref_potential,
		       &flags);

    } while (1);

    if (sys->MEM_var.flag_LOWMEM == TRUE) {	/* JFR- added 04.19.12: If in LOWMEM mode, override some other options */
	sys->MT_var.flag_Mcnt = FALSE;
    }

    check_input(*files, *sys, flags, *ref_potential);

    return 0;
}


/*****************************************************************************************
get_Bond_Inter_Types(): Copies the parameters for a intramolecular interaction in the 
sys data structure for a interaction type listed in the parameter file for a specific
group of site types.
*****************************************************************************************/
tW_Bonded_Inter *get_Bond_Inter_Types(FILE * fp_par, tW_line inp_line,
				      tW_system sys, int *N_Bond_Int_Types,
				      tW_Bonded_Inter * Bonded_Inter_Types,
				      tW_word keyword, int N_Int_sites,
				      int N_Inter_Types)
{
    int i;
    int i_inter = sys.N_Bond_Int_Types;
    tW_word test_keyword;
    tW_Bonded_Inter *Bond_ptr;

    strcpy(test_keyword, "\0");

    *N_Bond_Int_Types = sys.N_Bond_Int_Types + N_Inter_Types;
    Bonded_Inter_Types = (tW_Bonded_Inter *) erealloc(Bonded_Inter_Types,
						     (*N_Bond_Int_Types) *
						     sizeof
						     (tW_Bonded_Inter));

    for (i = 0; i < *N_Bond_Int_Types; i++) {
	Bonded_Inter_Types[i].Inter_List = NULL;
    }

    for (; i_inter < *N_Bond_Int_Types; i_inter++) {
	get_next_line(fp_par, inp_line);

	check_inp_line("par.txt", keyword, inp_line);

	Bond_ptr = &(Bonded_Inter_Types[i_inter]);

	strcpy(Bond_ptr->name, keyword);

	Bond_ptr->N_Int_Sites = N_Int_sites;

	Bond_ptr->N_instances = 0;

	Bond_ptr->Site_Types =
	    (tW_word *) emalloc(N_Int_sites * sizeof(tW_word));

	get_bond_inter_site_types(inp_line, Bond_ptr, N_Int_sites,
				  Bond_ptr->inter_name);

	get_bond_basis_parameters(Bond_ptr->inter_name, Bond_ptr,
				  N_Int_sites, inp_line, sys);
    }

    get_next_line(fp_par, inp_line);

    test_end_directive(inp_line, keyword, "par.txt");

    return Bonded_Inter_Types;
}


/*****************************************************************************************
get_structures(): Saves all the file names that contain structures from the list of 
files names read from the parameter file. Returns the combined list of file names that
contain structures. Also, makes sure that the topology distribution is normalized.
*****************************************************************************************/
tW_word *get_structures(tW_files * files, int *N_struct)
{
    int i_file;
    int N_files = 0;
    int N_files_1;
    int i_file_list = 0;
    int i_file_list_1;
    double pr_tot;
    int test_sscanf;
    FILE *fp_in;
    tW_line inp_line;
    tW_word *file_list = NULL;


    /* Loop over files that contain structures (e.g. inp.txt). */
    for (i_file = 0; i_file < files->N_struct_file; i_file++) {
	/* Open file i_file with list of structures. */
	if (file_exists(files->struct_file[i_file])) {
	    fp_in = (FILE *) fopen(files->struct_file[i_file], "r");
	} else {
	    printf("\nERROR: Could not open structure file: %s.\n",
		   files->struct_file[i_file]);
	    exit(EXIT_FAILURE);
	}

	/* Get the number of structure files in struct_file[i_file]. */
	get_next_line(fp_in, inp_line);
	delete_directive_from_line(inp_line);
	test_sscanf = sscanf(inp_line, " %d ", &N_files_1);
	if (test_sscanf != 1) {
	    printf("\nERROR: Check the number of files in %s.\n",
		   files->struct_file[i_file]);
	    printf("%s\n", inp_line);
	    exit(EXIT_FAILURE);
	}

	/* Update N_files: total number of structure files. */
	N_files = N_files + N_files_1;

	/* Resize array. */
	file_list =
	    (tW_word *) erealloc(file_list, N_files * sizeof(tW_word));
	files->p_struct =
	    (double *) erealloc(files->p_struct, N_files * sizeof(double));

	i_file_list_1 = 0;

	/* Read in the names of the structure files in current list. */
	do {
	    get_next_line(fp_in, inp_line);

	    test_sscanf =
		sscanf(inp_line, "%s %lf", file_list[i_file_list],
		       &(files->p_struct[i_file_list]));
	    if (!file_exists(file_list[i_file_list])) {
		printf
		    ("\nERROR: Structure file '%s' could not be opened.\n",
		     file_list[i_file_list]);
		exit(EXIT_FAILURE);
	    }
	    if (test_sscanf != 2) {
		printf("\nERROR: Problem reading file name in %s\n",
		       files->struct_file[i_file]);
		printf("%s\n", inp_line);
		exit(EXIT_FAILURE);
	    }

	    check_inp_line(files->struct_file[i_file], KEY_STRUCT_FILES,
			   inp_line);

	    i_file_list++;	/* Global counter.           */
	    i_file_list_1++;	/* Counter for current file. */

	} while (i_file_list_1 < N_files_1);

	get_next_line(fp_in, inp_line);
	test_end_directive(inp_line, KEY_STRUCT_FILES, "par.txt");

	fclose(fp_in);
    }

    /* Sum topology probabilities. */
    pr_tot = 0.;
    for (i_file = 0; i_file < N_files; i_file++) {
	pr_tot += files->p_struct[i_file];
    }

    /* Make sure probability is normalized. */
    if (fabs(pr_tot - 1.0) > FLOAT_EPS) {
	printf
	    ("WARNING!! input topology probability distributions unnormalized.  Renormalizing now.\n");
	for (i_file = 0; i_file < N_files; i_file++) {
	    files->p_struct[i_file] /= pr_tot;
	}
    }

    *N_struct = N_files;

    return file_list;
}


/*****************************************************************************************
get_Type_Inter2_List(): Reads the pair interaction from the parameter file. Returns
a pointer to the data structure Inter2_Type_List which contains the information read
from the parameter file.
*****************************************************************************************/
tW_type_inter2 *get_Type_Inter2_List(FILE * fp_par, tW_line inp_line,
				     tW_system sys, int N_Inter_Types,
				     tW_word keyword)
{
    int i_inter;
    tW_type_inter2 *Inter2_Type_List;

    Inter2_Type_List =
	(tW_type_inter2 *) emalloc(N_Inter_Types * sizeof(tW_type_inter2));

    for (i_inter = 0; i_inter < N_Inter_Types; i_inter++) {
	get_next_line(fp_par, inp_line);

	check_inp_line("par.txt", keyword, inp_line);

	read_Type_Inter2(inp_line, &(Inter2_Type_List[i_inter]), sys);
    }

    get_next_line(fp_par, inp_line);

    test_end_directive(inp_line, keyword, "par.txt");

    return Inter2_Type_List;
}


/*****************************************************************************************
read_Type_Inter2(): Finds the interaction type store in sys and copies that information
for a given pair of site types. Execution is stopped if the interaction name found for
a pair of sites listed in the parameter file is not a pair interaction (in sys) or the
interaction name is not found in the interaction types listed in sys.
*****************************************************************************************/
int read_Type_Inter2(tW_line inp_line, tW_type_inter2 * inter,
		     tW_system sys)
{
    int i;
    int flag = FALSE;

    sscanf(inp_line, "%s %s %s", inter->inter_name, inter->name1,
	   inter->name2);

    for (i = 0; i < sys.N_Inter_Types; i++) {
	if (strcmp(inter->inter_name, sys.Inter_Types[i].inter_name) == 0) {
	    strcpy(inter->basis, sys.Inter_Types[i].basis);
	    inter->N_inter = 0;
	    inter->N_pts = sys.Inter_Types[i].N_pts;
	    inter->N_coeff = sys.Inter_Types[i].N_coeff;
	    inter->i_basis = sys.Inter_Types[i].i_basis;
	    inter->i_0 = sys.Inter_Types[i].i_0;
	    inter->dr = sys.Inter_Types[i].dr;
	    inter->R_0 = sys.Inter_Types[i].R_min;
	    inter->R_max = sys.Inter_Types[i].R_max;
	    inter->n_smooth = sys.Inter_Types[i].n_smooth;
	    inter->N_powers = sys.Inter_Types[i].N_powers;
	    inter->powers = sys.Inter_Types[i].powers;
	    inter->kspline = sys.Inter_Types[i].kspline;	/* JFR - 07.22.12: Bspline order */

	    flag = TRUE;

	    if (strcmp(sys.Inter_Types[i].inter_type, NB_PAIR) != 0) {
		printf
		    ("\nERROR: %s is not of type %s.\n  Check par.txt.\n",
		     inter->inter_name, NB_PAIR);
		exit(EXIT_FAILURE);
	    }

	    break;
	}
    }

    if (flag == FALSE) {
	printf
	    ("\nERROR: Interaction Type \"%s\" not found.\n  Check par.txt.",
	     inter->inter_name);
	exit(EXIT_FAILURE);
    }

    return 0;
}


/*****************************************************************************************
setup_tables(): Allocates memories for the grids. Zeros the grids just in case.
*****************************************************************************************/
int setup_tables(tW_system * sys)
{
    int i;
    int i_start = 0;
    int N_coeff_tmp = 0;

    /* Determine the total number of coefficients. */
    for (i = 0; i < sys->N_Inter_Types; i++) {
	N_coeff_tmp += sys->Inter_Types[i].N_coeff;
    }
    sys->N_coeff = N_coeff_tmp;
    sys->N_pack = (N_coeff_tmp * N_coeff_tmp + N_coeff_tmp) / 2;

    /* Allocate memory for master/system wide grids. */
    sys->b = (double *) ecalloc(N_coeff_tmp, sizeof(double));
    sys->b_ref = (double *) ecalloc(N_coeff_tmp, sizeof(double));
    sys->g = (double *) ecalloc(N_coeff_tmp, sizeof(double));
    sys->g_cnt = (double *) ecalloc(N_coeff_tmp, sizeof(double));
    sys->L = (double *) ecalloc(N_coeff_tmp, sizeof(double));
    // JFR - added 04.11.12: put the matrix in packed form
    sys->M = (double *) ecalloc(sys->N_pack, sizeof(double));
    sys->M2 = (double *) ecalloc(sys->N_pack, sizeof(double));

    if (sys->MT_var.flag_Mcnt == TRUE) {
	sys->M_cnt = (double *) ecalloc(sys->N_pack, sizeof(double));
    } else {
	sys->M_cnt = NULL;
    }

    sys->d2b = (double *) ecalloc(N_coeff_tmp, sizeof(double));	/* JFR - 01.31.13 */
    sys->d2M = (double *) ecalloc(sys->N_pack, sizeof(double));	/* JFR - 01.31.13 */
    sys->rescale = (double *) ecalloc(N_coeff_tmp, sizeof(double));

    /* NJD - Here, we should conditionally allocate memory for framewise copies of the above arrays. */

    if (sys->FRAMEWEIGHT_var.flag_FRAMEWEIGHT == TRUE ) {
	    sys->b_wt = (double *) ecalloc(N_coeff_tmp, sizeof(double));		
	    sys->M_wt = (double *) ecalloc(sys->N_pack, sizeof(double));
	    sys->M2_wt = (double *) ecalloc(sys->N_pack, sizeof(double));
	    sys->wt_norm = 0;
    }

    /* These are not needed until end of calculation. */
    sys->x = NULL;
    sys->b_forces = NULL;
    sys->b_struct = NULL;
    sys->phi = NULL;
    sys->phi_forces = NULL;
    sys->phi_struct = NULL;

    /* Determine the initial grid index for each interaction type. */
    for (i = 0; i < sys->N_Inter_Types; i++) {
	sys->Inter_Types[i].i_0 = i_start;
	i_start += sys->Inter_Types[i].N_coeff;
    }

    /* Setup pointer to master grid for all interaction types. */
    for (i = 0; i < sys->N_Inter_Types; i++) {
	i_start = sys->Inter_Types[i].i_0;
	sys->Inter_Types[i].ptr_b = &(sys->b[i_start]);
	sys->Inter_Types[i].ptr_b_ref = &(sys->b_ref[i_start]);
	sys->Inter_Types[i].ptr_g = &(sys->g[i_start]);
	sys->Inter_Types[i].ptr_g_cnt = &(sys->g_cnt[i_start]);
	sys->Inter_Types[i].ptr_L = &(sys->L[i_start]);
    }

    /* Setup pointers to master grid for non-bonded and bonded interactions. */
    for (i = 0; i < sys->N_Inter2_Types; i++) {
	i_start =
	    get_interaction_i_0(sys->Inter2_Type_List[i].inter_name, sys);
	sys->Inter2_Type_List[i].i_0 = i_start;
	sys->Inter2_Type_List[i].ptr_b = &(sys->b[i_start]);
	sys->Inter2_Type_List[i].ptr_b_ref = &(sys->b_ref[i_start]);
	sys->Inter2_Type_List[i].ptr_g2 = &(sys->g[i_start]);
	sys->Inter2_Type_List[i].ptr_g2_cnt = &(sys->g_cnt[i_start]);
	sys->Inter2_Type_List[i].ptr_L = &(sys->L[i_start]);
    }

    for (i = 0; i < sys->N_Bond_Int_Types; i++) {
	i_start =
	    get_interaction_i_0(sys->Bonded_Inter_Types[i].inter_name,
				sys);
	sys->Bonded_Inter_Types[i].i_0 = i_start;
	sys->Bonded_Inter_Types[i].ptr_b = &(sys->b[i_start]);
	sys->Bonded_Inter_Types[i].ptr_b_ref = &(sys->b_ref[i_start]);
	sys->Bonded_Inter_Types[i].ptr_g = &(sys->g[i_start]);
	sys->Bonded_Inter_Types[i].ptr_g_cnt = &(sys->g_cnt[i_start]);
	sys->Bonded_Inter_Types[i].ptr_L = &(sys->L[i_start]);
    }

    /* Just in case. */
    clear_sys_arrays(sys);

    return N_coeff_tmp;
}


/*****************************************************************************************
summarize_input(): Summarizes the input read from the parameter files. 
*****************************************************************************************/
void summarize_input(tW_files files, tW_system sys,
		     tW_ref_potential ref_potential)
{
    int i;
    FILE *fp_log = fopen(SUM_INPUT_FNAME, "w");

    print_line_stars(fp_log);

    /* Print temperature, nrexcl, and mode. */
    fprintf(fp_log, "  Mode: %s\n", files.mode);
    fprintf(fp_log, "  Temperature: %f\n", sys.Temperature);
    if (strcmp(files.mode, PDB_MODE) == 0) {
	fprintf(fp_log, "  NREXCL: %d", sys.nrexcl);
	if (sys.nrexcl == -1) {
	    fprintf(fp_log, "  (undeclared).");
	}
    }

    /* MRD 03.04.2019 */
    if (sys.SKIP_TRIPLE_LOOP) { fprintf(fp_log,"  Evaluating the G matrix WITHOUT the triple loop\n"); }
    else { fprintf(fp_log,"  Evaluating the G matrix the old way (using the triple loop)\n"); }

    fprintf(fp_log, "\n\n");
    print_line_stars(fp_log);


    /* Print total number of coefficients. */
    fprintf(fp_log, "  N_coeff: %d\n\n", sys.N_coeff);
    print_line_stars(fp_log);

    /* List files that contain structures. */
    fprintf(fp_log, "  N_structures: %d\n", files.N_struct);
    for (i = 0; i < files.N_struct; i++) {
	fprintf(fp_log, "    Structure: %d\t Weight: %lf\t Name: %s\n",
		i + 1, files.p_struct[i], files.structures[i]);
    }
    fprintf(fp_log, "\n");
    print_line_stars(fp_log);

    /* List site types. */
    fprintf(fp_log, "  N_Site_Types: %d\n", sys.N_Site_Types);
    for (i = 0; i < sys.N_Site_Types; i++) {
	fprintf(fp_log, "    Site_Type: %d\t name: %s\n", i + 1,
		sys.Site_Types[i]);
    }
    fprintf(fp_log, "\n");
    print_line_stars(fp_log);

    /* List Interaction Types. */
    summarize_Inter_Types(fp_log, sys);
    print_line_stars(fp_log);

    summarize_input_pair_inter(fp_log, sys);
    print_line_stars(fp_log);

    summarize_input_bond_inter(fp_log, sys, B_BOND_STRETCH, 2, files.mode);
    print_line_stars(fp_log);

    summarize_input_bond_inter(fp_log, sys, B_ANGLE, 3, files.mode);
    print_line_stars(fp_log);

    summarize_input_bond_inter(fp_log, sys, B_DIHEDRAL, 4, files.mode);
    print_line_stars(fp_log);

    summarize_input_bond_inter(fp_log, sys, B_NB_PAIR_BOND, 2, files.mode);
    print_line_stars(fp_log);

    if (sys.flag_ref_potential == 1) {
	summarize_input_ref_potential(fp_log, sys, ref_potential);
    }

    if (sys.PT_var.flag_PT == TRUE) {
	summarize_input_PT(fp_log, sys);
    }

    if (sys.Eigen_var.flag_Eigen == TRUE) {
	summarize_input_Eigen(fp_log, sys);
    }

    if (sys.SVD_var.flag_SVD == TRUE) {
	summarize_input_SVD(fp_log, sys);
    }

    if (sys.TPR_var.flag_TPR == TRUE) {
	summarize_input_TPR(fp_log, sys);
    }

    if (sys.TPR_EXCL_var.flag_TPR_excl == TRUE) {
	summarize_input_TPR_EXCL(fp_log, sys);
    }

    if (sys.MT_var.flag_MT == TRUE) {
	summarize_input_MT(fp_log, sys);
    }

    if (sys.MFD_var.flag_MFD == TRUE) {
	summarize_input_MFD(fp_log, sys);
    }

    if (sys.CalcMODE_var.flag_CalcMODE == TRUE) {
	summarize_input_CalcMODE(fp_log, sys);
    }

    if (sys.PC_var.flag_PC == TRUE) {
	summarize_input_PC(fp_log, sys);
    }

    if (sys.MEM_var.flag_MEM == TRUE) {
	summarize_input_MEM(fp_log, sys);
    }

    if (sys.SOLN_var.flag_SOLN == TRUE) {
	summarize_input_SOLN(fp_log, sys);
    }

    if (sys.FRAMEWEIGHT_var.flag_FRAMEWEIGHT == TRUE ) {
	summarize_input_FRAMEWEIGHT(fp_log, sys);
    }

    if (sys.ERR_var.flag_ERR == TRUE) {
	summarize_input_ERR(fp_log, sys);
    }

    if (sys.REF_var.flag_REF == TRUE) {
	summarize_input_REF(fp_log, sys);
    }

    if (sys.TRIM_var.flag_TRIM == TRUE) {
	summarize_input_TRIM(fp_log, sys);
    }

    if (sys.CHISQD_var.flag_CHISQD == TRUE) {
	summarize_input_CHISQD(fp_log, sys);
    }

    if (sys.REG_var.flag_REG == TRUE) {
	summarize_input_REG(fp_log, sys);
    }

    if (sys.RESCALE_var.flag_RESCALE == TRUE) {
	summarize_input_RESCALE(fp_log, sys);
    }

    if (sys.CONSTRAIN_var.flag_CONSTRAIN == TRUE) {
	summarize_input_CONSTRAIN(fp_log, sys);
    }

    if (sys.ITER_var.flag_ITER == TRUE) {
	summarize_input_ITER(fp_log, sys);
    }

    fclose(fp_log);
}

/*****************************************************************************************
summarize_input_PT():
*****************************************************************************************/
void summarize_input_PT(FILE * fp_log, tW_system sys)
{

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing PT parameters ... ");

    fprintf(fp_log, "    %s: %d \n", "N_PT", sys.PT_var.N_PT);

    fprintf(fp_log, "    %s: %d \n", "dPT", sys.PT_var.dPT);

    fprintf(fp_log, "    %s: %d \n", "flag_MMOTF_SEP",
	    sys.PT_var.flag_MMOTF_SEP);

    fprintf(fp_log, "    %s: %d \n", "flag_eigen", sys.PT_var.flag_eigen);

    fprintf(fp_log, "    %s: %d \n", "N_eigen", sys.PT_var.N_eigen);

}

/*****************************************************************************************
summarize_input_Eigen():
*****************************************************************************************/
void summarize_input_Eigen(FILE * fp_log, tW_system sys)
{

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Eigen parameters ... ");

    fprintf(fp_log, "    %s: %d \n", "N_Eigen", sys.Eigen_var.N_Eigen);

    fprintf(fp_log, "    %s: %d \n", "flag_printn",
	    sys.Eigen_var.flag_printn);

    fprintf(fp_log, "    %s: %d \n", "DM", sys.Eigen_var.DM);

    fprintf(fp_log, "    %s: %d \n", "DL", sys.Eigen_var.DL);

    fprintf(fp_log, "    %s: %d \n", "flag_Gbar", sys.Eigen_var.flag_Gbar);

    fprintf(fp_log, "    %s: %d \n", "flag_norm", sys.Eigen_var.flag_norm);

}

/*****************************************************************************************
summarize_input_SVD():
*****************************************************************************************/
void summarize_input_SVD(FILE * fp_log, tW_system sys)
{
    tW_word tmp;

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing SVD parameters ... ");

    fprintf(fp_log, "    %s: %lf \n", "rcond", sys.SVD_var.rcond);

    if (sys.SVD_var.flag_printSV == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_printSV", tmp);

    if (sys.SVD_var.flag_printevecs == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_printevecs", tmp);

    if (sys.SVD_var.flag_solve == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_solve", tmp);

}

/*****************************************************************************************
summarize_input_TPR():
*****************************************************************************************/
void summarize_input_TPR(FILE * fp_log, tW_system sys)
{

    int i;

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing TPR filenames ... ");

    fprintf(fp_log, "    %s: %d \n", "N_TPR", sys.TPR_var.N_TPR);

    for (i = 0; i < sys.TPR_var.N_TPR; i++) {

	fprintf(fp_log, "    %s[%d]: %s \n", "TPR", i,
		sys.TPR_var.TPR_files[i]);

    }

}

/*****************************************************************************************
summarize_input_TPR_EXCL():
*****************************************************************************************/
void summarize_input_TPR_EXCL(FILE * fp_log, tW_system sys)
{

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing TPR_EXCL filename ... ");

    fprintf(fp_log, "    %s: %s \n", "TPR_EXCL",
	    sys.TPR_EXCL_var.TPR_excl);

}

/*****************************************************************************************
summarize_input_MT():
*****************************************************************************************/
void summarize_input_MT(FILE * fp_log, tW_system sys)
{
    tW_word tmp;

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Metric_Tensor parameters ... ");

    if (sys.MT_var.flag_print == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }

    if (sys.MT_var.flag_norm == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_norm", tmp);

    if (sys.MT_var.flag_Mcnt == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_Mcnt", tmp);

}

/*****************************************************************************************
summarize_input_MFD():
*****************************************************************************************/
void summarize_input_MFD(FILE * fp_log, tW_system sys)
{

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n",
	    "printing Mean_Force_Decomposition parameters ... ");

}

/*****************************************************************************************
summarize_input_CalcMODE():
*****************************************************************************************/
void summarize_input_CalcMODE(FILE * fp_log, tW_system sys)
{
    tW_word mode;

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n",
	    "printing Calculation_Mode parameters ... ");

    if (sys.CalcMODE_var.CalcMODE == IFULL) {
	strcpy(mode, FULL);
    } else if (sys.CalcMODE_var.CalcMODE == IFIRST_HALF) {
	strcpy(mode, FIRST_HALF);
    } else {
	strcpy(mode, SECOND_HALF);
    }

    fprintf(fp_log, "    %s: %s \n", "CalcMODE", mode);

}

/*****************************************************************************************
summarize_input_PC():
*****************************************************************************************/
void summarize_input_PC(FILE * fp_log, tW_system sys)
{
    tW_word tmp;
    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Preconditioning options ... ");

    fprintf(fp_log, "    %s: %s \n", "RPC_type", sys.PC_var.RPC);

    if (sys.PC_var.flag_normb == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_normb", tmp);

    fprintf(fp_log, "    %s: %s \n", "LPC_type", sys.PC_var.LPC);

    if (sys.PC_var.flag_normphi == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_normphi", tmp);

}

/*****************************************************************************************
summarize_input_MEM():
*****************************************************************************************/
void summarize_input_MEM(FILE * fp_log, tW_system sys)
{

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Memory parameters ... ");

    fprintf(fp_log, "    %s: %d \n", "flag_LOWMEM",
	    sys.MEM_var.flag_LOWMEM);

    fprintf(fp_log, "    %s: %s \n", "info", sys.MEM_var.info);

    fprintf(fp_log, "    %s: %d \n", "flag_mult_top",
	    sys.MEM_var.flag_mult_top);

}

/*****************************************************************************************
summarize_input_SOLN():
*****************************************************************************************/
void summarize_input_SOLN(FILE * fp_log, tW_system sys)
{
    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Solution Method options ... ");

    fprintf(fp_log, "    %s: %s \n", "SOLN_METH", sys.SOLN_var.SOLN_METH);

}

/*****************************************************************************************
summarize_input_EMSEMBLE():
*****************************************************************************************/
void summarize_input_FRAMEWEIGHT(FILE * fp_log, tW_system sys)
{
	print_line_stars(fp_log);

	fprintf(fp_log,"  %s \n", "printing Frame Weighting options ... ");

	fprintf(fp_log, "    %s: %s \n", "FRAMEWEIGHT", sys.FRAMEWEIGHT_var.FRAMEWEIGHT);
}


/*****************************************************************************************
summarize_input_ERR():
*****************************************************************************************/
void summarize_input_ERR(FILE * fp_log, tW_system sys)
{
    tW_word tmp;

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Error options ... ");

    if (sys.ERR_var.FACT == 'E') {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_Equil", tmp);

}

/*****************************************************************************************
summarize_input_REF():
*****************************************************************************************/
void summarize_input_REF(FILE * fp_log, tW_system sys)
{

    int i;

    tW_word tmp;

    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Reference options ... ");

    if (sys.REF_var.flag_calcbref == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_calcbref", tmp);

    if (sys.REF_var.flag_readbref == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_readbref", tmp);

    if (sys.REF_var.flag_splitfiles == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_splitfiles", tmp);

    if (sys.REF_var.flag_reftrr == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_reftrr", tmp);

    fprintf(fp_log, "    %s: %d \n", "N_fnm", sys.REF_var.N_fnm);

    for (i = 0; i < sys.REF_var.N_fnm; i++) {

	fprintf(fp_log, "    %s[%d]: %s \n", "REF", i,
		sys.REF_var.reftrr_fnm[i]);

    }


}

/*****************************************************************************************
summarize_input_TRIM():
*****************************************************************************************/
void summarize_input_TRIM(FILE * fp_log, tW_system sys)
{
    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Trim options ... ");

    fprintf(fp_log, "    %s: %lf \n", "FE", sys.TRIM_var.FE);

}

/*****************************************************************************************
summarize_input_CHISQD():
*****************************************************************************************/
void summarize_input_CHISQD(FILE * fp_log, tW_system sys)
{
    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing CHISQD options ... ");

    fprintf(fp_log, "    %s: %s \n", "force_fnm", sys.CHISQD_var.force_fnm);

}

/*****************************************************************************************
summarize_input_REG():
*****************************************************************************************/
void summarize_input_REG(FILE * fp_log, tW_system sys)
{
    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Regularization options ... ");

    fprintf(fp_log, "    %s: %s \n", "type", sys.REG_var.type);

    fprintf(fp_log, "    %s: %d \n", "Nmax", sys.REG_var.Nmax);

    fprintf(fp_log, "    %s: %lf \n", "tau_alpha", sys.REG_var.tau_alpha);

    fprintf(fp_log, "    %s: %lf \n", "tau_beta", sys.REG_var.tau_beta);

}

/*****************************************************************************************
summarize_input_RESCALE():
*****************************************************************************************/
void summarize_input_RESCALE(FILE * fp_log, tW_system sys)
{
    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Rescale Forces options ... ");

    fprintf(fp_log, "    %s: %d \n", "Nmax", sys.RESCALE_var.Nmax);

    fprintf(fp_log, "    %s: %lf \n", "tau_phi", sys.RESCALE_var.tau_phi);

}

/*****************************************************************************************
summarize_input_CONSTRAIN():
*****************************************************************************************/
void summarize_input_CONSTRAIN(FILE * fp_log, tW_system sys)
{
    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n",
	    "printing Constrain Dihedrals options ... ");

    fprintf(fp_log, "    %s: %d \n", "Nmax", sys.CONSTRAIN_var.Nmax);

    fprintf(fp_log, "    %s: %lf \n", "tau_phi",
	    sys.CONSTRAIN_var.tau_dih);

    fprintf(fp_log, "    %s: %lf \n", "lambda", sys.CONSTRAIN_var.lambda);

    fprintf(fp_log, "    %s: %lf \n", "dlambda",
	    sys.CONSTRAIN_var.dlambda);

}

/*****************************************************************************************
summarize_input_ITER():
*****************************************************************************************/
void summarize_input_ITER(FILE * fp_log, tW_system sys)
{
    tW_word tmp;
    print_line_stars(fp_log);

    fprintf(fp_log, "  %s \n", "printing Iter-gYBG options ... ");

    fprintf(fp_log, "    %s: %d \n", "Nmax", sys.ITER_var.flag_AAM2);

    fprintf(fp_log, "    %s: %s \n", "AAM2_fnm", sys.ITER_var.AAM2_fnm);

    if (sys.ITER_var.flag_bsolnerr == TRUE) {
	strcpy(tmp, "YES");
    } else {
	strcpy(tmp, "NO");
    }
    fprintf(fp_log, "    %s: %s \n", "flag_bsolnerr", tmp);

}

/*****************************************************************************************
summarize_input_bond_inter(): Summarizes the information stored for the intramolecular
interactions.
*****************************************************************************************/
void summarize_input_bond_inter(FILE * fp_log, tW_system sys,
				tW_word inter_name, int N_Int_Sites,
				tW_word mode)
{
    int i, j;
    int N_inter = 0;
    int *inter_list = NULL;
    tW_Bonded_Inter *tmp_Inter_Bond;

    /* Get indices for this bond interactions. */
    for (i = 0; i < sys.N_Bond_Int_Types; i++) {
	if (strcmp(inter_name, sys.Bonded_Inter_Types[i].name) == 0) {
	    N_inter++;
	    inter_list =
		(int *) erealloc(inter_list, (N_inter * sizeof(int)));
	    inter_list[N_inter - 1] = i;
	}
    }

    fprintf(fp_log, "  N_%s: %d \n", inter_name, N_inter);

    for (i = 0; i < N_inter; i++) {
	fprintf(fp_log, "    %s: %d \n", inter_name, i + 1);

	tmp_Inter_Bond = &(sys.Bonded_Inter_Types[inter_list[i]]);

	fprintf(fp_log, "    Name: %s", tmp_Inter_Bond->name);

	fprintf(fp_log, "    Inter_Name: %s", tmp_Inter_Bond->inter_name);

	fprintf(fp_log, "    Site Types: ");
	for (j = 0; j < N_Int_Sites; j++) {
	    fprintf(fp_log, "%s", tmp_Inter_Bond->Site_Types[j]);
	    if (j != (N_Int_Sites - 1)) {
		fprintf(fp_log, "-");
	    }
	}

	if ((strcmp(tmp_Inter_Bond->name, B_NB_PAIR_BOND) == 0) &&
	    (strcmp(mode, PDB_MODE) == 0)
	    ) {
	    fprintf(fp_log, "    Site Positions: ");
	    fprintf(fp_log, "1-%d",
		    1 + tmp_Inter_Bond->n_bonds_intramolec_pair_inter);
	}

	fprintf(fp_log, "\n");

	fprintf(fp_log, "    Basis: %s", tmp_Inter_Bond->basis);

	fprintf(fp_log, "    i_basis: %d    i_0: %d",
		tmp_Inter_Bond->i_basis, tmp_Inter_Bond->i_0);

	fprintf(fp_log, "    N_pts: %d    N_coeff: %d \n",
		tmp_Inter_Bond->N_pts, tmp_Inter_Bond->N_coeff);

	if (tmp_Inter_Bond->N_powers > 0) {
	    fprintf(fp_log, "    R_max: %f  N_powers: %d  Powers: ",
		    tmp_Inter_Bond->R_max, tmp_Inter_Bond->N_powers);
	    for (j = 0; j < tmp_Inter_Bond->N_powers; j++) {
		fprintf(fp_log, "%d", tmp_Inter_Bond->powers[j]);
		if (j < tmp_Inter_Bond->N_powers - 1) {
		    fprintf(fp_log, ", ");
		} else {
		    fprintf(fp_log, "\n");
		}
	    }
	} else {
	    fprintf(fp_log,
		    "    dr: %f     R_0: %f     R_max: %f     n_smooth: %d \n",
		    tmp_Inter_Bond->dr, tmp_Inter_Bond->R_0,
		    tmp_Inter_Bond->R_max, tmp_Inter_Bond->n_smooth);
	}

	if (tmp_Inter_Bond->i_basis == BSPLINE_BASIS_INDEX) {
	    fprintf(fp_log, "   kspline: %d\n", tmp_Inter_Bond->kspline);
	}
	/* JFR - 07.22.12: Bspline order */
	fprintf(fp_log, "\n");
    }

    if (N_inter == 0) {
	fprintf(fp_log, "\n");
    }

    free(inter_list);
}


/*****************************************************************************************
summarize_input_pair_inter(): Summarizes the information stored for each intermolecular
pair interaction.
*****************************************************************************************/
void summarize_input_pair_inter(FILE * fp_log, tW_system sys)
{
    int i, j;
    tW_type_inter2 *tmp_Inter2;

    /* List pair interactions. */
    fprintf(fp_log, "  N_Pair_Interactions: %d", sys.N_Inter2_Types);
    for (i = 0; i < sys.N_Inter2_Types; i++) {
	fprintf(fp_log, "\n    Pair_Interaction: %d\n", i + 1);
	tmp_Inter2 = &sys.Inter2_Type_List[i];

	fprintf(fp_log, "    Inter_Name: %s  ", tmp_Inter2->inter_name);
	fprintf(fp_log, "    Sites: %s-%s\n", tmp_Inter2->name1,
		tmp_Inter2->name2);

	fprintf(fp_log, "    Basis: %s", tmp_Inter2->basis);
	fprintf(fp_log, "    i_basis: %d    i_0: %d \n",
		tmp_Inter2->i_basis, tmp_Inter2->i_0);

	fprintf(fp_log, "    N_pts: %d \t N_coeff: %d \n",
		tmp_Inter2->N_pts, tmp_Inter2->N_coeff);

	if (tmp_Inter2->N_powers > 0) {
	    fprintf(fp_log, "    R_max: %f  N_powers: %d  Powers: ",
		    tmp_Inter2->R_max, tmp_Inter2->N_powers);
	    for (j = 0; j < tmp_Inter2->N_powers; j++) {
		fprintf(fp_log, "%d", tmp_Inter2->powers[j]);
		if (j < tmp_Inter2->N_powers - 1) {
		    fprintf(fp_log, ", ");
		} else {
		    fprintf(fp_log, "\n");
		}
	    }
	} else {
	    fprintf(fp_log,
		    "    dr: %lf   R_0: %lf   R_max: %lf   n_smooth: %d\n",
		    tmp_Inter2->dr, tmp_Inter2->R_0, tmp_Inter2->R_max,
		    tmp_Inter2->n_smooth);
	}

	if (tmp_Inter2->i_basis == BSPLINE_BASIS_INDEX) {
	    fprintf(fp_log, "    kspline: %d\n", tmp_Inter2->kspline);
	}			/* JFR - 07.22.12: Bspline order */
    }

    fprintf(fp_log, "\n");
}


/*****************************************************************************************
setup_sys_copy(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void setup_sys_copy(tW_system sys_orig, tW_system * sys_copy, int overwrite_arrays)
{
    /* Temperature. */
    sys_copy->Temperature = sys_orig.Temperature;

    /* nrexcl. */
    sys_copy->nrexcl = sys_orig.nrexcl;

    /* Flag for reference potential. */
    sys_copy->flag_ref_potential = sys_orig.flag_ref_potential;

    /* Copy site type information. */
    copy_site_types(sys_orig, sys_copy);

    if (overwrite_arrays)
    {
	    /* Allocates memory for arrays needed in calc_grids(). Sets other arrays to NULL. */
	    copy_arrays(sys_orig, sys_copy);
    }

    /* Allocates memory and copies information for all interaction types. */
    copy_Inter_Types(sys_orig, sys_copy);

    /* Allocates memory and copies information for all nb pair interactions. */
    copy_Inter2_Types(sys_orig, sys_copy);

    /* Allocates memory and copies information for all bond interactions. */
    copy_Bonded_Inter_Types(sys_orig, sys_copy);

    /* Copy PT variables */
    copy_PT(sys_orig, sys_copy);

    /* Copy Eigen variables */
    copy_Eigen(sys_orig, sys_copy);

    /* Copy SVD variables */
    copy_SVD(sys_orig, sys_copy);

    /* Copy TPR variables */
    copy_TPR(sys_orig, sys_copy);

    /* Copy TPR_EXCL variables */
    copy_TPR_EXCL(sys_orig, sys_copy);

    /* Copy MT variables */
    copy_MT(sys_orig, sys_copy);

    /* Copy MFD variables */
    copy_MFD(sys_orig, sys_copy);

    /* Copy CalcMODE variables */
    copy_CalcMODE(sys_orig, sys_copy);

    /* Copy PC variables */
    copy_PC(sys_orig, sys_copy);

    /* Copy MEM variables */
    copy_MEM(sys_orig, sys_copy);

    /* Copy SOLN variables */
    copy_SOLN(sys_orig, sys_copy);

    /* Copy FRAMEWEIGHT variables */
    copy_FRAMEWEIGHT(sys_orig, sys_copy);

    /* Copy ERR variables */
    copy_ERR(sys_orig, sys_copy);

    /* Copy REF variables */
    copy_REF(sys_orig, sys_copy);

    /* Copy TRIM variables */
    copy_TRIM(sys_orig, sys_copy);

    /* Copy CHISQD variables */
    copy_CHISQD(sys_orig, sys_copy);

    /* Copy REG variables */
    copy_REG(sys_orig, sys_copy);

    /* Copy RESCALE variables */
    copy_RESCALE(sys_orig, sys_copy);

    /* Copy CONSTRAIN variables */
    copy_CONSTRAIN(sys_orig, sys_copy);

    /* Copy ITER variables */
    copy_ITER(sys_orig, sys_copy);

    /* JFR - 06.27.12: Chi2 */
    sys_copy->Chi2 = sys_orig.Chi2;

    /* MRD 03.04.2019 */
    sys_copy->SKIP_TRIPLE_LOOP = sys_orig.SKIP_TRIPLE_LOOP;
    sys_copy->M_M2_proc = sys_orig.M_M2_proc;
}

/*****************************************************************************************
copy_PT(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_PT(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->PT_var.flag_PT = sys_orig.PT_var.flag_PT;

    sys_copy->PT_var.N_PT = sys_orig.PT_var.N_PT;

    sys_copy->PT_var.dPT = sys_orig.PT_var.dPT;

    sys_copy->PT_var.flag_MMOTF_SEP = sys_orig.PT_var.flag_MMOTF_SEP;

    sys_copy->PT_var.flag_eigen = sys_orig.PT_var.flag_eigen;

    sys_copy->PT_var.N_eigen = sys_orig.PT_var.N_eigen;

}

/*****************************************************************************************
copy_Eigen(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_Eigen(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->Eigen_var.flag_Eigen = sys_orig.Eigen_var.flag_Eigen;

    sys_copy->Eigen_var.N_Eigen = sys_orig.Eigen_var.N_Eigen;

    sys_copy->Eigen_var.flag_printn = sys_orig.Eigen_var.flag_printn;

    sys_copy->Eigen_var.DM = sys_orig.Eigen_var.DM;

    sys_copy->Eigen_var.DL = sys_orig.Eigen_var.DL;

    sys_copy->Eigen_var.flag_Gbar = sys_orig.Eigen_var.flag_Gbar;

    sys_copy->Eigen_var.flag_norm = sys_orig.Eigen_var.flag_norm;

}

/*****************************************************************************************
copy_SVD(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_SVD(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->SVD_var.flag_SVD = sys_orig.SVD_var.flag_SVD;

    sys_copy->SVD_var.rcond = sys_orig.SVD_var.rcond;

    sys_copy->SVD_var.flag_printSV = sys_orig.SVD_var.flag_printSV;

    sys_copy->SVD_var.flag_printevecs = sys_orig.SVD_var.flag_printevecs;

    sys_copy->SVD_var.flag_solve = sys_orig.SVD_var.flag_solve;

}

/*****************************************************************************************
copy_TPR(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_TPR(tW_system sys_orig, tW_system * sys_copy)
{
    int i;

    sys_copy->TPR_var.flag_TPR = sys_orig.TPR_var.flag_TPR;

    sys_copy->TPR_var.N_TPR = sys_orig.TPR_var.N_TPR;

    sys_copy->TPR_var.TPR_files = (tW_word *) ecalloc(sys_orig.TPR_var.N_TPR, sizeof(tW_word));	//JFR - added 04.06.12

    for (i = 0; i < sys_orig.TPR_var.N_TPR; i++) {
	strcpy(sys_copy->TPR_var.TPR_files[i],
	       sys_orig.TPR_var.TPR_files[i]);
    }

}

/*****************************************************************************************
copy_TPR_EXCL(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_TPR_EXCL(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->TPR_EXCL_var.flag_TPR_excl = sys_orig.TPR_EXCL_var.flag_TPR_excl;

    strcpy(sys_copy->TPR_EXCL_var.TPR_excl, sys_orig.TPR_EXCL_var.TPR_excl);

}

/*****************************************************************************************
copy_MT(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_MT(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->MT_var.flag_MT = sys_orig.MT_var.flag_MT;

    sys_copy->MT_var.flag_print = sys_orig.MT_var.flag_print;

    sys_copy->MT_var.flag_norm = sys_orig.MT_var.flag_norm;

    sys_copy->MT_var.flag_Mcnt = sys_orig.MT_var.flag_Mcnt;

}

/*****************************************************************************************
copy_MFD(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_MFD(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->MFD_var.flag_MFD = sys_orig.MFD_var.flag_MFD;

}

/*****************************************************************************************
copy_CalcMODE(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_CalcMODE(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->CalcMODE_var.flag_CalcMODE = sys_orig.CalcMODE_var.flag_CalcMODE;

    sys_copy->CalcMODE_var.CalcMODE = sys_orig.CalcMODE_var.CalcMODE;

}

/*****************************************************************************************
copy_PC(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_PC(tW_system sys_orig, tW_system * sys_copy)
{
    sys_copy->PC_var.flag_PC = sys_orig.PC_var.flag_PC;

    strcpy(sys_copy->PC_var.RPC, sys_orig.PC_var.RPC);

    sys_copy->PC_var.flag_normb = sys_orig.PC_var.flag_normb;

    strcpy(sys_copy->PC_var.LPC, sys_orig.PC_var.LPC);

    sys_copy->PC_var.flag_normphi = sys_orig.PC_var.flag_normphi;
}

/*****************************************************************************************
copy_MEM(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_MEM(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->MEM_var.flag_MEM = sys_orig.MEM_var.flag_MEM;

    sys_copy->MEM_var.flag_LOWMEM = sys_orig.MEM_var.flag_LOWMEM;

    strcpy(sys_copy->MEM_var.info, sys_orig.MEM_var.info);

    sys_copy->MEM_var.flag_mult_top = sys_orig.MEM_var.flag_mult_top;

}

/*****************************************************************************************
copy_SOLN(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_SOLN(tW_system sys_orig, tW_system * sys_copy)
{
    sys_copy->SOLN_var.flag_SOLN = sys_orig.SOLN_var.flag_SOLN;

    strcpy(sys_copy->SOLN_var.SOLN_METH, sys_orig.SOLN_var.SOLN_METH);
}

/*****************************************************************************************
copy_FRAMEWEIGHT(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_FRAMEWEIGHT(tW_system sys_orig, tW_system * sys_copy)
{
    sys_copy->FRAMEWEIGHT_var.flag_FRAMEWEIGHT = sys_orig.FRAMEWEIGHT_var.flag_FRAMEWEIGHT;

    sys_copy->wt_norm = sys_orig.wt_norm;

    strcpy(sys_copy->FRAMEWEIGHT_var.FRAMEWEIGHT, sys_orig.FRAMEWEIGHT_var.FRAMEWEIGHT);
}

/*****************************************************************************************
copy_ERR(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_ERR(tW_system sys_orig, tW_system * sys_copy)
{
    sys_copy->ERR_var.flag_ERR = sys_orig.ERR_var.flag_ERR;

    sys_copy->ERR_var.FACT = sys_orig.ERR_var.FACT;
}

/*****************************************************************************************
copy_REF(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_REF(tW_system sys_orig, tW_system * sys_copy)
{
    int i;

    sys_copy->REF_var.flag_REF = sys_orig.REF_var.flag_REF;

    sys_copy->REF_var.flag_calcbref = sys_orig.REF_var.flag_calcbref;

    sys_copy->REF_var.flag_readbref = sys_orig.REF_var.flag_readbref;

    sys_copy->REF_var.flag_splitfiles = sys_orig.REF_var.flag_splitfiles;

    sys_copy->REF_var.flag_reftrr = sys_orig.REF_var.flag_reftrr;

    sys_copy->REF_var.N_fnm = sys_orig.REF_var.N_fnm;

    sys_copy->REF_var.reftrr_fnm = emalloc(sys_orig.REF_var.N_fnm * sizeof(tW_word));

    for (i = 0; i < sys_orig.REF_var.N_fnm; i++) {
	strcpy(sys_copy->REF_var.reftrr_fnm[i],
	       sys_orig.REF_var.reftrr_fnm[i]);
    }

}

/*****************************************************************************************
copy_TRIM(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_TRIM(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->TRIM_var.flag_TRIM = sys_orig.TRIM_var.flag_TRIM;

    sys_copy->TRIM_var.FE = sys_orig.TRIM_var.FE;

}

/*****************************************************************************************
copy_CHISQD(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_CHISQD(tW_system sys_orig, tW_system * sys_copy)
{
    sys_copy->CHISQD_var.flag_CHISQD = sys_orig.CHISQD_var.flag_CHISQD;

    strcpy(sys_copy->CHISQD_var.force_fnm, sys_orig.CHISQD_var.force_fnm);

}

/*****************************************************************************************
copy_REG(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_REG(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->REG_var.flag_REG = sys_orig.REG_var.flag_REG;

    strcpy(sys_copy->REG_var.type, sys_orig.REG_var.type);

    sys_copy->REG_var.Nmax = sys_orig.REG_var.Nmax;

    sys_copy->REG_var.tau_alpha = sys_orig.REG_var.tau_alpha;

    sys_copy->REG_var.tau_beta = sys_orig.REG_var.tau_beta;

    sys_copy->REG_var.Nframes = sys_orig.REG_var.Nframes;

}

/*****************************************************************************************
copy_RESCALE(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_RESCALE(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->RESCALE_var.flag_RESCALE = sys_orig.RESCALE_var.flag_RESCALE;

    sys_copy->RESCALE_var.Nmax = sys_orig.RESCALE_var.Nmax;

    sys_copy->RESCALE_var.tau_phi = sys_orig.RESCALE_var.tau_phi;

}

/*****************************************************************************************
copy_CONSTRAIN(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_CONSTRAIN(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->CONSTRAIN_var.flag_CONSTRAIN =
	sys_orig.CONSTRAIN_var.flag_CONSTRAIN;

    sys_copy->CONSTRAIN_var.Nmax = sys_orig.CONSTRAIN_var.Nmax;

    sys_copy->CONSTRAIN_var.tau_dih = sys_orig.CONSTRAIN_var.tau_dih;

    sys_copy->CONSTRAIN_var.lambda = sys_orig.CONSTRAIN_var.lambda;

    sys_copy->CONSTRAIN_var.dlambda = sys_orig.CONSTRAIN_var.dlambda;

}

/*****************************************************************************************
copy_ITER(): Copies the information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_ITER(tW_system sys_orig, tW_system * sys_copy)
{

    sys_copy->ITER_var.flag_ITER = sys_orig.ITER_var.flag_ITER;

    sys_copy->ITER_var.flag_AAM2 = sys_orig.ITER_var.flag_AAM2;

    strcpy(sys_copy->ITER_var.AAM2_fnm, sys_orig.ITER_var.AAM2_fnm);

    sys_copy->ITER_var.flag_bsolnerr = sys_orig.ITER_var.flag_bsolnerr;
}

/*****************************************************************************************
copy_Bonded_Inter_Types(): Copies information from sys_orig to sys_copy for all
intramolecular interactions.
*****************************************************************************************/
void copy_Bonded_Inter_Types(tW_system sys_orig, tW_system * sys_copy)
{
    int i, j;
    int N_Bond_Types = sys_orig.N_Bond_Int_Types;
    int i_0;
    tW_Bonded_Inter *bond_ptr;

    sys_copy->N_Bond_Int_Types = N_Bond_Types;
    sys_copy->Bonded_Inter_Types = (tW_Bonded_Inter *) emalloc(N_Bond_Types * sizeof(tW_Bonded_Inter));

    /* Copy information for each interaction. */
    for (i = 0; i < N_Bond_Types; i++) {
	bond_ptr = &(sys_copy->Bonded_Inter_Types[i]);

	strcpy(bond_ptr->inter_name,
	       sys_orig.Bonded_Inter_Types[i].inter_name);
	strcpy(bond_ptr->name, sys_orig.Bonded_Inter_Types[i].name);

	bond_ptr->N_Int_Sites = sys_orig.Bonded_Inter_Types[i].N_Int_Sites;
	bond_ptr->Site_Types =
	    (tW_word *) ecalloc(bond_ptr->N_Int_Sites, sizeof(tW_word));
	for (j = 0; j < bond_ptr->N_Int_Sites; j++) {
	    strcpy(bond_ptr->Site_Types[j],
		   sys_orig.Bonded_Inter_Types[i].Site_Types[j]);
	}

	bond_ptr->i_basis = sys_orig.Bonded_Inter_Types[i].i_basis;
	bond_ptr->N_coeff = sys_orig.Bonded_Inter_Types[i].N_coeff;
	i_0 = sys_orig.Bonded_Inter_Types[i].i_0;
	bond_ptr->i_0 = i_0;

	bond_ptr->ptr_b = &(sys_copy->b[i_0]);
	bond_ptr->ptr_b_ref = &(sys_copy->b_ref[i_0]);
	bond_ptr->ptr_g = &(sys_copy->g[i_0]);
	bond_ptr->ptr_g_cnt = &(sys_copy->g_cnt[i_0]);
	bond_ptr->ptr_L = &(sys_copy->L[i_0]);
	bond_ptr->dr = sys_orig.Bonded_Inter_Types[i].dr;
	bond_ptr->R_0 = sys_orig.Bonded_Inter_Types[i].R_0;
	bond_ptr->R_max = sys_orig.Bonded_Inter_Types[i].R_max;
	bond_ptr->n_smooth = sys_orig.Bonded_Inter_Types[i].n_smooth;
	bond_ptr->N_pts = sys_orig.Bonded_Inter_Types[i].N_pts;
	bond_ptr->kspline = sys_orig.Bonded_Inter_Types[i].kspline;	/* JFR - 07.23.12: Bspline */

	bond_ptr->N_instances = 0;
	bond_ptr->Inter_List = NULL;

	bond_ptr->N_powers = sys_orig.Bonded_Inter_Types[i].N_powers;
	bond_ptr->powers = sys_orig.Bonded_Inter_Types[i].powers;

	bond_ptr->n_bonds_intramolec_pair_inter = sys_orig.Bonded_Inter_Types[i].n_bonds_intramolec_pair_inter;
    }
}


/*****************************************************************************************
copy_site_types(): Copies information from sys_orig to sys_copy.
*****************************************************************************************/
void copy_site_types(tW_system sys_orig, tW_system * sys_copy)
{
    int i, j;
    int N_Site_Types = sys_orig.N_Site_Types;

    sys_copy->N_Site_Types = N_Site_Types;

    sys_copy->Site_Types = (tW_word *) emalloc(N_Site_Types * sizeof(tW_word));

    sys_copy->Inter_Map = ecalloc(N_Site_Types, sizeof(tW_word *));
    sys_copy->Inter_Map_Len = ecalloc(N_Site_Types, sizeof(int));
    sys_copy->Inter_iMap = ecalloc(N_Site_Types, sizeof(int *));

    for (i = 0; i < N_Site_Types; i++) 
    {
	strcpy(sys_copy->Site_Types[i], sys_orig.Site_Types[i]);

	sys_copy->Inter_Map_Len[i] = sys_orig.Inter_Map_Len[i];

	sys_copy->Inter_iMap[i] = ecalloc(sys_orig.Inter_Map_Len[i], sizeof(int));
	sys_copy->Inter_Map[i] = ecalloc(sys_orig.Inter_Map_Len[i], sizeof(tW_word));

	for (j=0; j<sys_orig.Inter_Map_Len[i]; j++)
	{
		sys_copy->Inter_iMap[i][j] = sys_orig.Inter_iMap[i][j];
		strcpy(sys_copy->Inter_Map[i][j], sys_orig.Inter_Map[i][j]);
	}

    }
}


/*****************************************************************************************
copy_arrays(): Allocates memory for arrays in sys_copy and makes sure they are zeroed.
*****************************************************************************************/
void copy_arrays(tW_system sys_orig, tW_system * sys_copy)
{
    int N_coeff = sys_orig.N_coeff;
    int N_pack = sys_orig.N_pack;

    sys_copy->N_coeff = N_coeff;
    sys_copy->N_pack = N_pack;

    // Try free'ing the arrays here before reallocating them
    if (sys_copy->b != NULL)
    {
	free(sys_copy->b);
    }
    if (sys_copy->b_ref != NULL)
    {
	free(sys_copy->b_ref);
    }
    if (sys_copy->g != NULL)
    {
	free(sys_copy->g);
    }
    if (sys_copy->g_cnt != NULL)
    {
	free(sys_copy->g_cnt);
    }
    if (sys_copy->L != NULL)
    {
	free(sys_copy->L);
    }
    if (sys_copy->M != NULL)
    {
	free(sys_copy->M);
    }
    if (sys_copy->M2 != NULL)
    {
	free(sys_copy->M2);
    }
    if (sys_copy->M_cnt != NULL)
    {
	free(sys_copy->M_cnt);
    }
    if (sys_copy->d2b != NULL)
    {
	free(sys_copy->d2b);
    }
    if (sys_copy->d2M != NULL)
    {
	free(sys_copy->d2M);
    }
    if (sys_copy->rescale != NULL)
    {
	free(sys_copy->rescale);
    }


    /* These are needed for the calc_grids(). */
    sys_copy->b = (double *) ecalloc(N_coeff, sizeof(double));
    sys_copy->b_ref = (double *) ecalloc(N_coeff, sizeof(double));
    sys_copy->g = (double *) ecalloc(N_coeff, sizeof(double));
    sys_copy->g_cnt = (double *) ecalloc(N_coeff, sizeof(double));
    sys_copy->L = (double *) ecalloc(N_coeff, sizeof(double));
    // JFR - added 04.11.12: put the matrix in packed form
    sys_copy->M = (double *) ecalloc(N_pack, sizeof(double));
    sys_copy->M2 = (double *) ecalloc(N_pack, sizeof(double));
    if (sys_orig.MT_var.flag_Mcnt == TRUE) {
	sys_copy->M_cnt = (double *) ecalloc(N_pack, sizeof(double));
    } else {
	sys_copy->M_cnt = NULL;
    }

    if (sys_orig.FRAMEWEIGHT_var.flag_FRAMEWEIGHT == TRUE ) {
	    if (sys_copy->b_wt != NULL)
	    {
		free(sys_copy->b_wt);
	    }
	    if (sys_copy->M_wt != NULL)
	    {
		free(sys_copy->M_wt);
	    }
	    if (sys_copy->M2_wt != NULL)
	    {
		free(sys_copy->M2_wt);
	    }

	    sys_copy->b_wt = (double *) ecalloc(N_coeff, sizeof(double));		
	    sys_copy->M_wt = (double *) ecalloc(N_pack, sizeof(double));
	    sys_copy->M2_wt = (double *) ecalloc(N_pack, sizeof(double));
    }

    //if ( strcmp( sys_orig.PC_var.LPC, "bvar" ) == 0 )
    //{
    sys_copy->d2b = (double *) ecalloc(N_coeff, sizeof(double));	/* JFR - 01.31.13 */
    //}
    //else { sys_copy->d2b = NULL; }
    //if ( strcmp( sys_orig.PC_var.RPC, "MTvar" ) == 0 )
    //{
    sys_copy->d2M = (double *) ecalloc(N_pack, sizeof(double));	/* JFR - 01.31.13 */
    //}
    //else { sys_copy->d2M = NULL; }

    sys_copy->rescale = (double *) ecalloc(N_coeff, sizeof(double));




    if (sys_copy->b_forces != NULL)
    {
	free(sys_copy->b_forces);
    }
    if (sys_copy->b_struct != NULL)
    {
	free(sys_copy->b_struct);
    }
    if (sys_copy->phi != NULL)
    {
	free(sys_copy->phi);
    }
    if (sys_copy->phi_forces != NULL)
    {
	free(sys_copy->phi_forces);
    }
    if (sys_copy->phi_struct != NULL)
    {
	free(sys_copy->phi_struct);
    }


    /* These are not needed until the end. */
    sys_copy->x = NULL;
    sys_copy->b_forces = NULL;
    sys_copy->b_struct = NULL;
    sys_copy->phi = NULL;
    sys_copy->phi_forces = NULL;
    sys_copy->phi_struct = NULL;

    /* To be sure, I should probably zero out these arrays by calling clear arrays. */
    clear_sys_arrays(sys_copy);

}


/*****************************************************************************************
copy_Inter_Types(): Copies information from sys_orig to sys_copy for all interaction 
types.
*****************************************************************************************/
void copy_Inter_Types(tW_system sys_orig, tW_system * sys_copy)
{
    int i;
    int i_0;
    int N_Inter_Types = sys_orig.N_Inter_Types;
    tW_Inter_Types *inter_ptr;

    sys_copy->N_Inter_Types = N_Inter_Types;

    sys_copy->Inter_Types = (tW_Inter_Types *) emalloc(N_Inter_Types * sizeof(tW_Inter_Types));

    /* Copy information for each interaction. */
    for (i = 0; i < N_Inter_Types; i++) {
	inter_ptr = &(sys_copy->Inter_Types[i]);
	strcpy(inter_ptr->inter_name, sys_orig.Inter_Types[i].inter_name);
	strcpy(inter_ptr->inter_type, sys_orig.Inter_Types[i].inter_type);
	strcpy(inter_ptr->basis, sys_orig.Inter_Types[i].basis);
	inter_ptr->i_basis = sys_orig.Inter_Types[i].i_basis;
	inter_ptr->dr = sys_orig.Inter_Types[i].dr;
	inter_ptr->R_min = sys_orig.Inter_Types[i].R_min;
	inter_ptr->R_max = sys_orig.Inter_Types[i].R_max;
	inter_ptr->N_pts = sys_orig.Inter_Types[i].N_pts;
	inter_ptr->N_coeff = sys_orig.Inter_Types[i].N_coeff;
	inter_ptr->n_smooth = sys_orig.Inter_Types[i].n_smooth;
	inter_ptr->i_0 = sys_orig.Inter_Types[i].i_0;
	i_0 = sys_orig.Inter_Types[i].i_0;
	inter_ptr->ptr_b = &(sys_copy->b[i_0]);
	inter_ptr->ptr_b_ref = &(sys_copy->b_ref[i_0]);
	inter_ptr->ptr_g = &(sys_copy->g[i_0]);
	inter_ptr->ptr_g_cnt = &(sys_copy->g_cnt[i_0]);
	inter_ptr->ptr_L = &(sys_copy->L[i_0]);

	inter_ptr->N_powers = sys_orig.Inter_Types[i].N_powers;
	inter_ptr->powers = sys_orig.Inter_Types[i].powers;

	inter_ptr->kspline = sys_orig.Inter_Types[i].kspline;	/* JFR - 07.22.12: copy Bspline order */
    }
}


/*****************************************************************************************
copy_Inter2_Types(): Copies information from sys_orig to sys_copy for all pair 
intermolecular interactions.
*****************************************************************************************/
void copy_Inter2_Types(tW_system sys_orig, tW_system * sys_copy)
{
    int i;
    int N_Inter2_Types = sys_orig.N_Inter2_Types;
    int i_0;
    tW_type_inter2 *inter_ptr;

    sys_copy->N_Inter2_Types = N_Inter2_Types;

    sys_copy->Inter2_Type_List = (tW_type_inter2 *) emalloc(N_Inter2_Types * sizeof(tW_type_inter2));

    /* Copy information for each interaction. */
    for (i = 0; i < N_Inter2_Types; i++) {
	inter_ptr = &(sys_copy->Inter2_Type_List[i]);

	strcpy(inter_ptr->inter_name, sys_orig.Inter2_Type_List[i].inter_name);
	strcpy(inter_ptr->name1, sys_orig.Inter2_Type_List[i].name1);
	strcpy(inter_ptr->name2, sys_orig.Inter2_Type_List[i].name2);

	inter_ptr->N_inter = 0;
	inter_ptr->N_coeff = sys_orig.Inter2_Type_List[i].N_coeff;
	inter_ptr->i_basis = sys_orig.Inter2_Type_List[i].i_basis;
	inter_ptr->i_0 = sys_orig.Inter2_Type_List[i].i_0;
	inter_ptr->dr = sys_orig.Inter2_Type_List[i].dr;
	inter_ptr->R_0 = sys_orig.Inter2_Type_List[i].R_0;
	inter_ptr->R_max = sys_orig.Inter2_Type_List[i].R_max;
	inter_ptr->n_smooth = sys_orig.Inter2_Type_List[i].n_smooth;
	inter_ptr->N_pts = sys_orig.Inter2_Type_List[i].N_pts;

	i_0 = sys_orig.Inter2_Type_List[i].i_0;
	inter_ptr->ptr_b = &(sys_copy->b[i_0]);
	inter_ptr->ptr_b_ref = &(sys_copy->b_ref[i_0]);
	inter_ptr->ptr_g2 = &(sys_copy->g[i_0]);
	inter_ptr->ptr_g2_cnt = &(sys_copy->g_cnt[i_0]);
	inter_ptr->ptr_L = &(sys_copy->L[i_0]);

	inter_ptr->N_powers = sys_orig.Inter2_Type_List[i].N_powers;
	inter_ptr->powers = sys_orig.Inter2_Type_List[i].powers;

	inter_ptr->kspline = sys_orig.Inter2_Type_List[i].kspline;	/* JFR - 07.22.12: Bspline order */
    }
}


/*****************************************************************************************
clear_sys_arrays(): Zeros the arrays that are not pointing to NULL.
*****************************************************************************************/
void clear_sys_arrays(tW_system * sys)
{
    int i;

    for (i = 0; i < sys->N_coeff; i++) {
	/* These should never point to NULL. */
	sys->b[i] = 0.0;
	sys->b_ref[i] = 0.0;
	sys->g[i] = 0.0;
	sys->g_cnt[i] = 0.0;
	sys->L[i] = 0.0;

	/* These are initially set to NULL. */
	if (sys->x != NULL) {
	    sys->x[i] = 0.0;
	}
	if (sys->b_forces != NULL) {
	    sys->b_forces[i] = 0.0;
	}
	if (sys->b_struct != NULL) {
	    sys->b_struct[i] = 0.0;
	}
	if (sys->phi != NULL) {
	    sys->phi[i] = 0.0;
	}
	if (sys->phi_forces != NULL) {
	    sys->phi_forces[i] = 0.0;
	}
	if (sys->phi_struct != NULL) {
	    sys->phi_struct[i] = 0.0;
	}

	if (sys->d2b != NULL) {
	    sys->d2b[i] = 0.0;
	}
	sys->rescale[i] = 0.0;
    }

    // JFR - added 04.11.12: put the matrix in packed form
    for (i = 0; i < sys->N_pack; i++) {
	sys->M[i] = 0.0;
	sys->M2[i] = 0.0;
	if (sys->M_cnt != NULL) {
	    sys->M_cnt[i] = 0.0;
	}
	if (sys->d2M != NULL) {
	    sys->d2M[i] = 0.0;
	}
    }


}


/*****************************************************************************************
free_sys_copy(): Frees memory allocated for the data structure sys_top.
*****************************************************************************************/
int free_sys_copy(tW_system * sys_top)
{
    int i, j;

    if (sys_top->Site_Types != NULL) {
	free(sys_top->Site_Types);
    }

    if (sys_top->x != NULL) {
	free(sys_top->x);
    }
    if (sys_top->b != NULL) {
	free(sys_top->b);
    }
    if (sys_top->b_ref != NULL) {
	free(sys_top->b_ref);
    }
    if (sys_top->b_forces != NULL) {
	free(sys_top->b_forces);
    }
    if (sys_top->b_struct != NULL) {
	free(sys_top->b_struct);
    }
    if (sys_top->g != NULL) {
	free(sys_top->g);
    }
    if (sys_top->g_cnt != NULL) {
	free(sys_top->g_cnt);
    }
    if (sys_top->L != NULL) {
	free(sys_top->L);
    }
    if (sys_top->phi != NULL) {
	free(sys_top->phi);
    }
    if (sys_top->phi_forces != NULL) {
	free(sys_top->phi_forces);
    }
    if (sys_top->phi_struct != NULL) {
	free(sys_top->phi_struct);
    }
    // JFR - added 04.11.12: put the matrix in packed form
    if (sys_top->M != NULL) {
	free(sys_top->M);
    }
    if (sys_top->M2 != NULL) {
	free(sys_top->M2);
    }

    if (sys_top->b_wt != NULL){
	free(sys_top->b_wt);
    }

    if (sys_top->M_wt != NULL){
	free(sys_top->M_wt);
    }

    if (sys_top->M2_wt != NULL){
	free(sys_top->M2_wt);
    }

    if (sys_top->M_cnt != NULL) {
	free(sys_top->M_cnt);
    }

    if (sys_top->d2b != NULL) {
	free(sys_top->d2b);
    }				/* JFR - 01.31.13 */
    if (sys_top->d2M != NULL) {
	free(sys_top->d2M);
    }
    /* JFR - 01.31.13 */
    if (sys_top->rescale != NULL) {
	free(sys_top->rescale);
    }
    /* JFR - 01.31.13 */
    if (sys_top->TPR_var.TPR_files != NULL) {
	free(sys_top->TPR_var.TPR_files);
    }				//JFR - added 04.06.12

    if (sys_top->Inter_Map_Len != NULL) {
	free(sys_top->Inter_Map_Len);
    }

    if (sys_top->Inter_Map != NULL) {
        for (i=0; i<sys_top->N_Site_Types; i++) {
            if (sys_top->Inter_Map[i] != NULL) {
		free(sys_top->Inter_Map[i]);
	    }
	}
	free(sys_top->Inter_Map);
    }

    if (sys_top->Inter_iMap != NULL) {
        for (i=0; i<sys_top->N_Site_Types; i++) {
            if (sys_top->Inter_iMap[i] != NULL) {
		free(sys_top->Inter_iMap[i]);
	    }
	}
	free(sys_top->Inter_iMap);
    }

    if (sys_top->Inter_Types != NULL) {
	free(sys_top->Inter_Types);
    }

    if (sys_top->Inter2_Type_List != NULL) {
	free(sys_top->Inter2_Type_List);
    }

    if (sys_top->Bonded_Inter_Types != NULL) {
	for (i = 0; i < sys_top->N_Bond_Int_Types; i++) {
	    free(sys_top->Bonded_Inter_Types[i].Site_Types);

	    if (sys_top->Bonded_Inter_Types[i].Inter_List != NULL) {
		for (j = 0; j < sys_top->Bonded_Inter_Types[i].N_instances; j++) {
		    free(sys_top->Bonded_Inter_Types[i].Inter_List[j]);

		}
		free(sys_top->Bonded_Inter_Types[i].Inter_List);
	    }
	}
	free(sys_top->Bonded_Inter_Types);
    }

    return 0;
}


/*****************************************************************************************
initialize_sys(): Initializes the variables in the sys data structure. This ensures that
variables have known values for later use, particularly, the pointers are set to NULL.
*****************************************************************************************/
void initialize_sys(tW_system * sys)
{
    sys->N_Site_Types = 0;
    sys->Site_Types = NULL;
    sys->Inter_Map = NULL;
    sys->Inter_Map_Len = NULL;
    sys->Inter_iMap = NULL;
    sys->N_Inter_Types = 0;
    sys->Inter_Types = NULL;
    sys->N_Inter2_Types = 0;
    sys->Inter2_Type_List = NULL;
    sys->N_Bond_Int_Types = 0;
    sys->Bonded_Inter_Types = NULL;
    sys->N_coeff = 0;
    sys->N_pack = 0;		// JFR - added 04.11.12: put the matrix in packed form
    sys->Temperature = 0.0;
    sys->nrexcl = -1;
    sys->x = NULL;
    sys->b = NULL;
    sys->b_ref = NULL;
    sys->b_forces = NULL;
    sys->b_struct = NULL;
    sys->g = NULL;
    sys->g_cnt = NULL;
    sys->L = NULL;
    sys->phi = NULL;
    sys->phi_forces = NULL;
    sys->phi_struct = NULL;
    sys->M = NULL;
    sys->M2 = NULL;
    sys->M_cnt = NULL;
    sys->d2b = NULL;
    sys->d2M = NULL;
    sys->rescale = NULL;
    sys->b_wt = NULL;
    sys->M_wt = NULL;
    sys->M2_wt = NULL;

    /* MRD 03.04.2019 */
    sys->SKIP_TRIPLE_LOOP = FALSE;
    sys->half_matrix = NULL;
    sys->bm_half_mat = NULL;
    sys->M_M2_proc = FALSE;

    sys->flag_ref_potential = FALSE;

    /*JFR - 08.10.11 - PT variables */
    sys->PT_var.flag_PT = FALSE;
    sys->PT_var.N_PT = 0;
    sys->PT_var.dPT = 1;
    sys->PT_var.flag_MMOTF_SEP = FALSE;
    sys->PT_var.flag_eigen = FALSE;
    sys->PT_var.N_eigen = 0;

    /*JFR - 08.10.11 - Eigen variables */
    sys->Eigen_var.flag_Eigen = FALSE;
    sys->Eigen_var.N_Eigen = 0;
    sys->Eigen_var.flag_printn = FALSE;
    sys->Eigen_var.DM = 1000;
    sys->Eigen_var.DL = 1000;
    sys->Eigen_var.flag_Gbar = FALSE;
    sys->Eigen_var.flag_norm = FALSE;

    /*JFR - 04.06.12 - SVD variables */
    sys->SVD_var.flag_SVD = FALSE;
    sys->SVD_var.rcond = FLOAT_EPS;
    sys->SVD_var.flag_printSV = FALSE;
    sys->SVD_var.flag_printevecs = FALSE;
    sys->SVD_var.flag_solve = TRUE;

    /*JFR - 04.06.12 - TPR variables */
    sys->TPR_var.flag_TPR = FALSE;
    sys->TPR_var.N_TPR = 0;
    sys->TPR_var.TPR_files = NULL;

    /*JFR - 06.27.12 - TPR variables */
    sys->TPR_EXCL_var.flag_TPR_excl = FALSE;
    strcpy(sys->TPR_EXCL_var.TPR_excl, " ");

    /*JFR - 04.06.12 - MT variables */
    sys->MT_var.flag_MT = FALSE;
    sys->MT_var.flag_print = FALSE;
    sys->MT_var.flag_norm = FALSE;
    sys->MT_var.flag_Mcnt = FALSE;

    /*JFR - 04.06.12 - MFD variables */
    sys->MFD_var.flag_MFD = FALSE;

    /*JFR - 04.06.12 - CalcMODE variables */
    sys->CalcMODE_var.flag_CalcMODE = FALSE;
    sys->CalcMODE_var.CalcMODE = IFULL;

    /*JFR - 04.06.12 - PC variables */
    sys->PC_var.flag_PC = FALSE;
    strcpy(sys->PC_var.RPC, "colnorm");
    strcpy(sys->PC_var.LPC, "rowmax");
    sys->PC_var.flag_normb = FALSE;
    sys->PC_var.flag_normphi = FALSE;

    /*JFR - 04.13.12 - MEM variables */
    sys->MEM_var.flag_MEM = FALSE;
    sys->MEM_var.flag_LOWMEM = FALSE;	/* Run in normal memory mode unless the user says otherwise */
    strcpy(sys->MEM_var.info, "structures");	/* This default will never apply */
    sys->MEM_var.flag_mult_top = FALSE;	/* DON'T solve each topology separately, unless the user says to */

    /*JFR - 04.16.12 - SOLN variables */
    sys->SOLN_var.flag_SOLN = FALSE;
    strcpy(sys->SOLN_var.SOLN_METH, "SVD");	/* Use SVD by default to avoid numerical issues */

    /*NJD - 03.10.15 - FRAMEWEIGHT variables */
    sys->FRAMEWEIGHT_var.flag_FRAMEWEIGHT = FALSE;
    strcpy(sys->FRAMEWEIGHT_var.FRAMEWEIGHT, "NONE"); /* Default to no weighting as this is most commonly used*/
    sys->wt_norm = 0.00;


    /*JFR - 04.16.12 - ERR variables */
    sys->ERR_var.flag_ERR = FALSE;
    sys->ERR_var.FACT = 'N';

    /* JFR - 06.27.12: Chi2 */
    sys->Chi2 = 0.00;

    /* JFR - 07.16.12 - REF variables */
    sys->REF_var.flag_REF = FALSE;
    sys->REF_var.flag_calcbref = FALSE;
    sys->REF_var.flag_readbref = FALSE;
    sys->REF_var.flag_splitfiles = FALSE;
    sys->REF_var.flag_reftrr = FALSE;
    sys->REF_var.N_fnm = 0;
    sys->REF_var.reftrr_fnm = NULL;

    /*JFR - 01.29.13 - TRIM variables */
    sys->TRIM_var.flag_TRIM = FALSE;
    sys->TRIM_var.FE = 0.001;

    /* JFR - 01.29.13 - CHISQD variables */
    sys->CHISQD_var.flag_CHISQD = FALSE;
    strcpy(sys->CHISQD_var.force_fnm, "");

    /* JFR - 12.03.13 - REG variables */
    sys->REG_var.flag_REG = FALSE;
    strcpy(sys->REG_var.type, "");
    sys->REG_var.Nmax = 100;
    sys->REG_var.tau_alpha = 0.01;
    sys->REG_var.tau_beta = 0.01;
    sys->REG_var.Nframes = 0;

    /* JFR - 01.31.13 - RESCALE variables */
    sys->RESCALE_var.flag_RESCALE = FALSE;
    sys->RESCALE_var.Nmax = 100;
    sys->RESCALE_var.tau_phi = 0.01;

    /* JFR - 01.31.13 - CONSTRAIN variables */
    sys->CONSTRAIN_var.flag_CONSTRAIN = FALSE;
    sys->CONSTRAIN_var.Nmax = 100;
    sys->CONSTRAIN_var.tau_dih = 0.01;
    sys->CONSTRAIN_var.lambda = 0.01;
    sys->CONSTRAIN_var.dlambda = 0.01;

    /* JFR - 12.03.13 - ITER variables */
    sys->ITER_var.flag_ITER = FALSE;
    sys->ITER_var.flag_AAM2 = FALSE;
    strcpy(sys->ITER_var.AAM2_fnm, " ");
    sys->ITER_var.flag_bsolnerr = FALSE;

}


/*****************************************************************************************
initialize_files(): Initializes the variables in the files data structure. 
*****************************************************************************************/
void initialize_files(tW_files * files)
{
    files->N_struct_file = -1;
    files->struct_file = NULL;
    files->N_struct = -1;
    files->structures = NULL;
    files->p_struct = NULL;
    strcpy(files->mode, "VOID");
}


/*****************************************************************************************
get_bond_inter_site_types(): Reads in the site types for intramolecular interactions. 
If two sites are present, then it also looks for the number of bonds separating the
two interaction sites. This is only relevant for the PDB loop for intramolecular
interactions.
*****************************************************************************************/
void get_bond_inter_site_types(tW_line inp_line,
			       tW_Bonded_Inter * Bond_ptr, int N_Int_Sites,
			       tW_word inter_name)
{
    int test_sscanf;

    if (N_Int_Sites == 2) {
	test_sscanf = sscanf(inp_line, "%s %s %s",
			     inter_name, Bond_ptr->Site_Types[0],
			     Bond_ptr->Site_Types[1]);
	if (test_sscanf != 3) {
	    printf("\nERROR: Problem reading line %s.\n", inp_line);
	    exit(EXIT_FAILURE);
	}

	test_sscanf =
	    sscanf(inp_line, "%*s %*s %*s %d",
		   &Bond_ptr->n_bonds_intramolec_pair_inter);
	if (test_sscanf == EOF) {
	    Bond_ptr->n_bonds_intramolec_pair_inter = NOT_SET;
	}
    }

    if (N_Int_Sites == 3) {
	test_sscanf = sscanf(inp_line, "%s %s %s %s",
			     inter_name, Bond_ptr->Site_Types[0],
			     Bond_ptr->Site_Types[1],
			     Bond_ptr->Site_Types[2]);
	if (test_sscanf != 4) {
	    printf("\nERROR: Problem reading line %s.\n", inp_line);
	    exit(EXIT_FAILURE);
	}
    }

    if (N_Int_Sites == 4) {
	test_sscanf = sscanf(inp_line, "%s %s %s %s %s",
			     inter_name, Bond_ptr->Site_Types[0],
			     Bond_ptr->Site_Types[1],
			     Bond_ptr->Site_Types[2],
			     Bond_ptr->Site_Types[3]);
	if (test_sscanf != 5) {
	    printf("\nERROR: Problem reading line %s.\n", inp_line);
	    exit(EXIT_FAILURE);
	}
    }

}


/*****************************************************************************************
get_power_basis_param(): Stores the parameters for a 'powers' basis (e.g. for a LJ 
potential, the powers listed would be 6, 12). Returns index to this basis. The powers
must be separated by commas and listed in parentheses. 
*****************************************************************************************/
int get_power_basis_param(tW_line inp_line, tW_Inter_Types * inter)
{
    int i, j, n;
    int terms_cnt = 0;
    int test_sscanf;
    char *string_ptr;
    tW_line string;
    tW_word inter_name;

    char *saveptr;
    char *tok;

    test_sscanf = sscanf(inp_line, " %s %*s power %lf %d",
			 inter_name, &(inter->R_max), &(inter->N_powers));

    if (test_sscanf != 3) {
	printf("\nERROR: Did not read %s correctly.\n  Check par.txt.\n",
	       inter_name);
	exit(EXIT_FAILURE);
    }

    inter->powers = ecalloc(inter->N_powers, sizeof(int));

    inter->N_coeff = inter->N_powers;

    strcpy(string, inp_line);
    string_ptr = strstr(string, "(");
    strcpy(string, string_ptr);
    n = strspn(string, "(");
    string_ptr = &string[n];

    tok = strtok_r(string_ptr, ",", &saveptr);
    while (tok != NULL) {
	terms_cnt++;
	if (terms_cnt > inter->N_powers) {
	    fprintf(stderr,
		    "ERROR: Found more powers than indicated for the power basis\n");
	    fprintf(stderr,
		    "Check par.txt (are there any extra commas?)\n");
	    exit(EXIT_FAILURE);
	}

	inter->powers[terms_cnt - 1] = atoi(tok);
	tok = strtok_r(NULL, ",", &saveptr);
    }
    if (terms_cnt < inter->N_powers) {
	fprintf(stderr,
		"ERROR: Found less powers than indicated for the power basis\n");
	fprintf(stderr, "Check par.txt (are there any missing commas?)\n");
	exit(EXIT_FAILURE);
    }

    /* Make sure there are no zero powers */
    for (i = 0; i < inter->N_powers; i++) {
	if (inter->powers[i] == 0) {
	    printf("ERROR in get_powers(): powers[%d] = %d.\n", i,
		   inter->powers[i]);
	    printf("Check par.txt.\n");
	    exit(EXIT_FAILURE);
	}
    }

    /* Maker sure no power is repeated. */
    for (i = 0; i < inter->N_powers - 1; i++) {
	for (j = i + 1; j < inter->N_powers; j++) {
	    if (inter->powers[i] == inter->powers[j]) {
		printf
		    ("ERROR in get_powers(): %d repeated more than once.\n",
		     inter->powers[i]);
		printf("Check par.txt\n");
		exit(EXIT_FAILURE);
	    }
	}
    }

    /* Make sure the powers are listed in increasing order. This is not necessary. */
    for (i = 0; i < inter->N_powers - 1; i++) {
	if (inter->N_powers <= 1) {
	    break;
	}

	if (inter->powers[i] > inter->powers[i + 1]) {
	    printf
		("ERROR: Please list the powers for %s in increasing order.\n",
		 inter->inter_name);
	    printf("    Thank You.\n");
	    exit(EXIT_FAILURE);
	}
    }


    return POWER_INDEX;
}


/*****************************************************************************************
get_delta_basis_param(): Reads in the parameters for a interaction parameterized using
a disrete delta basis. Returns index for this basis.
*****************************************************************************************/
int get_delta_basis_param(tW_line inp_line, double *dr, double *R_0,
			  double *R_max, int *N_pts, int *N_coeff,
			  int *n_smooth, tW_word inter_type)
{
    int test_sscanf;
    tW_word inter_name;

    test_sscanf = sscanf(inp_line, " %s %*s delta %lf %lf %lf %d",
			 inter_name, dr, R_0, R_max, n_smooth);

    if (test_sscanf != 5) {
	printf("\nERROR: Did not read %s correctly.\n  Check par.txt.\n",
	       inter_name);
	exit(EXIT_FAILURE);
    }

    /* Perform calculations in radians. Results will be converted to degrees. */
    if ((strcmp(inter_type, B_ANGLE) == 0)
	|| (strcmp(inter_type, B_DIHEDRAL) == 0)) {
	*dr *= M_PI / 180.0;
	*R_0 *= M_PI / 180.0;
	*R_max *= M_PI / 180.0;
    }

    *N_coeff = 1 + rint((*R_max - *R_0 + FLOAT_EPS) / *dr);
    *N_pts = *N_coeff;
    *R_max = *R_0 + (*N_coeff - 1) * (*dr);

    /* For dihedral types, make sure that R_0 and R_max are -180 and 180. */
//  if ( strcmp( inter_type, B_DIHEDRAL ) == 0 )
//  {
//    if (
//         ( fabs( *R_0   + M_PI ) > FLOAT_EPS ) ||
//         ( fabs( *R_max - M_PI ) > FLOAT_EPS )
//       )
//    {
//      printf( "\nERROR: Check interaction '%s'.\n", inter_name );
//      printf( "    For dihedral interactions, the end-points must be -180 and 180.\n" );
//      printf( "    R_0: %f   R_max: %f\n", *R_0 * 180.0 / M_PI, *R_max * 180.0 / M_PI );
//      exit( EXIT_FAILURE );
//    }
// 
//  }

    return DELTA_BASIS_INDEX;
}

/*START JFR*/
/*****************************************************************************************
get_linear_basis_param(): Reads in the parameters for a interaction parameterized using
a linear spline basis. Returns index for this basis.  Identical to get_delta_basis_param()
except returns LINEAR_BASIS_INDEX.
*****************************************************************************************/
int get_linear_basis_param(tW_line inp_line, double *dr, double *R_0,
			   double *R_max, int *N_pts, int *N_coeff,
			   int *n_smooth, tW_word inter_type)
{
    int test_sscanf;
    tW_word inter_name;

    test_sscanf = sscanf(inp_line, " %s %*s linear %lf %lf %lf %d",
			 inter_name, dr, R_0, R_max, n_smooth);

    if (test_sscanf != 5) {
	printf("\nERROR: Did not read %s correctly.\n  Check par.txt.\n",
	       inter_name);
	exit(EXIT_FAILURE);
    }

    /* Perform calculations in radians. Results will be converted to degrees. */
    if ((strcmp(inter_type, B_ANGLE) == 0)
	|| (strcmp(inter_type, B_DIHEDRAL) == 0)) {
	*dr *= M_PI / 180.0;
	*R_0 *= M_PI / 180.0;
	*R_max *= M_PI / 180.0;
    }

    *N_coeff = 1 + rint((*R_max - *R_0 + FLOAT_EPS) / *dr);
    *N_pts = *N_coeff;
    *R_max = *R_0 + (*N_coeff - 1) * (*dr);

    /* For dihedral types, make sure that R_0 and R_max are -180 and 180. */
    if (strcmp(inter_type, B_DIHEDRAL) == 0) {
	if ((fabs(*R_0 + M_PI) > FLOAT_EPS) ||
	    (fabs(*R_max - M_PI) > FLOAT_EPS)
	    ) {
	    printf("\nERROR: Check interaction '%s'.\n", inter_name);
	    printf
		("    For dihedral interactions, the end-points must be -180 and 180.\n");
	    printf("    R_0: %f   R_max: %f\n", *R_0 * 180.0 / M_PI,
		   *R_max * 180.0 / M_PI);
	    exit(EXIT_FAILURE);
	}

    }

    return LINEAR_BASIS_INDEX;
}

/*END JFR*/

/*****************************************************************************************
get_Bspline_basis_param(): Reads in the parameters for a interaction parameterized using
a Bspline basis. Returns index for this basis.  Identical to get_delta_basis_param()
except returns BSPLINE_BASIS_INDEX and also reads in kspline = order of the Bspline.
JFR - 07.22.12
*****************************************************************************************/
int get_Bspline_basis_param(tW_line inp_line, double *dr, double *R_0,
			    double *R_max, int *N_pts, int *kspline,
			    int *N_coeff, int *n_smooth,
			    tW_word inter_type)
{
    int test_sscanf;
    tW_word inter_name;

    test_sscanf = sscanf(inp_line, " %s %*s Bspline %lf %lf %lf %d %d",
			 inter_name, dr, R_0, R_max, kspline, n_smooth);

    if (test_sscanf != 6) {
	printf("\nERROR: Did not read %s correctly.\n  Check par.txt.\n",
	       inter_name);
	exit(EXIT_FAILURE);
    }

    /* Perform calculations in radians. Results will be converted to degrees. */
    if ((strcmp(inter_type, B_ANGLE) == 0)
	|| (strcmp(inter_type, B_DIHEDRAL) == 0)) {
	*dr *= M_PI / 180.0;
	*R_0 *= M_PI / 180.0;
	*R_max *= M_PI / 180.0;
    }

    *N_coeff = 1 + rint((*R_max - *R_0 + FLOAT_EPS) / *dr);
    *N_pts = *N_coeff;
    *R_max = *R_0 + (*N_coeff - 1) * (*dr);

    /* MRD 02.05.2019 added extra padding for Bspline NB */
    if (strcmp(inter_type,NB_PAIR) == 0)
    {
	(*N_pts) += (*kspline);
        (*N_coeff) += (*kspline);
        (*R_max) += ((*dr) * (double)(*kspline));
        if ((*R_0) > 0.0)
        {
            (*R_0) -= ((*dr) * (double)(*kspline));
	    (*N_coeff) += (*kspline);
	    (*N_pts) += (*kspline);
        }
    }

    /* For dihedral types, make sure that R_0 and R_max are -180 and 180. */
    if (strcmp(inter_type, B_DIHEDRAL) == 0) {
	if ((fabs(*R_0 + M_PI) > FLOAT_EPS) ||
	    (fabs(*R_max - M_PI) > FLOAT_EPS)
	    ) {
	    printf("\nERROR: Check interaction '%s'.\n", inter_name);
	    printf
		("    For dihedral interactions, the end-points must be -180 and 180.\n");
	    printf("    R_0: %f   R_max: %f\n", *R_0 * 180.0 / M_PI,
		   *R_max * 180.0 / M_PI);
	    exit(EXIT_FAILURE);
	}

    }

    return BSPLINE_BASIS_INDEX;
}

/*****************************************************************************************
get_no_inter_types(): Returns the number of interactions listed after a directive in-
dicated by keyword. Stops execution if the number is not found. 
*****************************************************************************************/
int get_no_inter_types(tW_word keyword, tW_line inp_line)
{
    int N_Inter_Types;
    int test_sscanf;

    delete_directive_from_line(inp_line);

    test_sscanf = sscanf(inp_line, " %d ", &N_Inter_Types);

    if (test_sscanf != 1) {
	printf("\nERROR: Number of interactions for %s not found.\n",
	       keyword);
	exit(EXIT_FAILURE);
    }

    return N_Inter_Types;
}


/*****************************************************************************************
get_bond_basis_parameters(): Copies the parameters for a intramolecular interaction 
type for a specific group of site types. Stops execution if the interaction type is
not found or is of the incorrect type. 
*****************************************************************************************/
void get_bond_basis_parameters(tW_word inter_name,
			       tW_Bonded_Inter * Bond_ptr, int N_Int_sites,
			       tW_line inp_line, tW_system sys)
{
    int i;
    int flag = FALSE;

    for (i = 0; i < sys.N_Inter_Types; i++) {
	if (strcmp(inter_name, sys.Inter_Types[i].inter_name) == 0) {
	    strcpy(Bond_ptr->basis, sys.Inter_Types[i].basis);
	    Bond_ptr->N_instances = 0;
	    Bond_ptr->N_pts = sys.Inter_Types[i].N_pts;
	    Bond_ptr->N_coeff = sys.Inter_Types[i].N_coeff;
	    Bond_ptr->i_basis = sys.Inter_Types[i].i_basis;
	    Bond_ptr->i_0 = sys.Inter_Types[i].i_0;
	    Bond_ptr->dr = sys.Inter_Types[i].dr;
	    Bond_ptr->R_0 = sys.Inter_Types[i].R_min;
	    Bond_ptr->R_max = sys.Inter_Types[i].R_max;
	    Bond_ptr->n_smooth = sys.Inter_Types[i].n_smooth;
	    Bond_ptr->N_powers = sys.Inter_Types[i].N_powers;
	    Bond_ptr->powers = sys.Inter_Types[i].powers;
	    Bond_ptr->kspline = sys.Inter_Types[i].kspline;

	    flag = TRUE;

	    if (strcmp(sys.Inter_Types[i].inter_type, Bond_ptr->name) != 0) {
		printf
		    ("\nERROR: %s is not of type %s.\n  Check par.txt.\n",
		     inter_name, Bond_ptr->name);
		exit(EXIT_FAILURE);
	    }

	    break;
	}

    }

    if (flag == FALSE) {
	printf
	    ("\nERROR: Interaction Type \"%s\" not found.\n  Check par.txt.\n",
	     inter_name);
	exit(EXIT_FAILURE);
    }

}


/*****************************************************************************************
check_input(): Makes sure that the input read from the parameter file is consistent.
*****************************************************************************************/
void check_input(tW_files files, tW_system sys, Par_Flags flags,
		 tW_ref_potential ref_potential)
{
    /* Were structures found? */
    if (flags.b_Structures == FALSE) {
	printf("ERROR: No structure file listed. Check input.\n");
	exit(EXIT_FAILURE);
    }

    /* Were site types found? */
    if (flags.b_SiteTypes == FALSE) {
	printf("ERROR: No site types listed. Check input.\n");
	exit(EXIT_FAILURE);
    }

    /* Temperature found in par.txt? */
    if (flags.b_Temperature == FALSE) {
	printf("ERROR: Temperature not specified in par.txt.\n");
	exit(EXIT_FAILURE);
    }

    /* Inter_Types found in par.txt? */
    if (flags.b_Inter_Types == FALSE) {
	printf("ERROR: No Inter_Types found in par.txt.\n");
	exit(EXIT_FAILURE);
    }

    /* If PDB Mode, nrexcl must be declared. */
    if ((strcmp(files.mode, PDB_MODE) == 0) && (flags.b_nrexcl == FALSE)) {
	printf
	    ("ERROR: nrexcl must be specified in input file when using PDB mode.\n");
	exit(EXIT_FAILURE);
    }

    check_struct_files(files);

    check_site_types(sys);

    check_Inter_Types(sys);

    check_Type_Inter2(sys);

    check_Bonded_Inter_Types(sys, B_BOND_STRETCH, 2, files.mode);

    check_Bonded_Inter_Types(sys, B_ANGLE, 3, files.mode);

    check_Bonded_Inter_Types(sys, B_DIHEDRAL, 4, files.mode);

    check_Bonded_Inter_Types(sys, B_NB_PAIR_BOND, 2, files.mode);

    if (sys.flag_ref_potential == TRUE) {
	check_ref_input(ref_potential, sys);
    }
}


/*****************************************************************************************
check_Inter_Types(): Makes sure no interaction type is listed more than once.
*****************************************************************************************/
void check_Inter_Types(tW_system sys)
{
    int i, j;

    for (i = 0; i < (sys.N_Inter_Types - 1); i++) {
	for (j = i + 1; j < (sys.N_Inter_Types); j++) {
	    if (strcmp
		(sys.Inter_Types[i].inter_name,
		 sys.Inter_Types[j].inter_name) == 0) {
		printf
		    ("\nERROR: Inter_Type \"%s\" listed twice. Check par.txt.\n",
		     sys.Inter_Types[i].inter_name);
		exit(EXIT_FAILURE);
	    }
	}
    }

}


/*****************************************************************************************
check_Bonded_Inter_Types(): Makes sure information for the intramolecular interaction
is consistent.
*****************************************************************************************/
void check_Bonded_Inter_Types(tW_system sys, tW_word inter_name,
			      int N_Int_Sites, tW_word mode)
{
    int i, j;
    int N_inter = 0;
    int *inter_list = NULL;

    /* Make sure site types are listed under [Site_Types]. */
    for (i = 0; i < sys.N_Bond_Int_Types; i++) {
	if (strcmp(inter_name, sys.Bonded_Inter_Types[i].name) == 0) {
	    check_site_list(sys.Bonded_Inter_Types[i].Site_Types,
			    N_Int_Sites, sys, inter_name, "par.txt");
	}
    }

    /* Get indices for this bond interactions. */
    for (i = 0; i < sys.N_Bond_Int_Types; i++) {
	if (strcmp(inter_name, sys.Bonded_Inter_Types[i].name) == 0) {
	    N_inter++;
	    inter_list =
		(int *) erealloc(inter_list, (N_inter * sizeof(int)));
	    inter_list[N_inter - 1] = i;
	}
    }

    /* Check for repeated interactions. */
    for (i = 0; i < N_inter - 1; i++) {
	for (j = i + 1; j < N_inter; j++) {
	    check_inter_sites(sys.Bonded_Inter_Types[inter_list[i]].
			      Site_Types,
			      sys.Bonded_Inter_Types[inter_list[j]].
			      Site_Types, N_Int_Sites, inter_name);
	}
    }

    /* For InterMolec_NB_Pair, make sure that n_bonds_intramolec_pair_inter is specified for PDB mode. */
    if ((strcmp(mode, PDB_MODE) == 0)
	&& (strcmp(inter_name, B_NB_PAIR_BOND) == 0)) {
	for (i = 0; i < sys.N_Bond_Int_Types; i++) {
	    if (strcmp(sys.Bonded_Inter_Types[i].name, B_NB_PAIR_BOND) ==
		0) {
		if (sys.Bonded_Inter_Types[i].
		    n_bonds_intramolec_pair_inter == NOT_SET) {
		    printf
			("\nERROR: For %s, the number of bonds separating sites must be specified for PDB mode.\n",
			 B_NB_PAIR_BOND);
		    printf("    Interaction: %s\n",
			   sys.Bonded_Inter_Types[i].inter_name);
		    exit(EXIT_FAILURE);
		}
		if (sys.Bonded_Inter_Types[i].
		    n_bonds_intramolec_pair_inter < 2) {
		    printf
			("\nERROR: For %s, the number of bonds separating sites must be greater than 1.\n",
			 B_NB_PAIR_BOND);
		    printf("    Interaction: %s (%d)\n",
			   sys.Bonded_Inter_Types[i].inter_name,
			   sys.Bonded_Inter_Types[i].
			   n_bonds_intramolec_pair_inter);
		    exit(EXIT_FAILURE);
		}
	    }
	}
    }

    free(inter_list);
}


/*****************************************************************************************
check_Type_Inter2(): Makes sure that the pair intermolecular interaction are not listed
more than once and that they include sites listed in the site type list.
*****************************************************************************************/
void check_Type_Inter2(tW_system sys)
{
    int i, j;
    tW_word list1[2], list2[2];

    /* Make sure site types are listed under [Site_Types]. */
    for (i = 0; i < sys.N_Inter2_Types; i++) {
	strcpy(list1[0], sys.Inter2_Type_List[i].name1);
	strcpy(list1[1], sys.Inter2_Type_List[i].name2);
	check_site_list(list1, 2, sys, NB_PAIR, "par.txt");
    }


    /* Make sure interactions aren't listed twice. */
    for (i = 0; i < sys.N_Inter2_Types - 1; i++) {
	for (j = i + 1; j < sys.N_Inter2_Types; j++) {
	    strcpy(list1[0], sys.Inter2_Type_List[i].name1);
	    strcpy(list1[1], sys.Inter2_Type_List[i].name2);
	    strcpy(list2[0], sys.Inter2_Type_List[j].name1);
	    strcpy(list2[1], sys.Inter2_Type_List[j].name2);
	    check_inter_sites(list1, list2, 2, "Pair_Interactions");
	}
    }
}


/*****************************************************************************************
check_inter_sites(): Makes sure interaction was not listed more than once.
*****************************************************************************************/
void check_inter_sites(tW_word list1[], tW_word list2[], int N_words,
		       tW_word inter_name)
{
    int i;
    int cnt;

    cnt = 0;
    for (i = 0; i < N_words; i++) {
	if (strcmp(list1[i], list2[i]) != 0) {
	    break;
	}
	cnt++;
    }
    if (cnt == N_words) {
	printf("ERROR: Interaction (%s) found twice for sites: ",
	       inter_name);
	for (i = 0; i < N_words; i++) {
	    printf("%s ", list1[i]);
	}
	printf("\n");
	exit(EXIT_FAILURE);
    }

    cnt = 0;
    for (i = 0; i < N_words; i++) {
	if (strcmp(list1[i], list2[N_words - 1 - i]) != 0) {
	    break;
	}
	cnt++;
    }
    if (cnt == N_words) {
	printf
	    ("ERROR: Interaction (%s) found twice for sites (reverse order): ",
	     inter_name);
	for (i = 0; i < N_words; i++) {
	    printf("%s ", list1[i]);
	}
	printf("\n");
	exit(EXIT_FAILURE);
    }

}


/*****************************************************************************************
check_site_types(): Makes sure no site type is more than once.
*****************************************************************************************/
void check_site_types(tW_system sys)
{
    int i, j;
    int N_duplicates = 0;
    tW_word *duplicate_list = NULL;

    /* Look for site types listed more than once. */
    for (i = 0; i < sys.N_Site_Types - 1; i++) {
	for (j = i + 1; j < sys.N_Site_Types; j++) {
	    if (strcmp(sys.Site_Types[i], sys.Site_Types[j]) == 0) {
		N_duplicates++;
		duplicate_list =
		    (tW_word *) erealloc(duplicate_list,
					(N_duplicates * sizeof(tW_word)));
		strcpy(duplicate_list[N_duplicates - 1],
		       sys.Site_Types[i]);
	    }
	}
    }

    if (N_duplicates > 0) {
	/* Report repeated site types. */
	printf
	    ("\nERROR: The following site types were listed more than once.\n");
	for (i = 0; i < N_duplicates; i++) {
	    printf("  %d. %s\n", i + 1, duplicate_list[i]);
	}
	exit(EXIT_FAILURE);
    }

    free(duplicate_list);
}


/*****************************************************************************************
check_struct_files(): Makes sure no duplicate structure files listed. Also, makes sure
that the mode (the formate of the structure files) is supported.
*****************************************************************************************/
void check_struct_files(tW_files files)
{
    int i, j;
    int N_duplicates = 0;
    tW_word *duplicate_list = NULL;

    /* Look for files listed more than once. */
    for (i = 0; i < files.N_struct - 1; i++) {
	for (j = i + 1; j < files.N_struct; j++) {
	    if (strcmp(files.structures[i], files.structures[j]) == 0) {
		N_duplicates++;
		duplicate_list =
		    (tW_word *) erealloc(duplicate_list,
					(N_duplicates * sizeof(tW_word)));
		strcpy(duplicate_list[N_duplicates - 1],
		       files.structures[i]);
	    }
	}
    }

    if (N_duplicates > 0) {
	/* Report repeated structure files. */
	printf
	    ("\nERROR: The following structures files were listed more than once.\n");
	for (i = 0; i < N_duplicates; i++) {
	    printf("  %d. %s\n", i + 1, duplicate_list[i]);
	}
	exit(EXIT_FAILURE);
    }

    /* Make sure input was assigned. */
    if (strcmp(files.mode, "VOID") == 0) {
	printf
	    ("\nERROR: Input mode was not specified in parameter file.\n");
	exit(EXIT_FAILURE);
    }

    /* Make sure that mode is supported. */
    if ((strcmp(files.mode, GMX_MODE) != 0) &&
	(strcmp(files.mode, PDB_MODE) != 0)
	) {
	printf("\nERROR: Input mode '%s' not supported.\n", files.mode);
	exit(EXIT_FAILURE);
    }

    free(duplicate_list);
}


/*****************************************************************************************
get_ref_potential(): Reads and stores the information for a reference potential. 
Execution is stopped if this information has been read previously or expected information
is not found.
*****************************************************************************************/
int get_ref_potential(tW_line inp_line, tW_ref_potential * ref_potential,
		      int ref_flag, tW_word mode)
{
    FILE *fp_ref;
    Ref_Flags flags = { FALSE, FALSE, FALSE, FALSE, FALSE, FALSE };

    if (ref_flag == TRUE) {
	printf("\nERROR: Reading reference potential twice.\n");
	exit(EXIT_FAILURE);
    }

    initialize_ref_potential(ref_potential);

    initialize_ref_flags(&flags);

    fp_ref = open_ref_potential_file(inp_line);

    do {
	if (get_next_line(fp_ref, inp_line) == -1) {
	    break;
	}
	read_REF_information(fp_ref, inp_line, ref_potential, &flags,
			     mode);
    } while (1);

    if (flags.b_interpolation == FALSE) {
	printf("ERROR: Specify interpolation to use for ref. pot.\n");
	exit(EXIT_FAILURE);
    }

    return TRUE;
}

/*****************************************************************************************
get_PT():
*****************************************************************************************/
int get_PT(FILE * fp_par, tW_line inp_line, int PT_flag, tW_system * sys)
{
    int l_inp, test_sscanf;
    tW_word tmp;

    if (PT_flag == TRUE) {
	printf("\nERROR: Reading Perturbation variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PT, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->PT_var.N_PT);	/* Number of perturbation steps */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PT, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->PT_var.dPT);	/* spacing of output */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PT, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* True => calc the decomposition of b at each perturbation step */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->PT_var.flag_MMOTF_SEP = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->PT_var.flag_MMOTF_SEP = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_MMOTF under [Iterative_Inversion], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PT, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* True => calc the eigenspectrum at each perturbation step */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->PT_var.flag_eigen = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->PT_var.flag_eigen = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_eigen under [Iterative_Inversion], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PT, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->PT_var.N_eigen);	/* Number of eigenvectors to print out from each side of the spectrum (pos and neg) */

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_PT, "par.txt");

    return TRUE;		/* True => do the PT calculation */
}

/*****************************************************************************************
get_Eigen():
*****************************************************************************************/
int get_Eigen(FILE * fp_par, tW_line inp_line, int Eigen_flag,
	      tW_system * sys)
{
    int test_sscanf, l_inp;
    tW_word tmp;

    if (Eigen_flag == TRUE) {
	printf("\nERROR: Reading Eigen variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_EIGEN, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->Eigen_var.N_Eigen);	/* Number of eigenvectors to print out from each side of the spectrum (pos and neg) */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_EIGEN, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* True => print fn, bn, Mn, for various sets of eigenvectors */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->Eigen_var.flag_printn = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->Eigen_var.flag_printn = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_fnbnMn under [Eigen], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_EIGEN, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->Eigen_var.DM);	/* spacing for sets of eigenvectors from the top */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_EIGEN, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->Eigen_var.DL);	/* spacing for sets of eigenvectors from the bottom */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_EIGEN, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->Eigen_var.flag_Gbar);	/* True => calc the eigenspectrum of the Gbar matrix (i.e., take out two body terms) */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_EIGEN, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->Eigen_var.flag_norm);	/* True => normalize the nb parts of the matrix by r^2 before calc the eigenspectrum */


    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_EIGEN, "par.txt");

    return TRUE;		/* True => calc the eigenspectrum of the G matrix */
}

/*****************************************************************************************
get_SVD():
*****************************************************************************************/
int get_SVD(FILE * fp_par, tW_line inp_line, int SVD_flag, tW_system * sys)
{
    int test_sscanf, l_inp;
    tW_word tmp;

    if (SVD_flag == TRUE) {
	printf("\nERROR: Reading SVD variables twice.\n");
	exit(EXIT_FAILURE);
    }
    strcpy(sys->SOLN_var.SOLN_METH, SVD_SOLV);

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_SVD, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->SVD_var.rcond);	/* throw out all eigenvalues below this value in SVD calculation */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_SVD, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* explicitly calculate and print out the singular values */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->SVD_var.flag_printSV = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->SVD_var.flag_printSV = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_printSV under [SVD], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_SVD, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* explicitly calculate and print out the singular values */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->SVD_var.flag_printevecs = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->SVD_var.flag_printevecs = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_printevecs under [SVD], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_SVD, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* explicitly calculate and print out the singular values */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->SVD_var.flag_solve = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->SVD_var.flag_solve = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_solve under [SVD], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_SVD, "par.txt");

    return TRUE;		/* True => do the matrix inversion using SVD as well as LU decomp. Also, calculate the condition number with all possible preconditioning */
}


/*****************************************************************************************
get_TPR():
*****************************************************************************************/
int get_TPR(FILE * fp_par, tW_line inp_line, int TPR_flag, tW_system * sys)
{
    int test_sscanf;

    if (TPR_flag == TRUE) {
	printf("ERROR: Reading TPR variables twice.\n");
	exit(EXIT_FAILURE);
    }

    delete_directive_from_line(inp_line);

    test_sscanf = sscanf(inp_line, " %d ", &sys->TPR_var.N_TPR);

    if (test_sscanf != 1) {
	printf("ERROR: Number of TPR files not found in par.txt.\n");
	exit(EXIT_FAILURE);
    }

    sys->TPR_var.TPR_files = (tW_word *) ecalloc(sys->TPR_var.N_TPR, sizeof(tW_word));	/* allocate memory for the filenames */

    sys->TPR_var.TPR_files =
	get_word_list(fp_par, KEY_TPR, sys->TPR_var.N_TPR);

    return TRUE;		/* replace default tpr filenames */
}

/*****************************************************************************************
get_TPR_EXCL():
*****************************************************************************************/
int get_TPR_EXCL(FILE * fp_par, tW_line inp_line, int TPR_EXCL_flag,
		 tW_system * sys)
{
    int test_sscanf, l_inp;

    if (TPR_EXCL_flag == TRUE) {
	printf("ERROR: Reading TPR_EXCL variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_TPR_EXCL, inp_line);
    test_sscanf = sscanf(inp_line, "%s", sys->TPR_EXCL_var.TPR_excl);	/* The TPR filename used for exclusions */

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_TPR_EXCL, "par.txt");

    return TRUE;		/* use alternative tpr file for exclusions */
}

/*****************************************************************************************
get_MT():
*****************************************************************************************/
int get_MT(FILE * fp_par, tW_line inp_line, int MT_flag, tW_system * sys)
{
    int test_sscanf, l_inp;
    tW_word tmp;

    if (MT_flag == TRUE) {
	printf("\nERROR: Reading Metric_Tensor variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_MT, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* print out the metric tensor blocks */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->MT_var.flag_print = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->MT_var.flag_print = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_print under [Metric_Tensor], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_MT, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* normalize the metric tensor by r^2 for nonbonded interactions */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->MT_var.flag_norm = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->MT_var.flag_norm = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_norm under [Metric_Tensor], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_MT, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* calculate the unweighted version of the matrix, Mcnt */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->MT_var.flag_Mcnt = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->MT_var.flag_Mcnt = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_Mcnt under [Metric_Tensor], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_MT, "par.txt");

    return TRUE;		/* True => print out the metric tensor related files */
}

/*****************************************************************************************
get_MFD():
*****************************************************************************************/
int get_MFD(FILE * fp_par, tW_line inp_line, int MFD_flag, tW_system * sys)
{
    int /*test_sscanf, */ l_inp;

    if (MFD_flag == TRUE) {
	printf
	    ("\nERROR: Reading Mean_Force_Decomposition variables twice.\n");
	exit(EXIT_FAILURE);
    }

    /* currently no variables to read in */

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_MFD, "par.txt");

    return TRUE;		/* True => calculate the decomposition of b */
}

/*****************************************************************************************
get_CalcMODE():
*****************************************************************************************/
int get_CalcMODE(FILE * fp_par, tW_line inp_line, int CalcMODE_flag,
		 tW_system * sys)
{
    int test_sscanf, l_inp;
    tW_word mode;

    if (CalcMODE_flag == TRUE) {
	printf("\nERROR: Reading Calculation_Mode variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_CalcMODE, inp_line);
    test_sscanf = sscanf(inp_line, "%s", mode);	/* The Calculation Mode */
    if (strcmp(mode, FULL) == 0) {
	sys->CalcMODE_var.CalcMODE = IFULL;
    } else if (strcmp(mode, FIRST_HALF) == 0) {
	sys->CalcMODE_var.CalcMODE = IFIRST_HALF;
    } else if (strcmp(mode, SECOND_HALF) == 0) {
	sys->CalcMODE_var.CalcMODE = ISECOND_HALF;
    } else if (strcmp(mode, TEST_INP) == 0) {
	sys->CalcMODE_var.CalcMODE = ITEST_INP;
    } else {
	printf("ERROR: %s is not a valid CalcMODE, see par.txt \n", mode);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_CalcMODE, "par.txt");

    return TRUE;		/* True => specify the calculation mode */
}

/*****************************************************************************************
get_PC():
*****************************************************************************************/
int get_PC(FILE * fp_par, tW_line inp_line, int PC_flag, tW_system * sys)
{
    int test_sscanf, l_inp;
    tW_word PC_type, tmp;

    if (PC_flag == TRUE) {
	printf("\nERROR: Reading Preconditioning variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PC, inp_line);
    test_sscanf = sscanf(inp_line, "%s", PC_type);	/* The Right Preconditioning Type */
    if (strcmp(PC_type, "no") == 0) {
	strcpy(sys->PC_var.RPC, "NO");
    } else if (strcmp(PC_type, "No") == 0) {
	strcpy(sys->PC_var.RPC, "NO");
    } else if (strcmp(PC_type, "NO") == 0) {
	strcpy(sys->PC_var.RPC, "NO");
    } else if (strcmp(PC_type, "dimless") == 0) {
	strcpy(sys->PC_var.RPC, "dimless");
    } else if (strcmp(PC_type, "colnorm") == 0) {
	strcpy(sys->PC_var.RPC, "colnorm");
    } else if (strcmp(PC_type, "MTvar") == 0) {
	strcpy(sys->PC_var.RPC, "MTvar");
    } else {
	printf("ERROR: %s is not a valid PC_type, see par.txt \n",
	       PC_type);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PC, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* True => rescale so that norm phi = 1 */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->PC_var.flag_normphi = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->PC_var.flag_normphi = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_normphi under [Preconditiong], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PC, inp_line);
    test_sscanf = sscanf(inp_line, "%s", PC_type);	/* The Left Preconditioning Type */
    if (strcmp(PC_type, "no") == 0) {
	strcpy(sys->PC_var.LPC, "NO");
    } else if (strcmp(PC_type, "No") == 0) {
	strcpy(sys->PC_var.LPC, "NO");
    } else if (strcmp(PC_type, "NO") == 0) {
	strcpy(sys->PC_var.LPC, "NO");
    } else if (strcmp(PC_type, "dimless") == 0) {
	strcpy(sys->PC_var.LPC, "dimless");
    } else if (strcmp(PC_type, "rowmax") == 0) {
	strcpy(sys->PC_var.LPC, "rowmax");
    } else if (strcmp(PC_type, "bvar") == 0) {
	strcpy(sys->PC_var.LPC, "bvar");
    } else {
	printf("ERROR: %s is not a valid PC_type, see par.txt \n",
	       PC_type);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_PC, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* True => rescale so that norm phi = 1 */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->PC_var.flag_normb = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->PC_var.flag_normb = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_normb under [Preconditiong], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_PC, "par.txt");

    return TRUE;		/* True => specify the Preconditioning for the matrix inversion */
}

/*****************************************************************************************
get_MEM():
*****************************************************************************************/
int get_MEM(FILE * fp_par, tW_line inp_line, int MEM_flag, tW_system * sys)
{
    int test_sscanf, l_inp;
    int flagf = 0;
    int flags = 0;
    tW_word tmp;

    if (MEM_flag == TRUE) {
	printf("\nERROR: Reading Memory variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_MEM, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* TRUE => run the program in Low Memory mode */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->MEM_var.flag_LOWMEM = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->MEM_var.flag_LOWMEM = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_Mcnt under [Metric_Tensor], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_MEM, inp_line);
    test_sscanf = sscanf(inp_line, "%s", sys->MEM_var.info);	/* Use either structures or forces when in Low Memory mode */
    if (sys->MEM_var.flag_LOWMEM == TRUE) {
	if (strcmp(sys->MEM_var.info, "forces") != 0) {
	    flagf++;
	}
	if (strcmp(sys->MEM_var.info, "structures") != 0) {
	    flags++;
	}
	printf("%s\n", sys->MEM_var.info);
	if ((flagf != 0) && (flags != 0)) {
	    printf
		("ERROR: You are in Low Memory mode and sys->MEM_var.info does not match structures or forces. \n");
	    exit(0);
	}
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_MEM, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* FALSE => do not calculate the answer to each topology separately */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->MEM_var.flag_mult_top = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->MEM_var.flag_mult_top = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_Mcnt under [Metric_Tensor], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_MEM, "par.txt");

    return TRUE;		/* True => Using Memory options */
}

/*****************************************************************************************
get_SOLN():
*****************************************************************************************/
int get_SOLN(FILE * fp_par, tW_line inp_line, int SOLN_flag,
	     tW_system * sys)
{
    int test_sscanf, l_inp;

    if (SOLN_flag == TRUE) {
	printf("\nERROR: Reading Solution_Method variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_SOLN, inp_line);
    test_sscanf = sscanf(inp_line, "%s", sys->SOLN_var.SOLN_METH);	/* solve the matrix inversion with this technique */
    if ((strcmp(sys->SOLN_var.SOLN_METH, LU) != 0)
	&& (strcmp(sys->SOLN_var.SOLN_METH, UU) != 0)
	&& (strcmp(sys->SOLN_var.SOLN_METH, CHOLESKY) != 0)
	&& (strcmp(sys->SOLN_var.SOLN_METH, SVD_SOLV) != 0) ) {
	printf("ERROR: %s is not a valid SOLN_METH, see par.txt. \n",
	       sys->SOLN_var.SOLN_METH);
	exit(0);
    }

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_SOLN, "par.txt");

    return TRUE;		/* True => specify the Solution Method/s for the matrix inversion */
}

/*****************************************************************************************
get_FRAMEWEIGHT():
*****************************************************************************************/
int get_FRAMEWEIGHT(FILE * fp_par, tW_line inp_line, int FRAMEWEIGHT_flag,
	     tW_system * sys)
{
    int test_sscanf, l_inp;

    if (FRAMEWEIGHT_flag == TRUE) {
	printf("\nERROR: Reading Frame Weighting variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_FRAMEWEIGHT, inp_line);
    test_sscanf = sscanf(inp_line, "%s", sys->FRAMEWEIGHT_var.FRAMEWEIGHT); 	/* The input trajectory was simulated with this ensemble */

    if ( (strcmp(sys->FRAMEWEIGHT_var.FRAMEWEIGHT, NONE) != 0)
	&& (strcmp(sys->FRAMEWEIGHT_var.FRAMEWEIGHT, NPT) != 0)) {
	printf("\nERROR: %s is not a valid FRAMEWEIGHT, see par.txt. \n",
		sys->FRAMEWEIGHT_var.FRAMEWEIGHT);
	exit(0);	
    }
	
    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_FRAMEWEIGHT, "par.txt");

    return TRUE;		/* True => specify the ensemble of the input trr file*/
}

/*****************************************************************************************
get_ERR():
*****************************************************************************************/
int get_ERR(FILE * fp_par, tW_line inp_line, int ERR_flag, tW_system * sys)
{
    int test_sscanf, l_inp;
    tW_word tmp;

    if (ERR_flag == TRUE) {
	printf("\nERROR: Reading Calculate_Errors variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_ERR, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->ERR_var.FACT = 'E';
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->ERR_var.FACT = 'N';
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_Equil under [Error_Estimates], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }


    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_ERR, "par.txt");

    return TRUE;		/* True => calculate errors associated with the matrix inversion */
}

/*****************************************************************************************
get_REF():
*****************************************************************************************/
int get_REF(FILE * fp_par, tW_line inp_line, int REF_flag, tW_system * sys)
{
    int i;
    int test_sscanf, l_inp;
    tW_word tmp;

    if (REF_flag == TRUE) {
	printf("\nERROR: Reading Reference_Options variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_REF, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->REF_var.flag_calcbref = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->REF_var.flag_calcbref = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_calcbref under [Reference_Options], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_REF, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->REF_var.flag_readbref = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->REF_var.flag_readbref = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_readbref under [Reference_Options], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_REF, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->REF_var.flag_splitfiles = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->REF_var.flag_splitfiles = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_splitfiles under [Reference_Options], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_REF, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->REF_var.flag_reftrr = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->REF_var.flag_reftrr = FALSE;
    } else {
	printf ("ERROR: %s is not a valid option for flag_reftrr under [Reference_Options], see par.txt \n", tmp);
	exit(EXIT_FAILURE);
    }

    // This is for a single file
    //l_inp = get_next_line( fp_par, inp_line );
    //check_inp_line( "par.txt", KEY_REF, inp_line );
    //test_sscanf = sscanf( inp_line, "%s", sys->REF_var.reftrr_fnm );

    sys->REF_var.reftrr_fnm = (tW_word *) ecalloc(sys->REF_var.N_fnm, sizeof(tW_word));	/* allocate memory for the filenames */

    sys->REF_var.reftrr_fnm = get_word_list(fp_par, KEY_REF, sys->REF_var.N_fnm);

    if (sys->REF_var.flag_reftrr) /*NJD 5.25.15 - we only need to check for these files if we plan to use them */
    {
        for (i = 0; i < sys->REF_var.N_fnm; i++) {
    	    if (!file_exists(sys->REF_var.reftrr_fnm[i])) {
    	        printf ("ERROR: Unable to open reference file '%s' for reading.\n", sys->REF_var.reftrr_fnm[i]);
    	        exit(EXIT_FAILURE);
    	    }
        }
    }

    return TRUE;		/* True => use reference options specified in this directive */
}

/*****************************************************************************************
get_TRIM():
*****************************************************************************************/
int get_TRIM(FILE * fp_par, tW_line inp_line, int TRIM_flag,
	     tW_system * sys)
{
    int test_sscanf, l_inp;

    if (TRIM_flag == TRUE) {
	printf("\nERROR: Reading TRIM variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_TRIM, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->TRIM_var.FE);	/* Trim out basis functions which are sampled less than FE */

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_TRIM, "par.txt");

    return TRUE;		/* True => Replace default FLOAT_EPS with FE */
}

/*****************************************************************************************
get_CHISQD():
*****************************************************************************************/
int get_CHISQD(FILE * fp_par, tW_line inp_line, int CHISQD_flag,
	       tW_system * sys)
{
    int test_sscanf, l_inp;

    if (CHISQD_flag == TRUE) {
	printf("\nERROR: Reading CHISQD variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_CHISQD, inp_line);
    test_sscanf = sscanf(inp_line, "%s", sys->CHISQD_var.force_fnm);

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_CHISQD, "par.txt");

    return TRUE;		/* True => read in a force file to use with the Chi2 calculation */
}

/*****************************************************************************************
get_REG():
*****************************************************************************************/
int get_REG(FILE * fp_par, tW_line inp_line, int REG_flag, tW_system * sys)
{
    int test_sscanf, l_inp;

    if (REG_flag == TRUE) {
	printf("\nERROR: Reading Regularization variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_REG, inp_line);
    test_sscanf = sscanf(inp_line, "%s", sys->REG_var.type);

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_REG, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->REG_var.Nmax);	/* max iterations when optimizing regularization parameters */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_REG, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->REG_var.tau_alpha);	/* tolerance in the convergence of the regularization matrix parameters */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_REG, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->REG_var.tau_beta);	/* tolerance in the convergence of the precision parameter */

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_REG, "par.txt");

    return TRUE;		/* True => Use Bayesian Inference Method */
}

/*****************************************************************************************
get_RESCALE():
*****************************************************************************************/
int get_RESCALE(FILE * fp_par, tW_line inp_line, int RESCALE_flag,
		tW_system * sys)
{
    int test_sscanf, l_inp;

    if (RESCALE_flag == TRUE) {
	printf("\nERROR: Reading RESCALE variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_RESCALE, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->RESCALE_var.Nmax);	/* max iterations when performing force rescale right preconditioning */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_RESCALE, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->RESCALE_var.tau_phi);	/* tolerance in the convergence of the force rescale parameters */

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_RESCALE, "par.txt");

    return TRUE;		/* True => Do force rescale right preconditioning */
}

/*****************************************************************************************
get_CONSTRAIN():
*****************************************************************************************/
int get_CONSTRAIN(FILE * fp_par, tW_line inp_line, int CONSTRAIN_flag,
		  tW_system * sys)
{
    int test_sscanf, l_inp;

    if (CONSTRAIN_flag == TRUE) {
	printf("\nERROR: Reading CONSTRAIN variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_CONSTRAIN, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->CONSTRAIN_var.Nmax);	/* max iterations when constraining dihedrals */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_CONSTRAIN, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->CONSTRAIN_var.tau_dih);	/* tolerance in the convergence of the contraint scaling parameter */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_CONSTRAIN, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->CONSTRAIN_var.lambda);	/* initial value of the constraint scaling parameter */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_CONSTRAIN, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->CONSTRAIN_var.dlambda);	/* amount to increase the constraint parameter each iteraction */

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_CONSTRAIN, "par.txt");

    return TRUE;		/* True => Constrain the dihedral force coefficients of each dih interaction to sum to zero */
}

/*****************************************************************************************
get_ITER():
*****************************************************************************************/
int get_ITER(FILE * fp_par, tW_line inp_line, int ITER_flag,
	     tW_system * sys)
{
    tW_word tmp;
    int test_sscanf, l_inp;

    if (ITER_flag == TRUE) {
	printf("\nERROR: Reading Iter-gYBG variables twice.\n");
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_ITER, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->ITER_var.flag_AAM2);	/* TRUE => use AA 2-body terms for bonds and angles */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_ITER, inp_line);
    test_sscanf = sscanf(inp_line, "%s", sys->ITER_var.AAM2_fnm);	/* filename to read AA 2-body terms */

    l_inp = get_next_line(fp_par, inp_line);
    check_inp_line("par.txt", KEY_ITER, inp_line);
    test_sscanf = sscanf(inp_line, "%s", tmp);	/* TRUE => calc diff in projected forces from current and previous step  */
    if ((strcmp(tmp, "YES") == 0) || (strcmp(tmp, "Yes") == 0)
	|| (strcmp(tmp, "yes") == 0)) {
	sys->ITER_var.flag_bsolnerr = TRUE;
    } else if ((strcmp(tmp, "NO") == 0) || (strcmp(tmp, "No") == 0)
	       || (strcmp(tmp, "no") == 0)) {
	sys->ITER_var.flag_bsolnerr = FALSE;
    } else {
	printf
	    ("ERROR: %s is not a valid option for flag_bsolnerr under [Iterative_gYBG], see par.txt \n",
	     tmp);
	exit(EXIT_FAILURE);
    }

    l_inp = get_next_line(fp_par, inp_line);
    test_end_directive(inp_line, KEY_ITER, "par.txt");

    return TRUE;		/* True => Use options for iter-gYBG procedure */
}

/*****************************************************************************************
initialize_ref_potential(): Initializes the variables in the ref_potential data 
structure.
*****************************************************************************************/
void initialize_ref_potential(tW_ref_potential * ref_potential)
{
    ref_potential->N_nb_pair_inter = 0;
    ref_potential->nb_pair_inter = NULL;
    ref_potential->N_bondstretch_inter = 0;
    ref_potential->bondstretch_inter = NULL;
    ref_potential->N_angle_inter = 0;
    ref_potential->angle_inter = NULL;
    ref_potential->N_dihedral_inter = 0;
    ref_potential->dihedral_inter = NULL;
    ref_potential->N_IntraMolecPairs_inter = 0;
    ref_potential->IntraMolecPairs_inter = NULL;
    ref_potential->interpolation_index = -1;
}


/*****************************************************************************************
get_ref_inter(): Reads information for a interaction type in the reference potential.
*****************************************************************************************/
tW_ref_inter *get_ref_inter(tW_word inter_name, int N_Inter_Sites,
			    int N_Inter_Types, FILE * fp,
			    tW_ref_inter * inter_list,
			    tW_ref_potential * ref_potential, tW_word mode)
{
    int i_inter;
    tW_line inp_line;

    inter_list =
	(tW_ref_inter *) emalloc(N_Inter_Types * sizeof(tW_ref_inter));

    for (i_inter = 0; i_inter < N_Inter_Types; i_inter++) {
	get_next_line(fp, inp_line);

	check_inp_line("reference potential input", inter_name, inp_line);

	get_ref_potential_info(&(inter_list[i_inter]), N_Inter_Sites,
			       inter_name, inp_line, ref_potential, mode);
    }

    get_next_line(fp, inp_line);

    test_end_directive(inp_line, inter_name, "reference potential input");

    return inter_list;
}


/*****************************************************************************************
get_ref_potential_info(): Stores information for a interaction in the reference
potential.
*****************************************************************************************/
void get_ref_potential_info(tW_ref_inter * inter_ptr, int N_Inter_Sites,
			    tW_word inter_name, tW_line inp_line,
			    tW_ref_potential * ref_potential, tW_word mode)
{
    int test_sscanf;
    tW_word fname;

    strcpy(inter_ptr->name, inter_name);

    inter_ptr->N_inter_sites = N_Inter_Sites;

    inter_ptr->site_types =
	(tW_word *) emalloc(N_Inter_Sites * sizeof(tW_word));

    if (N_Inter_Sites == 2) {
	sscanf(inp_line, "%s %s %s", inter_ptr->site_types[0],
	       inter_ptr->site_types[1], fname);
	inter_ptr->bonds_IntraMolecPairs = 0;
	if ((strcmp(mode, PDB_MODE) == 0) &&
	    (strcmp(inter_name, B_NB_PAIR_BOND) == 0)
	    ) {
	    test_sscanf =
		sscanf(inp_line, "%*s %*s %*s %d",
		       &(inter_ptr->bonds_IntraMolecPairs));
	    if (test_sscanf != 1) {
		printf
		    ("\nERROR: Must specify the number of bonds separating sites for\n");
		printf
		    ("PDB mode for intramolecular pairs for reference potential.\n");
		exit(EXIT_FAILURE);
	    }
	}
    } else if (N_Inter_Sites == 3) {
	inter_ptr->bonds_IntraMolecPairs = 0;
	sscanf(inp_line, "%s %s %s %s", inter_ptr->site_types[0],
	       inter_ptr->site_types[1], inter_ptr->site_types[2], fname);
    } else if (N_Inter_Sites == 4) {
	inter_ptr->bonds_IntraMolecPairs = 0;
	sscanf(inp_line, "%s %s %s %s %s", inter_ptr->site_types[0],
	       inter_ptr->site_types[1], inter_ptr->site_types[2],
	       inter_ptr->site_types[3], fname);
    } else {
	printf
	    ("\nERROR: Number of sites are incorrect for ref. interaction: %s\n",
	     inter_name);
	exit(EXIT_FAILURE);
    }

    get_ref_potential_ptr(fname, inter_ptr, ref_potential);
}


/*****************************************************************************************
get_ref_potential_ptr(): Stores information read from ref. force file and sets pointers.
*****************************************************************************************/
void get_ref_potential_ptr(tW_word fname, tW_ref_inter * inter_ptr,
			   tW_ref_potential * ref_potential)
{
    int f_index;
    tW_ref_force *ref_force;

    f_index =
	match_word(ref_potential->N_files, fname, ref_potential->fname);

    if (f_index == -1) {
	printf("\nERROR: File '%s' not found for ref. forces.\n", fname);
	exit(EXIT_FAILURE);
    }

    ref_force = &(ref_potential->forces[f_index]);

    inter_ptr->N_pts = ref_force->N_pts;
    inter_ptr->x = ref_force->x;
    inter_ptr->f = ref_force->f;
    inter_ptr->dx = ref_force->dx;
    inter_ptr->x_0 = ref_force->x_0;
    inter_ptr->x_max = ref_force->x_max;
}


/*****************************************************************************************
read_ref_potential(): Stores information for a interaction in the reference potential.
*****************************************************************************************/
void read_ref_potential(tW_word fname, tW_ref_force * inter_ptr,
			tW_word inter_type)
{
    FILE *fp_ref = open_file(fname, 'r');

    inter_ptr->N_pts = get_no_points(fp_ref);

    inter_ptr->x = ecalloc(inter_ptr->N_pts, sizeof(double));

    inter_ptr->f = ecalloc(inter_ptr->N_pts, sizeof(double));

    get_force(fp_ref, inter_ptr, inter_type);

}


/*****************************************************************************************
get_interpolation_index(): Returns the interpolation index or stops execution if the
specified interpolation type is not supported.
*****************************************************************************************/
int get_interpolation_index(tW_line inp_line)
{
    tW_word word;

    sscanf(inp_line, " [ Interpolation ] %s", word);

    if (strcmp(word, REF_NO_INTERPOL_N) == 0) {
	return REF_NO_INTERPOL_I;
    }

    else if (strcmp(word, REF_LINEAR_INTERPOL_N) == 0) {
	return REF_LINEAR_INTERPOL_I;
    }

    else {
	printf("\nERROR: Interpolation scheme not supported yet.\n");
	exit(EXIT_FAILURE);
    }

}


/*****************************************************************************************
get_no_points(): Returns the number of lines in a file before the first blank line.
*****************************************************************************************/
int get_no_points(FILE * fp)
{
    int n_pts = 0;
    tW_line line;

    do {

	if (get_next_line(fp, line) == -1) {
	    break;
	}

	n_pts++;

    } while (1);

    return n_pts;
}


/*****************************************************************************************
get_force(): Reads the force for a given interaction for the reference potential. 
Converts angles and forces to radians for angle bending and dihedrals interactions.
*****************************************************************************************/
int get_force(FILE * fp, tW_ref_force * inter_ptr, tW_word inter_type)
{
    int i = 0;
    int l_inp;
    double dx_test;
    tW_line line;

    rewind(fp);

    for (i = 0; i < inter_ptr->N_pts; i++) {
	l_inp = get_next_line(fp, line);

	if (l_inp == -1) {
	    printf("i: %d   N_pts: %d\n", i, inter_ptr->N_pts);
	    printf("\nERROR: Did not read all lines for a ref. pot.\n");
	    exit(EXIT_FAILURE);
	}

	sscanf(line, "%lf %lf", &(inter_ptr->x[i]), &(inter_ptr->f[i]));
    }

    inter_ptr->x_0 = inter_ptr->x[0];

    inter_ptr->x_max = inter_ptr->x[inter_ptr->N_pts - 1];

    inter_ptr->dx = inter_ptr->x[1] - inter_ptr->x[0];

    for (i = 0; i < inter_ptr->N_pts - 1; i++) {
	dx_test = inter_ptr->x[i + 1] - inter_ptr->x[i];

	if (fabs(dx_test - inter_ptr->dx) > FLOAT_EPS) {
	    printf("\ndx_test: %.8e  dx: %.8e  diff: %.12e\n",
		   dx_test, inter_ptr->dx, fabs(dx_test - inter_ptr->dx));
	    printf
		("\nERROR: Check ref. pot. Only supports even grid spacings.\n");
	    exit(EXIT_FAILURE);
	}
    }

    if ((strcmp(inter_type, B_ANGLE) == 0) ||
	(strcmp(inter_type, B_DIHEDRAL) == 0)
	) {
	inter_ptr->x_0 *= M_PI / 180.0;
	inter_ptr->x_max *= M_PI / 180.0;
	inter_ptr->dx *= M_PI / 180.0;
	for (i = 0; i < inter_ptr->N_pts; i++) {
	    inter_ptr->x[i] *= M_PI / 180.0;
	    inter_ptr->f[i] *= 180.0 / M_PI;
	}
    }

    return 0;
}


/*****************************************************************************************
summarize_input_ref_potential(): Summarizes the information store for the reference 
potential.
*****************************************************************************************/
void summarize_input_ref_potential(FILE * fp_log, tW_system sys,
				   tW_ref_potential ref_potential)
{
    int n_types =
	ref_potential.N_nb_pair_inter + ref_potential.N_bondstretch_inter +
	ref_potential.N_angle_inter + ref_potential.N_dihedral_inter +
	ref_potential.N_IntraMolecPairs_inter;
    int ctr = 0;

    fprintf(fp_log, "\n\n  Summarize reference potential.\n");
    print_interpolation_type(ref_potential.interpolation_index, fp_log);
    fprintf(fp_log,
	    "    There are %d interaction types for the reference potential.\n",
	    n_types);

    if (ref_potential.N_nb_pair_inter > 0) {
	print_ref_inter(ref_potential.nb_pair_inter, fp_log, &ctr,
			ref_potential.N_nb_pair_inter);
    }

    if (ref_potential.N_bondstretch_inter > 0) {
	print_ref_inter(ref_potential.bondstretch_inter, fp_log, &ctr,
			ref_potential.N_bondstretch_inter);
    }

    if (ref_potential.N_angle_inter > 0) {
	print_ref_inter(ref_potential.angle_inter, fp_log, &ctr,
			ref_potential.N_angle_inter);
    }

    if (ref_potential.N_dihedral_inter > 0) {
	print_ref_inter(ref_potential.dihedral_inter, fp_log, &ctr,
			ref_potential.N_dihedral_inter);
    }

    if (ref_potential.N_IntraMolecPairs_inter > 0) {
	print_ref_inter(ref_potential.IntraMolecPairs_inter, fp_log, &ctr,
			ref_potential.N_IntraMolecPairs_inter);
    }

}


/*****************************************************************************************
print_ref_inter(): Prints the information stored for a interaction in the reference 
potential.
*****************************************************************************************/
void print_ref_inter(tW_ref_inter * ref_inter, FILE * fp_log, int *ctr,
		     int N_inter)
{
    int i;
    int inter;

    for (inter = 0; inter < N_inter; inter++) {
	*ctr = *ctr + 1;
	fprintf(fp_log, "      %d. Reference Interaction Type: %s\n", *ctr,
		ref_inter[inter].name);
	fprintf(fp_log, "        N_Inter_Sites: %-4d    Sites: ",
		ref_inter[inter].N_inter_sites);
	for (i = 0; i < ref_inter[inter].N_inter_sites; i++) {
	    fprintf(fp_log, "%s", ref_inter[inter].site_types[i]);
	    if (i < ref_inter[inter].N_inter_sites - 1) {
		fprintf(fp_log, "-");
	    }
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log,
		"        N_grid_pnts: %-5d  dx: %-12.6f  x_min: %-12.6f  x_max: %-12.6f",
		ref_inter[inter].N_pts, ref_inter[inter].dx,
		ref_inter[inter].x_0, ref_inter[inter].x_max);
//    if ( ref_inter[inter].bonds_IntraMolecPairs > 0 )
//    {
	fprintf(fp_log, "  nr_bonds: %d",
		ref_inter[inter].bonds_IntraMolecPairs);
//    }

	fprintf(fp_log, "\n\n");
    }
}


/*****************************************************************************************
print_interpolation_type(): Prints the interpolation type for the reference potential.
*****************************************************************************************/
void print_interpolation_type(int interpolation_index, FILE * fp_log)
{
    if (interpolation_index == REF_NO_INTERPOL_I) {
	fprintf(fp_log, "  Interopolation: %s\n", REF_NO_INTERPOL_N);
    }

    else if (interpolation_index == REF_LINEAR_INTERPOL_I) {
	fprintf(fp_log, "  Interopolation: %s\n", REF_LINEAR_INTERPOL_N);
    }
}


/*****************************************************************************************
read_mode(): Reads and stores the mode for the calculation. Currently, only structures
read from GROMACS trajectories or PDB files are supproted.
*****************************************************************************************/
int read_mode(int b_mode, tW_line inp_line, tW_files * files)
{
    int test_sscanf;

    if (b_mode == TRUE) {
	printf("ERROR: Reading Mode twice.\n");
	exit(EXIT_FAILURE);
    }

    delete_directive_from_line(inp_line);

    test_sscanf = sscanf(inp_line, " %s ", files->mode);

    if (test_sscanf != 1) {
	printf("ERROR: Mode not found.\n");
	exit(EXIT_FAILURE);
    }

    return TRUE;
}


/*****************************************************************************************
read_nrexcl(): Reads the number of bonds for determining exclusion lists. Not used in the
GROMACS loop, but used in building topology for the PDB loop.
*****************************************************************************************/
int read_nrexcl(int b_nrexcl, tW_line inp_line, tW_system * sys)
{
    int test_sscanf;

    if (b_nrexcl == TRUE) {
	printf("ERROR: Reading nrexcl twice.\n");
	exit(EXIT_FAILURE);
    }

    delete_directive_from_line(inp_line);

    test_sscanf = sscanf(inp_line, " %d ", &(sys->nrexcl));

    if (test_sscanf != 1) {
	printf("ERROR: nrexcl not found.\n");
	exit(EXIT_FAILURE);
    }

    return TRUE;
}


/*****************************************************************************************
read_temperature(): Reads the temperature from the parameter file. If the temperature has
been previously read or if unable to read temperature, execution is stopped. Otherwise
TRUE is returned.
*****************************************************************************************/
int read_temperature(int b_Temperature, tW_line inp_line, tW_system * sys)
{
    int test_sscanf;

    if (b_Temperature == TRUE) {
	printf("ERROR: Reading temperature twice.\n");
	exit(EXIT_FAILURE);
    }

    delete_directive_from_line(inp_line);

    test_sscanf = sscanf(inp_line, " %lf ", &(sys->Temperature));

    if (test_sscanf != 1) {
	printf("ERROR: Temperature not found.\n");
	exit(EXIT_FAILURE);
    }

    return TRUE;
}


/*****************************************************************************************
read_structures(): Reads in the file names that contain structures. Stops execution if 
the file names have been read previously or the number of file names is not found. 
Otherwise, calls get_word_list() to read the file names that contain lists of file names 
of files that contain structures (confusing?). Then, calls get_structures() to read
the lists of files that contain structures. The idea is that there can be several file 
listed in the parameter file under [Structures] and these files can contain different
lists of files that actually contain the coordinates (and forces). Returns TRUE.
*****************************************************************************************/
int read_structures(int b_Structures, tW_line inp_line, tW_files * files,
		    FILE * fp_par)
{
    int test_sscanf;

    if (b_Structures == TRUE) {
	printf("ERROR: Reading structures twice.\n");
	exit(EXIT_FAILURE);
    }

    delete_directive_from_line(inp_line);

    test_sscanf = sscanf(inp_line, " %d ", &files->N_struct_file);

    if (test_sscanf != 1) {
	printf("ERROR: Number of structure files not found in par.txt.\n");
	exit(EXIT_FAILURE);
    }

    files->struct_file =
	get_word_list(fp_par, KEY_STRUCT, files->N_struct_file);
    files->structures = get_structures(files, &(files->N_struct));

    return TRUE;
}


/*****************************************************************************************
read_site_types(): Reads the site types listed in the parameter file. Stops execution
if the site types has been read previously or the number of site types is not found.
Returns TRUE. Calls get_word_list to retrieve the list from the parameter file.
*****************************************************************************************/
int read_site_types(int b_SiteTypes, tW_line inp_line, tW_system * sys,
		    FILE * fp_par)
{
    int test_sscanf;

    if (b_SiteTypes == TRUE) {
	printf("ERROR: Reading SiteTypes twice.\n");
	exit(EXIT_FAILURE);
    }

    delete_directive_from_line(inp_line);

    test_sscanf = sscanf(inp_line, " %d ", &sys->N_Site_Types);

    if (test_sscanf != 1) {
	printf("ERROR: Number of site types not found in par.txt.\n");
	exit(EXIT_FAILURE);
    }

    sys->Site_Types = get_word_list(fp_par, KEY_SITES, sys->N_Site_Types);


    return TRUE;
}


/*****************************************************************************************
read_nb_pair_inter(): Reads in the non-bonded intermolecular pair interactions. Returns
TRUE. Stops execution if interactions have been read previously.
*****************************************************************************************/
int read_nb_pair_inter(int b_Pairs, tW_line inp_line, tW_system * sys,
		       FILE * fp_par)
{

    if (b_Pairs == TRUE) {
	printf("\nERROR: Reading Pair_Interactions twice.\n");
	exit(EXIT_FAILURE);
    }

    sys->N_Inter2_Types = get_no_inter_types(NB_PAIR, inp_line);

    sys->Inter2_Type_List =
	get_Type_Inter2_List(fp_par, inp_line, *sys, sys->N_Inter2_Types,
			     NB_PAIR);

    return TRUE;
}


/*****************************************************************************************
read_BOND_stretch_inter(): Reads the bond stretch intramolecular interactions from the 
parameter file. Stops execution if they have been read previously. Returns TRUE.
*****************************************************************************************/
int read_BOND_stretch_inter(int b_Bonds, tW_line inp_line, tW_system * sys,
			    FILE * fp_par)
{
    int N_Inter_Types;

    if (b_Bonds == TRUE) {
	printf("\nERROR: Reading Bonds twice.\n");
	exit(EXIT_FAILURE);
    }

    N_Inter_Types = get_no_inter_types(B_BOND_STRETCH, inp_line);

    sys->Bonded_Inter_Types =
	get_Bond_Inter_Types(fp_par, inp_line, *sys,
			     &(sys->N_Bond_Int_Types),
			     sys->Bonded_Inter_Types, B_BOND_STRETCH, 2,
			     N_Inter_Types);
    return TRUE;
}


/*****************************************************************************************
read_BOND_angle_inter(): Reads the angle bending intramolecular interactions from the
parameter file. Stops execution if they have been read previously. Returns TRUE.
*****************************************************************************************/
int read_BOND_angle_inter(int b_Angles, tW_line inp_line, tW_system * sys,
			  FILE * fp_par)
{
    int N_Inter_Types;

    if (b_Angles == 1) {
	printf("ERROR: Reading Angles twice.\n");
	exit(EXIT_FAILURE);
    }

    N_Inter_Types = get_no_inter_types(B_ANGLE, inp_line);

    sys->Bonded_Inter_Types =
	get_Bond_Inter_Types(fp_par, inp_line, *sys,
			     &(sys->N_Bond_Int_Types),
			     sys->Bonded_Inter_Types, B_ANGLE, 3,
			     N_Inter_Types);

    return TRUE;
}


/*****************************************************************************************
read_Bond_dihedral_inter(): Reads the dihedral intramolecular interactions from the
parameter file. Stops execution if they have been read previously. Returns TRUE.
*****************************************************************************************/
int read_Bond_dihedral_inter(int b_Dihedrals, tW_line inp_line,
			     tW_system * sys, FILE * fp_par)
{
    int N_Inter_Types;

    if (b_Dihedrals == 1) {
	printf("ERROR: Reading Dihedrals twice.\n");
	exit(EXIT_FAILURE);
    }

    N_Inter_Types = get_no_inter_types(B_DIHEDRAL, inp_line);

    sys->Bonded_Inter_Types =
	get_Bond_Inter_Types(fp_par, inp_line, *sys,
			     &(sys->N_Bond_Int_Types),
			     sys->Bonded_Inter_Types, B_DIHEDRAL, 4,
			     N_Inter_Types);

    return TRUE;
}


/*****************************************************************************************
read_Bond_nb_pair_inter(): Reads the non-bonded intramolecular interactions from the
parameter file. Stops execution if they have been read previously. Returns TRUE.
*****************************************************************************************/
int read_Bond_nb_pair_inter(int b_Pair_bond, tW_line inp_line,
			    tW_system * sys, FILE * fp_par)
{
    int N_Inter_Types;

    if (b_Pair_bond == 1) {
	printf("ERROR: Reading InterMolec_Pairs twice.\n");
	exit(EXIT_FAILURE);
    }

    N_Inter_Types = get_no_inter_types(B_NB_PAIR_BOND, inp_line);

    sys->Bonded_Inter_Types =
	get_Bond_Inter_Types(fp_par, inp_line, *sys,
			     &(sys->N_Bond_Int_Types),
			     sys->Bonded_Inter_Types, B_NB_PAIR_BOND, 2,
			     N_Inter_Types);

    return TRUE;
}


/*****************************************************************************************
open_ref_potential_file(): Opens the file that contains the reference potential 
information. Execulation is stopped if the file name is not found or the file could
not be opened. Returns a pointer to this file.
*****************************************************************************************/
FILE *open_ref_potential_file(tW_line inp_line)
{
    int test_sscanf;
    FILE *fp_ref;
    tW_word fname;

    delete_directive_from_line(inp_line);

    test_sscanf = sscanf(inp_line, " %s ", fname);

    if (test_sscanf != 1) {
	printf("\nERROR: Did not read reference potential file name.\n");
	exit(EXIT_FAILURE);
    }


    if (!file_exists(fname)) {
	printf("\nERROR: Cannot open %s.\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_ref = fopen(fname, "r");



    return fp_ref;
}


/*****************************************************************************************
read_REF_nb_inter(): Reads information for the intermolecular interactions. Stops 
execution if this information has been read previously. Returns TRUE.
*****************************************************************************************/
int read_REF_nb_inter(int b_Pairs, int b_Forces, tW_line inp_line,
		      tW_ref_potential * ref_potential, FILE * fp_ref,
		      tW_word mode)
{
    if (b_Pairs == TRUE) {
	printf("ERROR: Reading Pair_Interactions twice (ref).\n");
	exit(EXIT_FAILURE);
    }

    test_read_ref_forces(b_Forces, NB_PAIR);

    ref_potential->N_nb_pair_inter = get_no_inter_types(NB_PAIR, inp_line);

    ref_potential->nb_pair_inter =
	get_ref_inter(NB_PAIR, 2, ref_potential->N_nb_pair_inter, fp_ref,
		      ref_potential->nb_pair_inter, ref_potential, mode);

    return TRUE;
}


/*****************************************************************************************
test_read_ref_forces(): Tests to see if the reference forces have been read in yet. 
Stops execution if they have not.
*****************************************************************************************/
void test_read_ref_forces(int b_Forces, tW_word keyword)
{
    if (b_Forces == FALSE) {
	printf
	    ("\nERROR: Must read files under [%s] before reading [%s].\n",
	     REF_FORCE, keyword);
	exit(EXIT_FAILURE);
    }
}



/*****************************************************************************************
read_REF_BOND_stretch_inter(): Reads information for the bond stretch interactions. Stops
execution if this information has been read previously. Returns TRUE.
*****************************************************************************************/
int read_REF_BOND_stretch_inter(int b_Bonds, int b_Forces,
				tW_line inp_line,
				tW_ref_potential * ref_potential,
				FILE * fp_ref, tW_word mode)
{
    if (b_Bonds == TRUE) {
	printf("ERROR: Reading Bonds twice (ref).\n");
	exit(EXIT_FAILURE);
    }

    test_read_ref_forces(b_Forces, B_BOND_STRETCH);

    ref_potential->N_bondstretch_inter =
	get_no_inter_types(B_BOND_STRETCH, inp_line);

    ref_potential->bondstretch_inter = get_ref_inter(B_BOND_STRETCH, 2,
						     ref_potential->
						     N_bondstretch_inter,
						     fp_ref,
						     ref_potential->
						     bondstretch_inter,
						     ref_potential, mode);

    return TRUE;
}


/*****************************************************************************************
read_REF_BOND_angle_inter(): Reads information for the angle bending interactions. Stops
execution if this information has been read previously. Returns TRUE.
*****************************************************************************************/
int read_REF_BOND_angle_inter(int b_Angles, int b_Forces, tW_line inp_line,
			      tW_ref_potential * ref_potential,
			      FILE * fp_ref, tW_word mode)
{
    if (b_Angles == TRUE) {
	printf("ERROR: Reading Angles twice (ref).\n");
	exit(EXIT_FAILURE);
    }

    test_read_ref_forces(b_Forces, B_ANGLE);

    ref_potential->N_angle_inter = get_no_inter_types(B_ANGLE, inp_line);

    ref_potential->angle_inter =
	get_ref_inter(B_ANGLE, 3, ref_potential->N_angle_inter, fp_ref,
		      ref_potential->angle_inter, ref_potential, mode);

    return TRUE;
}


/*****************************************************************************************
read_REF_Bond_dihedral_inter(): Reads information for the dihedral interactions. Stops
execution if this information has been read previously. Returns TRUE.
*****************************************************************************************/
int read_REF_Bond_dihedral_inter(int b_Dihedrals, int b_Forces,
				 tW_line inp_line,
				 tW_ref_potential * ref_potential,
				 FILE * fp_ref, tW_word mode)
{
    if (b_Dihedrals == TRUE) {
	printf("ERROR: Reading Dihedrals twice (ref).\n");
	exit(EXIT_FAILURE);
    }

    test_read_ref_forces(b_Forces, B_DIHEDRAL);

    ref_potential->N_dihedral_inter =
	get_no_inter_types(B_DIHEDRAL, inp_line);

    ref_potential->dihedral_inter = get_ref_inter(B_DIHEDRAL, 4,
						  ref_potential->
						  N_dihedral_inter, fp_ref,
						  ref_potential->
						  dihedral_inter,
						  ref_potential, mode);

    return TRUE;
}


/*****************************************************************************************
read_REF_IntraMolecPair_inter(): Reads information for the intramolecular pair inter-
actions. Stops execution if this information has been read previously. Returns TRUE.
*****************************************************************************************/
int read_REF_IntraMolecPair_inter(int b_IntraMolecPair, int b_Forces,
				  tW_line inp_line,
				  tW_ref_potential * ref_potential,
				  FILE * fp_ref, tW_word mode)
{
    if (b_IntraMolecPair == TRUE) {
	printf("ERROR: Reading intramolecular pairs twice (ref).\n");
	exit(EXIT_FAILURE);
    }

    test_read_ref_forces(b_Forces, B_NB_PAIR_BOND);

    ref_potential->N_IntraMolecPairs_inter =
	get_no_inter_types(B_NB_PAIR_BOND, inp_line);

    ref_potential->IntraMolecPairs_inter = get_ref_inter(B_NB_PAIR_BOND, 2,
							 ref_potential->
							 N_IntraMolecPairs_inter,
							 fp_ref,
							 ref_potential->
							 IntraMolecPairs_inter,
							 ref_potential,
							 mode);

    return TRUE;
}



/*****************************************************************************************
read_REF_interpol_type(): Reads the interpolation type used for the reference potential.
*****************************************************************************************/
int read_REF_interpol_type(int b_interpolation, tW_line inp_line,
			   tW_ref_potential * ref_potential)
{
    if (b_interpolation == 1) {
	printf("ERROR: Specified interpolation twice (ref).\n");
	exit(EXIT_FAILURE);
    }

    ref_potential->interpolation_index = get_interpolation_index(inp_line);

    return TRUE;
}


/*****************************************************************************************
read_REF_information(): Reads the information for a reference potential. Calls the 
relevant function for each keyword that is found. 
*****************************************************************************************/
void read_REF_information(FILE * fp_ref, tW_line inp_line,
			  tW_ref_potential * ref_potential,
			  Ref_Flags * flags, tW_word mode)
{
    if (strstr(inp_line, NB_PAIR) != NULL) {
	flags->b_Pairs =
	    read_REF_nb_inter(flags->b_Pairs, flags->b_Forces, inp_line,
			      ref_potential, fp_ref, mode);
    }

    else if (strstr(inp_line, B_BOND_STRETCH) != NULL) {
	flags->b_Bonds =
	    read_REF_BOND_stretch_inter(flags->b_Bonds, flags->b_Forces,
					inp_line, ref_potential, fp_ref,
					mode);
    }

    else if (strstr(inp_line, B_ANGLE) != NULL) {
	flags->b_Angles =
	    read_REF_BOND_angle_inter(flags->b_Angles, flags->b_Forces,
				      inp_line, ref_potential, fp_ref,
				      mode);
    }

    else if (strstr(inp_line, B_DIHEDRAL) != NULL) {
	flags->b_Dihedrals =
	    read_REF_Bond_dihedral_inter(flags->b_Dihedrals,
					 flags->b_Forces, inp_line,
					 ref_potential, fp_ref, mode);
    }

    else if (strstr(inp_line, B_NB_PAIR_BOND) != NULL) {
	flags->b_IntraMolecPair =
	    read_REF_IntraMolecPair_inter(flags->b_IntraMolecPair,
					  flags->b_Forces, inp_line,
					  ref_potential, fp_ref, mode);
    }

    else if (strstr(inp_line, "Interpolation") != NULL) {
	flags->b_interpolation =
	    read_REF_interpol_type(flags->b_interpolation, inp_line,
				   ref_potential);
    }

    else if (strstr(inp_line, "RefForces") != NULL) {
	flags->b_Forces =
	    read_REF_forces(flags->b_Forces, inp_line, ref_potential,
			    fp_ref);
    }


    else {
	printf
	    ("\nERROR: Did non understand line in the reference potential file.\n");
	printf("    line: '%s'\n", inp_line);
	exit(EXIT_FAILURE);
    }
}


/*****************************************************************************************
read_REF_forces(): Read in the reference forces listed in input.
*****************************************************************************************/
int read_REF_forces(int b_Forces, tW_line inp_line,
		    tW_ref_potential * ref_potential, FILE * fp)
{
    int i;
    int test_sscanf;
    tW_word inter_type;

    if (b_Forces == TRUE) {
	printf("\nERROR: Reading reference forces twice.\n");
	exit(EXIT_FAILURE);
    }

    ref_potential->N_files = get_no_inter_types(REF_FORCE, inp_line);
    ref_potential->fname =
	emalloc(ref_potential->N_files * sizeof(tW_word));
    ref_potential->forces =
	emalloc(ref_potential->N_files * sizeof(tW_ref_force));

    for (i = 0; i < ref_potential->N_files; i++) {
	get_next_line(fp, inp_line);
	test_sscanf =
	    sscanf(inp_line, " %s %s", ref_potential->fname[i],
		   inter_type);

	if (test_sscanf != 2) {
	    printf
		("\nERROR: Problem reading reference force file name and type.\n");
	    printf("****'%s'\n", inp_line);
	    exit(EXIT_FAILURE);
	}

	read_ref_potential(ref_potential->fname[i],
			   &(ref_potential->forces[i]), inter_type);
    }

    get_next_line(fp, inp_line);
    test_end_directive(inp_line, REF_FORCE, "ref.txt");

    return TRUE;
}


/*****************************************************************************************
get_parameters(): Reads a line from the parameter file passed from read_par(). This line
should start with a directive (i.e. a keyword enclosed by [...]). This function reads
that keyword and calls the appropriate function to read associated line(s) from the 
parameter file. If the line is not understood, (i.e. a known keyword is not found), 
execution is stopped. After a successful call to this function, the next line read in 
read_par() will be a line with a different directive. The MACROS for the keywords are 
defined in read_parameters.h and wnoid_types.h. 
*****************************************************************************************/
void get_parameters(FILE * fp_par, tW_line inp_line, tW_system * sys,
		    tW_files * files, tW_ref_potential * ref_potential,
		    Par_Flags * flags)
{
    /* [Temperature] */
    if (strstr(inp_line, KEY_TEMP) != NULL) {
	flags->b_Temperature =
	    read_temperature(flags->b_Temperature, inp_line, sys);
    }

    /* [Structures] */
    else if (strstr(inp_line, KEY_STRUCT) != NULL) {
	flags->b_Structures =
	    read_structures(flags->b_Structures, inp_line, files, fp_par);
    }

    /* [Site_Types] */
    else if (strstr(inp_line, KEY_SITES) != NULL) {
	flags->b_SiteTypes =
	    read_site_types(flags->b_SiteTypes, inp_line, sys, fp_par);
    }

    /* [Inter_Types] */
    else if (strstr(inp_line, INTER_TYPES) != NULL) {
	flags->b_Inter_Types =
	    read_Inter_Types(flags->b_Inter_Types, inp_line, sys, fp_par);
    }

    /* [Pair_Interaction] */
    else if (strstr(inp_line, NB_PAIR) != NULL) {
	flags->b_Pairs =
	    read_nb_pair_inter(flags->b_Pairs, inp_line, sys, fp_par);
    }

    /* [BondStretch] */
    else if (strstr(inp_line, B_BOND_STRETCH) != NULL) {
	flags->b_Bonds =
	    read_BOND_stretch_inter(flags->b_Bonds, inp_line, sys, fp_par);
    }

    /* [Angle] */
    else if (strstr(inp_line, B_ANGLE) != NULL) {
	flags->b_Angles =
	    read_BOND_angle_inter(flags->b_Angles, inp_line, sys, fp_par);
    }

    /* [Dihedral] */
    else if (strstr(inp_line, B_DIHEDRAL) != NULL) {
	flags->b_Dihedrals =
	    read_Bond_dihedral_inter(flags->b_Dihedrals, inp_line, sys,
				     fp_par);
    }

    /* [InterMolec_NB_Pair] */
    else if (strstr(inp_line, B_NB_PAIR_BOND) != NULL) {
	flags->b_Pair_bond =
	    read_Bond_nb_pair_inter(flags->b_Pair_bond, inp_line, sys,
				    fp_par);
    }

    /* [Reference_Potential] */
    else if (strstr(inp_line, KEY_REF_POT) != NULL) {
	sys->flag_ref_potential =
	    get_ref_potential(inp_line, ref_potential,
			      sys->flag_ref_potential, files->mode);
    }

    /* [Perturbation_Theory] *//* JFR - 08.10.11 */
    else if (strstr(inp_line, KEY_PT) != NULL) {
	sys->PT_var.flag_PT =
	    get_PT(fp_par, inp_line, sys->PT_var.flag_PT, sys);
    }

    /* [Eigen] *//* JFR - 08.11.11 */
    else if (strstr(inp_line, KEY_EIGEN) != NULL) {
	sys->Eigen_var.flag_Eigen =
	    get_Eigen(fp_par, inp_line, sys->Eigen_var.flag_Eigen, sys);
    }

    /* [SVD] *//* JFR - 04.06.12 */
    else if (strstr(inp_line, KEY_SVD) != NULL) {
	sys->SVD_var.flag_SVD =
	    get_SVD(fp_par, inp_line, sys->SVD_var.flag_SVD, sys);
    }

    /* [TPR] *//* JFR - 04.06.12 */
    else if (strstr(inp_line, KEY_TPR) != NULL) {
	sys->TPR_var.flag_TPR =
	    get_TPR(fp_par, inp_line, sys->TPR_var.flag_TPR, sys);
    }

    /* [TOP_EXCL] *//* JFR - 06.27.12 */
    else if (strstr(inp_line, KEY_TPR_EXCL) != NULL) {
	sys->TPR_EXCL_var.flag_TPR_excl =
	    get_TPR_EXCL(fp_par, inp_line, sys->TPR_EXCL_var.flag_TPR_excl,
			 sys);
    }

    /* [Metric_Tensor] *//* JFR - 04.06.12 */
    else if (strstr(inp_line, KEY_MT) != NULL) {
	sys->MT_var.flag_MT =
	    get_MT(fp_par, inp_line, sys->MT_var.flag_MT, sys);
    }

    /* [Mean_Force_Decomposition] *//* JFR - 04.06.12 */
    else if (strstr(inp_line, KEY_MFD) != NULL) {
	sys->MFD_var.flag_MFD =
	    get_MFD(fp_par, inp_line, sys->MFD_var.flag_MFD, sys);
    }

    /* [Calculation_Mode] *//* JFR - 04.06.12 */
    else if (strstr(inp_line, KEY_CalcMODE) != NULL) {
	sys->CalcMODE_var.flag_CalcMODE =
	    get_CalcMODE(fp_par, inp_line, sys->CalcMODE_var.flag_CalcMODE,
			 sys);
    }

    /* [Preconditioning] *//* JFR - 04.06.12 */
    else if (strstr(inp_line, KEY_PC) != NULL) {
	sys->PC_var.flag_PC =
	    get_PC(fp_par, inp_line, sys->PC_var.flag_PC, sys);
    }

    /* [Memory] *//* JFR - 04.06.12 */
    else if (strstr(inp_line, KEY_MEM) != NULL) {
	sys->MEM_var.flag_MEM =
	    get_MEM(fp_par, inp_line, sys->MEM_var.flag_MEM, sys);
    }

    /* [Solution_Method] *//* JFR - 04.16.12 */
    else if (strstr(inp_line, KEY_SOLN) != NULL) {
	sys->SOLN_var.flag_SOLN =
	    get_SOLN(fp_par, inp_line, sys->SOLN_var.flag_SOLN, sys);
    }

    /* [Frame_Weighting] *//* NJD - 03.10.15 */
    else if (strstr(inp_line, KEY_FRAMEWEIGHT) != NULL) {
	sys->FRAMEWEIGHT_var.flag_FRAMEWEIGHT =
	    get_FRAMEWEIGHT(fp_par, inp_line, sys->FRAMEWEIGHT_var.flag_FRAMEWEIGHT, sys);
    }

    /* [Calculate_Errors] *//* JFR - 04.16.12 */
    else if (strstr(inp_line, KEY_ERR) != NULL) {
	sys->ERR_var.flag_ERR =
	    get_ERR(fp_par, inp_line, sys->ERR_var.flag_ERR, sys);
    }

    /* [Reference_Options] *//* JFR - 07.16.12 */
    else if (strstr(inp_line, KEY_REF) != NULL) {
	sys->REF_var.N_fnm = files->N_struct;
	sys->REF_var.flag_REF =
	    get_REF(fp_par, inp_line, sys->REF_var.flag_REF, sys);
    }

    /* [TRIM] *//* JFR - 01.29.13 */
    else if (strstr(inp_line, KEY_TRIM) != NULL) {
	sys->TRIM_var.flag_TRIM =
	    get_TRIM(fp_par, inp_line, sys->TRIM_var.flag_TRIM, sys);
    }

    /* [CHISQD] *//* JFR - 01.29.13 */
    else if (strstr(inp_line, KEY_CHISQD) != NULL) {
	sys->CHISQD_var.flag_CHISQD = get_CHISQD(fp_par, inp_line, sys->CHISQD_var.flag_CHISQD, sys);
    }

    /* [Regularization] *//* JFR - 12.03.13 */
    else if (strstr(inp_line, KEY_REG) != NULL) {
	sys->REG_var.flag_REG =
	    get_REG(fp_par, inp_line, sys->REG_var.flag_REG, sys);
    }

    /* [Rescale_Forces] *//* JFR - 01.31.13 */
    else if (strstr(inp_line, KEY_RESCALE) != NULL) {
	sys->RESCALE_var.flag_RESCALE =
	    get_RESCALE(fp_par, inp_line, sys->RESCALE_var.flag_RESCALE,
			sys);
    }

    /* [Constrain_Dihedrals] *//* JFR - 01.29.13 */
    else if (strstr(inp_line, KEY_CONSTRAIN) != NULL) {
	sys->CONSTRAIN_var.flag_CONSTRAIN =
	    get_CONSTRAIN(fp_par, inp_line,
			  sys->CONSTRAIN_var.flag_CONSTRAIN, sys);
    }

    /* [Iterative_gYBG] *//* JFR - 12.03.13 */
    else if (strstr(inp_line, KEY_ITER) != NULL) {
	sys->ITER_var.flag_ITER =
	    get_ITER(fp_par, inp_line, sys->ITER_var.flag_ITER, sys);
    }

    /* [nrexcl] */
    else if (strstr(inp_line, KEY_NREXCL) != NULL) {
	flags->b_nrexcl = read_nrexcl(flags->b_nrexcl, inp_line, sys);
    }

    /* [Mode] */
    else if (strstr(inp_line, KEY_MODE) != NULL) {
	flags->b_mode = read_mode(flags->b_mode, inp_line, files);
    }

    /* [Skip_Triple_Loop] */
    else if (strstr(inp_line, KEY_SKIPTL) != NULL) {
        sys->SKIP_TRIPLE_LOOP = TRUE;
    }
    /* Stop execution if line is not understood. */
    else {
	printf
	    ("\nERROR: The following line from par.txt is not understood:\n%s\n",
	     inp_line);
	exit(EXIT_FAILURE);
    }

}


/*****************************************************************************************
read_Inter_Types(): Reads the interaction types listed in the parameter file. Stops
execution if interaction types has been previously read. Calls get_no_inter_types() to 
retrieve the number of interaction types listed. Calls, get_Inter_Type_List() to get
the list of interaction types. Returns TRUE.
*****************************************************************************************/
int read_Inter_Types(int b_Inter_Types, tW_line inp_line, tW_system * sys,
		     FILE * fp_par)
{

    if (b_Inter_Types == TRUE) {
	printf("\nERROR: Reading Inter_Types twice.\n");
	exit(EXIT_FAILURE);
    }

    sys->N_Inter_Types = get_no_inter_types(INTER_TYPES, inp_line);

    sys->Inter_Types =
	get_Inter_Type_List(fp_par, inp_line, *sys, sys->N_Inter_Types,
			    INTER_TYPES);

    return TRUE;
}


/*****************************************************************************************
get_Inter_Type_List(): Reads and stores the information for each interaction type listed
in the parameter file into the data Inter_Types. A pointer to this data structure is
returned. Ultimately, this structure is part of the sys data structure. Calls to 
get_Inter_Type() retrieves the information for a single interaction type. 
*****************************************************************************************/
tW_Inter_Types *get_Inter_Type_List(FILE * fp_par, tW_line inp_line,
				    tW_system sys, int N_Inter_Types,
				    tW_word keyword)
{
    int i_inter;
    tW_Inter_Types *Inter_Types;

    Inter_Types =
	(tW_Inter_Types *) emalloc(N_Inter_Types * sizeof(tW_Inter_Types));

    initialize_Inter_Types(N_Inter_Types, Inter_Types);

    for (i_inter = 0; i_inter < N_Inter_Types; i_inter++) {
	get_next_line(fp_par, inp_line);

	check_inp_line("par.txt", keyword, inp_line);

	get_Inter_Type(inp_line, &(Inter_Types[i_inter]));
    }

    get_next_line(fp_par, inp_line);

    test_end_directive(inp_line, keyword, "par.txt");

    return Inter_Types;

}


/*****************************************************************************************
get_Inter_Type(): Reads and stores the relevant information for a single interaction type
listed in the parameter file. The MACROS are defined in wnoid_types.h. Execution is 
stopped if an unsupported interaction type is listed.
*****************************************************************************************/
void get_Inter_Type(tW_line inp_line, tW_Inter_Types * Inter_Types)
{
    int test_sscanf;

    test_sscanf = sscanf(inp_line, "%s %s %s", Inter_Types->inter_name,
			 Inter_Types->inter_type, Inter_Types->basis);

    if (test_sscanf != 3) {
	printf
	    ("\nERROR: Not able to read interaction in parameter file.\n");
	printf("    line: '%s'\n", inp_line);
	exit(EXIT_FAILURE);
    }

    /* Delta basis. */
    if (strcmp(Inter_Types->basis, DELTA_NAME) == 0) {

	Inter_Types->i_basis =
	    get_delta_basis_param(inp_line, &(Inter_Types->dr),
				  &(Inter_Types->R_min),
				  &(Inter_Types->R_max),
				  &(Inter_Types->N_pts),
				  &(Inter_Types->N_coeff),
				  &(Inter_Types->n_smooth),
				  Inter_Types->inter_type);
    }

    /* START JFR */
    /* Linear basis. */
    else if (strcmp(Inter_Types->basis, LINEAR_NAME) == 0) {

	Inter_Types->i_basis =
	    get_linear_basis_param(inp_line, &(Inter_Types->dr),
				   &(Inter_Types->R_min),
				   &(Inter_Types->R_max),
				   &(Inter_Types->N_pts),
				   &(Inter_Types->N_coeff),
				   &(Inter_Types->n_smooth),
				   Inter_Types->inter_type);
    }
    /* END JFR */

    /* Bspline basis. */
    else if (strcmp(Inter_Types->basis, BSPLINE_NAME) == 0) {	/* JFR - 07.22.12 */

	Inter_Types->i_basis =
	    get_Bspline_basis_param(inp_line, &(Inter_Types->dr),
				    &(Inter_Types->R_min),
				    &(Inter_Types->R_max),
				    &(Inter_Types->N_pts),
				    &(Inter_Types->kspline),
				    &(Inter_Types->N_coeff),
				    &(Inter_Types->n_smooth),
				    Inter_Types->inter_type);
    }

    /* Harmonic basis. */
    else if (strcmp(Inter_Types->basis, HARMONIC_NAME) == 0) {
	Inter_Types->i_basis = HARMONIC_BASIS_INDEX;
	Inter_Types->N_coeff = 2;
    }

    /* Ryckaert-Bellemans basis for dihedrals. */
    else if (strcmp(Inter_Types->basis, RYCKAERT_BELLEMANS_NAME) == 0) {
	Inter_Types->i_basis = RYCKAERT_BELLEMANS_BASIS_INDEX;
	Inter_Types->N_coeff = 5;
    }

    /* Dihedral basis for the TOY PDB. */
    else if (strcmp(Inter_Types->basis, TOY_DIHED_NAME) == 0) {
	Inter_Types->i_basis = TOY_DIHED_INDEX;
	Inter_Types->N_coeff = 3;
    }

    /* Power basis for non-bonded inter or intra molecular pair interactions. */
    else if (strcmp(Inter_Types->basis, POWER_NAME) == 0) {
	Inter_Types->i_basis =
	    get_power_basis_param(inp_line, Inter_Types);
    }

    /* Stop execution if unknown basis is found. */
    else {
	printf("\nERROR: Unsupported basis type: %s\n",
	       Inter_Types->basis);
	printf("   inp_line: %s\n", inp_line);
	exit(EXIT_FAILURE);
    }

}


/*****************************************************************************************
initialize_Inter_Types(): Initializes the variables in the data structure Inter_Types.
*****************************************************************************************/
void initialize_Inter_Types(int N_Inter_Types,
			    tW_Inter_Types * Inter_Types)
{
    int i;

    for (i = 0; i < N_Inter_Types; i++) {
	Inter_Types[i].i_basis = -1;
	Inter_Types[i].dr = 0.0;
	Inter_Types[i].R_min = 0.0;
	Inter_Types[i].R_max = 0.0;
	Inter_Types[i].N_pts = 0;
	Inter_Types[i].N_coeff = 0;
	Inter_Types[i].n_smooth = 0;
	Inter_Types[i].i_0 = -1;
	Inter_Types[i].N_powers = 0;
	Inter_Types[i].kspline = 0;	/* JFR - 07.22.12: Bspline order */
	Inter_Types[i].powers = NULL;
	Inter_Types[i].ptr_x = NULL;
	Inter_Types[i].ptr_g = NULL;
	Inter_Types[i].ptr_g_cnt = NULL;
	Inter_Types[i].ptr_L = NULL;
	Inter_Types[i].ptr_b = NULL;
	Inter_Types[i].ptr_b_ref = NULL;
	Inter_Types[i].ptr_b_forces = NULL;
	Inter_Types[i].ptr_b_struct = NULL;
	Inter_Types[i].ptr_phi = NULL;
	Inter_Types[i].ptr_phi_forces = NULL;
	Inter_Types[i].ptr_phi_struct = NULL;
    }
}


/*****************************************************************************************
summarize_Inter_Types(): Summarizes the input stored for each interaction type.
*****************************************************************************************/
void summarize_Inter_Types(FILE * fp_log, tW_system sys)
{
    int i, j;
    tW_Inter_Types *tmp_Inter;

    /* List Interaction Types. */
    fprintf(fp_log, "  N_Inter_Types: %d", sys.N_Inter_Types);
    for (i = 0; i < sys.N_Inter_Types; i++) {
	fprintf(fp_log, "\n    Inter_Type: %d\n", i + 1);
	tmp_Inter = &sys.Inter_Types[i];
	fprintf(fp_log, "    Inter_Name: %s", tmp_Inter->inter_name);
	fprintf(fp_log, "    Inter_Type: %s\n", tmp_Inter->inter_type);
	fprintf(fp_log, "    Basis: %s    i_basis: %d\n", tmp_Inter->basis,
		tmp_Inter->i_basis);

	if (tmp_Inter->N_powers > 0) {
	    fprintf(fp_log, "    R_max: %f  N_powers: %d  Powers: ",
		    tmp_Inter->R_max, tmp_Inter->N_powers);
	    for (j = 0; j < tmp_Inter->N_powers; j++) {
		fprintf(fp_log, "%d", tmp_Inter->powers[j]);
		if (j < tmp_Inter->N_powers - 1) {
		    fprintf(fp_log, ", ");
		} else {
		    fprintf(fp_log, "\n");
		}
	    }
	} else {
	    fprintf(fp_log, "    dr: %f    R_min: %f    R_max: %f\n",
		    tmp_Inter->dr, tmp_Inter->R_min, tmp_Inter->R_max);
	}

	fprintf(fp_log,
		"    N_pts: %d    N_coeff: %d    n_smooth: %d    i_0: %d\n",
		tmp_Inter->N_pts, tmp_Inter->N_coeff, tmp_Inter->n_smooth,
		tmp_Inter->i_0);

	if (tmp_Inter->i_basis == BSPLINE_BASIS_INDEX) {
	    fprintf(fp_log, "    kspline: %d\n", tmp_Inter->kspline);
	}			/* JFR - 07.22.12: Bspline order */
    }

    fprintf(fp_log, "\n");

}


/*****************************************************************************************
get_interaction_i_0(): Finds the index to the grids for an interaction. 
*****************************************************************************************/
int get_interaction_i_0(tW_word inter_name, tW_system * sys)
{
    int i;
    int i_0 = -1;

    for (i = 0; i < sys->N_Inter_Types; i++) {
	if (strcmp(inter_name, sys->Inter_Types[i].inter_name) == 0) {
	    i_0 = sys->Inter_Types[i].i_0;
	    break;
	}
    }

    if (i_0 < 0) {
	printf("ERROR: Could not find i_0 for Inter_Type: %s.\n",
	       inter_name);
	exit(EXIT_FAILURE);
    }

    return i_0;
}


/*****************************************************************************************
check_site_list(): Makes sure that sites were listed under [Site_Types].
*****************************************************************************************/
void check_site_list(tW_word * trial_list, int N_inter_sites,
		     tW_system sys, tW_word keyword, tW_word input_file)
{
    int i, j;

    for (i = 0; i < N_inter_sites; i++) {
	for (j = 0; j < sys.N_Site_Types; j++) {
	    if (strcmp(trial_list[i], sys.Site_Types[j]) == 0) {
		break;
	    }
	    if (j == (sys.N_Site_Types - 1)) {
		printf
		    ("ERROR: Site Type \"%s\" listed under [%s] was not listed under[Site_Types].\n",
		     trial_list[i], keyword);
		printf("    Check %s.\n", input_file);
		exit(EXIT_FAILURE);
	    }
	}
    }

}


/*****************************************************************************************
check_ref_input(): Makes sure that the information for the reference potential is 
consisten with the rest of the input parameters.
*****************************************************************************************/
void check_ref_input(tW_ref_potential ref_potential, tW_system sys)
{
    /* Make sure site types listed for each reference force is listed under [Site_Types]. */
    check_sites_ref_potential(ref_potential.N_nb_pair_inter,
			      ref_potential.nb_pair_inter, sys);
    check_sites_ref_potential(ref_potential.N_bondstretch_inter,
			      ref_potential.bondstretch_inter, sys);
    check_sites_ref_potential(ref_potential.N_angle_inter,
			      ref_potential.angle_inter, sys);
    check_sites_ref_potential(ref_potential.N_dihedral_inter,
			      ref_potential.dihedral_inter, sys);

    /* At some point, it might be good to write a check to see if reference interactions are listed 
     * listed int par.txt under the relevant sections since I am relying on the same GROMACS topology.
     */

}


/*****************************************************************************************
check_sites_ref_potential(): Makes sure the interaction sites for the reference potential
are consisten with the site type list.
*****************************************************************************************/
void check_sites_ref_potential(int N_interactions, tW_ref_inter * inter,
			       tW_system sys)
{
    int i;

    /* Make sure site types are listed under [Site_Types]. */
    for (i = 0; i < N_interactions; i++) {
	check_site_list(inter[i].site_types, inter[i].N_inter_sites, sys,
			inter[i].name, "reference input");
    }

}


/*****************************************************************************************
initialize_par_flags(): Initializes the variables in the flags data structure.
*****************************************************************************************/
void initialize_par_flags(Par_Flags * flags)
{
    flags->b_Structures = FALSE;
    flags->b_SiteTypes = FALSE;
    flags->b_Temperature = FALSE;
    flags->b_Inter_Types = FALSE;
    flags->b_Pairs = FALSE;
    flags->b_Bonds = FALSE;
    flags->b_Angles = FALSE;
    flags->b_Dihedrals = FALSE;
    flags->b_Pair_bond = FALSE;
    flags->b_nrexcl = FALSE;
    flags->b_mode = FALSE;
}


/*****************************************************************************************
initialize_ref_flags(): Initializes the variables in the flags data structure.
*****************************************************************************************/
void initialize_ref_flags(Ref_Flags * flags)
{
    flags->b_Forces = FALSE;
    flags->b_Pairs = FALSE;
    flags->b_Bonds = FALSE;
    flags->b_Angles = FALSE;
    flags->b_Dihedrals = FALSE;
    flags->b_IntraMolecPair = FALSE;
    flags->b_interpolation = FALSE;
}

/*************************************************************************************************
read_save_state(): for resuming a calculation after the loop over particles.  JFR - added 04.11.12 
**************************************************************************************************/
int read_save_state(FILE * fp_log, tW_system * sys)
{
    int i;
    tW_word fname;
    FILE *fp_b_struct, *fp_b_forces, *fp_g, *fp_g_cnt, *fp_L, *fp_M,
	*fp_M2, *fp_M_cnt, *fp_d2b, *fp_d2M, *fp_FF, *fp_Nf, *fp_rescale;

    int l_inp, test_sscanf;
    tW_line inp_line;

    sprintf(fname, "%s", "save.b_struct.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_b_struct = fopen(fname, "r");


    sprintf(fname, "%s", "save.b_forces.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_b_forces = fopen(fname, "r");


    sprintf(fname, "%s", "save.g.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_g = fopen(fname, "r");

    sprintf(fname, "%s", "save.g_cnt.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_g_cnt = fopen(fname, "r");


    sprintf(fname, "%s", "save.L.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_L = fopen(fname, "r");


    sprintf(fname, "%s", "save.M.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_M = fopen(fname, "r");


    sprintf(fname, "%s", "save.M2.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_M2 = fopen(fname, "r");


    if (sys->MT_var.flag_Mcnt == TRUE) {
	sprintf(fname, "%s", "save.M_cnt.dat");
	if (!file_exists(fname)) {
	    printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	    exit(EXIT_FAILURE);
	}
	fp_M_cnt = fopen(fname, "r");

    }
    //if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 )
    //{
    sprintf(fname, "%s", "save.d2b.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_d2b = fopen(fname, "r");

    //}

    //if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 )
    //{
    sprintf(fname, "%s", "save.d2M.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_d2M = fopen(fname, "r");

    //}

    sprintf(fname, "%s", "save.FF.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_FF = fopen(fname, "r");


    sprintf(fname, "%s", "save.Nf.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_Nf = fopen(fname, "r");


    sprintf(fname, "%s", "save.rescale.dat");
    if (!file_exists(fname)) {
	printf("ERROR: Unable to find checkpoint file '%s'\n", fname);
	exit(EXIT_FAILURE);
    }
    fp_rescale = fopen(fname, "r");


    /* Allocate memory for the b's, since we won't be accessing get_b_X() */
    sys->b_forces = (double *) ecalloc(sys->N_coeff, sizeof(double));
    sys->b_struct = (double *) ecalloc(sys->N_coeff, sizeof(double));


    for (i = 0; i < sys->N_coeff; i++) {
	l_inp = get_next_line(fp_b_struct, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->b_struct[i]);

	l_inp = get_next_line(fp_b_forces, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->b_forces[i]);

	l_inp = get_next_line(fp_g, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->g[i]);

	l_inp = get_next_line(fp_g_cnt, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->g_cnt[i]);

	l_inp = get_next_line(fp_L, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->L[i]);

	//if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 )
	//{
	l_inp = get_next_line(fp_d2b, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->d2b[i]);
	//}

	l_inp = get_next_line(fp_rescale, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->rescale[i]);
    }

    for (i = 0; i < sys->N_pack; i++) {
	l_inp = get_next_line(fp_M, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->M[i]);

	l_inp = get_next_line(fp_M2, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->M2[i]);

	if (sys->MT_var.flag_Mcnt == TRUE) {
	    l_inp = get_next_line(fp_M_cnt, inp_line);
	    test_sscanf = sscanf(inp_line, "%lf", &sys->M_cnt[i]);
	}
	//if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 )
	//{
	l_inp = get_next_line(fp_d2M, inp_line);
	test_sscanf = sscanf(inp_line, "%lf", &sys->d2M[i]);
	//}
    }

    l_inp = get_next_line(fp_FF, inp_line);
    test_sscanf = sscanf(inp_line, "%lf", &sys->Chi2);

    l_inp = get_next_line(fp_Nf, inp_line);
    test_sscanf = sscanf(inp_line, "%d", &sys->REG_var.Nframes);

    fclose(fp_b_struct);
    fclose(fp_b_forces);
    fclose(fp_g);
    fclose(fp_g_cnt);
    fclose(fp_L);
    fclose(fp_M);
    fclose(fp_M2);
    if (sys->MT_var.flag_Mcnt == TRUE) {
	fclose(fp_M_cnt);
    }
    fclose(fp_d2b);
    fclose(fp_d2M);
    ///*if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) {*/ fclose(fp_d2b); /*}*/
    ///*if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 ) {*/ fclose(fp_d2M); /*}*/
    fclose(fp_FF);
    fclose(fp_Nf);
    fclose(fp_rescale);

    return 0;
}

/*************************************************************************************************
estimate_memory_usage(): estimates the memory that will be used for the calculation based on the 
number of coefficients and some other input options.  JFR - added 04.19.12
**************************************************************************************************/
void estimate_memory_usage(tW_files files, tW_system * sys)
{
    int i;
    int N_coeff = 0;
    double coeffN2, coeffN;
    long double mem;

    fprintf(files.fp_log, "Estimating Memory Usage. . . \n");

    /* Determine the total number of coefficients. */
    for (i = 0; i < sys->N_Inter_Types; i++) {
	N_coeff += sys->Inter_Types[i].N_coeff;
    }
    fprintf(files.fp_log, "There are %d CGFF parameters. \n", N_coeff);

    /* Case 1 - serial, first half */
    coeffN2 = 2.0;
    coeffN = 12.0;
    if ((files.N_struct == 1) || (sys->MEM_var.flag_mult_top == FALSE) || (sys->REF_var.flag_splitfiles == TRUE)) {	/* sys_top and sys_top_global not made */
	coeffN2 = 1.0;
	coeffN = 6.0;
	if (sys->MT_var.flag_Mcnt == FALSE) {
	    coeffN2 = 0.5;
	    coeffN = 5.5;
	}
    } else if (sys->MT_var.flag_Mcnt == FALSE) {
	coeffN2 = 1.0;
	coeffN = 6.0;
    }
    mem =
	((long double) (coeffN2 * ((double) (N_coeff * N_coeff)))) +
	((long double) coeffN * ((double) N_coeff));
    mem = mem / 134217728.0;
    fprintf(files.fp_log,
	    "For the first half of the calculation, you will need about: \n");
    fprintf(files.fp_log, "Serial -> %Lf GB \n", mem);

    /* Case 2 - parallel, first half */
    coeffN2 = 4.0;
    coeffN = 24.0;
    if ((files.N_struct == 1) || (sys->MEM_var.flag_mult_top == FALSE) || (sys->REF_var.flag_splitfiles == TRUE)) {	/* sys_top and sys_top_global not made */
	coeffN2 = 2.0;
	coeffN = 12.0;
    }
    if (sys->MT_var.flag_Mcnt == FALSE) {
	coeffN2 = 1.0;
	coeffN = 11.0;
    }
    mem =
	coeffN2 * N_coeff * ((double) N_coeff) +
	coeffN * ((double) N_coeff);
    mem = mem / 134217728.0;
    fprintf(files.fp_log, "Parallel -> %Lf GB \n", mem);

    /* Case 3 - serial, second half */
    coeffN2 = 1.0;
    coeffN = 16.0;
    if (sys->MEM_var.flag_LOWMEM == TRUE) {
	coeffN2 = 0.5;
	coeffN = 15.5;
    } else if (sys->SVD_var.flag_SVD == TRUE) {
	coeffN2 = 2.0;
	coeffN = 24.5;
    } else if (strcmp(sys->SOLN_var.SOLN_METH, LU) == 0) {
	coeffN2 = 2.0;
	coeffN = 15.0;
    }
    mem =
	((long double) (coeffN2 * ((double) (N_coeff * N_coeff)))) +
	((long double) coeffN * ((double) N_coeff));
    mem = mem / 134217728.0;
    fprintf(files.fp_log,
	    "For the second half of the calculation, you will need about: \n");
    fprintf(files.fp_log, "Serial -> %Lf GB \n", mem);

    /* Case 4 - parallel, second half */
    coeffN2 = 1.5;
    coeffN = 21.5;
    if (sys->MEM_var.flag_LOWMEM == TRUE) {
	coeffN2 = 1.0;
	coeffN = 21.0;
    } else if (sys->SVD_var.flag_SVD == TRUE) {
	coeffN2 = 3.0;
	coeffN = 29.5;
    } else if (strcmp(sys->SOLN_var.SOLN_METH, LU) == 0) {
	coeffN2 = 3.0;
	coeffN = 20.0;
    }
    mem =
	coeffN2 * N_coeff * ((double) N_coeff) +
	coeffN * ((double) N_coeff);
    mem = mem / 134217728.0;
    fprintf(files.fp_log, "Parallel -> %Lf GB \n", mem);

}
