/**
@file read_parameters.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef READ_PARAM
#define READ_PARAM

#include "cgff_types.h"


/* Names for reading input */
#define KEY_TEMP         "Temperature"
#define KEY_STRUCT       "Structures"
#define KEY_SITES        "Site_Types"
#define INTER_TYPES      "Inter_Types"
#define KEY_NREXCL       "nrexcl"
#define KEY_REF_POT      "Reference_Potential"
#define KEY_MODE         "Mode"
#define KEY_STRUCT_FILES "Struct_Files"
#define KEY_PT           "Iterative_Inversion"	//JFR - 08.10.11: Note that I used to call this "Perturbation_Theory" so all the
					       // corresponding functions are labeled with "PT"
#define KEY_EIGEN        "Eigen"	//JFR - 08.11.11
#define KEY_SVD          "SVD"	//JFR - 04.06.12
#define KEY_TPR          "TPR"	//JFR - 04.06.12
#define KEY_TPR_EXCL     "TPR_EXCL"	//JFR - 06.27.12
#define KEY_MT           "Metric_Tensor"	//JFR - 04.06.12
#define KEY_MFD          "Mean_Force_Decomposition"	//JFR - 04.06.12
#define KEY_CalcMODE     "Calculation_Mode"	//JFR - 04.06.12
#define KEY_PC           "Preconditioning"	//JFR - 04.06.12
#define KEY_MEM          "Memory"	//JFR - 04.13.12
#define KEY_SOLN         "Solution_Method"	//JFR - 04.16.12
#define KEY_FRAMEWEIGHT     "Frame_Weighting"	//NJD - 03.10.15
#define KEY_ERR          "Error_Estimates"	//JFR - 04.16.12
#define KEY_REF          "Reference_Options"	//JFR - 07.16.12
#define KEY_TRIM         "TRIM"	//JFR - 01.29.13
#define KEY_CHISQD       "CHISQD"	//JFR - 01.29.13
#define KEY_REG          "Regularization"	//JFR - 12.03.13
#define KEY_RESCALE      "Rescale_Forces"	//JFR - 01.31.13
#define KEY_CONSTRAIN    "Constrain_Dih"	//JFR - 01.31.13
#define KEY_ITER         "Iterative_gYBG"	//JFR - 12.03.13
#define KEY_SKIPTL       "Skip_Triple_Loop"     //MRD - 03.04.19

/* File name for summarizing input. */
#define SUM_INPUT_FNAME  "summarize_input.txt"


/****************************************************************************************/
/*Types defined just for read_parameters.c***********************************************/
/****************************************************************************************/

typedef struct {
    int b_Structures;
    int b_SiteTypes;
    int b_Temperature;
    int b_Inter_Types;
    int b_Pairs;
    int b_Bonds;
    int b_Angles;
    int b_Dihedrals;
    int b_Pair_bond;
    int b_nrexcl;
    int b_mode;
} Par_Flags;

typedef struct {
    int b_Forces;
    int b_Pairs;
    int b_Bonds;
    int b_Angles;
    int b_Dihedrals;
    int b_IntraMolecPair;
    int b_interpolation;
} Ref_Flags;

/****************************************************************************************/
/*Function called outside of read_parameters.c*******************************************/
/****************************************************************************************/

int read_par(tW_files * files, tW_system * sys,
	     tW_ref_potential * ref_potential);

int setup_tables(tW_system * sys);

void summarize_input(tW_files files, tW_system sys,
		     tW_ref_potential ref_potential);

void setup_sys_copy(tW_system sys_orig, tW_system * sys_copy, int overwrite_arrays);

int free_sys_copy(tW_system * sys_top);

/****************************************************************************************/
/*Functions associated with read_par()***************************************************/
/****************************************************************************************/

void initialize_sys(tW_system * sys);

void initialize_files(tW_files * files);

void initialize_par_flags(Par_Flags * flags);

void get_parameters(FILE * fp_par, tW_line inp_line, tW_system * sys,
		    tW_files * files, tW_ref_potential * ref_potential,
		    Par_Flags * flags);

int read_temperature(int b_Temperature, tW_line inp_line, tW_system * sys);

int read_structures(int b_Structures, tW_line inp_line, tW_files * files,
		    FILE * fp_par);

tW_word *get_structures(tW_files * files, int *N_struct);

int read_site_types(int b_SiteTypes, tW_line inp_line, tW_system * sys,
		    FILE * fp_par);

int read_Inter_Types(int b_Inter_Types, tW_line inp_line, tW_system * sys,
		     FILE * fp_par);

int get_no_inter_types(tW_word keyword, tW_line inp_line);

tW_Inter_Types *get_Inter_Type_List(FILE * fp_par, tW_line inp_line,
				    tW_system sys, int N_Inter_Types,
				    tW_word keyword);

void initialize_Inter_Types(int N_Inter_Types,
			    tW_Inter_Types * Inter_Types);

void get_Inter_Type(tW_line inp_line, tW_Inter_Types * Inter_Types);

int get_delta_basis_param(tW_line inp_line, double *dr, double *R_0,
			  double *R_max, int *N_pts, int *N_coeff,
			  int *n_smooth, tW_word inter_type);

/* START JFR */
int get_linear_basis_param(tW_line inp_line, double *dr, double *R_0,
			   double *R_max, int *N_pts, int *N_coeff,
			   int *n_smooth, tW_word inter_type);
/* END JFR */

/* JFR - 07.22.12 */
int get_Bspline_basis_param(tW_line inp_line, double *dr, double *R_0,
			    double *R_max, int *N_pts, int *kspline,
			    int *N_coeff, int *n_smooth,
			    tW_word inter_type);

int get_power_basis_param(tW_line inp_line, tW_Inter_Types * inter);

tW_type_inter2 *get_Type_Inter2_List(FILE * fp_par, tW_line inp_line,
				     tW_system sys, int N_Inter_Types,
				     tW_word keyword);

int read_Type_Inter2(tW_line inp_line, tW_type_inter2 * inter,
		     tW_system sys);

int read_nb_pair_inter(int b_Pairs, tW_line inp_line, tW_system * sys,
		       FILE * fp_par);

int read_BOND_stretch_inter(int b_Bonds, tW_line inp_line, tW_system * sys,
			    FILE * fp_par);

int read_BOND_angle_inter(int b_Angles, tW_line inp_line, tW_system * sys,
			  FILE * fp_par);

int read_Bond_dihedral_inter(int b_Dihedrals, tW_line inp_line,
			     tW_system * sys, FILE * fp_par);

int read_Bond_nb_pair_inter(int b_Pair_bond, tW_line inp_line,
			    tW_system * sys, FILE * fp_par);

tW_Bonded_Inter *get_Bond_Inter_Types(FILE * fp_par, tW_line inp_line,
				      tW_system sys, int *N_Bond_Int_Types,
				      tW_Bonded_Inter * Bonded_Inter_Types,
				      tW_word keyword, int N_Int_sites,
				      int N_Inter_Types);

void get_bond_inter_site_types(tW_line inp_line,
			       tW_Bonded_Inter * Bond_ptr, int N_Int_Sites,
			       tW_word inter_name);

void get_bond_basis_parameters(tW_word inter_name,
			       tW_Bonded_Inter * Bond_ptr, int N_Int_sites,
			       tW_line inp_line, tW_system sys);

int get_ref_potential(tW_line inp_line, tW_ref_potential * ref_potential,
		      int ref_flag, tW_word mode);

/* START JFR */
int get_PT(FILE * fp_par, tW_line inp_line, int PT_flag, tW_system * sys);

int get_Eigen(FILE * fp_par, tW_line inp_line, int Eigen_flag,
	      tW_system * sys);

int get_SVD(FILE * fp_par, tW_line inp_line, int SVD_flag,
	    tW_system * sys);

int get_TPR(FILE * fp_par, tW_line inp_line, int TPR_flag,
	    tW_system * sys);

int get_TPR_EXCL(FILE * fp_par, tW_line inp_line, int TPR_EXCL_flag,
		 tW_system * sys);

int get_MT(FILE * fp_par, tW_line inp_line, int MT_flag, tW_system * sys);

int get_MFD(FILE * fp_par, tW_line inp_line, int MFD_flag,
	    tW_system * sys);

int get_CalcMODE(FILE * fp_par, tW_line inp_line, int CalcMODE_flag,
		 tW_system * sys);

int get_PC(FILE * fp_par, tW_line inp_line, int PC_flag, tW_system * sys);

int get_MEM(FILE * fp_par, tW_line inp_line, int MEM_flag,
	    tW_system * sys);

int get_SOLN(FILE * fp_par, tW_line inp_line, int SOLN_flag,
	     tW_system * sys);

int get_FRAMEWEIGHT(FILE * fp_par, tW_line inp_line, int FRAMEWEIGHT_flag,
	     tW_system * sys);

int get_ERR(FILE * fp_par, tW_line inp_line, int ERR_flag,
	    tW_system * sys);

int get_REF(FILE * fp_par, tW_line inp_line, int REF_flag,
	    tW_system * sys);

int get_TRIM(FILE * fp_par, tW_line inp_line, int TRIM_flag,
	     tW_system * sys);

int get_CHISQD(FILE * fp_par, tW_line inp_line, int CHISQD_flag,
	       tW_system * sys);

int get_REG(FILE * fp_par, tW_line inp_line, int REG_flag,
	    tW_system * sys);

int get_RESCALE(FILE * fp_par, tW_line inp_line, int RESCALE_flag,
		tW_system * sys);

int get_CONSTRAIN(FILE * fp_par, tW_line inp_line, int CONSTRAIN_flag,
		  tW_system * sys);

int get_ITER(FILE * fp_par, tW_line inp_line, int ITER_flag,
	     tW_system * sys);
/* END JFR */

void initialize_ref_potential(tW_ref_potential * ref_potential);

void initialize_ref_flags(Ref_Flags * flags);

FILE *open_ref_potential_file(tW_line inp_line);

void read_REF_information(FILE * fp_ref, tW_line inp_line,
			  tW_ref_potential * ref_potential,
			  Ref_Flags * flags, tW_word mode);

int read_REF_nb_inter(int b_Pairs, int b_Forces, tW_line inp_line,
		      tW_ref_potential * ref_potential, FILE * fp_ref,
		      tW_word mode);

int read_REF_BOND_stretch_inter(int b_Bonds, int b_Forces,
				tW_line inp_line,
				tW_ref_potential * ref_potential,
				FILE * fp_ref, tW_word mode);

int read_REF_BOND_angle_inter(int b_Angles, int b_Forces, tW_line inp_line,
			      tW_ref_potential * ref_potential,
			      FILE * fp_ref, tW_word mode);

int read_REF_Bond_dihedral_inter(int b_Dihedrals, int b_Forces,
				 tW_line inp_line,
				 tW_ref_potential * ref_potential,
				 FILE * fp_ref, tW_word mode);

int read_REF_interpol_type(int b_interpolation, tW_line inp_line,
			   tW_ref_potential * ref_potential);

int get_interpolation_index(tW_line inp_line);

tW_ref_inter *get_ref_inter(tW_word inter_name, int N_Inter_Sites,
			    int N_Inter_Types, FILE * fp,
			    tW_ref_inter * inter_list,
			    tW_ref_potential * ref_potential,
			    tW_word mode);

void get_ref_potential_info(tW_ref_inter * inter_ptr, int N_Inter_Sites,
			    tW_word inter_name, tW_line inp_line,
			    tW_ref_potential * ref_potential,
			    tW_word mode);

void read_ref_potential(tW_word fname, tW_ref_force * inter_ptr,
			tW_word inter_type);

int get_no_points(FILE * fp);

int get_force(FILE * fp, tW_ref_force * inter_ptr, tW_word inter_type);

int read_nrexcl(int b_nrexcl, tW_line inp_line, tW_system * sys);

int read_mode(int b_mode, tW_line inp_line, tW_files * files);

void check_input(tW_files files, tW_system sys, Par_Flags flags,
		 tW_ref_potential ref_potential);



void check_struct_files(tW_files files);

void check_site_types(tW_system sys);

void check_Inter_Types(tW_system sys);

void check_Type_Inter2(tW_system sys);

void check_Bonded_Inter_Types(tW_system, tW_word inter_name,
			      int N_Int_Sites, tW_word mode);

void check_site_list(tW_word * trial_list, int N_inter_sites,
		     tW_system sys, tW_word keyword, tW_word input_file);

void check_inter_sites(tW_word list1[], tW_word list2[], int N_words,
		       tW_word inter_name);

void check_ref_input(tW_ref_potential ref_potential, tW_system sys);

void check_sites_ref_potential(int N_interactions, tW_ref_inter * inter,
			       tW_system sys);

int read_REF_forces(int b_Forces, tW_line inp_line,
		    tW_ref_potential * ref_potential, FILE * fp);

void get_ref_potential_ptr(tW_word fname, tW_ref_inter * inter_ptr,
			   tW_ref_potential * ref_potential);

void test_read_ref_forces(int b_Forces, tW_word keyword);

int read_REF_IntraMolecPair_inter(int b_IntraMolecPair, int b_Forces,
				  tW_line inp_line,
				  tW_ref_potential * ref_potential,
				  FILE * fp_ref, tW_word mode);


/****************************************************************************************/
/*Functions associated with summarize_input()********************************************/
/****************************************************************************************/

void summarize_Inter_Types(FILE * fp_log, tW_system sys);

void summarize_input_pair_inter(FILE * fp, tW_system sys);

void summarize_input_bond_inter(FILE * fp, tW_system sys,
				tW_word inter_name, int N_Int_Sites,
				tW_word mode);

void summarize_input_ref_potential(FILE * fp_log, tW_system sys,
				   tW_ref_potential ref_potential);

void summarize_input_PT(FILE * fp_log, tW_system sys);	//JFR - added 09.27.11

void summarize_input_Eigen(FILE * fp_log, tW_system sys);	//JFR - added 09.27.11

void summarize_input_SVD(FILE * fp_log, tW_system sys);	//JFR - added 04.06.12

void summarize_input_TPR(FILE * fp_log, tW_system sys);	//JFR - added 04.06.12

void summarize_input_TPR_EXCL(FILE * fp_log, tW_system sys);	//JFR - added 06.27.12

void summarize_input_MT(FILE * fp_log, tW_system sys);	//JFR - added 04.06.12

void summarize_input_MFD(FILE * fp_log, tW_system sys);	//JFR - added 04.06.12

void summarize_input_CalcMODE(FILE * fp_log, tW_system sys);	//JFR - added 04.06.12

void summarize_input_PC(FILE * fp_log, tW_system sys);	//JFR - added 04.06.12

void summarize_input_MEM(FILE * fp_log, tW_system sys);	//JFR - added 04.13.12

void summarize_input_SOLN(FILE * fp_log, tW_system sys);	//JFR - added 04.16.12

void summarize_input_FRAMEWEIGHT(FILE * fp_log, tW_system sys);	//NJD - added 03.10.15

void summarize_input_ERR(FILE * fp_log, tW_system sys);	//JFR - added 04.16.12

void summarize_input_REF(FILE * fp_log, tW_system sys);	//JFR - added 07.16.12

void summarize_input_TRIM(FILE * fp_log, tW_system sys);	//JFR - added 01.29.13

void summarize_input_CHISQD(FILE * fp_log, tW_system sys);	//JFR - added 01.29.13

void summarize_input_REG(FILE * fp_log, tW_system sys);	//JFR - added 12.03.13

void summarize_input_RESCALE(FILE * fp_log, tW_system sys);	//JFR - added 01.31.13

void summarize_input_CONSTRAIN(FILE * fp_log, tW_system sys);	//JFR - added 01.31.13

void summarize_input_ITER(FILE * fp_log, tW_system sys);	//JFR - added 12.03.13

void print_interpolation_type(int interpolation_index, FILE * fp_log);

void print_ref_inter(tW_ref_inter * ref_inter, FILE * fp_log, int *ctr,
		     int N_inter);

/****************************************************************************************/
/*Functions associated with setup_sys_copy()*********************************************/
/****************************************************************************************/

void copy_site_types(tW_system sys_orig, tW_system * sys_copy);

void copy_arrays(tW_system sys_orig, tW_system * sys_copy);

void copy_Inter_Types(tW_system sys_orig, tW_system * sys_copy);

void copy_Inter2_Types(tW_system sys_orig, tW_system * sys_copy);

void copy_Bonded_Inter_Types(tW_system sys_orig, tW_system * sys_copy);

void copy_PT(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 09.27.11

void copy_Eigen(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 09.27.11

void copy_SVD(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.06.12

void copy_TPR(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.06.12

void copy_TPR_EXCL(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 06.27.12

void copy_MT(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.06.12

void copy_MFD(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.06.12

void copy_CalcMODE(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.06.12

void copy_PC(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.06.12

void copy_MEM(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.13.12

void copy_SOLN(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.16.12

void copy_FRAMEWEIGHT(tW_system sys_orig, tW_system * sys_copy);	//NJD - added 03.10.15

void copy_ERR(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 04.16.12

void copy_REF(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 07.16.12

void copy_TRIM(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 01.29.13

void copy_CHISQD(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 01.29.13

void copy_REG(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 12.03.13

void copy_RESCALE(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 01.31.13

void copy_CONSTRAIN(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 01.31.13

void copy_ITER(tW_system sys_orig, tW_system * sys_copy);	//JFR - added 01.31.13

/****************************************************************************************/
/*Functions associated with setup_tables()***********************************************/
/****************************************************************************************/

int get_interaction_i_0(tW_word inter_name, tW_system * sys);

/****************************************************************************************/
/*Functions called by more than one fuction in this file*********************************/
/****************************************************************************************/

void clear_sys_arrays(tW_system * sys);

/****************************************************************************************/
/****************************************************************************************/
int read_save_state(FILE * fp_log, tW_system * sys);

void estimate_memory_usage(tW_files files, tW_system * sys);

#endif
