/**
@file calc_grids.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef CALC_GRIDS
#define CALC_CRIDS


#define Max_Num_Prot_Inter1 100

//#include "gmx-interface.h"
#include "gromacs_topology.h"

/****************************************************************************************/
/*Functions called outside of calc_grids-gmx.c*******************************************/
/****************************************************************************************/

//int calc_grids( FILE *fp, tW_gmx_info info, int N_sites, tW_CG_site *CG_struct, 
//                tW_system *sys );

int calc_grids2(FILE * fp, tW_gmx_info info, int N_sites,
		tW_CG_site * CG_struct, tW_system * sys);

int normalize_arrays_top(int N_sites, int N_frames, tW_system * sys_top);

int update_N_instances_nb(tW_system sys_top_global,
			  tW_system * sys_global);

int update_total_arrays(tW_system sys_top, double pr_top, tW_system * sys);

int get_results(FILE * fp_log, tW_system * sys, int forces, tW_word tag);

int weight_local_top(tW_system * sys_top, double w_local, int N_coeff);

/****************************************************************************************/
/*Functions associated with calc_grids().************************************************/
/****************************************************************************************/

bool setup_box_dimensions(dvec box, double matrix[DIM][DIM]);

int det_min_image(dvec box, dvec Origin, dvec x);

void subtract_ref_forces(int N_sites, tW_CG_site * CG_struct);

int eval_bond_basis_vectors(FILE * fp, tW_CG_site * CG_struct,
			    tW_Bonded_Inter * Bonded_Inter_Types,
			    int i_site, int *D_b, dvec * calG_b,
			    bool i_flag,
			    int *D_b_inter_index
			    /* JFR - 02.26.13: for sep M2 */ );

int get_grid_index_for_delta_basis(double r, int i_0, double dr,
				   double R_0);

/* START JFR */
int get_grid_index_for_linear_basis(double r, int i_0, double dr,
				    double R_0);
/* END JFR */

int eval_bonds_M_b(tW_system * sys, int nr_bond_coeff, int D_b[],
		   dvec calG_b[], bool b_F, dvec f_i, dvec f_i_ref,
		   int *D_b_inter_index /* JFR - 02.26.13: for sep M2 */ );

int get_iList(tW_word name, int N_Inter_Types, tW_type_inter2 InterList[],
	      int *i1_list, tW_word * n1_list);

void get_all_iList(tW_system *sys);

bool skip_excl(int nr_excl, int *excl_list, int j);

tW_type_inter2 *get_nb_pair_inter_ptr(tW_CG_site * j_site, tW_system * sys,
				      int N_i1, tW_word *nList1,
				      int *iList1);

void get_nb_pair_info(double *r_ij, dvec x_i, dvec x_j, dvec x_ij,
		      int b_PBC, dvec box);

int check_inter_range(double r_ij, double R_0, double R_max, double dr);

/* START JFR */
int check_linear_inter_range(double r_ij, double R_0, double R_max,
			     double dr);
/* END JFR */

/* JFR - 07.22.12 */
int check_Bspline_inter_range(double r_ij, double R_0, double R_max,
			      double dr);

/* JFR - 07.22.12 */
int eval_delta_basis_vectors(tW_type_inter2 * ij_inter, double n_basis_ij,
			     dvec u_ij, dvec * ij_basis, double r_ij,
			     int *ij_index, tW_system * sys,
			     tW_CG_site site_i, bool b_F);

/* START JFR */
int eval_linear_basis_vectors(tW_type_inter2 * ij_inter, double n_basis_ij,
			      dvec u_ij, dvec * ij_basis, double r_ij,
			      int *ij_index, tW_system * sys,
			      tW_CG_site site_i, bool b_F, int flag_grids);
/* END JFR */

/* JFR - 07.22.12 */
int eval_Bspline_basis_vectors(tW_type_inter2 * ij_inter, dvec u_ij,
			       dvec * ij_basis, double r_ij, int *ij_index,
			       tW_system * sys, tW_CG_site site_i,
			       bool b_F, int flag_grids);

/* START JFR */
double calc_linear_spline_A(int bond_coeff_ctr, int *D_b, double dr,
			    double R_0, int i_0, int N_pts, double R,
			    double *Ap, double *Bp, int periodic);
/* END JFR */

/* JFR - 07.16.12: Bspline basis set */
double calc_Bspline(double dr, double R_0, int i_0, int N_pts, double R,
		    int i, int k, double *Nik, double *Npik);

double calc_Bspline_deriv(double dr, double R_0, int i_0, int N_pts,
			  double R, int i, int j, int k, double *Nik,
			  double *Npik);

int eval_power_basis_vectors(tW_type_inter2 * ij_inter, dvec u_ij,
			     dvec * power_ij_basis, double r_ij,
			     int *ij_index, tW_system * sys,
			     tW_CG_site site_i, bool b_F, int flag_grids);

void update_nb_pair_grids(int n_basis_ij, int ij_index[], dvec ij_basis[],
			  bool b_F, double r_ij, tW_system * sys,
			  tW_CG_site site_i, tW_type_inter2 * ij_inter);

int eval_M_nbPair_bonds(tW_system * sys, int n_basis_ij, int ij_index[],
			dvec ij_basis[], int nr_i_bond_coeff, int D_b[],
			dvec calG_b[]);

void eval_M_nb_pair_inter(int n_basis_ij, int ij_index[], dvec ij_basis[],
			  int n_basis_ik, int ik_index[], dvec ik_basis[],
			  tW_system * sys, int flag_M2);

/****************************************************************************************/
/*Functions associated with eval_bond_basis_vectors().***********************************/
/****************************************************************************************/

int get_BondStretch_basis_vectors(int i_site, int i_inter,
				  int bond_coeff_ctr,
				  tW_CG_site * i_site_ptr,
				  tW_Bonded_Inter * Bonded_Inter_Types,
				  int *D_b, dvec * calG_b,
				  tW_CG_site * CG_struct, bool i_flag);

double get_BondStretch_info(int i_site, int *sites, tW_CG_site * CG_struct,
			    dvec u_ij);

int get_harmonic_bond_basis_vector(int bond_coeff_ctr, int i_bond_type,
				   int *D_b, dvec * calG_b,
				   tW_Bonded_Inter * Bonded_Inter_Types,
				   double norm_ij, dvec u_ij, bool i_flag);

int get_delta_bond_basis_vector(int bond_coeff_ctr, int i_bond_type,
				int *D_b, dvec * calG_b,
				tW_Bonded_Inter * Bonded_Inter_Types,
				double norm_ij, dvec u_ij, bool i_flag);

/* START JFR */
int get_linear_bond_basis_vector(int bond_coeff_ctr, int i_bond_type,
				 int *D_b, dvec * calG_b,
				 tW_Bonded_Inter * Bonded_Inter_Types,
				 double norm_ij, dvec u_ij, bool i_flag);
/* END JFR */

/* JFR - 07.23.12: Bspline */
int get_Bspline_bond_basis_vector(int bond_coeff_ctr, int i_bond_type,
				  int *D_b, dvec * calG_b,
				  tW_Bonded_Inter * Bonded_Inter_Types,
				  double norm_ij, dvec u_ij, bool i_flag);

/****************************************************************************************/

int get_Angle_basis_vectors(int i_site, int i_inter, int bond_coeff_ctr,
			    tW_CG_site * i_site_ptr,
			    tW_Bonded_Inter * Bonded_Inter_Types, int *D_b,
			    dvec * calG_b, tW_CG_site * CG_struct,
			    bool i_flag);

int get_Angle_info(int i_site, int *sites, tW_CG_site * CG_struct,
		   dvec L_2, dvec L_3, double *cos_theta,
		   double *sin_theta, double *theta, double *norm_12,
		   double *norm_13);

int get_atom_position_angle(int angle_triple[], int i_site);

int get_L2_L3_angles(dvec u_12, dvec u_13, double cos_theta,
		     double norm_12, double norm_13, dvec L_2, dvec L_3);

double get_g_theta_angle(double cos_theta, double r_12, double r_13);

double get_laplacian_theta_angle(double cos_theta, double sin_theta,
				 double r_12, double r_13);

int get_harmonic_angle_basis_vector(int i_site, int bond_coeff_ctr,
				    int i_bond_type, int *sites, int *D_b,
				    dvec * calG_b,
				    tW_Bonded_Inter * Bonded_Inter_Types,
				    tW_CG_site * CG_struct, bool i_flag);

int get_delta_angle_basis_vector(int i_site, int bond_coeff_ctr,
				 int i_bond_type, int *sites, int *D_b,
				 dvec * calG_b,
				 tW_Bonded_Inter * Bonded_Inter_Types,
				 tW_CG_site * CG_struct, bool i_flag);

/* START JFR */
int get_linear_angle_basis_vector(int i_site, int bond_coeff_ctr,
				  int i_bond_type, int *sites, int *D_b,
				  dvec * calG_b,
				  tW_Bonded_Inter * Bonded_Inter_Types,
				  tW_CG_site * CG_struct, bool i_flag);
/* END JFR */

/* JFR - 07.23.12: Bspline */
int get_Bspline_angle_basis_vector(int i_site, int bond_coeff_ctr,
				   int i_bond_type, int *sites, int *D_b,
				   dvec * calG_b,
				   tW_Bonded_Inter * Bonded_Inter_Types,
				   tW_CG_site * CG_struct, bool i_flag);

/****************************************************************************************/

int get_Dihedral_basis_vectors(int i_site, int i_inter, int bond_coeff_ctr,
			       tW_CG_site * i_site_ptr,
			       tW_Bonded_Inter * Bonded_Inter_Types,
			       int *D_b, dvec * calG_b,
			       tW_CG_site * CG_struct, bool i_flag);

int get_atom_position_dihedral(int dihedral_quartet[], int i_site);

void eval_grad_B(tW_dihedral * dihedral, tW_CG_site * CG_struct);

int get_RB_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				 int *D_b, dvec * calG_b,
				 tW_Bonded_Inter * Bonded_Inter_Types,
				 tW_dihedral dihedral, bool i_flag);

int get_delta_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				    int *D_b, dvec * calG_b,
				    tW_Bonded_Inter * Bonded_Inter_Types,
				    tW_dihedral dihedral, bool i_flag);

/* START JFR */
int get_linear_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				     int *D_b, dvec * calG_b,
				     tW_Bonded_Inter * Bonded_Inter_Types,
				     tW_dihedral dihedral, bool i_flag);
/* END JFR */

/* JFR - 07.23.12: Bspline */
int get_Bspline_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				      int *D_b, dvec * calG_b,
				      tW_Bonded_Inter * Bonded_Inter_Types,
				      tW_dihedral dihedral, bool i_flag);

int get_TOY_dihedral_basis_vector(int bond_coeff_ctr, int i_bond_type,
				  int *D_b, dvec * calG_b,
				  tW_Bonded_Inter * Bonded_Inter_Types,
				  tW_dihedral dihedral, bool i_flag);

double get_g_dihedral_angle(dvec basis);

/****************************************************************************************/

int get_IntraMolec_NB_Pair_basis_vectors(int i_site, int i_inter,
					 int bond_coeff_ctr,
					 tW_CG_site * i_site_ptr,
					 tW_Bonded_Inter *
					 Bonded_Inter_Types, int *D_b,
					 dvec * calG_b,
					 tW_CG_site * CG_struct,
					 bool i_flag);

double get_IntraMolec_NB_Pair_info(int i_site, int *sites,
				   tW_CG_site * CG_struct, dvec u_ij);

int get_delta_IntraMolec_NB_Pair_basis_vector(int bond_coeff_ctr,
					      int i_bond_type, int *D_b,
					      dvec * calG_b,
					      tW_Bonded_Inter *
					      Bonded_Inter_Types,
					      double norm_ij, dvec u_ij,
					      bool i_flag);

/* START JFR */
int get_linear_IntraMolec_NB_Pair_basis_vector(int bond_coeff_ctr,
					       int i_bond_type, int *D_b,
					       dvec * calG_b,
					       tW_Bonded_Inter *
					       Bonded_Inter_Types,
					       double norm_ij, dvec u_ij,
					       bool i_flag);
/* END JFR */

/* JFR - 07.23.12: Bspline */
int get_Bspline_IntraMolec_NB_Pair_basis_vector(int bond_coeff_ctr,
						int i_bond_type, int *D_b,
						dvec * calG_b,
						tW_Bonded_Inter *
						Bonded_Inter_Types,
						double norm_ij, dvec u_ij,
						bool i_flag);

void print_debug_eval_bond_basis(tW_CG_site * i_site_ptr,
				 tW_Bonded_Inter * Bonded_Inter_Types,
				 int i, FILE * fp);

/****************************************************************************************/
/*Functions associated with get_results()************************************************/
/****************************************************************************************/

int get_b_struct_and_b_forces(FILE * fp_log, tW_system * sys, tW_word tag);

int trim_Mb(FILE * fp, tW_system * sys);

int get_phi_struct_and_phi_forces(FILE * fp_log, tW_system * sys,
				  int forces);

int get_b_soln(FILE * fp_log, tW_system * sys);	//JFR

int get_b_soln_err(FILE * fp_log, tW_system * sys);	//JFR

int get_b_soln_errAA(FILE * fp_log, tW_system * sys);	//JFR

int calc_d2b(tW_system * sys);	/* JFR - 01.31.13 */

int calc_d2M(tW_system * sys);	/* JFR - 01.31.13 */

int get_rescale(tW_system * sys);	/* JFR - 01.31.13 */

/****************************************************************************************/
/*Functions associated with get_b_struct_and_b_forces()**********************************/
/****************************************************************************************/

int calc_b_struct(tW_system * sys, double temperature, tW_word tag);

void calc_dg_dz(tW_Inter_Types * inter, double *dg_dz, double *sm_dg_dz,
		tW_word tag);

int get_b_harmonic_bond(tW_Inter_Types * bond_inter, double kT);

int get_b_rb_dihedral(tW_Inter_Types * bond_inter, double kT);

void get_b_TOY_dihedral(tW_Inter_Types * inter, double kT);

void get_b_power(tW_Inter_Types * inter, double kT);

/****************************************************************************************/
/*Functions associated with trim_Mb()****************************************************/
/****************************************************************************************/

void initialize_Zero_list(int N_coeff0, int *Zero_list);

int get_Zero_list(int *Zero_list, int *N_coeff, double *g_cnt,
		  /*JFR*/ tW_system * sys);

//void setup_arrays_trim_Mb(tW_system * sys, int N_coeff);

//int remove_zero_rows(FILE * fp, int N_coeff0, int N_zero, int *Zero_list,
//		     tW_system * sys, tW_system * sys_temp);

int remove_zero_rows_mem(FILE * fp, int N_coeff0, int N_zero,
			 int *Zero_list, tW_system * sys);

// JFR - added 04.12.12:  resize the arrays instead of making a copy
int resize_arrays_trim_Mb(tW_system * sys, int N_coeff);

//void free_arrays_trim_Mb(tW_system * sys, int N_coeff);

//int copy_arrays_trim_Mb(tW_system * sys, tW_system sys_temp, int N_coeff);

int set_pointers_trim_Mb(tW_system * sys);

/****************************************************************************************/
/*Functions associated with get_phi_struct_and_phi_forces()******************************/
/****************************************************************************************/

int get_phi_forces(FILE * fp_log, tW_system * sys);

int get_phi_struct(FILE * fp_log, tW_system * sys);

int direct_solve_lin_eqns(FILE * fp_log, tW_system * sys, tW_word info);

/****************************************************************************************/

void eval_grad_B2(tW_dihedral * dihedral, tW_CG_site * CG_struct);

#endif
