/**
@file ref_potential.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef REF_POT
#define REF_POT

/****************************************************************************************/
/*Functions called outside of calc_ref_potential.c***************************************/
/****************************************************************************************/

void get_ref_forces(FILE * fp, int N_sites, tW_CG_site * CG_struct,
		    tW_gmx_info info, tW_gmx_topology *top,
		    tW_ref_potential ref_potential);

void print_bref(int N_coeff, double *bref, tW_word tag);	/* JFR - 07.16.12: for printing bref */

/****************************************************************************************/
/*Functions associated with get_ref_forces()*********************************************/
/****************************************************************************************/

void initialize_ref_forces(int N_sites, tW_CG_site * CG_struct);

void calc_ref_nb_pair_forces(int N_sites, tW_CG_site * CG_struct,
			     tW_gmx_info info, tW_ref_potential ref_potential, FILE * fp);

void calc_ref_BondStretch_forces(tW_gmx_topology *top, tW_CG_site * CG_struct,
				 tW_ref_potential ref_potential,
				 FILE * fp);

void calc_ref_Angle_forces(tW_gmx_topology *top, tW_CG_site * CG_struct,
			   tW_ref_potential ref_potential, FILE * fp);

void calc_ref_Dihedral_forces(tW_gmx_topology *top,
			      tW_CG_site * CG_struct,
			      tW_ref_potential ref_potential, FILE * fp);

void calc_ref_IntraMolecPairs_forces(tW_gmx_topology *top,
				     tW_CG_site * CG_struct,
				     tW_ref_potential ref_potential,
				     FILE * fp);


void print_CG_ref_f(int N_sites, tW_CG_site * CG_struct, FILE * fp);

/****************************************************************************************/


/****************************************************************************************/
/*Misc. functions that could probably be moved or changed********************************/
/****************************************************************************************/

double get_difference_vector_PBC(bool b_PBC, dvec box, dvec x_i, dvec x_j,
				 dvec x_ij);

void copy_2_site_names(tW_word * site_names, tW_word name1, tW_word name2);

void copy_3_site_names(tW_word * site_names, tW_word name1, tW_word name2,
		       tW_word name3);

void copy_4_site_names(tW_word * site_names, tW_word name1, tW_word name2,
		       tW_word name3, tW_word name4);

double get_scalar_f_from_table(double x, tW_ref_inter * inter,
			       int interpolation_index);

double get_no_interpolation_f(tW_ref_inter * inter, int coeff, double x);

double get_linear_interpolation_f(tW_ref_inter * inter, int coeff,
				  double x);

void get_central_pair_vector_forces(dvec x_ij, double r_ij, double f_r,
				    tW_CG_site * i_site,
				    tW_CG_site * j_site);

tW_ref_inter *get_inter_ptr(tW_word * site_names, tW_ref_inter * inter,
			    int N_inter);

tW_ref_inter *get_angle_inter_ptr(tW_word * site_names,
				  tW_ref_inter * inter, int N_inter);

tW_ref_inter *get_dihedral_inter_ptr(tW_word * site_names,
				     tW_ref_inter * inter, int N_inter);

tW_ref_inter *get_IntraMolecPair_inter_ptr(tW_word * site_names,
					   tW_ref_inter * inter,
					   int N_inter);

double get_ref_Angle_info(dvec x_1, dvec x_2, dvec x_3, dvec L_2, dvec L_3,
			  FILE * fp);

void get_angle_vector_forces(double f_theta, dvec L_2, dvec L_3,
			     double theta, tW_CG_site * site1,
			     tW_CG_site * site2, tW_CG_site * site3);

double get_ref_Dihedral_info(tW_dihedral * dihedral,
			     tW_CG_site * CG_struct, dvec * B);

void get_dihedral_vector_forces(double f_phi, double phi, dvec * B,
				tW_CG_site * site1, tW_CG_site * site2,
				tW_CG_site * site3, tW_CG_site * site4);

/****************************************************************************************/

#endif
