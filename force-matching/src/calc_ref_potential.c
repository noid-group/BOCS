/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file ref_potential.c 
@authors Joseph Rudzinski, Nicholas Dunn
@brief Functions related to the use of the reference potential functionality. 
*/

//c libary includes
#include <stdio.h>
#include <string.h>

//local includes
#include "safe_mem.h"
#include "gmx-interface.h"
#include "calc_grids.h"
#include "wnoid_math.h"
#include "calc_ref_potential.h"
#include "read_parameters.h"
#include "io_read.h"

#define DEBUG_ref_forces   FALSE
#define DEBUG_ref_nb       FALSE
#define DEBUG_ref_bs       FALSE
#define DEBUG_ref_ab       FALSE
#define DEBUG_ref_dh       FALSE
#define DEBUG_ref_ip       FALSE


/*****************************************************************************************
get_ref_forces(): Calculate the forces due to the reference potential.
*****************************************************************************************/
void get_ref_forces(FILE * fp, int N_sites, tW_CG_site * CG_struct,
		    tW_gmx_info info, tW_gmx_topology *top,
		    tW_ref_potential ref_potential)
{
    if (DEBUG_ref_forces) {
	fprintf(fp, "In get_ref_forces().\n");
    }

    initialize_ref_forces(N_sites, CG_struct);

    calc_ref_nb_pair_forces(N_sites, CG_struct, info, ref_potential,
			    fp);

    calc_ref_BondStretch_forces(top, CG_struct,
				ref_potential, fp);

    calc_ref_Angle_forces(top, CG_struct, ref_potential,
			  fp);

    calc_ref_Dihedral_forces(top, CG_struct,
			     ref_potential, fp);

    calc_ref_IntraMolecPairs_forces(top, CG_struct,
				    ref_potential, fp);

    if (DEBUG_ref_forces) {
	print_CG_ref_f(N_sites, CG_struct, fp);
    }

}


/*****************************************************************************************
initialize_ref_forces(): Clear out the reference force variables.
*****************************************************************************************/
void initialize_ref_forces(int N_sites, tW_CG_site * CG_struct)
{
    int i, j;

    for (i = 0; i < N_sites; i++) {
	for (j = 0; j < 3; j++) {
	    CG_struct[i].ref_f[j] = 0.0;
	}
    }
}


/*****************************************************************************************
calc_ref_nb_pair_forces(): Evaluate the pair forces over all pairs due to the reference
potential.
*****************************************************************************************/
void calc_ref_nb_pair_forces(int N_sites, 
			     tW_CG_site * CG_struct,
			     tW_gmx_info info, 
			     tW_ref_potential ref_potential, 
			     FILE * fp)
{
    int i, j;
    double r_ij;
    double f_r;
    dvec box;
    dvec x_ij;
    tW_word site_names[2];
    tW_ref_inter *ij_inter;

    if (DEBUG_ref_nb) {
	fprintf(fp, "\n  In calc_ref_nb_pair_forces().\n");
    }

    if (info.b_PBC) {
	if (!setup_box_dimensions(box, info.box)) {
	    exit(EXIT_FAILURE);
	}
    }

    for (i = 0; i < N_sites - 1; i++) {
	for (j = i + 1; j < N_sites; j++) {
	    if (skip_excl(CG_struct[i].nr_excl, CG_struct[i].excl_list, j)) {
		continue;
	    }

	    r_ij =
		get_difference_vector_PBC(info.b_PBC, box, CG_struct[i].r,
					  CG_struct[j].r, x_ij);

	    if (DEBUG_ref_nb) {
		fprintf(fp, "    r_%d%d: %f\n", i, j, r_ij);
	    }

	    copy_2_site_names(site_names, CG_struct[i].name,
			      CG_struct[j].name);

	    ij_inter =
		get_inter_ptr(site_names, ref_potential.nb_pair_inter,
			      ref_potential.N_nb_pair_inter);

	    if (ij_inter == NULL) {
		continue;
	    }

	    if (r_ij >= ij_inter->x_max) {
		continue;
	    }

	    f_r =
		get_scalar_f_from_table(r_ij, ij_inter,
					ref_potential.interpolation_index);

	    if (DEBUG_ref_nb) {
		fprintf(fp, "    f_%d%d: %f\n", i, j, f_r);
	    }

	    get_central_pair_vector_forces(x_ij, r_ij, f_r,
					   &(CG_struct[i]),
					   &(CG_struct[j]));

	}
    }

    if (DEBUG_ref_nb) {
	fprintf(fp, "  Leaving calc_ref_nb_pair_forces().\n");
    }

}


/*****************************************************************************************
get_difference_vector_PBC(): If periodic boundary conditions are present, then applies
the minimum image convention. Only works for cubic boxes. 
*****************************************************************************************/
double get_difference_vector_PBC(bool b_PBC, dvec box, dvec x_i, dvec x_j,
				 dvec x_ij)
{
    double r_ij;
    dvec x_i_copy;
    dvec x_j_copy;

    copy_vector(x_i, x_i_copy);
    copy_vector(x_j, x_j_copy);

    if (b_PBC) {
	det_min_image(box, x_i_copy, x_j_copy);
    }
    vect_diff(x_i_copy, x_j_copy, x_ij);
    r_ij = calc_norm(x_ij);

    return r_ij;
}


/*****************************************************************************************
copy_2_site_names():
*****************************************************************************************/
void copy_2_site_names(tW_word * site_names, tW_word name1, tW_word name2)
{
    strcpy(site_names[0], name1);
    strcpy(site_names[1], name2);
}


/*****************************************************************************************
copy_3_site_names():
*****************************************************************************************/
void copy_3_site_names(tW_word * site_names, tW_word name1, tW_word name2,
		       tW_word name3)
{
    strcpy(site_names[0], name1);
    strcpy(site_names[1], name2);
    strcpy(site_names[2], name3);
}


/*****************************************************************************************
copy_4_site_names():
*****************************************************************************************/
void copy_4_site_names(tW_word * site_names, tW_word name1, tW_word name2,
		       tW_word name3, tW_word name4)
{
    strcpy(site_names[0], name1);
    strcpy(site_names[1], name2);
    strcpy(site_names[2], name3);
    strcpy(site_names[3], name4);
}


/*****************************************************************************************
get_IntraMolecPair_inter_ptr(): Returns the pointer to a intramolecular interaction 
depending on the sequence in site_names. Returns NULL if no interaction for the pair.
I need to be able to specify the nuber of bonds away so that I could have for example a
B-B 1-4 interaction and a B-L 1-4 interaction.
*****************************************************************************************/
tW_ref_inter *get_IntraMolecPair_inter_ptr(tW_word * site_names,
					   tW_ref_inter * inter,
					   int N_inter)
{
    int i;

    for (i = 0; i < N_inter; i++) {
	if (match_word_list(2, site_names, inter[i].site_types)) {
	    return &(inter[i]);
	}
    }

    return NULL;
}


/*****************************************************************************************
get_dihedral_inter_ptr(): Returns the pointer to a dihedral interaction depending on the
squence found in site_names. Returns NULL if no interaction for the quartet.
*****************************************************************************************/
tW_ref_inter *get_dihedral_inter_ptr(tW_word * site_names,
				     tW_ref_inter * inter, int N_inter)
{
    int i;

    for (i = 0; i < N_inter; i++) {
	if (match_word_list(4, site_names, inter[i].site_types)) {
	    return &(inter[i]);
	}
    }

    return NULL;
}


/*****************************************************************************************
get_angle_inter_ptr(): Returns the pointer to the angle interaction depending on the
sequence found in site_names. Returns NULL if no interaction for the triplet.
*****************************************************************************************/
tW_ref_inter *get_angle_inter_ptr(tW_word * site_names,
				  tW_ref_inter * inter, int N_inter)
{
    int i;

    for (i = 0; i < N_inter; i++) {
//I could just use match_word_list.
	if (strcmp(site_names[1], inter[i].site_types[1]) == 0) {
	    if (check_word_list
		(inter[i].N_inter_sites, site_names,
		 inter[i].site_types)) {
		return &(inter[i]);
	    }
	}
    }

    return NULL;
}


/*****************************************************************************************
get_inter_ptr(): Returns an pointer to the interaction based on the sequence of sites
in site_names. Returns NULL if no interaction is present for that sequence of sites.
*****************************************************************************************/
tW_ref_inter *get_inter_ptr(tW_word * site_names, tW_ref_inter * inter,
			    int N_inter)
{
    int i;

    for (i = 0; i < N_inter; i++) {
	if (check_word_list
	    (inter[i].N_inter_sites, site_names, inter[i].site_types)) {
	    return &(inter[i]);
	}
    }

    return NULL;
}


/*****************************************************************************************
get_scalar_f_from_table(): Returns f(x) which were read in from input.
*****************************************************************************************/
double get_scalar_f_from_table(double x, tW_ref_inter * inter,
			       int interpolation_index)
{
    int coeff;
    double f = 0;

    coeff = (int) floor((x - inter->x_0) / inter->dx);

    if (interpolation_index == REF_NO_INTERPOL_I) {
	f = get_no_interpolation_f(inter, coeff, x);
    }

    else if (interpolation_index == REF_LINEAR_INTERPOL_I) {
	f = get_linear_interpolation_f(inter, coeff, x);
    }

    else {
	printf("\nERROR: Interpolation method not accepted.\n");
	exit(EXIT_FAILURE);
    }

    return f;
}


/*****************************************************************************************
get_no_interpolation_f(): Returns f(x) without interpolating between grid points.
*****************************************************************************************/
double get_no_interpolation_f(tW_ref_inter * inter, int coeff, double x)
{
    double resid;

    resid = x - (inter->x_0 + coeff * inter->dx);
    if (resid > 0.5 * inter->dx) {
	coeff += 1;
    }

    return inter->f[coeff];
}


/*****************************************************************************************
get_linear_interpolation_f(): Returns f(x) by linearly interpolating between grid points
on either side of x.
*****************************************************************************************/
double get_linear_interpolation_f(tW_ref_inter * inter, int coeff,
				  double x)
{
    double f, L_0, L_1;

//I need to double check this.
    if (coeff == 0) {
	coeff++;
    }
    if (coeff == inter->N_pts - 1) {
	coeff--;
    }

    L_0 =
	(x - inter->x[coeff + 1]) / (inter->x[coeff] -
				     inter->x[coeff + 1]);
    L_1 = (x - inter->x[coeff]) / (inter->x[coeff + 1] - inter->x[coeff]);

    f = (L_0 * inter->f[coeff]) + (L_1 * inter->f[coeff + 1]);

    return f;
}


/*****************************************************************************************
get_central_pair_vector_forces(): Calculates the force between two pairs of site.
*****************************************************************************************/
void get_central_pair_vector_forces(dvec x_ij, double r_ij, double f_r,
				    tW_CG_site * i_site,
				    tW_CG_site * j_site)
{
    dvec f_ij;

    scal_times_vect(f_r / r_ij, x_ij, f_ij);

    vect_sum(i_site->ref_f, f_ij, i_site->ref_f);

    vect_diff(j_site->ref_f, f_ij, j_site->ref_f);
}


/*****************************************************************************************
calc_ref_BondStretch_forces(): Calculates the forces due to bond stretching.
*****************************************************************************************/
void calc_ref_BondStretch_forces(tW_gmx_topology *top, 
                                 tW_CG_site * CG_struct,
				 tW_ref_potential ref_potential, 
                                 FILE * fp)
{
    int i = 0;
    double r_ij;
    double f_r;
    dvec x_ij;
    tW_word site_names[2];
    tW_ref_inter *bond_inter;
    int *bond_list;
    int n_bonds;

    bond_list = top->get_bond_list(top);
    n_bonds = top->get_nbonds(top);


    if (DEBUG_ref_bs) {
	fprintf(fp, "  In calc_ref_BondStretch_forces().");
    }


    /* Loop over the number of bonds. */
    while (i < n_bonds) {

	copy_2_site_names(site_names,
			  CG_struct[bond_list[i + 1]].name,
			  CG_struct[bond_list[i + 2]].name);

	bond_inter =
	    get_inter_ptr(site_names, ref_potential.bondstretch_inter,
			  ref_potential.N_bondstretch_inter);
	if (bond_inter == NULL) {
	    i += 3;
	    continue;
	}

	vect_diff(CG_struct[bond_list[i + 1]].r,
		  CG_struct[bond_list[i + 2]].r, x_ij);

	r_ij = calc_norm(x_ij);

	if (check_inter_range
	    (r_ij, bond_inter->x_0, bond_inter->x_max, bond_inter->dx)
	    != 0) {
	    printf
		("\nERROR: Ref. Potent. out of range for BondStretch.\n");
	    exit(EXIT_FAILURE);
	}

	f_r =
	    get_scalar_f_from_table(r_ij, bond_inter,
				    ref_potential.interpolation_index);

	if (DEBUG_ref_bs) {
	    fprintf(fp, "    r_%d%d: %.8e  f_r: %.8e\n",
		    bond_list[i + 1], bond_list[i + 2], r_ij,
		    f_r);
	}

	get_central_pair_vector_forces(x_ij, r_ij, f_r,
				       &(CG_struct
					 [bond_list[i + 1]]),
				       &(CG_struct
					 [bond_list[i + 2]]));

	i += 3;
    }

    if (DEBUG_ref_bs) {
	fprintf(fp, "  Leaving calc_ref_BondStretch_forces().\n");
    }

    free(bond_list);
}


/*****************************************************************************************
calc_ref_Angle_forces():
*****************************************************************************************/
void calc_ref_Angle_forces(tW_gmx_topology *top, tW_CG_site * CG_struct,
			   tW_ref_potential ref_potential, FILE * fp)
{
    int i = 0;
    double theta, f_theta;
    dvec L_2, L_3;
    tW_word site_names[3];
    tW_ref_inter *angle_inter;
    int *angle_list;
    int n_angles;

    angle_list = top->get_angle_list(top);
    n_angles = top->get_nangles(top);


    if (DEBUG_ref_ab) {
	fprintf(fp, "  In calc_ref_Angle_forces().\n");
    }

    while (i < n_angles) {
	copy_3_site_names(site_names,
			  CG_struct[angle_list[i + 1]].name,
			  CG_struct[angle_list[i + 2]].name,
			  CG_struct[angle_list[i + 3]].name);

	angle_inter =
	    get_angle_inter_ptr(site_names, ref_potential.angle_inter,
				ref_potential.N_angle_inter);
	if (angle_inter == NULL) {
	    i += 4;
	    continue;
	}

	if (DEBUG_ref_ab) {
	    fprintf(fp, "    Did not skip inter. for %s-%s-%s.\n",
		    site_names[0], site_names[1], site_names[2]);
	    fprintf(fp, "    Matched with: %s-%s-%s.\n",
		    angle_inter->site_types[0], angle_inter->site_types[1],
		    angle_inter->site_types[2]);
	}

	theta = get_ref_Angle_info(CG_struct[angle_list[i + 2]].r,
				   CG_struct[angle_list[i + 1]].r,
				   CG_struct[angle_list[i + 3]].r,
				   L_2, L_3, fp);

	if (check_inter_range
	    (theta, angle_inter->x_0, angle_inter->x_max, angle_inter->dx)
	    != 0) {
	    printf("\nERROR: Ref. Potent. out of range for Angle: %f.\n",
		   theta);
	    printf("****R_min: %f  R_max: %f\n", angle_inter->x_0,
		   angle_inter->x_max);
	    exit(EXIT_FAILURE);
	}

	f_theta = get_scalar_f_from_table(theta, angle_inter,
					  ref_potential.
					  interpolation_index);

	if (DEBUG_ref_ab) {
	    fprintf(fp, "    theta: %f  f_theta: %f\n", theta, f_theta);
	}

	get_angle_vector_forces(f_theta, L_2, L_3, theta,
				&(CG_struct[angle_list[i + 2]]),
				&(CG_struct[angle_list[i + 1]]),
				&(CG_struct[angle_list[i + 3]]));

	i += 4;
    }

    if (DEBUG_ref_ab) {
	fprintf(fp, "  Leaving calc_ref_Angle_forces().\n");
    }

    free(angle_list);
}


/*****************************************************************************************
get_ref_Angle_info(): Calculates the forces on the sites for the angle defined as <213.
*****************************************************************************************/
double get_ref_Angle_info(dvec x_1, dvec x_2, dvec x_3, dvec L_2, dvec L_3,
			  FILE * fp)
{
    double theta, cos_theta, norm_12, norm_13;
    dvec r_12, r_13, u_12, u_13;

    get_difference_unit_vector(x_1, x_2, r_12, u_12, &norm_12);
    get_difference_unit_vector(x_1, x_3, r_13, u_13, &norm_13);
    cos_theta = dot_prod(u_12, u_13);
    theta = acos(cos_theta);
    get_L2_L3_angles(u_12, u_13, cos_theta, norm_12, norm_13, L_2, L_3);

    return theta;
}


/*****************************************************************************************
get_angle_vector_forces(): Calculate the forces due to the angle beding reference 
potential.
*****************************************************************************************/
void get_angle_vector_forces(double f_theta, dvec L_2, dvec L_3,
			     double theta, tW_CG_site * site1,
			     tW_CG_site * site2, tW_CG_site * site3)
{
    int i;
    double sin_theta = sin(theta);
    dvec f_1, f_2, f_3;

    if (fabs(sin_theta) < FLOAT_EPS) {
	sin_theta = FLOAT_EPS;
    }

    scal_times_vect(f_theta / sin_theta, L_2, f_2);
    scal_times_vect(f_theta / sin_theta, L_3, f_3);
    vect_sum(f_2, f_3, f_1);

    vect_diff(site1->ref_f, f_1, site1->ref_f);
    vect_sum(site2->ref_f, f_2, site2->ref_f);
    vect_sum(site3->ref_f, f_3, site3->ref_f);
}


/*****************************************************************************************
calc_ref_Dihedral_forces(): Calculate the forces due to the reference dihedral potential.
*****************************************************************************************/
void calc_ref_Dihedral_forces(tW_gmx_topology *top,
			      tW_CG_site * CG_struct,
			      tW_ref_potential ref_potential, 
			      FILE * fp)
{
    int i = 0;
    double phi, f_phi;
    dvec B[4];
    tW_word site_names[4];
    tW_dihedral dihedral;
    tW_ref_inter *dihedral_inter;
    int *dihedral_list;
    int n_dihs;

    dihedral_list = top->get_dih_list(top);
    n_dihs = top->get_ndihs(top);


    if (DEBUG_ref_dh) {
	fprintf(fp, "  In calc_ref_Dihedral_forces().\n");
    }

    while (i < n_dihs) {
	copy_4_site_names(site_names,
			  CG_struct[dihedral_list[i + 1]].name,
			  CG_struct[dihedral_list[i + 2]].name,
			  CG_struct[dihedral_list[i + 3]].name,
			  CG_struct[dihedral_list[i + 4]].name);


	dihedral_inter =
	    get_dihedral_inter_ptr(site_names,
				   ref_potential.dihedral_inter,
				   ref_potential.N_dihedral_inter);
	if (dihedral_inter == NULL) {
	    i += 5;
	    continue;
	}


	dihedral.dihedral_quartet = &(dihedral_list[i + 1]);
	phi = get_ref_Dihedral_info(&dihedral, CG_struct, B);

	if (check_inter_range
	    (phi, dihedral_inter->x_0, dihedral_inter->x_max,
	     dihedral_inter->dx)
	    != 0) {
	    printf
		("\nERROR: Ref. Potent. out of range for Dihedral: %f.\n",
		 phi);
	    exit(EXIT_FAILURE);
	}

	f_phi = get_scalar_f_from_table(phi, dihedral_inter,
					ref_potential.interpolation_index);

	if (DEBUG_ref_dh) {
	    fprintf(fp, "    phi: %f  f_phi: %f  sin_phi: %f\n",
		    phi, f_phi, sin(phi * M_PI / 180.0));
	    fprintf(fp, "    Sites: %s-%s-%s-%s\n", site_names[0],
		    site_names[1], site_names[2], site_names[3]);
	    fprintf(fp, "    Matched: %s-%s-%s-%s\n",
		    dihedral_inter->site_types[0],
		    dihedral_inter->site_types[1],
		    dihedral_inter->site_types[2],
		    dihedral_inter->site_types[3]);
	}

	get_dihedral_vector_forces(f_phi, phi, B,
				   &(CG_struct
				     [dihedral_list[i + 1]]),
				   &(CG_struct
				     [dihedral_list[i + 2]]),
				   &(CG_struct
				     [dihedral_list[i + 3]]),
				   &(CG_struct
				     [dihedral_list[i + 4]]));

	i += 5;
    }

    if (DEBUG_ref_dh) {
	fprintf(fp, "  Leaving calc_ref_Dihedral_forces().\n");
    }


    free(dihedral_list);
}


/*****************************************************************************************
calc_ref_IntraMolecPairs_forces(): Calculate forces due to intramolecular pair forces.
*****************************************************************************************/
void calc_ref_IntraMolecPairs_forces(tW_gmx_topology *top,
				     tW_CG_site * CG_struct,
				     tW_ref_potential ref_potential,
				     FILE * fp)
{
    int i = 0;
    double r_ij, f_r;
    tW_word site_names[2];
    tW_ref_inter *pair_inter;
    dvec x_ij;
    int *pair_list;
    int n_pairs;

    pair_list = top->get_pair_list(top);
    n_pairs = top->get_npairs(top);

    if (DEBUG_ref_ip) {
	fprintf(fp, "  In calc_ref_IntraMolecPairs_forces().\n");
    }

    while (i < n_pairs) {
	copy_2_site_names(site_names,
			  CG_struct[pair_list[i + 1]].name,
			  CG_struct[pair_list[i + 2]].name);

	pair_inter =
	    get_IntraMolecPair_inter_ptr(site_names,
					 ref_potential.
					 IntraMolecPairs_inter,
					 ref_potential.
					 N_IntraMolecPairs_inter);
	if (pair_inter == NULL) {
	    i += 3;
	    continue;
	}

	vect_diff(CG_struct[pair_list[i + 1]].r,
		  CG_struct[pair_list[i + 2]].r, x_ij);

	r_ij = calc_norm(x_ij);

	if (check_inter_range
	    (r_ij, pair_inter->x_0, pair_inter->x_max, pair_inter->dx)
	    != 0) {
	    printf
		("\nERROR: Ref. Potent. out of range for BondStretch.\n");
	    exit(EXIT_FAILURE);
	}

	f_r =
	    get_scalar_f_from_table(r_ij, pair_inter,
				    ref_potential.interpolation_index);

	get_central_pair_vector_forces(x_ij, r_ij, f_r,
				       &(CG_struct
					 [pair_list[i + 1]]),
				       &(CG_struct
					 [pair_list[i + 2]]));

	i += 3;
    }

    if (DEBUG_ref_ip) {
	fprintf(fp, "  Leaving calc_ref_IntraMolecPairs_forces().\n");
    }

    free(pair_list);
}


/*****************************************************************************************
get_ref_Dihedral_info(): Copies information for a dihedral interactions into the 
variable dihedral.
*****************************************************************************************/
double get_ref_Dihedral_info(tW_dihedral * dihedral,
			     tW_CG_site * CG_struct, dvec * B)
{
    double phi;

    eval_grad_B(dihedral, CG_struct);
    copy_vector(dihedral->grad_B[0], B[0]);
    copy_vector(dihedral->grad_B[1], B[1]);
    copy_vector(dihedral->grad_B[2], B[2]);
    copy_vector(dihedral->grad_B[3], B[3]);

    phi = dihedral->angle;

    if (phi > M_PI) {
	phi = M_PI;
    }
    if (phi < (-1.0 * M_PI)) {
	phi = (-1.0 * M_PI);
    }

    return phi;
}


/*****************************************************************************************
get_dihedral_vector_forces(): Calculates the forces due to the dihedral reference
potential.
*****************************************************************************************/
void get_dihedral_vector_forces(double f_phi, double phi, dvec * B,
				tW_CG_site * site1, tW_CG_site * site2,
				tW_CG_site * site3, tW_CG_site * site4)
{
    double sin_phi = sin(phi);
    dvec f_1, f_2, f_3, f_4;

    if (fabs(sin_phi) < FLOAT_EPS) {
	scal_times_vect(0.0, B[0], f_1);
	scal_times_vect(0.0, B[1], f_2);
	scal_times_vect(0.0, B[2], f_3);
	scal_times_vect(0.0, B[3], f_4);

	vect_sum(site1->ref_f, f_1, site1->ref_f);
	vect_sum(site2->ref_f, f_2, site2->ref_f);
	vect_sum(site3->ref_f, f_3, site3->ref_f);
	vect_sum(site4->ref_f, f_4, site4->ref_f);
    } else {
	scal_times_vect(f_phi / sin_phi, B[0], f_1);
	scal_times_vect(f_phi / sin_phi, B[1], f_2);
	scal_times_vect(f_phi / sin_phi, B[2], f_3);
	scal_times_vect(f_phi / sin_phi, B[3], f_4);

	vect_sum(site1->ref_f, f_1, site1->ref_f);
	vect_sum(site2->ref_f, f_2, site2->ref_f);
	vect_sum(site3->ref_f, f_3, site3->ref_f);
	vect_sum(site4->ref_f, f_4, site4->ref_f);
    }
}


/*****************************************************************************************
print_CG_ref_f(): Prints the forces due to the reference potential.
*****************************************************************************************/
void print_CG_ref_f(int N_sites, tW_CG_site * CG_struct, FILE * fp)
{
    int i, j;

    fprintf(fp, "\n\n  In print_CG_ref_f().\n");

    for (i = 0; i < N_sites; i++) {
	for (j = 0; j < 3; j++) {
	    fprintf(fp,
		    "   Site %-4d: f_gmx[%d]: %16.6f f_ref[%d]: %16.6f diff: %16.6f percent: %16.6f\n",
		    i, j, CG_struct[i].f[j], j, CG_struct[i].ref_f[j],
		    (CG_struct[i].f[j] - CG_struct[i].ref_f[j]),
		    (CG_struct[i].f[j] -
		     CG_struct[i].ref_f[j]) * 100.0 /
		    CG_struct[i].ref_f[j]);
	}
    }
    fprintf(fp, "  Leaving print_CG_ref_f().\n\n");
}

/*****************************************************************************************
print_bref(): JFR - 07.16.12: Prints bref.
*****************************************************************************************/
void print_bref(int N_coeff, double *bref, tW_word tag)
{
    int i;
    FILE *fp;
    tW_word fnm;

    sprintf(fnm, "bref.%s.dat", tag);
    printf("ref filename = %s \n", fnm);
    fp = fopen(fnm, "w");

    for (i = 0; i < N_coeff; i++) {
	fprintf(fp, "%.15lf \n", bref[i]);
    }

    fclose(fp);
}
