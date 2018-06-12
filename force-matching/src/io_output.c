/**
@file io_output.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Functions related to writing cgff output
*/

//c library includes
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//local includes
#include "cgff_types.h"
#include "io_read.h"
#include "io_output.h"
//#include "gmx-interface.h"
#include "gromacs_topology.h"
#include "wnoid_math.h"
#include "calc_grids.h"


/*****************************************************************************************
get_top_tag(): Saves the name of the GROMACS trajectory file minus the extension into
the variable tag. This is used for naming files for each topology.
NJD 5.24.16 - Now also strips the folder / file location preceding the file name. Currently
will not support windows file locations.
*****************************************************************************************/
int get_top_tag(tW_word coord_fnm, tW_word tag)
{
    char *marker;
    char *name_start;

    // Remove file location
    name_start = strrchr(coord_fnm, '/');

    if (name_start != NULL)
    {
        strcpy(tag, name_start+1);
    } else
    {
        strcpy(tag, coord_fnm);
    }

    // Remove file extension
    marker = strstr(tag, ".trr\0");
    strcpy(marker, "\0");

    return 0;
}


/*****************************************************************************************
open_output_file1(): Opens a file to print general information about the results for 
the topology labeled with tag.
*****************************************************************************************/
FILE *open_output_file1(tW_word tag)
{
    FILE *fp;
    tW_word fname;

    sprintf(fname, "out.%s.txt", tag);
    if ((fp = fopen(fname, "w"))) {
	return fp;
    }

    printf("\nERROR: Could not open %s for writing output.\n", fname);
    exit(EXIT_FAILURE);

    return fp;
}


/*****************************************************************************************
print_output(): Print the results. General information and results for functional forms
are printed to the file pointed to by fp.
*****************************************************************************************/
int print_output(int forces, tW_system sys, tW_word tag)
{
    int i;
    FILE *fp;
    tW_Inter_Types *inter;


    fp = open_output_file1(tag);

    print_line_stars(fp);

    /* JFR - 06.16.12: Chi2 */
    print_Chi2(fp, &sys);

    print_line_stars(fp);

    /* JFR - 09.10.12: decomp Chi2 */
    print_decomp_Chi2(fp, &sys);

    print_line_stars(fp);

    for (i = 0; i < sys.N_Inter_Types; i++) {
	inter = &(sys.Inter_Types[i]);

	if (was_inter_present(inter->ptr_g_cnt, inter->N_coeff) == FALSE) {
	    continue;
	}

	if (inter->i_basis == DELTA_BASIS_INDEX) {
	    if (sys.CalcMODE_var.CalcMODE != ISECOND_HALF) {
		print_delta_basis_output(fp, forces, sys, tag, inter, i);
	    }

	    if (sys.CalcMODE_var.CalcMODE != IFIRST_HALF) {
		print_delta_basis_forces(fp, forces, sys, tag, inter, i);
	    }
	}
	/* START JFR */
	else if (inter->i_basis == LINEAR_BASIS_INDEX) {
	    if (sys.CalcMODE_var.CalcMODE != ISECOND_HALF) {
		print_linear_basis_output(fp, forces, sys, tag, inter, i);
	    }

	    if (sys.CalcMODE_var.CalcMODE != IFIRST_HALF) {
		print_linear_basis_forces(fp, forces, sys, tag, inter, i);
	    }
	}
	/* END JFR */
	else if (inter->i_basis == BSPLINE_BASIS_INDEX) {	/* JFR - 07.22.12 */
	    if (sys.CalcMODE_var.CalcMODE != ISECOND_HALF) {
		print_Bspline_basis_output(fp, forces, sys, tag, inter, i);
	    }

	    if (sys.CalcMODE_var.CalcMODE != IFIRST_HALF) {
		print_Bspline_basis_forces(fp, forces, sys, tag, inter, i);
	    }
	} else if (inter->i_basis == HARMONIC_BASIS_INDEX) {
	    if (sys.CalcMODE_var.CalcMODE != ISECOND_HALF) {
		print_harmonic_basis_output(fp, forces, sys, tag, inter,
					    i);
	    }

	    if (sys.CalcMODE_var.CalcMODE != IFIRST_HALF) {
		print_harmonic_basis_forces(fp, forces, sys, tag, inter,
					    i);
	    }
	} else if (inter->i_basis == RYCKAERT_BELLEMANS_BASIS_INDEX) {
	    if (sys.CalcMODE_var.CalcMODE != ISECOND_HALF) {
		print_RB_basis_output(fp, forces, sys, tag, inter, i);
	    }

	    if (sys.CalcMODE_var.CalcMODE != IFIRST_HALF) {
		print_RB_basis_forces(fp, forces, sys, tag, inter, i);
	    }
	} else if (inter->i_basis == TOY_DIHED_INDEX) {
	    if (sys.CalcMODE_var.CalcMODE != ISECOND_HALF) {
		print_TOY_basis_output(fp, forces, sys, tag, inter, i);
	    }

	    if (sys.CalcMODE_var.CalcMODE != IFIRST_HALF) {
		print_TOY_basis_forces(fp, forces, sys, tag, inter, i);
	    }
	} else if (inter->i_basis == POWER_INDEX) {
	    if (sys.CalcMODE_var.CalcMODE != ISECOND_HALF) {
		print_power_basis_output(fp, forces, sys, tag, inter, i);
	    }

	    if (sys.CalcMODE_var.CalcMODE != IFIRST_HALF) {
		print_power_basis_forces(fp, forces, sys, tag, inter, i);
	    }
	} else {
	    printf
		("\nERROR: Basis type '%d' not supported in print_output().\n",
		 inter->i_basis);
	    exit(EXIT_FAILURE);
	}

    }

    return 0;
}


/*****************************************************************************************
open_output_file2(): Opens files for printing tabulated results.
*****************************************************************************************/
FILE *open_output_file2(tW_word function, tW_word inter_type,
			tW_word inter_name, tW_word tag)
{
    FILE *fp;
    tW_word fname;

    sprintf(fname, "%s.%s.%s.%s.dat", function, inter_type, tag,
	    inter_name);
    if ((fp = fopen(fname, "w"))) {
	return fp;
    }

    printf("\nERROR: Could not open %s for writing output.\n", fname);
    exit(EXIT_FAILURE);

}


/*****************************************************************************************
print_delta_basis_output(): Prints the results for delta basis functions. Prints results
in degrees for angles and dihedrals. Adds the 180 degree grid point for dihedrals.
*****************************************************************************************/
void print_delta_basis_output(FILE * fp, int forces, tW_system sys,
			      tW_word tag, tW_Inter_Types * inter, int i)
{
    int j;
    double dx = inter->dr;
    FILE *fp_g, *fp_g_cnt, *fp_L, *fp_b_forces, *fp_b_struct;

    print_general_inter_output(fp, tag, inter, i);

    /* Open files for writing output. */
    fp_g = open_output_file2("g",
			     inter->inter_type, inter->inter_name, tag);
    fp_g_cnt = open_output_file2("g_cnt",
				 inter->inter_type, inter->inter_name,
				 tag);
    fp_L =
	open_output_file2("L", inter->inter_type, inter->inter_name, tag);
    fp_b_struct =
	open_output_file2("b_struct", inter->inter_type, inter->inter_name,
			  tag);
    if (forces) {
	fp_b_forces = open_output_file2("b_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Account for angle and dihedral forces in degrees. */
    if ((strcmp(B_ANGLE, inter->inter_type) == 0) ||
	(strcmp(B_DIHEDRAL, inter->inter_type) == 0)
	) {
	for (j = 0; j < inter->N_pts; j++) {
	    inter->ptr_x[j] *= 180.0 / M_PI;
	}
    }

    /* Print results. */
    for (j = 0; j < inter->N_pts; j++) {
	fprintf(fp_g, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_g[j] / dx);
	fprintf(fp_g_cnt, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_g_cnt[j] / dx);
	fprintf(fp_L, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_L[j] / dx);
	fprintf(fp_b_struct, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_b_struct[j] / dx);
	if (forces) {
	    fprintf(fp_b_forces, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		    inter->ptr_b_forces[j] / dx);
	}
    }

    /* Add the 180 grid point. Periodicity. */
    if ((strcmp(B_DIHEDRAL, inter->inter_type) == 0) &&
	(fabs(inter->ptr_x[0] + 180.0) < FLOAT_EPS)
	) {
	fprintf(fp_g, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_g[0] / dx);
	fprintf(fp_g_cnt, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_g_cnt[0] / dx);
	fprintf(fp_L, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_L[0] / dx);
	fprintf(fp_b_struct, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_b_struct[0] / dx);
	if (forces) {
	    fprintf(fp_b_forces, "%10.5lf    %15.5lf \n", 180.0,
		    inter->ptr_b_forces[0] / dx);
	}
    }

    /* Close files for writing output. */
    fclose(fp_g);
    fclose(fp_g_cnt);
    fclose(fp_L);
    fclose(fp_b_struct);
    if (forces) {
	fclose(fp_b_forces);
    }

    fprintf(fp, "\n");
    print_line_stars(fp);
}

/*****************************************************************************************
print_delta_basis_forces(): JFR - added 04.06.12:  This function is the same as
print_delta_basis_output().  I just separated the printing of the forces for the different
Calculation Modes.
*****************************************************************************************/
void print_delta_basis_forces(FILE * fp, int forces, tW_system sys,
			      tW_word tag, tW_Inter_Types * inter, int i)
{
    int j;
    FILE *fp_f_forces, *fp_f_struct;

    print_general_inter_output(fp, tag, inter, i);

    /* Open files for writing output. */
    fp_f_struct = open_output_file2("f_struct",
				    inter->inter_type, inter->inter_name,
				    tag);

    if (forces) {
	fp_f_forces = open_output_file2("f_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Account for angle and dihedral forces in degrees. */
    if ((strcmp(B_ANGLE, inter->inter_type) == 0) ||
	(strcmp(B_DIHEDRAL, inter->inter_type) == 0)
	) {
	for (j = 0; j < inter->N_pts; j++) {
	    if (sys.CalcMODE_var.CalcMODE != IFULL) {	/* JFR - prevent this from being done twice */
		inter->ptr_x[j] *= 180.0 / M_PI;
	    }
	    inter->ptr_phi_struct[j] *= M_PI / 180.0;
	    if (forces) {
		inter->ptr_phi_forces[j] *= M_PI / 180.0;
	    }
	}
    }

/* Print results. */
    for (j = 0; j < inter->N_pts; j++) {
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_phi_struct[j]);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		    inter->ptr_phi_forces[j]);
	}
    }

    /* Add the 180 grid point. Periodicity. */
    if ((strcmp(B_DIHEDRAL, inter->inter_type) == 0) &&
	(fabs(inter->ptr_x[0] + 180.0) < FLOAT_EPS)
	) {
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_phi_struct[0]);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", 180.0,
		    inter->ptr_phi_forces[0]);
	}
    }

    /* Close files for writing output. */
    fclose(fp_f_struct);
    if (forces) {
	fclose(fp_f_forces);
    }

    fprintf(fp, "\n");
    print_line_stars(fp);
}

/* START JFR */
/*****************************************************************************************
print_linear_basis_output(): Prints the results for delta basis functions. Prints results
in degrees for angles and dihedrals. Adds the 180 degree grid point for dihedrals.
*****************************************************************************************/
void print_linear_basis_output(FILE * fp, int forces, tW_system sys,
			       tW_word tag, tW_Inter_Types * inter, int i)
{
    int j;
    double dx = inter->dr;
    FILE *fp_g, *fp_g_cnt, *fp_L, *fp_b_forces, *fp_b_struct;

    print_general_inter_output(fp, tag, inter, i);

    /* Open files for writing output. */
    fp_g = open_output_file2("g",
			     inter->inter_type, inter->inter_name, tag);
    fp_g_cnt = open_output_file2("g_cnt",
				 inter->inter_type, inter->inter_name,
				 tag);
    fp_L =
	open_output_file2("L", inter->inter_type, inter->inter_name, tag);
    fp_b_struct =
	open_output_file2("b_struct", inter->inter_type, inter->inter_name,
			  tag);
    if (forces) {
	fp_b_forces = open_output_file2("b_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Account for angle and dihedral forces in degrees. */
    if ((strcmp(B_ANGLE, inter->inter_type) == 0) ||
	(strcmp(B_DIHEDRAL, inter->inter_type) == 0)
	) {
	for (j = 0; j < inter->N_pts; j++) {
	    inter->ptr_x[j] *= 180.0 / M_PI;
	}
    }

    /* Print results. */
    for (j = 0; j < inter->N_pts; j++) {
	fprintf(fp_g, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_g[j] / dx);
	fprintf(fp_g_cnt, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_g_cnt[j] / dx);
	fprintf(fp_L, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_L[j] / dx);
	fprintf(fp_b_struct, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_b_struct[j] / dx);
	if (forces) {
	    fprintf(fp_b_forces, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		    inter->ptr_b_forces[j] / dx);
	}
    }

    /* Add the 180 grid point. Periodicity. */
    if ((strcmp(B_DIHEDRAL, inter->inter_type) == 0) &&
	(fabs(inter->ptr_x[0] + 180.0) < FLOAT_EPS)
	) {
	fprintf(fp_g, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_g[0] / dx);
	fprintf(fp_g_cnt, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_g_cnt[0] / dx);
	fprintf(fp_L, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_L[0] / dx);
	fprintf(fp_b_struct, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_b_struct[0] / dx);
	if (forces) {
	    fprintf(fp_b_forces, "%10.5lf    %15.5lf \n", 180.0,
		    inter->ptr_b_forces[0] / dx);
	}
    }

    /* Close files for writing output. */
    fclose(fp_g);
    fclose(fp_g_cnt);
    fclose(fp_L);
    fclose(fp_b_struct);
    if (forces) {
	fclose(fp_b_forces);
    }

    fprintf(fp, "\n");
    print_line_stars(fp);
}

/* END JFR */

/*****************************************************************************************
print_linear_basis_forces(): JFR - added 04.06.12:  This function is the same as
print_linear_basis_output().  I just separated the printing of the forces for the different
Calculation Modes.
*****************************************************************************************/
void print_linear_basis_forces(FILE * fp, int forces, tW_system sys,
			       tW_word tag, tW_Inter_Types * inter, int i)
{
    int j;
    FILE *fp_f_forces, *fp_f_struct;

    print_general_inter_output(fp, tag, inter, i);

    /* Open files for writing output. */
    fp_f_struct = open_output_file2("f_struct",
				    inter->inter_type, inter->inter_name,
				    tag);
    if (forces) {
	fp_f_forces = open_output_file2("f_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Account for angle and dihedral forces in degrees. */
    if ((strcmp(B_ANGLE, inter->inter_type) == 0) ||
	(strcmp(B_DIHEDRAL, inter->inter_type) == 0)
	) {
	for (j = 0; j < inter->N_pts; j++) {
	    if (sys.CalcMODE_var.CalcMODE != IFULL) {	/* JFR - prevent this from being done twice */
		inter->ptr_x[j] *= 180.0 / M_PI;
	    }
	    inter->ptr_phi_struct[j] *= M_PI / 180.0;
	    if (forces) {
		inter->ptr_phi_forces[j] *= M_PI / 180.0;
	    }
	}
    }

    /* Print results. */
    for (j = 0; j < inter->N_pts; j++) {
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		inter->ptr_phi_struct[j]);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", inter->ptr_x[j],
		    inter->ptr_phi_forces[j]);
	}
    }

    /* Add the 180 grid point. Periodicity. */
    if ((strcmp(B_DIHEDRAL, inter->inter_type) == 0) &&
	(fabs(inter->ptr_x[0] + 180.0) < FLOAT_EPS)
	) {
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_phi_struct[0]);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", 180.0,
		    inter->ptr_phi_forces[0]);
	}
    }

    /* Close files for writing output. */
    fclose(fp_f_struct);
    if (forces) {
	fclose(fp_f_forces);
    }

    fprintf(fp, "\n");
    print_line_stars(fp);
}

/*****************************************************************************************
print_Bspline_basis_output(): Prints the results for Bspline basis functions. Prints results
in degrees for angles and dihedrals. Adds the 180 degree grid point for dihedrals.
JFR - 07.22.12: this is the same as the other spline output
*****************************************************************************************/
void print_Bspline_basis_output(FILE * fp, int forces, tW_system sys,
				tW_word tag, tW_Inter_Types * inter, int m)
{
    int j, i, l, grid, index;
    double dx = inter->dr;
    FILE *fp_g, *fp_g_cnt, *fp_L, *fp_b_forces, *fp_b_struct, *fp_g_coeff,
	*fp_g_cnt_coeff, *fp_L_coeff, *fp_b_forces_coeff,
	*fp_b_struct_coeff;
    double R, B;
    int N_finepts = 10;		/* JFR - 07.22.12: For now, print out a grid 10 times finer than calcualted */
    int k = inter->kspline;
    double L, g, g_cnt, b_forces, b_struct;
    double Nik[MAX_BSPLINE_COEFF];
    double Npik[MAX_BSPLINE_COEFF];

    double ptr_x;

    print_general_inter_output(fp, tag, inter, m);

    /* Open files for writing output. */
    fp_g_coeff = open_output_file2("g_coeff",
				   inter->inter_type, inter->inter_name,
				   tag);
    fp_g_cnt_coeff =
	open_output_file2("g_cnt_coeff", inter->inter_type,
			  inter->inter_name, tag);
    fp_L_coeff =
	open_output_file2("L_coeff", inter->inter_type, inter->inter_name,
			  tag);
    fp_b_struct_coeff =
	open_output_file2("b_struct_coeff", inter->inter_type,
			  inter->inter_name, tag);
    if (forces) {
	fp_b_forces_coeff = open_output_file2("b_forces_coeff",
					      inter->inter_type,
					      inter->inter_name, tag);
    }

    fp_g = open_output_file2("g",
			     inter->inter_type, inter->inter_name, tag);
    fp_g_cnt = open_output_file2("g_cnt",
				 inter->inter_type, inter->inter_name,
				 tag);
    fp_L =
	open_output_file2("L", inter->inter_type, inter->inter_name, tag);
    fp_b_struct =
	open_output_file2("b_struct", inter->inter_type, inter->inter_name,
			  tag);
    if (forces) {
	fp_b_forces = open_output_file2("b_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Print results. */
    for (j = 0; j < inter->N_pts; j++) {
	ptr_x = inter->ptr_x[j];
	/* Account for angle and dihedral forces in degrees. */
	if ((strcmp(B_ANGLE, inter->inter_type) == 0)
	    || (strcmp(B_DIHEDRAL, inter->inter_type) == 0)) {
	    ptr_x *= 180.0 / M_PI;
	}

	fprintf(fp_g_coeff, "%10.5lf    %15.5lf \n", ptr_x,
		inter->ptr_g[j] / dx);
	fprintf(fp_g_cnt_coeff, "%10.5lf    %15.5lf \n", ptr_x,
		inter->ptr_g_cnt[j] / dx);
	fprintf(fp_L_coeff, "%10.5lf    %15.5lf \n", ptr_x,
		inter->ptr_L[j] / dx);
	fprintf(fp_b_struct_coeff, "%10.5lf    %15.5lf \n", ptr_x,
		inter->ptr_b_struct[j] / dx);
	if (forces) {
	    fprintf(fp_b_forces_coeff, "%10.5lf    %15.5lf \n", ptr_x,
		    inter->ptr_b_forces[j] / dx);
	}
    }

    /* Add the 180 grid point. Periodicity. */
    if ((strcmp(B_DIHEDRAL, inter->inter_type) == 0)
	&& (fabs(inter->ptr_x[0] + 180.0) < FLOAT_EPS)) {
	fprintf(fp_g_coeff, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_g[0] / dx);
	fprintf(fp_g_cnt_coeff, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_g_cnt[0] / dx);
	fprintf(fp_L_coeff, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_L[0] / dx);
	fprintf(fp_b_struct_coeff, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_b_struct[0] / dx);
	if (forces) {
	    fprintf(fp_b_forces_coeff, "%10.5lf    %15.5lf \n", 180.0,
		    inter->ptr_b_forces[0] / dx);
	}
    }

    /* print output on a finer grid than the calculation */
    for (j = 0; j < inter->N_pts; j++) {
	index =
	    get_grid_index_for_linear_basis(inter->ptr_x[j] + FLOAT_EPS,
					    inter->i_0, inter->dr,
					    inter->R_min) - inter->i_0;

	for (l = 0; l < N_finepts; l++) {
	    b_forces = 0.00;
	    b_struct = 0.00;
	    g = 0.00;
	    g_cnt = 0.00;
	    L = 0.00;

	    R = inter->ptr_x[j] + l * (inter->dr * (1.0 / N_finepts));

	    if ((R > inter->ptr_x[inter->N_pts - 1])
		&& (strcmp(B_DIHEDRAL, inter->inter_type) != 0)) {
		continue;
	    }

	    B = calc_Bspline(inter->dr, inter->R_min, inter->i_0,
			     inter->N_pts, R, index, k, Nik, Npik);
	    for (i = 0; i < k; i++) {
		grid = j - i + k / 2;
		if (strcmp(B_DIHEDRAL, inter->inter_type) == 0) {
		    /* Account for periodicity */
		    if (grid < 0) {
			grid += inter->N_pts;
		    } else if (grid >= inter->N_pts) {
			grid -= (inter->N_pts);
		    }
		}

		B = Nik[k - 1 - i];

		if ((grid >= 0) && (grid < inter->N_pts)) {
		    g += B * inter->ptr_g[grid];
		    g_cnt += B * inter->ptr_g_cnt[grid];
		    L += B * inter->ptr_L[grid];
		    b_struct += B * inter->ptr_b_struct[grid];
		    if (forces) {
			b_forces += B * inter->ptr_b_forces[grid];
		    }
		}
	    }

	    if ((strcmp(B_ANGLE, inter->inter_type) == 0)
		|| (strcmp(B_DIHEDRAL, inter->inter_type) == 0)) {
		R *= 180.0 / M_PI;
	    }
	    fprintf(fp_g, "%10.5lf    %15.5lf \n", R, g / dx);
	    fprintf(fp_g_cnt, "%10.5lf    %15.5lf \n", R, g_cnt / dx);
	    fprintf(fp_L, "%10.5lf    %15.5lf \n", R, L / dx);
	    fprintf(fp_b_struct, "%10.5lf    %15.5lf \n", R,
		    b_struct / dx);
	    if (forces) {
		fprintf(fp_b_forces, "%10.5lf    %15.5lf \n", R,
			b_forces / dx);
	    }

	}
    }

    /* Close files for writing output. */
    fclose(fp_g_coeff);
    fclose(fp_g_cnt_coeff);
    fclose(fp_L_coeff);
    fclose(fp_b_struct_coeff);
    if (forces) {
	fclose(fp_b_forces_coeff);
    }
    fclose(fp_g);
    fclose(fp_g_cnt);
    fclose(fp_L);
    fclose(fp_b_struct);
    if (forces) {
	fclose(fp_b_forces);
    }


    fprintf(fp, "\n");
    print_line_stars(fp);
}

/*****************************************************************************************
print_Bspline_basis_forces(): JFR - 07.22.12:  This function is the same as
print_Bspline_basis_output().  I just separated the printing of the forces for the different
Calculation Modes.
For the Bspline, we'll print out the calculated coefficients, but also need to print a
finer tabulated version of the actual force to use in calculations.
*****************************************************************************************/
void print_Bspline_basis_forces(FILE * fp, int forces, tW_system sys,
				tW_word tag, tW_Inter_Types * inter, int m)
{
    int j, i, l, grid, index;
    FILE *fp_f_forces, *fp_f_struct, *fp_f_forces_coeff,
	*fp_f_struct_coeff;
    double f_forces, f_struct;
    double R, B;
    int N_finepts = 10;		/* JFR - 07.22.12: For now, print out a grid 10 times finer than calculated */
    int k = inter->kspline;
    double Nik[MAX_BSPLINE_COEFF];
    double Npik[MAX_BSPLINE_COEFF];

    double ptr_x;
    double grid180_forces, grid180_struct;

    print_general_inter_output(fp, tag, inter, m);

    /* Open files for writing output. */
    fp_f_struct = open_output_file2("f_struct",
				    inter->inter_type, inter->inter_name,
				    tag);
    fp_f_struct_coeff =
	open_output_file2("f_struct_coeff", inter->inter_type,
			  inter->inter_name, tag);
    if (forces) {
	fp_f_forces = open_output_file2("f_forces",
					inter->inter_type,
					inter->inter_name, tag);
	fp_f_forces_coeff =
	    open_output_file2("f_forces_coeff", inter->inter_type,
			      inter->inter_name, tag);
    }

    /* Print results. */

    /* print the coefficients calculated */
    for (j = 0; j < inter->N_pts; j++) {
	ptr_x = inter->ptr_x[j];
	/* Account for angle and dihedral forces in degrees. */
	if ((strcmp(B_ANGLE, inter->inter_type) == 0)
	    || (strcmp(B_DIHEDRAL, inter->inter_type) == 0)) {
	    ptr_x *= 180.0 / M_PI;
	}

	fprintf(fp_f_struct_coeff, "%10.5lf    %15.5lf \n", ptr_x,
		inter->ptr_phi_struct[j] * (M_PI / 180.0));
	if (forces) {
	    fprintf(fp_f_forces_coeff, "%10.5lf    %15.5lf \n", ptr_x,
		    inter->ptr_phi_forces[j] * (M_PI / 180.0));
	}
    }

    /* Add the 180 grid point. Periodicity. */
    if ((strcmp(B_DIHEDRAL, inter->inter_type) == 0)
	&& (fabs(inter->ptr_x[0] * (180.0 / M_PI) + 180.0) < FLOAT_EPS)) {
	fprintf(fp_f_struct_coeff, "%10.5lf    %15.5lf \n", 180.0,
		inter->ptr_phi_struct[0] * (M_PI / 180.0));
	if (forces) {
	    fprintf(fp_f_forces_coeff, "%10.5lf    %15.5lf \n", 180.0,
		    inter->ptr_phi_forces[0] * (M_PI / 180.0));
	}
    }

    /* print the forces tabulated on a finer grid than the calculation */
    for (j = 0; j < inter->N_pts; j++) {

	index =
	    get_grid_index_for_linear_basis(inter->ptr_x[j] + FLOAT_EPS,
					    inter->i_0, inter->dr,
					    inter->R_min) - inter->i_0;

	for (l = 0; l < N_finepts; l++) {
	    f_forces = 0.00;
	    f_struct = 0.00;

	    R = inter->ptr_x[j] + l * (inter->dr * (1.0 / N_finepts));

	    if ((R > inter->ptr_x[inter->N_pts - 1])
		&& (strcmp(B_DIHEDRAL, inter->inter_type) != 0)) {
		continue;
	    }

	    B = calc_Bspline(inter->dr, inter->R_min, inter->i_0,
			     inter->N_pts, R, index, k, Nik, Npik);
	    for (i = 0; i < k; i++) {

		grid = j - i + k / 2;
		if (strcmp(B_DIHEDRAL, inter->inter_type) == 0) {
		    /* Account for periodicity */
		    if (grid < 0) {
			grid += inter->N_pts;
		    } else if (grid >= inter->N_pts) {
			grid -= inter->N_pts;
		    }
		}

		B = Nik[k - 1 - i];

		if ((grid >= 0) && ((grid) < inter->N_pts)) {
		    f_struct += B * inter->ptr_phi_struct[grid];
		    if (forces) {
			f_forces += B * inter->ptr_phi_forces[grid];
		    }
		}
	    }

	    if ((strcmp(B_ANGLE, inter->inter_type) == 0)
		|| (strcmp(B_DIHEDRAL, inter->inter_type) == 0)) {
		R *= 180.0 / M_PI;
		f_struct *= M_PI / 180.0;
		if (forces) {
		    f_forces *= M_PI / 180.0;
		}
	    }
	    fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", R, f_struct);
	    if (forces) {
		fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", R, f_forces);
	    }

	    if (strcmp(B_DIHEDRAL, inter->inter_type) == 0) {
		if (j == 0 && l == 0) {
		    grid180_struct = f_struct;
		    if (forces) {
			grid180_forces = f_forces;
		    }
		}
	    }

	}
    }

    /* Add the 180 grid point. Periodicity. */
    if ((strcmp(B_DIHEDRAL, inter->inter_type) == 0)
	&& (fabs(inter->ptr_x[0] * (180.0 / M_PI) + 180.0) < FLOAT_EPS)) {
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", 180.0,
		grid180_struct);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", 180.0,
		    grid180_forces);
	}
    }

    /* Close files for writing output. */
    fclose(fp_f_struct);
    fclose(fp_f_struct_coeff);
    if (forces) {
	fclose(fp_f_forces);
	fclose(fp_f_forces_coeff);
    }

    fprintf(fp, "\n");
    print_line_stars(fp);
}

/*****************************************************************************************
print_harmonic_basis_output(): Print the results for a harmonic basis.
*****************************************************************************************/
void print_harmonic_basis_output(FILE * fp, int forces, tW_system sys,
				 tW_word tag, tW_Inter_Types * inter,
				 int i)
{
    print_general_inter_output(fp, tag, inter, i);

    if (forces) {
	fprintf(fp, "        b[0]_forces: %lf\n", inter->ptr_b_forces[0]);
	fprintf(fp, "        b[1]_forces: %lf\n", inter->ptr_b_forces[1]);
    }
    fprintf(fp, "        b[0]_struct: %lf\n", inter->ptr_b_struct[0]);
    fprintf(fp, "        b[1]_struct: %lf\n", inter->ptr_b_struct[1]);

    fprintf(fp, "               g[0]: %f\n", inter->ptr_g[0]);
    fprintf(fp, "               g[1]: %f\n", inter->ptr_g[1]);

    fprintf(fp, "           g_cnt[0]: %f\n", inter->ptr_g_cnt[0]);
    fprintf(fp, "           g_cnt[1]: %f\n", inter->ptr_g_cnt[1]);

    fprintf(fp, "               L[0]: %f\n", inter->ptr_L[0]);
    fprintf(fp, "               L[1]: %f\n", inter->ptr_L[1]);

    fprintf(fp, "\n");
    print_line_stars(fp);

}

/*****************************************************************************************
print_harmonic_basis_forces():  JFR - added 04.06.12:  This function is the same as
print_harmonic_basis_output().  I just separated the printing of the forces for the different
Calculation Modes.
*****************************************************************************************/
void print_harmonic_basis_forces(FILE * fp, int forces, tW_system sys,
				 tW_word tag, tW_Inter_Types * inter,
				 int i)
{
    double k, r0;

    print_general_inter_output(fp, tag, inter, i);

    if (forces) {
	k = inter->ptr_phi_forces[0];
	r0 = inter->ptr_phi_forces[1] / k;
	if (strcmp(inter->inter_type, B_ANGLE) == 0) {
	    r0 *= 180.0 / M_PI;
	}

	fprintf(fp, "           k_forces: %lf \n", k);
	fprintf(fp, "          r0_forces: %lf \n", r0);
    }

    k = inter->ptr_phi_struct[0];
    r0 = inter->ptr_phi_struct[1] / k;
    if (strcmp(inter->inter_type, B_ANGLE) == 0) {
	r0 *= 180.0 / M_PI;
    }
    fprintf(fp, "           k_struct: %lf \n", k);
    fprintf(fp, "          r0_struct: %lf \n", r0);

    fprintf(fp, "\n");
    print_line_stars(fp);

    /* Also print a tabulated version of the force */
    int j;
    FILE *fp_f_forces, *fp_f_struct;
    double R_max;
    double R_min;
    double dr;
    double k_forces = inter->ptr_phi_forces[0];
    double r0_forces = inter->ptr_phi_forces[1];
    double k_struct = inter->ptr_phi_struct[0];
    double r0_struct = inter->ptr_phi_struct[1];
    double pi = M_PI;
    double conv = pi / 180.00;

    double r, F_struct, F_forces;
    int N;

    if (strcmp(inter->inter_type, B_ANGLE) == 0) {
	R_max = 180.0;
	R_min = 0.00;
	dr = 0.1;
	r0_forces /= conv;
	k_forces *= (conv * conv);
	r0_struct /= conv;
	k_struct *= (conv * conv);
    } else {
	R_max = 0.80;
	R_min = 0.00;
	dr = 0.0001;
    }

    N = (int) ceil((R_max - R_min) / dr);
    r = R_min;

    /* Open files for writing output. */
    fp_f_struct = open_output_file2("f_struct",
				    inter->inter_type, inter->inter_name,
				    tag);
    if (forces) {
	fp_f_forces = open_output_file2("f_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Print results. */
    for (j = 0; j < N + 1; j++) {
	F_struct = (-1.0) * (k_struct) * (r - r0_struct);
	F_forces = (-1.0) * (k_forces) * (r - r0_forces);
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", r, F_struct);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", r, F_forces);
	}
	r += dr;
    }

    fclose(fp_f_forces);
    fclose(fp_f_struct);

}

/*****************************************************************************************
print_RB_basis_output(): Print results for a RB basis.
*****************************************************************************************/
void print_RB_basis_output(FILE * fp, int forces, tW_system sys,
			   tW_word tag, tW_Inter_Types * inter, int i)
{
    print_general_inter_output(fp, tag, inter, i);

    if (forces) {
	fprintf(fp, "        b[0]_forces: %lf \n", inter->ptr_b_forces[0]);
	fprintf(fp, "        b[1]_forces: %lf \n", inter->ptr_b_forces[1]);
	fprintf(fp, "        b[2]_forces: %lf \n", inter->ptr_b_forces[2]);
	fprintf(fp, "        b[3]_forces: %lf \n", inter->ptr_b_forces[3]);
	fprintf(fp, "        b[4]_forces: %lf \n", inter->ptr_b_forces[4]);
    }

    fprintf(fp, "        b[0]_struct: %lf \n", inter->ptr_b_struct[0]);
    fprintf(fp, "        b[1]_struct: %lf \n", inter->ptr_b_struct[1]);
    fprintf(fp, "        b[2]_struct: %lf \n", inter->ptr_b_struct[2]);
    fprintf(fp, "        b[3]_struct: %lf \n", inter->ptr_b_struct[3]);
    fprintf(fp, "        b[4]_struct: %lf \n", inter->ptr_b_struct[4]);

    fprintf(fp, "               g[0]: %f\n", inter->ptr_g[0]);
    fprintf(fp, "               g[1]: %f\n", inter->ptr_g[1]);
    fprintf(fp, "               g[2]: %f\n", inter->ptr_g[2]);
    fprintf(fp, "               g[3]: %f\n", inter->ptr_g[3]);
    fprintf(fp, "               g[4]: %f\n", inter->ptr_g[4]);

    fprintf(fp, "           g_cnt[0]: %f\n", inter->ptr_g_cnt[0]);
    fprintf(fp, "           g_cnt[1]: %f\n", inter->ptr_g_cnt[1]);
    fprintf(fp, "           g_cnt[2]: %f\n", inter->ptr_g_cnt[2]);
    fprintf(fp, "           g_cnt[3]: %f\n", inter->ptr_g_cnt[3]);
    fprintf(fp, "           g_cnt[4]: %f\n", inter->ptr_g_cnt[4]);

    fprintf(fp, "               L[0]: %f\n", inter->ptr_L[0]);
    fprintf(fp, "               L[1]: %f\n", inter->ptr_L[1]);
    fprintf(fp, "               L[2]: %f\n", inter->ptr_L[2]);
    fprintf(fp, "               L[3]: %f\n", inter->ptr_L[3]);
    fprintf(fp, "               L[4]: %f\n", inter->ptr_L[4]);

    fprintf(fp, "\n");
    print_line_stars(fp);
}

/*****************************************************************************************
print_RB_basis_forces(): JFR - added 04.06.12:  This function is the same as
print_RB_basis_output().  I just separated the printing of the forces for the different
Calculation Modes.
*****************************************************************************************/
void print_RB_basis_forces(FILE * fp, int forces, tW_system sys,
			   tW_word tag, tW_Inter_Types * inter, int i)
{
    print_general_inter_output(fp, tag, inter, i);

    if (forces) {
	fprintf(fp, "          c1_forces: %lf \n",
		inter->ptr_phi_forces[0]);
	fprintf(fp, "          c2_forces: %lf \n",
		inter->ptr_phi_forces[1]);
	fprintf(fp, "          c3_forces: %lf \n",
		inter->ptr_phi_forces[2]);
	fprintf(fp, "          c4_forces: %lf \n",
		inter->ptr_phi_forces[3]);
	fprintf(fp, "          c5_forces: %lf \n",
		inter->ptr_phi_forces[4]);
    }

    fprintf(fp, "          c1_struct: %lf \n", inter->ptr_phi_struct[0]);
    fprintf(fp, "          c2_struct: %lf \n", inter->ptr_phi_struct[1]);
    fprintf(fp, "          c3_struct: %lf \n", inter->ptr_phi_struct[2]);
    fprintf(fp, "          c4_struct: %lf \n", inter->ptr_phi_struct[3]);
    fprintf(fp, "          c5_struct: %lf \n", inter->ptr_phi_struct[4]);

    fprintf(fp, "\n");
    print_line_stars(fp);

    /* Also print a tabulated version of the force */
    int j, k;
    FILE *fp_f_forces, *fp_f_struct;
    double R_max = 180.0;
    double R_min = -180.0;
    double dr = 0.5;
    int NRB = 5;
    double C_struct[NRB];
    double C_forces[NRB];
    for (k = 0; k < NRB; k++) {
	C_struct[k] = inter->ptr_phi_struct[k];
	C_forces[k] = inter->ptr_phi_forces[k];
    }
    double pi = M_PI;
    double conv = pi / 180.00;

    double psi, F_struct, F_forces;
    int N;

    N = (int) ceil((R_max - R_min) / dr);
    psi = R_min;

    /* Open files for writing output. */
    fp_f_struct = open_output_file2("f_struct",
				    inter->inter_type, inter->inter_name,
				    tag);
    if (forces) {
	fp_f_forces = open_output_file2("f_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Print results. */
    for (j = 0; j < N + 1; j++) {
	F_struct = 0.00;
	F_forces = 0.00;
	for (k = 0; k < NRB; k++) {
	    F_struct +=
		(conv) * (pow(-1.0, k + 1) * C_struct[k] *
			  ((double) k + 1.0) * pow(cos(psi * conv),
						   (double) k) * sin(psi *
								     conv));
	    F_forces +=
		(conv) * (pow(-1.0, k + 1) * C_forces[k] *
			  ((double) k + 1.0) * pow(cos(psi * conv),
						   (double) k) * sin(psi *
								     conv));
	}
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", psi, F_struct);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", psi, F_forces);
	}
	psi += dr;
    }

    fclose(fp_f_forces);
    fclose(fp_f_struct);

}

/*****************************************************************************************
print_TOY_basis_output(): Prints resutls for a TOY dihedral basis.
*****************************************************************************************/
void print_TOY_basis_output(FILE * fp, int forces, tW_system sys,
			    tW_word tag, tW_Inter_Types * inter, int i)
{
    print_general_inter_output(fp, tag, inter, i);

    if (forces) {
	fprintf(fp, "        b[0]_forces: %lf \n", inter->ptr_b_forces[0]);
	fprintf(fp, "        b[1]_forces: %lf \n", inter->ptr_b_forces[1]);
	fprintf(fp, "        b[2]_forces: %lf \n", inter->ptr_b_forces[2]);
    }

    fprintf(fp, "        b[0]_struct: %lf \n", inter->ptr_b_struct[0]);
    fprintf(fp, "        b[1]_struct: %lf \n", inter->ptr_b_struct[1]);
    fprintf(fp, "        b[2]_struct: %lf \n", inter->ptr_b_struct[2]);

    fprintf(fp, "               g[0]: %f\n", inter->ptr_g[0]);
    fprintf(fp, "               g[1]: %f\n", inter->ptr_g[1]);
    fprintf(fp, "               g[2]: %f\n", inter->ptr_g[2]);

    fprintf(fp, "           g_cnt[0]: %f\n", inter->ptr_g_cnt[0]);
    fprintf(fp, "           g_cnt[1]: %f\n", inter->ptr_g_cnt[1]);
    fprintf(fp, "           g_cnt[2]: %f\n", inter->ptr_g_cnt[2]);

    fprintf(fp, "               L[0]: %f\n", inter->ptr_L[0]);
    fprintf(fp, "               L[1]: %f\n", inter->ptr_L[1]);
    fprintf(fp, "               L[2]: %f\n", inter->ptr_L[2]);

    fprintf(fp, "\n");
    print_line_stars(fp);

}

/*****************************************************************************************
print_TOY_basis_forces(): JFR - added 04.06.12:  This function is the same as
print_TOY_basis_output().  I just separated the printing of the forces for the different
Calculation Modes.
*****************************************************************************************/
void print_TOY_basis_forces(FILE * fp, int forces, tW_system sys,
			    tW_word tag, tW_Inter_Types * inter, int i)
{
    print_general_inter_output(fp, tag, inter, i);

    if (forces) {
	fprintf(fp, "         A-B_forces: %lf \n",
		inter->ptr_phi_forces[0]);
	fprintf(fp, "           C_forces: %lf \n",
		inter->ptr_phi_forces[1]);
	fprintf(fp, "           D_forces: %lf \n",
		inter->ptr_phi_forces[2]);
    }

    fprintf(fp, "         A-B_struct: %lf \n", inter->ptr_phi_struct[0]);
    fprintf(fp, "           C_struct: %lf \n", inter->ptr_phi_struct[1]);
    fprintf(fp, "           D_struct: %lf \n", inter->ptr_phi_struct[2]);

    fprintf(fp, "\n");
    print_line_stars(fp);

    /* Also print a tabulated version of the force */
    int j;
    FILE *fp_f_forces, *fp_f_struct;
    double R_max = 180.0;
    double R_min = -180.0;
    double dr = 0.5;
    double AB_forces = inter->ptr_phi_forces[0];
    double C_forces = inter->ptr_phi_forces[1];
    double D_forces = inter->ptr_phi_forces[2];
    double AB_struct = inter->ptr_phi_struct[0];
    double C_struct = inter->ptr_phi_struct[1];
    double D_struct = inter->ptr_phi_struct[2];
    double pi = M_PI;
    double conv = pi / 180.00;

    double psi, F_struct, F_forces;
    int N;

    N = (int) ceil((R_max - R_min) / dr);
    psi = R_min;

    /* Open files for writing output. */
    fp_f_struct = open_output_file2("f_struct",
				    inter->inter_type, inter->inter_name,
				    tag);
    if (forces) {
	fp_f_forces = open_output_file2("f_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Print results. */
    for (j = 0; j < N + 1; j++) {
	F_struct =
	    (conv) * (AB_struct * sin(psi * conv) +
		      C_struct * sin(3 * psi * conv) +
		      D_struct * sin(psi * conv + pi / 4));
	F_forces =
	    (conv) * (AB_forces * sin(psi * conv) +
		      C_forces * sin(3 * psi * conv) +
		      D_forces * sin(psi * conv + pi / 4));
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", psi, F_struct);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", psi, F_forces);
	}
	psi += dr;
    }

    fclose(fp_f_forces);
    fclose(fp_f_struct);

}

/*****************************************************************************************
print_general_inter_output(): Prints general information about a particular interaction.
*****************************************************************************************/
void print_general_inter_output(FILE * fp, tW_word tag,
				tW_Inter_Types * inter, int i)
{

    fprintf(fp,
	    "  %3d. Inter_Name: %-10s    Inter_Type: %-20s    Basis: %-20s\n",
	    i + 1, inter->inter_name, inter->inter_type, inter->basis);

    fprintf(fp, "        System: %-20s\n", tag);

}


/*****************************************************************************************
print_power_basis_output(): Print results for a power basis.
*****************************************************************************************/
void print_power_basis_output(FILE * fp, int forces, tW_system sys,
			      tW_word tag, tW_Inter_Types * inter, int i)
{
    int i_coeff;
    int N_coeff = inter->N_powers;
    tW_word label;

    print_general_inter_output(fp, tag, inter, i);

    fprintf(fp, "    R_max: %f\n", inter->R_max);

    if (forces) {
	for (i_coeff = 0; i_coeff < N_coeff; i_coeff++) {
	    sprintf(label, "b_forces_C%d", inter->powers[i_coeff]);
	    fprintf(fp, "%19s: %e\n", label, inter->ptr_b_forces[i_coeff]);
	}
    }

    for (i_coeff = 0; i_coeff < N_coeff; i_coeff++) {
	sprintf(label, "b_struct_C%d", inter->powers[i_coeff]);
	fprintf(fp, "%19s: %e\n", label, inter->ptr_b_struct[i_coeff]);
    }

    for (i_coeff = 0; i_coeff < N_coeff; i_coeff++) {
	sprintf(label, "g_C%d", inter->powers[i_coeff]);
	fprintf(fp, "%19s: %f\n", label, inter->ptr_g[i_coeff]);
    }

    for (i_coeff = 0; i_coeff < N_coeff; i_coeff++) {
	sprintf(label, "g_cnt_C%d", inter->powers[i_coeff]);
	fprintf(fp, "%19s: %f\n", label, inter->ptr_g_cnt[i_coeff]);
    }

    for (i_coeff = 0; i_coeff < N_coeff; i_coeff++) {
	sprintf(label, "L_C%d", inter->powers[i_coeff]);
	fprintf(fp, "%19s: %f\n", label, inter->ptr_L[i_coeff]);
    }

    fprintf(fp, "\n");
    print_line_stars(fp);
}

/*****************************************************************************************
print_power_basis_forces():  JFR - added 04.06.12:  This function is the same as
print_power_basis_output().  I just separated the printing of the forces for the different
Calculation Modes.
*****************************************************************************************/
void print_power_basis_forces(FILE * fp, int forces, tW_system sys,
			      tW_word tag, tW_Inter_Types * inter, int i)
{
    int i_coeff;
    int N_coeff = inter->N_powers;
    tW_word label;

    print_general_inter_output(fp, tag, inter, i);

    fprintf(fp, "    R_max: %f\n", inter->R_max);

    if (forces) {
	for (i_coeff = 0; i_coeff < N_coeff; i_coeff++) {
	    sprintf(label, "C%d_forces", inter->powers[i_coeff]);
	    fprintf(fp, "%19s: %e\n", label,
		    inter->ptr_phi_forces[i_coeff]);
	}
    }

    for (i_coeff = 0; i_coeff < N_coeff; i_coeff++) {
	sprintf(label, "C%d_struct", inter->powers[i_coeff]);
	fprintf(fp, "%19s: %e\n", label, inter->ptr_phi_struct[i_coeff]);
    }

    fprintf(fp, "\n");
    print_line_stars(fp);

    /* Also print a tabulated version of the force */
    int j, k;
    FILE *fp_f_forces, *fp_f_struct;
    //double R_max = inter->R_max;
    //double R_min = inter->R_min;
    double R_max = 3.00;
    double R_min = 0.05;
    //double dr    = inter->dr;
    double dr = 0.002;
    double C_struct[N_coeff];
    double C_forces[N_coeff];
    for (j = 0; j < N_coeff; j++) {
	C_struct[j] = inter->ptr_phi_struct[j];
	C_forces[j] = inter->ptr_phi_forces[j];
    }

    double r, F_struct, F_forces;
    int N;

    N = (int) ceil((R_max - R_min) / dr);
    r = R_min;

    /* Open files for writing output. */
    fp_f_struct = open_output_file2("f_struct",
				    inter->inter_type, inter->inter_name,
				    tag);
    if (forces) {
	fp_f_forces = open_output_file2("f_forces",
					inter->inter_type,
					inter->inter_name, tag);
    }

    /* Print results. */
    for (k = 0; k < N + 1; k++) {
	F_struct = 0.00;
	F_forces = 0.00;
	for (j = 0; j < N_coeff; j++) {
	    F_struct +=
		C_struct[j] * (inter->powers[j] + 1.0) * pow(r,
							     -1.0 *
							     inter->
							     powers[j]);
	    F_forces +=
		C_forces[j] * (inter->powers[j] + 1.0) * pow(r,
							     -1.0 *
							     inter->
							     powers[j]);
	}
	fprintf(fp_f_struct, "%10.5lf    %15.5lf \n", r, F_struct);
	if (forces) {
	    fprintf(fp_f_forces, "%10.5lf    %15.5lf \n", r, F_forces);
	}
	r += dr;
    }

    fclose(fp_f_forces);
    fclose(fp_f_struct);

}

/*****************************************************************************************
was_inter_sampled(): Returns FALSE if the interaction was not present in calculation and
TRUE otherwise.
*****************************************************************************************/
int was_inter_present(double *g_cnt, int N_coeff)
{
    int i;
    int flag = FALSE;

    for (i = 0; i < N_coeff; i++) {
	if (fabs(g_cnt[i]) > FLOAT_EPS) {
	    flag = TRUE;
	    break;
	}
    }

    return flag;
}

/*****************************************************************************************
print_M_matrix2(): Print elements of the M matrix.  JFR - added 09.01.10 
*****************************************************************************************/
void print_M_matrix2(tW_system * sys)
{
    int i, j, ki, kj;
    int row_i, col_j;
    int i_flag, j_flag;
    double ri, rj;
    tW_Inter_Types *inter_i, *inter_j;
    tW_word fname;
    FILE *fp_M, *fp_Mcnt, *fp_M2, *fp_M2cnt, *fp_M3, *fp_T, *fp_T2, *fp_T3,
	*fp_Tcnt;
    int diag_flag = FALSE;
    int index;
    double normalized;

    fp_T = fopen("T.dat", "w");	//Whole matrix without diagonals

    fp_T2 = fopen("T2.dat", "w");	//whole matrix with diagonal elements

    fp_T3 = fopen("T3.dat", "w");	//whole matrix with diagonal elements, but NOT the pair correlation

    fp_Tcnt = fopen("Tcnt.dat", "w");	//whole matrix unweighted (i.e., pure distributions)

    for (i = 0; i < sys->N_Inter_Types; i++) {
	inter_i = &(sys->Inter_Types[i]);

	i_flag = strcmp(NB_PAIR, inter_i->inter_type);

	//for ( j=0; j<sys->N_Inter_Types; j++ )
	for (j = 0; j <= i; j++) {	/* Only print out the lower diagonal of the matrix */
	    inter_j = &(sys->Inter_Types[j]);

	    j_flag = strcmp(NB_PAIR, inter_j->inter_type);

	    sprintf(fname, "M.%s.%s.dat", inter_i->inter_name,
		    inter_j->inter_name);
	    fp_M = fopen(fname, "w");

	    sprintf(fname, "Mcnt.%s.%s.dat", inter_i->inter_name,
		    inter_j->inter_name);
	    fp_Mcnt = fopen(fname, "w");

	    sprintf(fname, "M2.%s.%s.dat", inter_i->inter_name,
		    inter_j->inter_name);
	    fp_M2 = fopen(fname, "w");	//M2 will be the matrix with diagonal elements

	    sprintf(fname, "M2cnt.%s.%s.dat", inter_i->inter_name,
		    inter_j->inter_name);
	    fp_M2cnt = fopen(fname, "w");	//M2 will be the matrix with diagonal elements

	    sprintf(fname, "M3.%s.%s.dat", inter_i->inter_name,
		    inter_j->inter_name);
	    fp_M3 = fopen(fname, "w");	//M3 will be the matrix with diagonal elements, but NOT the pair correlation

	    ri = inter_i->R_min;	// if before trim_Mb
//      ri = *inter_i->ptr_x; // if after trim_Mb
	    for (ki = 0; ki < inter_i->N_coeff; ki++) {
		if (fabs(ri) < FLOAT_EPS) {
		    ri += inter_i->dr;
		    continue;
		}
		/* JFR - avoid dividing through by zero */
		row_i = inter_i->i_0 + ki;

		rj = inter_j->R_min;	// if before trim_Mb
//        rj = *inter_j->ptr_x; // if after trim_Mb
		for (kj = 0; kj < inter_j->N_coeff; kj++) {
		    if (fabs(rj) < FLOAT_EPS) {
			rj += inter_j->dr;
			continue;
		    }
		    /* JFR - avoid dividing through by zero */
		    col_j = inter_j->i_0 + kj;

		    if (i_flag == 0) {
			if (j_flag == 0) {
			    normalized = ri * ri * rj * rj;
			} else {
			    normalized = ri * ri;
			}
		    } else if (j_flag == 0) {
			normalized = rj * rj;
		    } else {
			normalized = 1.0;
		    }

		    if (sys->MT_var.flag_norm == FALSE) {
			normalized = 1.0;
		    }		// JFR - added 04.06.12:  don't normalize if the user says not to

		    index = index_Lpacked(row_i, col_j, sys->N_coeff);

		    if ((inter_i->i_basis == DELTA_BASIS_INDEX)
			&& (inter_j->i_basis == DELTA_BASIS_INDEX)) {

			if (fabs(sys->M[index]) > FLOAT_EPS) {

			    if (row_i == col_j) {
				fprintf(fp_M2, "%.6e  %.6e  %.6e\n", ri,
					rj, sys->M[index] / normalized);
				if (sys->M_cnt != NULL)	//JFR - added 04.13.12: check if in LOWMEM mode
				{
				    fprintf(fp_M2cnt, "%.6e  %.6e  %.6e\n",
					    ri, rj,
					    sys->M_cnt[index] /
					    normalized);
				}
				fprintf(fp_T2, "%d  %d  %.6e\n", row_i,
					col_j, sys->M[index] / normalized);
			    } else {
				fprintf(fp_M, "%.6e  %.6e  %.6e\n", ri, rj,
					sys->M[index] / normalized);
				if (sys->M_cnt != NULL)	//JFR - added 04.13.12: check if in LOWMEM mode
				{
				    fprintf(fp_Mcnt, "%.6e  %.6e  %.6e\n",
					    ri, rj,
					    sys->M_cnt[index] /
					    normalized);
				    fprintf(fp_Tcnt, "%d  %d  %.6e\n",
					    row_i, col_j,
					    sys->M_cnt[index] /
					    normalized);
				    fprintf(fp_M2cnt, "%.6e  %.6e  %.6e\n",
					    ri, rj,
					    sys->M_cnt[index] /
					    normalized);
				}
				fprintf(fp_T, "%d  %d  %.6e\n", row_i,
					col_j, sys->M[index] / normalized);
				fprintf(fp_M2, "%.6e  %.6e  %.6e\n", ri,
					rj, sys->M[index] / normalized);
				fprintf(fp_T2, "%d  %d  %.6e\n", row_i,
					col_j, sys->M[index] / normalized);
			    }
			}

			if (row_i == col_j) {
			    if (fabs(sys->M[index] - sys->g[row_i]) >
				FLOAT_EPS) {
				fprintf(fp_M3, "%.6e  %.6e  %.6e\n", ri,
					rj,
					(sys->M[index] -
					 sys->g[row_i]) / normalized);
				fprintf(fp_T3, "%d  %d  %.6e\n", row_i,
					col_j,
					(sys->M[index] -
					 sys->g[row_i]) / normalized);
			    }
			} else {
			    if (fabs(sys->M[index]) > FLOAT_EPS) {
				fprintf(fp_M3, "%.6e  %.6e  %.6e\n", ri,
					rj, sys->M[index] / normalized);
				fprintf(fp_T3, "%d  %d  %.6e\n", row_i,
					col_j, sys->M[index] / normalized);
			    }
			}
		    } else if (inter_i->i_basis == LINEAR_BASIS_INDEX
			       && inter_j->i_basis == LINEAR_BASIS_INDEX) {

			if ((row_i == col_j) || (row_i == col_j + 1)
			    || (row_i == col_j - 1)) {
			    diag_flag = TRUE;
			}

			if (fabs(sys->M[index]) > FLOAT_EPS) {

			    if (diag_flag == TRUE) {
				fprintf(fp_M2, "%.6e  %.6e  %.6e\n", ri,
					rj, sys->M[index] / normalized);
				if (sys->M_cnt != NULL) {
				    sys->M_cnt[index] += 1.0;
				}	//JFR - added 04.13.12: check if in LOWMEM mode
				{
				    fprintf(fp_M2cnt, "%.6e  %.6e  %.6e\n",
					    ri, rj,
					    sys->M_cnt[index] /
					    normalized);
				}
				fprintf(fp_T2, "%d  %d  %.6e\n", row_i,
					col_j, sys->M[index] / normalized);
			    } else {
				fprintf(fp_M, "%.6e  %.6e  %.6e\n", ri, rj,
					sys->M[index] / normalized);
				if (sys->M_cnt != NULL) {
				    sys->M_cnt[index] += 1.0;
				}	//JFR - added 04.13.12: check if in LOWMEM mode
				{
				    fprintf(fp_Mcnt, "%.6e  %.6e  %.6e\n",
					    ri, rj,
					    sys->M_cnt[index] /
					    normalized);
				    fprintf(fp_Tcnt, "%d  %d  %.6e\n",
					    row_i, col_j,
					    sys->M_cnt[index] /
					    normalized);
				    fprintf(fp_M2cnt, "%.6e  %.6e  %.6e\n",
					    ri, rj,
					    sys->M_cnt[index] /
					    normalized);
				}
				fprintf(fp_T, "%d  %d  %.6e\n", row_i,
					col_j, sys->M[index] / normalized);
				fprintf(fp_M2, "%.6e  %.6e  %.6e\n", ri,
					rj, sys->M[index] / normalized);
				fprintf(fp_T2, "%d  %d  %.6e\n", row_i,
					col_j, sys->M[index] / normalized);
			    }
			}

			if (row_i == col_j) {
			    if (fabs(sys->M[index] - sys->g[row_i]) >
				FLOAT_EPS) {
				fprintf(fp_M3, "%.6e  %.6e  %.6e\n", ri,
					rj,
					(sys->M[index] -
					 sys->g[row_i]) / normalized);
				fprintf(fp_T3, "%d  %d  %.6e\n", row_i,
					col_j,
					(sys->M[index] -
					 sys->g[row_i]) / normalized);
			    }
			} else if ((row_i == col_j + 1)
				   || (row_i == col_j - 1)) {
			    if (fabs(sys->M[index] - sys->g[row_i]) >
				FLOAT_EPS) {
				fprintf(fp_M3, "%.6e  %.6e  %.6e\n", ri,
					rj,
					(sys->M[index] -
					 0.5 * sys->g[row_i] -
					 0.5 * sys->g[col_j]) /
					normalized);
				fprintf(fp_T3, "%d  %d  %.6e\n", row_i,
					col_j,
					(sys->M[index] -
					 0.5 * sys->g[row_i] -
					 0.5 * sys->g[col_j]) /
					normalized);
			    }
			} else {
			    if (fabs(sys->M[index]) > FLOAT_EPS) {
				fprintf(fp_M3, "%.6e  %.6e  %.6e\n", ri,
					rj, sys->M[index] / normalized);
				fprintf(fp_T3, "%d  %d  %.6e\n", row_i,
					col_j, sys->M[index] / normalized);
			    }
			}

			diag_flag = FALSE;
		    }

		    rj += inter_j->dr;
		}
		ri += inter_i->dr;
	    }

	    fclose(fp_M);
	    fclose(fp_Mcnt);
	    fclose(fp_M2);
	    fclose(fp_M2cnt);
	    fclose(fp_M3);
	}
    }

    fclose(fp_T);
    fclose(fp_T2);
    fclose(fp_T3);
    fclose(fp_Tcnt);
}

/*****************************************************************************************
print_save_state(): Prints high precision values of M,b,... for resuming a calculation
after the loop over particles.  JFR - added 04.11.12 
*****************************************************************************************/
void print_save_state(FILE * fp_log, tW_system * sys)
{
    int i;
    tW_word fname;
    FILE *fp_b_struct, *fp_b_forces, *fp_g, *fp_g_cnt, *fp_L, *fp_M,
	*fp_M2, *fp_M_cnt, *fp_d2b, *fp_d2M, *fp_FF, *fp_Nf, *fp_rescale;

    sprintf(fname, "%s", "save.b_struct.dat");
    fp_b_struct = fopen(fname, "w");

    sprintf(fname, "%s", "save.b_forces.dat");
    fp_b_forces = fopen(fname, "w");

    sprintf(fname, "%s", "save.g.dat");
    fp_g = fopen(fname, "w");

    sprintf(fname, "%s", "save.g_cnt.dat");
    fp_g_cnt = fopen(fname, "w");

    sprintf(fname, "%s", "save.L.dat");
    fp_L = fopen(fname, "w");

    sprintf(fname, "%s", "save.M.dat");
    fp_M = fopen(fname, "w");

    sprintf(fname, "%s", "save.M2.dat");
    fp_M2 = fopen(fname, "w");

    if (sys->MT_var.flag_Mcnt == TRUE) {
	sprintf(fname, "%s", "save.M_cnt.dat");
	fp_M_cnt = fopen(fname, "w");
    }
    //if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 )
    //{
    sprintf(fname, "%s", "save.d2b.dat");
    fp_d2b = fopen(fname, "w");
    //}

    //if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 )
    //{
    sprintf(fname, "%s", "save.d2M.dat");
    fp_d2M = fopen(fname, "w");
    //}

    sprintf(fname, "%s", "save.FF.dat");
    fp_FF = fopen(fname, "w");

    sprintf(fname, "%s", "save.Nf.dat");
    fp_Nf = fopen(fname, "w");

    sprintf(fname, "%s", "save.rescale.dat");
    fp_rescale = fopen(fname, "w");

    for (i = 0; i < sys->N_coeff; i++) {
	fprintf(fp_b_struct, "%.15lf \n", sys->b_struct[i]);
	fprintf(fp_b_forces, "%.15lf \n", sys->b_forces[i]);
	fprintf(fp_g, "%.15lf \n", sys->g[i]);
	fprintf(fp_g_cnt, "%.15lf \n", sys->g_cnt[i]);
	fprintf(fp_L, "%.15lf \n", sys->L[i]);
	fprintf(fp_d2b, "%.15lf \n", sys->d2b[i]);
	//*if ( strcmp( sys->PC_var.LPC, "bvar" ) == 0 ) {*/ fprintf( fp_d2b, "%.15lf \n", sys->d2b[i] ); /*}*/
	fprintf(fp_rescale, "%.15lf \n", sys->rescale[i]);
    }

    for (i = 0; i < sys->N_pack; i++) {
	fprintf(fp_M, "%.15lf \n", sys->M[i]);
	fprintf(fp_M2, "%.15lf \n", sys->M2[i]);
	if (sys->MT_var.flag_Mcnt == TRUE) {
	    fprintf(fp_M_cnt, "%.15lf \n", sys->M_cnt[i]);
	}
	fprintf(fp_d2M, "%.15lf \n", sys->d2M[i]);
	//*if ( strcmp( sys->PC_var.RPC, "MTvar" ) == 0 ) {*/ fprintf( fp_d2M, "%.15lf \n", sys->d2M[i] ); /*}*/
    }

    fprintf(fp_FF, "%.15lf \n", sys->Chi2);

    fprintf(fp_Nf, "%d \n", sys->REG_var.Nframes);

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

}

/*****************************************************************************************
print_Chi2(): Calc the CG-dependent part of Chi2 and then print the final result.  
JFR - added 06.16.12 
*****************************************************************************************/
void print_Chi2(FILE * fp, tW_system * sys)
{
    int i, j;
    double b[sys->N_coeff];
    double ff = sys->Chi2;
    double neg2bphi;
    int index;

    /* Get force vector */
    int l_inp, test_sscanf;
    tW_line inp_line;
    double force[sys->N_coeff];

    if (sys->CHISQD_var.flag_CHISQD == TRUE) {
	FILE *fp_force;

	if (!file_exists(sys->CHISQD_var.force_fnm)) {
	    printf
		("ERROR: Cannot open file '%s' in print_decomp_Chi2().\n",
		 sys->CHISQD_var.force_fnm);
	    exit(EXIT_FAILURE);
	}
	fp_force = fopen(sys->CHISQD_var.force_fnm, "r");

	for (i = 0; i < sys->N_coeff; i++) {
	    l_inp = get_next_line(fp_force, inp_line);
	    test_sscanf = sscanf(inp_line, "%lf", &force[i]);
	}

	fclose(fp_force);
    } else {
	for (i = 0; i < sys->N_coeff; i++) {
	    force[i] = sys->phi_forces[i];
	}
    }
    /* End Get force vector */

    fprintf(fp, "ff = %lf \n", ff);

    for (i = 0; i < sys->N_coeff; i++) {
	b[i] = 0.00;
    }

    for (i = 0; i < sys->N_coeff; i++) {
	sys->Chi2 -= 2.0 * sys->b_forces[i] * force[i];

	for (j = 0; j < sys->N_coeff; j++) {
	    index = index_Lpacked(i, j, sys->N_coeff);
	    b[i] += sys->M[index] * force[j];
	}
    }

    neg2bphi = sys->Chi2 - ff;
    fprintf(fp, "-2bphi = %lf \n", neg2bphi);

    for (i = 0; i < sys->N_coeff; i++) {
	sys->Chi2 += force[i] * b[i];
    }

    fprintf(fp, "phiGphi = %lf \n", sys->Chi2 - neg2bphi - ff);

    fprintf(fp, "Chi2 = %lf \n", sys->Chi2);

}

/*****************************************************************************************
print_decomp_Chi2(): Calc the CG-dependent part of Chi2 and then print the final result.  
Decompose into various contributions JFR - added 09.10.12
*****************************************************************************************/
void print_decomp_Chi2(FILE * fp, tW_system * sys)
{
    int i, j, ki, kj;
    int row_i, col_j;
    double ri, rj;
    tW_Inter_Types *inter_i, *inter_j;
    int index;

    /* Get force vector */
    int l_inp, test_sscanf;
    tW_line inp_line;
    double force[sys->N_coeff];
    if (sys->CHISQD_var.flag_CHISQD == TRUE) {
	FILE *fp_force;

	if (!file_exists(sys->CHISQD_var.force_fnm)) {
	    printf
		("ERROR: Cannot open file '%s' in print_decomp_Chi2().\n",
		 sys->CHISQD_var.force_fnm);
	    exit(EXIT_FAILURE);
	}
	fp_force = fopen(sys->CHISQD_var.force_fnm, "r");


	for (i = 0; i < sys->N_coeff; i++) {
	    l_inp = get_next_line(fp_force, inp_line);
	    test_sscanf = sscanf(inp_line, "%lf", &force[i]);
	}

	fclose(fp_force);
    } else {
	for (i = 0; i < sys->N_coeff; i++) {
	    force[i] = sys->phi_forces[i];
	}
    }
    /* End Get force vector */

    double neg2bphi = 0.00;
    double phiGphi = 0.00;

    fprintf(fp, "Printing Chi2 decomposition . . . \n");

    for (i = 0; i < sys->N_Inter_Types; i++) {
	inter_i = &(sys->Inter_Types[i]);

	neg2bphi = 0.00;
	for (ki = 0; ki < inter_i->N_coeff; ki++) {
	    row_i = inter_i->i_0 + ki;
	    neg2bphi -= 2.0 * sys->b_forces[row_i] * force[row_i];
	}
	fprintf(fp, "-2bphi(%s) = %lf \n", inter_i->inter_name, neg2bphi);

	for (j = 0; j <= i; j++) {	/* Only print out the lower diagonal of the matrix */
	    inter_j = &(sys->Inter_Types[j]);

	    phiGphi = 0.00;
//      ri = inter_i->R_min; // if before trim_Mb
	    ri = *inter_i->ptr_x;	// if after trim_Mb
	    for (ki = 0; ki < inter_i->N_coeff; ki++) {
		row_i = inter_i->i_0 + ki;

//        rj = inter_j->R_min; // if before trim_Mb
		rj = *inter_j->ptr_x;	// if after trim_Mb
		for (kj = 0; kj < inter_j->N_coeff; kj++) {
		    col_j = inter_j->i_0 + kj;

		    index = index_Lpacked(row_i, col_j, sys->N_coeff);

		    phiGphi += force[row_i] * sys->M[index] * force[col_j];

		    rj += inter_j->dr;
		}
		ri += inter_i->dr;
	    }
	    fprintf(fp, "phiGphi(%s,%s) = %lf \n", inter_i->inter_name,
		    inter_j->inter_name, phiGphi);
	}
    }

}

/*****************************************************************************************
print_sep_forces(): Separate the force vector into different interactions before printing
 JFR - added 10.05.12
*****************************************************************************************/

void print_sep_forces(tW_word fname, double *phi, tW_system * sys)
{
    int i, j;

    FILE *fp;
    tW_word fname2;

    //JFR - Separate the interactions for output
    for (i = 0; i < sys->N_Inter_Types; i++) {
	sprintf(fname2, "%s.%s.%s.%s.dat", fname,
		sys->Inter_Types[i].inter_type, "total",
		sys->Inter_Types[i].inter_name);
	fp = fopen(fname2, "w");

	for (j = 0; j < sys->Inter_Types[i].N_pts; j++) {
	    fprintf(fp, "%10.5lf    %15.5lf \n",
		    sys->Inter_Types[i].ptr_x[j],
		    phi[sys->Inter_Types[i].i_0 + j]);
	}

	fclose(fp);
    }

}

/*****************************************************************************************
print_sep_M(): Separate the metric tensor into different interaction blocks before printing
 JFR - added 10.05.12
*****************************************************************************************/

void print_sep_M(tW_word fname, double **M, tW_system * sys)
{
    int i, j, ki, kj;
    int i_flag, j_flag;
    int row_i, col_j;
    double ri, rj;
    double normalized;
    FILE *fp;
    tW_word fname2;

    //JFR - Separate the interactions for output
    for (i = 0; i < sys->N_Inter_Types; i++) {
	i_flag = strcmp(NB_PAIR, sys->Inter_Types[i].inter_type);
	for (j = 0; j <= i; j++) {
	    j_flag = strcmp(NB_PAIR, sys->Inter_Types[j].inter_type);
	    sprintf(fname2, "%s.%s.%s.dat", fname,
		    sys->Inter_Types[i].inter_name,
		    sys->Inter_Types[j].inter_name);
	    fp = fopen(fname2, "w");

	    for (ki = 0; ki < sys->Inter_Types[i].N_pts; ki++) {
		ri = sys->Inter_Types[i].ptr_x[ki];
		row_i = sys->Inter_Types[i].i_0 + ki;
		for (kj = 0; kj < sys->Inter_Types[j].N_pts; kj++) {
		    rj = sys->Inter_Types[j].ptr_x[kj];
		    col_j = sys->Inter_Types[j].i_0 + kj;

		    if (i_flag == 0) {
			if (j_flag == 0) {
			    normalized = ri * ri * rj * rj;
			} else {
			    normalized = ri * ri;
			}
		    } else if (j_flag == 0) {
			normalized = rj * rj;
		    } else {
			normalized = 1.0;
		    }

		    if (sys->MT_var.flag_norm == FALSE) {
			normalized = 1.0;
		    }

		    if ( /* (fabs(M[i][j]) > FLOAT_EPS) && */
			(row_i != col_j)) {
			fprintf(fp, "%10.5lf    %10.5lf    %15.5lf \n", ri,
				rj, M[row_i][col_j]);
		    }
		}
	    }

	    fclose(fp);
	}
    }

}
