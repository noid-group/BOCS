/**
@file io_output.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef IO_OUT
#define IO_OUT

#ifdef __cplusplus
extern "C"
{
#endif

/****************************************************************************************/
/*Functions called outside of io_output.c************************************************/
/****************************************************************************************/

int get_top_tag(tW_word coord_fnm, tW_word tag);

int print_output(int forces, tW_system sys, tW_word tag);

// JFR - added 04.11.12:  I got rid of this function, since I have been using the second
//version for some time now and I don't see any reason to update this so that the matrix
//is in packed form.  See older versions if you are interested.
//void print_M_matrix( tW_system *sys );

void print_M_matrix2(tW_system * sys);	/* JFR - added 09.01.10 */

void print_Chi2(FILE * fp, tW_system * sys);	/* JFR - added 06.16.12: Chi2 */

void print_decomp_Chi2(FILE * fp, tW_system * sys);	/* JFR - added 09.10.12: decomp Chi2 */

void print_sep_forces(tW_word fname, double *phi, tW_system * sys);	/* JFR - added 10.05.12: for solve_PT */

void print_sep_M(tW_word fname, double **M, tW_system * sys);	/* JFR - added 10.05.12: for solve_PT */

/****************************************************************************************/
/* Functions associated with print_output()**********************************************/
/****************************************************************************************/

FILE *open_output_file1(tW_word tag);

FILE *open_output_file2(tW_word function, tW_word inter_type,
			tW_word inter_name, tW_word tag);

void print_general_inter_output(FILE * fp, tW_word tag,
				tW_Inter_Types * inter, int i);

void print_delta_basis_output(FILE * fp, int forces, tW_system sys,
			      tW_word tag, tW_Inter_Types * inter, int i);

void print_delta_basis_forces(FILE * fp, int forces, tW_system sys, tW_word tag, tW_Inter_Types * inter, int i);	// JFR - added 04.06.12: Separate force printing for different CalcMODES

/* START JFR */
void print_linear_basis_output(FILE * fp, int forces, tW_system sys,
			       tW_word tag, tW_Inter_Types * inter, int i);
/* END JFR */

void print_linear_basis_forces(FILE * fp, int forces, tW_system sys, tW_word tag, tW_Inter_Types * inter, int i);	// JFR - added 04.06.12: Separate force printing for different CalcMODES

/* JFR - 07.22.12 */
void print_Bspline_basis_output(FILE * fp, int forces, tW_system sys,
				tW_word tag, tW_Inter_Types * inter,
				int m);

/* JFR - 07.22.12 */
void print_Bspline_basis_forces(FILE * fp, int forces, tW_system sys,
				tW_word tag, tW_Inter_Types * inter,
				int m);

void print_harmonic_basis_output(FILE * fp, int forces, tW_system sys,
				 tW_word tag, tW_Inter_Types * inter,
				 int i);

void print_harmonic_basis_forces(FILE * fp, int forces, tW_system sys, tW_word tag, tW_Inter_Types * inter, int i);	// JFR - added 04.06.12: Separate force printing for different CalcMODES

void print_RB_basis_output(FILE * fp, int forces, tW_system sys,
			   tW_word tag, tW_Inter_Types * inter, int i);

void print_RB_basis_forces(FILE * fp, int forces, tW_system sys, tW_word tag, tW_Inter_Types * inter, int i);	// JFR - added 04.06.12: Separate force printing for different CalcMODES

void print_TOY_basis_output(FILE * fp, int forces, tW_system sys,
			    tW_word tag, tW_Inter_Types * inter, int i);

void print_TOY_basis_forces(FILE * fp, int forces, tW_system sys, tW_word tag, tW_Inter_Types * inter, int i);	// JFR - added 04.06.12: Separate force printing for different CalcMODES

void print_power_basis_output(FILE * fp, int forces, tW_system sys,
			      tW_word tag, tW_Inter_Types * inter, int i);

void print_power_basis_forces(FILE * fp, int forces, tW_system sys, tW_word tag, tW_Inter_Types * inter, int i);	// JFR - added 04.06.12: Separate force printing for different CalcMODES

int was_inter_present(double *g_cnt, int N_coeff);


/****************************************************************************************/
/* Functions associated with save state**********************************************/
/****************************************************************************************/
void print_save_state(FILE * fp_log, tW_system * sys);

#ifdef __cplusplus
}
#endif

#endif
