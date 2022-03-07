/**
@file table_functions.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

# ifndef table_functions
# define table_functions

#ifdef __cplusplus
extern "C"
{
#endif

# include <math.h>

int get_inter_type( char *argv[], int argc );

int get_basis_type( char *argv[] );

void get_table_parameters( Input *input, char *argv[] );

void read_input_forces( int n_pts, Arrays *arrays, FILE *fp_in );

int get_potential( Input *input, Arrays *arrays );

int get_nb_potential( Input *input, Arrays *arrays );

int get_bond_potential( Input *input, Arrays *arrays );

int get_angle_potential( Input *input, Arrays *arrays );

int get_dihedral_potential( Input *input, Arrays *arrays );

int get_total_n_pts( int n_pts, double *r, int *total_n_pts, double *r_min );

void copy_array( int n_pts, double *array, double *copy );

void fill_r_nb( int total_n_pts, int flag, double *r, double *r_temp, double dr, double r_min );

void fill_f_nb( int total_n_pts, int flag, double *f, double *f_temp, double slope, double dr );

void scal_vect_mult( int total_n_pts, double scalar, double original[], double *result );

int integrate_trapezoid (double x[], double *int_y, double y[], int no_of_points);

void shift_potential( int total_n_pts, double *u );

void fill_r_bond( int total_n_pts, double *r, double *r_temp, Input *input );

int get_total_n_pts_bond( int n_pts, double *r_min, double *r_max, double dr, double *r );

void fill_f_bond( int total_n_pts, double *f, double *f_temp, double *r_temp, double dr, Input *input );

void shift_potential_bond( int total_n_pts, double *u );

void check_angle_min_max( double r_min, double r_max );

void print_table( Input input, Arrays arrays, FILE *fp_out );

void print_nb_table( Input input, Arrays arrays, FILE *fp_out );

void print_bond_table( Input input, Arrays arrays, FILE *fp_out );

void check_dihedral_endpoints( int n_pts, double *r );

/***************************************************************************************************/
/*Error Messages************************************************************************************/
/***************************************************************************************************/
void print_type_error_message( char *inter_name );

void print_type_index_error( int type_index );

void print_read_input_error1( int i );

void print_read_input_error2( int i );

void print_error_fill_r( double r_test, double r_temp );

void print_fill_r_bond_error( double r_test, double r_max, int total_n_pts );

void print_check_angle_min_max_error( double r_test, double r_exact, char *r_test_name );

void print_check_dihedral_endpoints_error( double r_test, double r_exact, char *r_test_name );

/************************ JFR - 02.27.13: new functions ******************************************/
int edit_forces( Input *input, Arrays *arrays, Arrays *smarrays );

int edit_nb_force( Input *input, Arrays *arrays, Arrays *smarrays );

int edit_bond_force( Input *input, Arrays *arrays, Arrays *smarrays );

int edit_angle_force( Input *input, Arrays *arrays, Arrays *tmparrays );

int smooth_forces( int n_pts, Input *input, Arrays *arrays );

int smooth_forces_periodic( int n_pts, Input *input, Arrays *arrays );

int interpolate_forces( Input *input, Arrays *arrays, Arrays *tmparrays );

int edit_intra_nb_force( Input *input, Arrays *arrays, Arrays *tmparrays );

#ifdef __cplusplus
}
#endif

# endif
