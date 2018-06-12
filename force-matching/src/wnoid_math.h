/**
@file wnoid_math.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef NOID_MATH
#define NOID_MATH

#include "cgff_types.h"

/****************************************************************************************/

int copy_vector(dvec source, dvec dest);

int copy_matrix(matrix source, matrix dest);

double dot_prod(dvec v1, dvec v2);

double calc_norm(dvec v);

int vect_diff(dvec v1, dvec v2, dvec diff);

int vect_sum(dvec v1, dvec v2, dvec sum);

void vect_inc( dvec v1, dvec v2 );

int scal_times_vect(double sc, dvec vect, dvec result);

int cross_product(dvec v1, dvec v2, dvec v_cross);

int calc_centered_diff(int Nmax, double dx, double *X, double *dX_dx);

int calc_running_avg(int Nmax, int Nlength, double *X, double *smX);

int get_difference_unit_vector(dvec x_i, dvec x_j, dvec x_ij, dvec u_ij, double *r_ij);

int find_max_int(int i, int j);

int compare_ints(const void *p, const void *q);

int calc_running_avg_periodic(int Nmax, int Nlength, double *X, double *smX);

int index_Lpacked(int i, int j, int N);

/****************************************************************************************/

#endif
