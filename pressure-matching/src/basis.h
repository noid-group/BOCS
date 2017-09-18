/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file basis.h
@author N. Dunn
@date Apr 30, 2014
@brief Functions for pressure matching basis representations

This file contains headers for portable basis functions used to represent functions on a grid
or as a sum of Taylor-expansion terms. These basis functions are suitable for use representing
general functions, particularly the linear and delta bases.  The Das_Andersen basis type is
essentially a specialized Taylor expansion, and could be used as such.

 */

#ifndef BASIS_FN
#define BASIS_FN

/*! Finds the delta basis function given a grid, an index, and a volume */
int get_delta_basis( int index, double vmin, double vmax, double dv, double v );

/*! Finds the linear basis function given a grid, an index, and the grid parameters.  
The [index+1]th basis function is 1.0-A.  Returns -1 if the v doesn't fall within
the specified gridpoint, or if it falls outside of the grid overall */
double get_linear_basis(int index, double vmin, double vmax, double dv, double v);

/*! Finds the specified pressure basis function in the force representation */
double get_pres_basis_force(int b, int N, double v_avg, double v_CG);

/*! Finds the specified pressure basis function in the force representation */
double get_pres_basis_force(int b, int N, double v_avg, double v_CG);

/*! Finds the specified pressure basis function in the potential representation */
double get_pres_basis_potential(int b, int N, double v_avg, double v_CG);

/*! Given a set of coefficients psi, returns the total pressure correction for the volume v_cg */
double get_das_andersen_corr(double *psi, int N_coeff, int N, double v_avg, double v_cg);


#endif
