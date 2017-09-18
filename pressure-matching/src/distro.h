/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file distro.h
@author N. Dunn
@date Apr 30, 2014
@brief Functions related to dataset statistics and distributions

These functions are related to statistical properties of a dataset,
and should be completely portable to any other application.

 */


#ifndef DISTRO_FN
#define DISTRO_FN


#include <stdio.h>


/*! Accepts an array of data and produces a distribution plot outfile for it */
void make_distro(double *data, int len_data, int bin_count, FILE *fout);

/*! Finds the average of a data set of length len_data */
double get_avg(double *data, int len_data);

/*! Finds and returns the maximum value in the array data */
double get_max(double *data, int len_data);

/*! Finds and returns the minimum value in the array data */
double get_min(double *data, int len_data);

/*! Finds the variance of a data set given the avg */
double get_variance(double *data, int len_data, double avg);

/*! Checks that the spacing in an array of doubles is uniform and returns the spacing.
If not uniform, the program will exit here */
double check_spacing(double *array, int N);

/*! Finds and returns the minimum value in the array data */
double get_min(double *data, int len_data);

/*! Finds the bin index for quantity X in a grid with spacing dx and starting at minX */
int  get_bin_index( double X, double dx, double minX );

/*! Inverse of get_bin_index.  Given the index value, returns the value 
 X that corresponds to that bin. */
double get_bin_val(int index, double dx, double minX);

/*! Returns an array containing the grid specified by the arguments */
double * get_grid(double xmin, double dx, int N);

/*! bins quanty along quantx with a bin width dx.  returns a pointer to a 2d array that contains
the binned data along element [0] and the number of samples in that bin along element [1]. (allocates memory) */
double ** bin_data_avg( double* quantX, double* quantY, double dx, int n_frames, double xmin, int n_bins );

#endif
