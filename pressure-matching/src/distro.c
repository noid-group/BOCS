/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**


@file distro.c 
@author N. Dunn
@date Apr 24, 2014
@brief Functions related to calculating distributions
*/

#include <string.h>
#include <math.h>

#include "safe_mem.h"
#include "distro.h"
#include "basis.h"
#include "solv_lin_eqns.h"

/*! Accepts an array of data and prints a normalized sampling 
    distribution to the file pointed to by \e fout */
void make_distro(double *data, int len_data, int bin_count, FILE *fout)
{
/**
@param data The array containing the data to be distributed
@param len_data The number of data points in data
@param bin_count The number of bins to split the data into
@param fout The pointer to the output file

*/

	int i, bindex;
	double max, min, avg;
	double bin_width;
	double variance;
	double norm_const;

	double *distro_bin = (double *) ecalloc(bin_count, sizeof(double));

	avg = get_avg(data, len_data);
	max = get_max(data, len_data);
	min = get_min(data, len_data);

	variance = get_variance(data, len_data, avg);


	fprintf(stdout, "Average is: %lg\n", avg);
	fprintf(stdout, "Max is: %lg\nMin is: %lg\n", max, min);
	fprintf(stdout, "Variance is: %lg\n", variance);

	//run through volumes and increment each bin appropriately
	bin_width = (max - min) / bin_count;

	fprintf(stdout, "Bin width is: %lg\n", bin_width);

	for (i=0; i < len_data; i++)
	{
		bindex = get_bin_index(data[i], bin_width, min);

		distro_bin[bindex]++;
	}


	norm_const = 0;
	for (i=0; i < bin_count; i++)
	{
		norm_const += distro_bin[i] * bin_width;
	}

	for (i=0; i < bin_count; i++)
	{
		distro_bin[i] = distro_bin[i] / norm_const;

		fprintf(fout, "%lg %lg\n", min + (i*bin_width), distro_bin[i]);

	}

	free(distro_bin);

}

/*! Finds the average (mean) of a data set of length len_data */
double get_avg(double *data, int len_data)
{
/**
@param data The array containing the data to be averaged
@param len_data The number of data points in data

@return The average of the numbers in \e data.
*/

	int i;
	double avg;

	avg = 0;

	for (i=0; i < len_data; i++)
	{
		avg += data[i];
	}

	avg = avg / len_data;

	return avg;
}


/*! Finds and returns the maximum value in the array data */
double get_max(double *data, int len_data)
{
/**
@param data The array containing the data of which to find the min
@param len_data The number of data points in data

@return The max of the numbers in \e data.
*/
	int i;
	double max;

	max = data[0];

	for (i=0; i < len_data; i++)
	{
		if (data[i] > max)
		{
			max = data[i] + (data[i] - data[i-1]);
		}
	}

	return max;
}




/*! Finds and returns the minimum value in the array data */
double get_min(double *data, int len_data)
{
/**
@param data The array containing the data of which to find the max
@param len_data The number of data points in data

@return The minimum of the numbers in \e data.
*/
	int i;
	double min;

	min = data[0];

	for (i=0; i < len_data; i++)
	{
		if (data[i] < min)
		{
			min = data[i];
		}
	}

	return min;
}


/*! Finds the variance of a data set given the avg */
double get_variance(double *data, int len_data, double avg)
{
/**
@param data The array containing the data for which to calculate the variance
@param len_data The number of data points in data
@param avg The average of the numbers in \e data.

@return The variance of the numbers in \e data.
*/
	int i;
	double variance;

	variance = 0;
	for (i=0; i < len_data; i++)
	{
		variance += pow(data[i] - avg, 2);
	}

	variance = variance / len_data;


	return variance;
}


/*! Checks that the spacing in an array of doubles is uniform and returns the spacing.
    If not uniform, exits the program and throws an error. */
double check_spacing(double *array, int N)
{
/**
@param array The array whose grid spacing is to be checked
@param N The length of array

@return The grid spacing, if it is found to be uniform. 
*/

	int i;
	double dx;
	double test_dx;
	double error;

	dx = fabs(array[1] - array[0]);

	for (i=2; i<N; i++)
	{
		test_dx = fabs(array[i] - array[i-1] );
		error = fabs(dx - test_dx);
		if (error > dx / 2.0 )
		{
			fprintf(stderr, "ERROR: Spacing in table file is uneven between elements %d and %d\n", i-1, i);
			exit(1);
		}
	}

	return dx;
}


/*! Finds the bin index for quantity X in a grid with spacing dx and starting at minX */
int  get_bin_index( double X, double dx, double minX )
{
/**
@param X The value to find the bin for
@param dx The grid spacing
@param minX The minimum value of the grid

@return The index of the bin where \e X falls
*/
	int bin_index;

	bin_index = (int) floor(((X - minX) / dx));

	return bin_index;
}


/*! Inverse of get_bin_index.  Given the bin index, returns the value 
 X that corresponds to that bin. */
double get_bin_val(int index, double dx, double minX)
{
/**
@param index The index within the grid
@param dx The grid spacing
@param minX The minimum value of the grid

@return The lower bound of the bin corresponding to the provided \e index.
*/
	double X;

	X = minX + index*dx;

	return X;
}


/*! Returns an array containing the grid specified by the arguments */
double * get_grid(double xmin, double dx, int N)
{
/**
@param xmin The minimum (starting) value of the grid
@param dx The grid spacing
@param N The number of grid points

@return An array containing the lower-bound values of each bin
*/

	int i;
	double *grid;

	grid = (double*) ecalloc(N, sizeof(double));

	for (i=0; i<N; i++)
	{
		grid[i] = xmin + (dx * i);
	}

	return grid;
}


/*! bins quanty along quantx with a bin width dx.  returns a pointer to a 2d array that contains
the binned data along element [0] and the number of samples in that bin along element [1]. (allocates memory) */
double ** bin_data_avg( double* quantX, double* quantY, double dx, int n_frames, double xmin, int n_bins )
{
/**
@param quantX The independent variable. Defines the grid bin that a point occupies
@param quantY The dependent variable.  This value is averaged per bin. 
@param dx The grid spacing
@param n_frames The number of entries in quantX and quantY
@param xmin The minimum grid value
@param n_bins The number of gridpoints to use

@return Pointer to a 2d array, in which element [0] is an array of the average value of quantY in that bin, and element \
[1] is an array containing the number of samples per bin.
*/
	double **binned_data;
	int i;
	int bin_index;


	binned_data = (double**) ecalloc(2, sizeof(double*));
	binned_data[0] = (double*) ecalloc(n_bins, sizeof(double));
	binned_data[1] = (double*) ecalloc(n_bins, sizeof(double));

	for (i=0; i<n_frames; i++)
	{
		bin_index = get_bin_index( quantX[i], dx, xmin );

		if (bin_index < n_bins && bin_index >= 0)
		{
			binned_data[1][bin_index] += 1;
			binned_data[0][bin_index] += quantY[i];
		}
	}

	for (i=0; i<n_bins; i++)
	{
		if (binned_data[1][i] != 0)
		{
			binned_data[0][i] = binned_data[0][i]/binned_data[1][i];
		} else
		{
			binned_data[0][i] = 0;
		}

		fprintf(stdout, "Data[%d] = %f  %f\n", i, binned_data[0][i], binned_data[1][i]);
	}

	return binned_data;
}
