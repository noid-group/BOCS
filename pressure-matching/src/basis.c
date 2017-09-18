/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**

@file basis.c 
@author N. Dunn
@date Apr 24, 2014
@brief Functions related to basis representations of pressure correction
*/

#include <stdlib.h>
#include <math.h>
#include "basis.h"
#include "distro.h"



/*! Finds the delta basis function given a grid, an index, and a volume.
  Returns 1 if v falls into the grid specified by index, 0 otherwise. */
int get_delta_basis( int index, double vmin, double vmax, double dv, double v )
{
/**
@param index The bin index to calculate A for
@param vmin The minimum (starting) value of the grid
@param vmax The maximum value of the grid
@param dv The grid spacing
@param v The volume for which the basis is being calculated

@return The contribution from the delta basis to bin i.
*/
	int delta = 0;
	int bin_id;

	if ( (v >= vmin) && (v <= vmax))
	{
		bin_id = get_bin_index(v, dv, vmin);

		if (bin_id == index) { delta = 1; } 
	}

	return delta;
}




/*! Returns a term in the pressure form of the das_andersen
    basis function. */
double get_pres_basis_force(int b, int N, double v_avg, double v_CG)
{
/**
@param b The order of the basis term to be calculated.  (b>0 for the force)
@param N The number of sites in the CG model
@param v_avg The average volume of the reference atomistic simulation
@param v_CG The volume for which the value of the basis function is being calculated


@return The <em>b</em>th term in the Das-Andersen correction to the pressure
*/
	double p_b;

	if (b < 1)
	{
		fprintf(stderr, "ERROR: invalid basis index %d\n", b);
		exit(1.0);

	} else if (b == 1)
	{
		p_b = N / v_avg;
	} else
	{
		p_b = (N * b / v_avg) * pow( (1 / v_avg) * (v_CG - v_avg), b-1);
	}

	return (-1) * p_b;
}


/*! Returns a term in the potential-pressure form of the das_andersen
    basis function. */
double get_pres_basis_potential(int b, int N, double v_avg, double v_CG)
{
/**
@param b The order of the basis term to be calculated.  (b>=0 for the potential)
@param N The number of sites in the CG model
@param v_avg The average volume of the reference atomistic simulation
@param v_CG The volume for which the value of the basis function is being calculated

@return The <em>b</em>th term in the Das-Andersen correction to the PV potential
*/

	double p_b;

	if (b < 1)
	{
		fprintf(stderr, "ERROR: invalid basis index %d\n", b);
		exit(1.0);

	} else if (b == 1)
	{
		p_b = N * v_CG / v_avg;
	} else
	{
		p_b = N * pow( ((v_CG - v_avg) / v_avg), b );
	}

	return p_b;
}

/*! Given a set of coefficients psi, returns the total das_andersen
    pressure correction for the volume v_cg */
double get_das_andersen_corr(double *psi, int N_coeff, int N, double v_avg, double v_cg)
{
/**
@param psi The solution coefficients of the pressure correction
@param N_coeff The number of terms in the Das_Andersen basis 
@param N The number of sites in the CG model
@param v_avg The average volume of the reference atomistic simulation
@param v_cg The volume for which the value of the basis function is being calculated

@return The total Das-Andersen pressure correction for the given volume

*/
	int i;
	double corr = 0;

	for (i=1; i<N_coeff+1; i++)
	{
		corr += psi[i] * get_pres_basis_force(i, N, v_avg, v_cg);
	}

	return corr;
}





