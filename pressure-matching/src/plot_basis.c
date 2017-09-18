/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */


/* plot_basis.c - creates tables of pressure correction values for visualization */

/**
@file plot_basis.c
@author N. Dunn
@date Apr 24, 2014
@brief Plots a Das_Andersen-type pressure correction given the psi.dat solution
file

*/

#include <math.h>
#include <string.h>

#include "io.h"
#include "pressure_matching_fn.h"

/*! Driver function for plotting a Das Andersen basis representation of Fv as a function of the volume. */
int main( int argc,char *argv[] )
{

	int i, j;
	int n_bins;
	int nargs = 5; //the number of arguments expected by this program
	double v_max, v_min, v_curr, dv;
	double press_corr, potential_press_corr;
	FILE *psi_table, *outfile;
	pres_match_sys sys;


	//checks that there are the expected number of arguments
	if (nargs + 1 != argc)
	{
		fprintf(stderr, "Error: unexpected input format\n");
		fprintf(stderr, "Use plot_basis [psi_table_filename] [outfile] [vmin] [vmax] [n_bins]\n");
		fprintf(stderr, "There's usually a script paired with this program that's all set up to go\n");
		exit(1);
	}

	//collects needed info from command line input
	psi_table = fopen( argv[1], "r" );
	outfile = fopen( argv[2], "w" );
	sscanf( argv[3], "%lg", &v_min );
	sscanf( argv[4], "%lg", &v_max );
	sscanf( argv[5], "%d", &n_bins);

	//read in psi table values
	read_psi_table(&sys, psi_table);

	//check input psi values
	fprintf(stdout, "%d coefficients\n%d sites\n%lg is the avg_vol\n",\
	 sys.n_basis, sys.N_sites, sys.avg_AA_vol);

	fprintf(stdout, "The coefficients are:\n");
	for (i=0; i<sys.n_basis; i++)
	{
		fprintf(stdout, "psi %d: %lg\n", i, sys.psi[i]);
	}


	//cut vmax-vmin into even increments based on n_bins
	dv = (v_max - v_min) / n_bins;
	v_curr = v_min;

	
	for (i=0; i<n_bins; i++) //loop over the range of volumes
	{
		v_curr += dv;
		press_corr = 0;
		potential_press_corr = 0;

		for (j=1; j <= sys.n_basis; j++) // loop over the pressure basis functions
		{
			press_corr += find_DA_pressure(&sys, j, v_curr);
			potential_press_corr += find_DA_potential(&sys, j, v_curr);
		}

		fprintf(outfile, "%lg, %lg, %lg\n", v_curr, press_corr, potential_press_corr);

	}



	fclose(psi_table);
	fclose(outfile);

	return 0;
}
