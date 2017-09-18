/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file io.h 
@author N. Dunn
@date July 2, 2013
@brief Structure and definitions for i/o operations with csv and xvg files
*/

#ifndef IO_FN
#define IO_FN

//#include <stdio.h>
#include <stdlib.h>
#include "pressure_matching_types.h"

typedef struct
{
/**
@struct xvg_cont
@brief The contents of an xvg file
*/
	int has_header; ///< 1 if the xvg_cont contains a header, 0 otherwise
	char * header; ///< A string containing the header of the file

	int N_lines; ///< The number of rows present in the file
	int N_columns; ///< The number of columns present in the file

	double **contents; ///< An array of dimenstion N_linesxN_columns containing actual data

} xvg_cont;

/*! Counts the lines in the supplied file and rewinds it */
int count_lines(FILE *file);

/*! Reads the contents of an xvg file to memory - [allocates memory]*/
void read_xvg(FILE *xvg, xvg_cont *file_cont);

/*! Counts the number of rows containing data in the provided xvg file (public) */
int get_xvg_rows(FILE *xvg);

/*! Prints the xvg contents to a duplicate file to check for read errors */
void print_xvg_contents(xvg_cont *file_cont, FILE *output);

/*! Frees the memory associated with an xvg_cont's data array */
void close_xvg_cont(xvg_cont *file_cont);

/*! Frees all the memory associated with an xvg_cont structure */
void close_xvg(xvg_cont *file_cont);

/*! Makes a new xvg object from memory */
xvg_cont * make_new_xvg(int n_lines, int n_columns, double **data);

/*! Gets a double setting from the configfile*/
double get_setting_d(char setting[], FILE *config);

/*! Gets an int setting from the configfile*/
int get_setting_i(char setting[], FILE *config);

/*! Gets a string from the file config  */
void get_string_setting(char setting[], FILE *config, char *string);

/*! Returns a pointer to the desired file */
FILE *safe_open_file(char *filename, char *mode);

/*! Read in the specified column of an xvg file to the provided array, which must already be allocated to the appropriate size */
void read_xvg_to_array(FILE *source_file, double *target_array, int column );

/*! Checks the provided xvg files to make sure they have the same # of frames */
int check_xvg_lengths(FILE **files, int num_files);

/*! Reads a Das_Andersen psi output table back into a pres_match_sys */
void read_psi_table(pres_match_sys *pres_sys, FILE* psi_table);

/*! Reads a tabulated Fv into the provided pres_sys. Will throw an error
and halt the program if the grid spacing is uneven. */
void read_tabulated_basis(pres_match_sys *pres_sys, FILE *table_file);

/*! Prints the correlation matrix Q to an ascii file in a direct row-by-column fashion */
void print_Q(pres_match_sys *pres_sys);

/*! Prints a generic array of doubles of length len to a file  */
void print_array_d(double *array, int len, FILE *file);

/*! Prints psi output to a "table" file for use by a simulation package with knowledge of the basis forms. */
void print_psi_table(pres_match_sys *pres_sys);




#endif
