/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**


@file io.c 
@author N. Dunn
@date Apr 24, 2014
@brief Functions related to file reading and writing
*/

#include <errno.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "io.h"
#include "safe_mem.h"
#include "pressure_matching_fn.h"
#include "distro.h"


/*! Returns the first non-whitespace character from a string. Returns an
empty character for strings that are only whitespace */
static char first_non_whitespace(char *str)
{
/**
@param str The string to find the first non-whitespace character of

@return The first non-whitespace character
*/
	int i=0;

	// Skip leading whitespace
	while (isspace(str[i]))
	{
		i++;
	}

	return str[i];

}

/*! Counts the lines in the supplied file and rewinds it */
int count_lines(FILE *file)
{
/**
@param file The file whose lines are to be counted

@return The number of lines in \e file.
*/
	int count = 0;
	size_t nbytes = 100;
	char *line;

	line = ecalloc(nbytes, sizeof(char));
	
	while ( getline(&line, &nbytes, file) != -1 && isdigit(first_non_whitespace(line)) )
	{
		count++;
	}

	rewind(file);

	free(line);

	return count;
}


/*! Allocates the memory for a NxM array for use in xvg_cont structs */
static double ** allocate_data_mem(double **contents, int N_lines, int N_columns, char *data_type)
{
/**
@param contents The pointer to the memory to be allocated (Will be N_columns by N_rows)
@param N_lines The number of lines in the *_cont structure
@param N_columns The number of columns in the *_cont structure
@param data_type A string that identifies the type of file struct being dealt with

@return Pointer to the allocated memory
*/
	int i;

	if (strcmp(data_type, "xvg") == 0)
	{
		contents = (double**) ecalloc(N_columns, sizeof(double*));

		for (i=0; i<N_columns; i++)
		{
			contents[i] =  (double *) ecalloc(N_lines, sizeof(double));
		}
	} else
	{
		fprintf(stderr, "ERROR: cannot allocate memory for unknown data type \'%s\'", data_type);
		exit(1);
	}

	return contents;
} 



/*! Reads over the header of an xvg file and returns the number of lines skipped */
static int eat_header(FILE *xvg)
{
/**
@param xvg The file whose header you wish to skip

@return The number of lines skipped by this function
*/
	int i;
	int n_header, bytes_read; //number of lines in the header
	int is_header = 1;
	char *line;
	size_t len_line;

	len_line = 256;
	line = ecalloc(len_line, sizeof(char));
	n_header = 0;

	//eats the header of the file and counts the number of lines
	bytes_read = getline(&line, &len_line, xvg);

	while ( is_header )
	{
//		if ((line[0] != '#') && (line[0] != '@'))
		if ( isdigit(first_non_whitespace(line)) )		
		{
			is_header = 0;
		} else
		{
			n_header++;

			bytes_read = getline(&line, &len_line, xvg);
		}
	}
	
	rewind(xvg);
	
	for (i=0; i<n_header; i++)
	{
		bytes_read = getline(&line, &len_line, xvg);
	}

	free(line);

	return n_header;
}







/*! Counts the number of rows containing data in the provided xvg file, skipping
any comment or header lines */
int get_xvg_rows(FILE *xvg)
{
/**
@param xvg Pointer to file for which this function will count the non-header rows

@return The number of rows of data in the xvg file
*/
	int rows;

	eat_header(xvg);

	rows = count_lines(xvg); 

	rewind(xvg);

	return rows;
}

/*! Counts the number of columns in an xvg file */
static int count_xvg_columns(FILE *xvg)
{
/**
@param xvg Pointer to file for which the space-separated columns will be counted

@return The number of columns in the xvg file
*/
	int i, j;
	int columns = 1; // we assume the file has *some* data in it
	int bytes_read;
	size_t len_line = 256;
	char *line;

	line = ecalloc(len_line, sizeof(char));

	eat_header(xvg);
	bytes_read = getline(&line, &len_line, xvg);

	//seek first # after whitespace
	/* |  99950.007812  5648.928223| */

	j = 0;
	while ( isspace(line[j]) ) //skip blanks at the beginning of the line, if any
	{
		j++;
	}

	for (i=j; i<bytes_read; i++)
	{
		if ( isspace(line[i]) )
		{
			//if next is within range && if next is a digit
			if ( i < (bytes_read -1) && !isspace(line[i+1])   )
			{
				columns++;
			}
		}
	}

	rewind(xvg);

	free(line);

	return columns;
}

/*! Reads a line of an xvg file into memory as doubles */
static void read_xvg_line(char *line, int len_line, int line_id, xvg_cont *file_cont)
{
/**
@param line The line to be read into the structure
@param len_line The number of characters in line
@param line_id The index of the line within the file
@param file_cont The structure to which the line will be read
*/
	int i = 0;
	int j = 0;
	int column = 0;
	char *phrase;

	phrase = ecalloc(len_line, sizeof(char));

	while ( line[i] != '\n' )
	{
		while ( isspace(line[i]) || line[i] == ','  ) //skip blanks at the beginning of the line, if any
		{
			i++;
		}

		if (line[i] == '\n' ){ break;}


		while ( !isspace(line[i]) && line[i] !=  ',' )
		{
			if ( isdigit(line[i]) || line[i] == '.' || line[i] == 'e' \
			  || line[i] == 'E' || line[i] == '-' || line[i] == '+' )
			{
				phrase[j] = line[i];
				j++;

			}

			i++;

			
			if  ( isspace( line[i] ) || line[i] == ',' )
			{
				sscanf(phrase, "%lg", &file_cont->contents[column][line_id]);
				memset(phrase, '\0', j+1); 

				j=0;
				column++;
			}
		}
	}

	free(phrase);
}






/*! Reads the contents of an xvg file to memory */
void read_xvg(FILE *xvg, xvg_cont *file_cont)
{
/**
@param xvg Pointer to the file to be read into memory
@param file_cont The structure to which the file will be read
*/
	char *line;
	int n_header, bytes_read, i;
	size_t len_line;

	len_line = 100; /* Suggested line length, but getline() will realloc if needed */
	n_header = 0;
	bytes_read = 0;

	line = ecalloc(len_line, sizeof(char));

	file_cont->N_lines = get_xvg_rows(xvg);

	file_cont->N_columns = count_xvg_columns(xvg);

	file_cont->contents = allocate_data_mem(file_cont->contents, file_cont->N_lines, file_cont->N_columns, "xvg");

	file_cont->has_header = 0;
	n_header = eat_header(xvg);

	for (i = 0; i < file_cont->N_lines; i++)
	{
		bytes_read = getline(&line, &len_line, xvg);

		if ( isdigit(first_non_whitespace(line)) )
		{
			read_xvg_line(line, bytes_read, i, file_cont);
		}
	}

	free(line);

	rewind(xvg);

}


/*! Prints the xvg contents to a duplicate file to check for read errors */
void print_xvg_contents(xvg_cont *file_cont, FILE *dupe)
{
/**
@param file_cont The xvg contents to be dumped to a file
@param dupe Pointer to the file where file_cont will be written
*/
	int i, j, buff_size;
	char *outline, *buffer;

	buff_size = 50;

	buffer  = ecalloc(buff_size, sizeof(char));
	outline = ecalloc(file_cont->N_columns * buff_size, sizeof(char));
	

	for (i=0; i<file_cont->N_lines; i++)
	{
		
		for (j=0; j<file_cont->N_columns; j++)
		{
			if(j==0)
			{
				sprintf(outline, "%f ", file_cont->contents[j][i]);

			} else if (j < (file_cont->N_columns - 1) )
			{
				sprintf(buffer, "%f ", file_cont->contents[j][i]);

				strcat(outline, buffer);
			} else
			{
				sprintf(buffer, "%f\n", file_cont->contents[j][i]);

				strcat(outline, buffer);
			}
		}

		fprintf(dupe, outline);

	}

	free( buffer );
	free( outline );
}




/*! Frees all the memory associated with an xvg_cont structure */
void close_xvg(xvg_cont *file_cont)
{
/**
@param file_cont The structure whose memory is to be freed
*/
	int i;

	for (i = 0; i < file_cont->N_columns; i++)
	{
		free(file_cont->contents[i]);
	}

	free(file_cont->contents);

	if (file_cont->has_header)
	{
		free(file_cont->header);
	}

	//free(file_cont);
}





/*! Allocates a new xvg_cont object in memory */
xvg_cont * make_new_xvg(int n_lines, int n_columns, double **data)
{
/**
@param n_lines The number of lines (rows) to allocate in memory
@param n_columns The number of columns to allocate in memory
@param data The array containing the data to put into the xvg struct

@return Pointer to the newly allocated xvg_cont object
*/
	int i, j;
	xvg_cont *new = (xvg_cont *) emalloc(sizeof(xvg_cont));

	new->N_lines = n_lines;
	new->N_columns = n_columns;

	new->contents = allocate_data_mem(new->contents, n_lines, n_columns, "xvg");

	for (i=0; i<n_columns; i++)
	{
		for (j=0; j<n_lines; j++)
		{
			new->contents[i][j] = data[i][j];
		}
	}

	return new;
}


/*! Gets a double setting from the config file*/
double get_setting_d(char setting[], FILE *config)
{
/**
@param setting The name of the setting to find in config
@param config The configuration file containing various settings in 'Setting Value' format

@return The value of the setting pulled from the config file
*/
	double value;
	char candidate[128];

	fscanf(config, "%s", candidate);

	while (strcmp(setting, candidate) != 0)
	{	
		fscanf(config, "%s", candidate);

		if (getc(config) == EOF)
		{
			rewind(config);
			printf("Setting for %s not found.\n", setting);
			exit(-1);
		}
	}
	
	fscanf(config, "%lg", &value);

	rewind(config);
	return value;
}


/*! Gets an int setting from the config file*/
int get_setting_i(char setting[], FILE *config)
{
/**
@param setting The name of the setting to find in config
@param config The configuration file containing various settings in 'Setting: Value' format

@return The value of the setting pulled from the config file
*/
	int value;
	char candidate[128];

	fscanf(config, "%s", candidate);

	while (strcmp(setting, candidate) != 0)
	{	
		fscanf(config, "%s", candidate);

		if (getc(config) == EOF)
		{
			rewind(config);
			printf("Setting for %s not found.\n", setting);
			exit(-1);
		}
	}
	
	fscanf(config, "%d", &value);

	rewind(config);
	return value;
}



/*! Gets a string setting from the config and stores it in the memory pointed to by \e string  */
void get_string_setting(char setting[], FILE *config, char *string)
{
/**
@param setting The name of the setting to find in config
@param config The configuration file containing various settings in 'Setting: Value' format
@param string The address where the value will be held

*/
	char candidate[128];

	fscanf(config, "%s", candidate);	

	while (strcmp(setting, candidate) != 0)
	{
		fscanf(config, "%s", candidate);

		if (getc(config) == EOF)
		{
			rewind(config);
			printf("Setting for %s not found.\n", setting);
			exit(-1);
		}
	}

	fscanf(config, "%s", string);
	
	rewind(config);
}


/*! Returns a pointer to the desired file in the designated mode. The program will
 exit if the file can't be opened. */
FILE *safe_open_file(char *filename, char *mode)
{
/**
@param filename The name of the file to open
@param mode The mode in which to open the file (w, r, a, etc.)

@return A pointer to the file that was generated
*/
	FILE *file;

	file = fopen(filename, mode);
	if (file == NULL)
	{
		fprintf(stderr, "The following error occurred when opening file \"%s\": %s\n", filename, strerror(errno));
		exit(1);
	}

	return file;
}

/*! Read in the specified column of an xvg file to the provided array, which must already be allocated to the appropriate size */
void read_xvg_to_array(FILE *source_file, double *target_array, int column ) 
{
/**
@param source_file The file fomr which to read the data
@param target_array The (preallocated) array to which the xvg_data will be read
@param column  The column (starting from 0) in the xvg file to read
*/
	int i;
	xvg_cont data;

	read_xvg(source_file, &data);

	for (i=0; i<data.N_lines; i++)
	{
		target_array[i] = data.contents[column][i];
	}

	close_xvg(&data);

}





/*! Checks the provided xvg files to make sure they have the same # of frames
and returns the number of frames if they match. If there is a mismatch, -1 is returned. */
int check_xvg_lengths(FILE **files, int num_files)
{
/**
@param files An array containing the array of files to check
@param num_files The number of files to check

@return The number of frames in each xvg file if they all match; -1 otherwise
*/

	int i, j;
	int *frame_counts, num_frames;
	int mismatch = 0;

	frame_counts = ecalloc(num_files, sizeof(int));

	for (i=0; i<num_files; i++)
	{
		frame_counts[i] = get_xvg_rows(files[i]);
	}

	for (i=0; i<num_files; i++)
	{
		for (j=i+1; j<num_files; j++)
		{
			if (frame_counts[i] != frame_counts[j])
			{
				mismatch = 1;
			}
		}
	} 

	if (mismatch == 1)
	{
		return -1;
	} else
	{
		num_frames = frame_counts[0];
	}
	free(frame_counts);

	return num_frames;
}




/*! Reads a Das_Andersen psi output table back into a pres_match_sys */
void read_psi_table(pres_match_sys *pres_sys, FILE* psi_table)
{
/**
@param pres_sys The psi table will be read into this data structure
@param psi_table The file containing the Das_Andersen basis parameters 
*/

	int i, j;
	bool avg_vol_read = false;
	bool psi_coeff_read = false;
	bool n_sites_read = false;
	int N_section = 3; //number of sections in the psi table file
	
	int bytes_read;
	size_t len_line = 100;
	char *line;

	line = ecalloc(len_line, sizeof(char));

	//loops over expected parameters
	for (i=0; i < N_section; i++)
	{
		bytes_read = getline(&line, &len_line, psi_table);

		if ( strncmp(line, "[AVG_AA_VOLUME]", 15) == 0  && avg_vol_read == false)
		{
			bytes_read = getline(&line, &len_line, psi_table);

			sscanf( line, "%lg", &(pres_sys->avg_AA_vol) );

			avg_vol_read = true;

		} else if ( strncmp(line, "[PSI_COEFF]", 11) == 0 && psi_coeff_read == false)
		{
			sscanf( line, "[PSI_COEFF] %d", &(pres_sys->n_basis) );

			pres_sys->psi =  ecalloc( pres_sys->n_basis, sizeof(double) );

			for (j = 0; j < pres_sys->n_basis; j++)
			{
				bytes_read = getline(&line, &len_line, psi_table);

				sscanf( line, "%lg", &(pres_sys->psi[j]) );
			}

			psi_coeff_read = true;

		} else if ( strncmp(line, "[N_SITES]", 9) == 0 && n_sites_read == false )
		{

			bytes_read = getline(&line, &len_line, psi_table);

			sscanf( line, "%d", &(pres_sys->N_sites) );

			n_sites_read = true;

		} else
		{
			fprintf(stderr, "ERROR: Unexpected file format for psi table\n");
			exit(1);
		}
	}	

	free(line);
}


/*! Reads a tabulated Fv into the provided pres_sys. Will throw an error
and halt the program if the grid spacing is uneven. */
void read_tabulated_basis(pres_match_sys *pres_sys, FILE *table_file)
{
/**
@param pres_sys The main structure for the calculation
@param table_file The file containing the tabulated basis function
*/

	int i;
	xvg_cont table;

	read_xvg(table_file, &table);

	pres_sys->dv = check_spacing(table.contents[0], table.N_lines);
	pres_sys->vmin = table.contents[0][0];
	pres_sys->vmax = table.contents[0][table.N_lines-1];
	pres_sys->n_basis = table.N_lines;

	pres_sys->psi = ecalloc(table.N_lines, sizeof(double));

	for (i=0; i<table.N_lines; i++ )
	{
		pres_sys->psi[i] = table.contents[1][i];
	}

	close_xvg(&table);
}

/*! Prints the correlation matrix Q to an ascii file in a direct row-by-column fashion */
void print_Q(pres_match_sys *pres_sys)
{
/**
@param pres_sys The main structure for the calculation
*/
	int i, j;

	for (i=0; i<pres_sys->n_basis; i++)
	{
		for (j=0; j<pres_sys->n_basis; j++)
		{
			fprintf(pres_sys->files->Q_out, " %f ", pres_sys->Q[i][j]);
		}
		fprintf(pres_sys->files->Q_out, "\n");
	}
}


/*! Prints a generic array of doubles of length len to a file  */
void print_array_d(double *array, int len, FILE *file)
{
/**
@param array The array containing data to be printed
@param len The number of entries in \e array
@param file The file to which the array will be printed
*/
	int i;

	for (i=0; i<len; i++)
	{
		fprintf(file, "%f\n", array[i]);
	}
}


/*! Prints psi output to a "table" file for use by a simulation package with knowledge of the basis forms. */
void print_psi_table(pres_match_sys *pres_sys)
{
/**
@param pres_sys The main structure for the calculation
*/
	int i;
	int index = 0;
	double vol;
	double ref_psi;

	if ( strcmp(DAS_ANDERSEN, pres_sys->basis_type) == 0 )
	{
		fprintf(pres_sys->files->psi_out, "[AVG_AA_VOLUME]\n%f\n", pres_sys->avg_AA_vol);

		fprintf(pres_sys->files->psi_out, "[PSI_COEFF] %d\n", pres_sys->n_basis);
	
		if (pres_sys->ref_Fv_flag == 1)
		{
			for (i=0; i<pres_sys->n_basis; i++)
			{
				fprintf(pres_sys->files->psi_out, "%f\n", pres_sys->psi[i] + pres_sys->ref_psi[i]);
			}
		} else
		{
			for (i=0; i<pres_sys->n_basis; i++)
			{
				fprintf(pres_sys->files->psi_out, "%f\n", pres_sys->psi[i]);
			}
		}

		fprintf(pres_sys->files->psi_out, "[N_SITES]\n%d\n", pres_sys->N_sites);

	} else if ( strcmp(DELTA, pres_sys->basis_type) == 0 )
	{
		for (i=0; i < (pres_sys->n_basis + pres_sys->n_zeroes); i++)
		{
			if ( !pres_sys->zero_list[i] )
			{
				vol = get_bin_val(i, pres_sys->dv, pres_sys->vmin);

				if (pres_sys->ref_Fv_flag == 1)
				{
					ref_psi = get_ref_psi_value(pres_sys, vol);
					fprintf(pres_sys->files->psi_out, "%f, %f\n", vol, pres_sys->psi[index] + ref_psi);
				} else
				{
					fprintf(pres_sys->files->psi_out, "%f, %f\n", vol, pres_sys->psi[index]);
				}
				index++;
			}
		}
	} 
}

