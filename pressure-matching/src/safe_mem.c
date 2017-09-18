/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file safe_mem.c 
@author N. Dunn
@date May 12, 2014
@brief Functions related to safe memory management
*/

#include "safe_mem.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>

/*! This function replaces calloc in code, and causes the program to 
exit when memory cannot be allocated rather than permit undefined behavior.*/
void * ecalloc(size_t number, size_t size)
{
/**
@param number The number of blocks to allocate in memory
@param size The size of each block in bytes

@return Pointer to a block of memory of the requested size, if it can be allocated
*/
	void *v = calloc(number, size);

	if (v == NULL  && number * size != 0)
	{
		fprintf(stderr, "ERROR: unable to calloc block of size %zu bytes: %s\n", number*size, strerror(errno));
		exit(1);
	}

	return v;
}



/*! This function replaces malloc in code, and causes the program to
exit when memory cannot be allocated rather than permit undefined behavior */
void * emalloc(size_t size)
{
	void *v = malloc(size);

	if (v == NULL && size != 0)
	{
		fprintf(stderr, "ERROR: unable to malloc block of size %zu bytes: %s\n", size, strerror(errno));
		exit(1);
	}

	return v;
}
