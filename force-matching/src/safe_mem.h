/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file safe_mem.h 
@author Nicholas Dunn
@brief Functions related to safe memory management
*/

#ifndef SAFE_MEM
#define SAFE_MEM

#include <stdlib.h>


/** This function replaces calloc in code, and causes the program to 
exit when memory cannot be allocated rather than permit undefined behavior.*/
void *ecalloc(size_t number, size_t size);

/*! This function replaces malloc in code, and causes the program to
exit when memory cannot be allocated rather than permit undefined behavior */
void *emalloc(size_t size);

/*! This function replaces realloc in code, and causes the program to
 * exit with an error when memory cannoy be allocated rather than 
 * permit undefined behavior. */
void *erealloc(void *old_ptr, size_t size);

#endif
