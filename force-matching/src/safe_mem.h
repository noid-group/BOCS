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

/*! This function replaces free in code, and checks
 * if a pointer is NULL before trying to free it. */
void efree(void *ptr);

#endif
