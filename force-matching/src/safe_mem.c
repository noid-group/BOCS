/**
@file safe_mem.c 
@author Nicholas Dunn
@brief Functions related to safe memory management
*/

//c library includes
#include <errno.h>
#include <stdio.h>
#include <string.h>

//local includes
#include "safe_mem.h"


/*! This function replaces calloc in code, and causes the program to 
exit when memory cannot be allocated rather than permit undefined behavior.*/
void *ecalloc(size_t number, size_t size)
{
/**
@param number The number of blocks to allocate in memory
@param size The size of each block in bytes

@return Pointer to a block of memory of the requested size, if it can be allocated
*/
    void *v = calloc(number, size);

    if (v == NULL && number * size != 0) {
	fprintf(stderr,
		"ERROR: unable to calloc block of size %zu bytes: %s\n",
		number * size, strerror(errno));
	exit(1);
    }

    return v;
}



/*! This function replaces malloc in code, and causes the program to
exit when memory cannot be allocated rather than permit undefined behavior */
void *emalloc(size_t size)
{
/**
@param size The size of each block in bytes

@return Pointer to a block of memory of the requested size, if it can be allocated
*/

    void *v = malloc(size);

    if (v == NULL && size != 0) {
	fprintf(stderr,
		"ERROR: unable to malloc block of size %zu bytes: %s\n",
		size, strerror(errno));
	exit(1);
    }

    return v;
}


/*! This function replaces realloc in code, and causes the program to
  exit with an error when memory cannoy be allocated rather than 
  permit undefined behavior. */
void *erealloc(void *old_ptr, size_t size)
{
/**
@param old_ptr The location of the pointer to be reallocated
@param size The size of each block in bytes

@return Pointer to a block of memory of the requested size, if it can be allocated
*/

    void *v = realloc(old_ptr, size);

    if (v == NULL && size != 0) {
	fprintf(stderr,
		"ERROR: unable to realloc block of size %zu bytes: %s\n",
		size, strerror(errno));
	exit(1);
    }


    return v;
}


/*! This function replaces free in code, and checks
 * if a pointer is NULL before trying to free it. */
void efree(void *ptr)
{
/**
@param ptr The location of the pointer to be free'd
*/

    if (ptr != NULL ) {
	free(ptr);
    }

}
