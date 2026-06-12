/**
@file safe_mem.c 
@author Nicholas Dunn, Maria Lesniewski
@brief Functions related to safe memory management
*/

//c library includes
#include <errno.h>
#include <stdio.h>
#include <string.h>

//local includes
#include "safe_mem.h"

// These Functions are called with just e.g. number/size args in the code, without the _impl appended name
// The macros in safe_mem.h link those calls to these implementations

/*! This function replaces calloc in code, and causes the program to 
exit when memory cannot be allocated rather than permit undefined behavior.*/
void *ecalloc_impl(size_t number, size_t size, const char* file, int line, const char* func)
{
/**
@param number The number of blocks to allocate in memory
@param size The size of each block in bytes
@param file The filename of the function calling ecalloc
@param line The line number ecalloc was called
@param func The function calling ecalloc

@return Pointer to a block of memory of the requested size, if it can be allocated
*/
    void *v = calloc(number, size);

    if (v == NULL && number * size != 0) {
	fprintf(stderr,
		"ERROR: unable to calloc block of size %zu bytes: %s\n",
		number * size, strerror(errno));
	fprintf(stderr, "File %s, Line: %d, Function: %s()\n", file, line, func);
	exit(1);
    }

    return v;
}



/*! This function replaces malloc in code, and causes the program to
exit when memory cannot be allocated rather than permit undefined behavior */
void *emalloc_impl(size_t size, const char* file, int line, const char* func)
{
/**
@param size The size of each block in bytes
@param file The filename of the function calling emalloc
@param line The line number emalloc was called
@param func The function calling emalloc

@return Pointer to a block of memory of the requested size, if it can be allocated
*/

    void *v = malloc(size);

    if (v == NULL && size != 0) {
	fprintf(stderr,
		"ERROR: unable to malloc block of size %zu bytes: %s\n",
		size, strerror(errno));
	fprintf(stderr, "File %s, Line: %d, Function: %s()\n", file, line, func);

	exit(1);
    }

    return v;
}


/*! This function replaces realloc in code, and causes the program to
  exit with an error when memory cannoy be allocated rather than 
  permit undefined behavior. */
void *erealloc_impl(void *old_ptr, size_t size, const char* file, int line, const char* func)
{
/**
@param old_ptr The location of the pointer to be reallocated
@param size The size of each block in bytes
@param file The filename of the function calling erealloc
@param line The line number erealloc was called
@param func The function calling erealloc


@return Pointer to a block of memory of the requested size, if it can be allocated
*/

    void *v = realloc(old_ptr, size);

    if (v == NULL && size != 0) {
	fprintf(stderr,
		"ERROR: unable to realloc block of size %zu bytes: %s\n",
		size, strerror(errno));
	fprintf(stderr, "File %s, Line: %d, Function: %s()\n", file, line, func);
	exit(1);
    }


    return v;
}


/*! This function replaces free in code, and checks
 * if a pointer is NULL before trying to free it. */
void efree_impl(void **ptr, const char* file, int line, const char* func)
{
/**
@param ptr The location of the pointer to be free'd
@param file The filename of the function calling efree
@param line The line number efree was called
@param func The function calling efree

*/
    if (ptr == NULL || *ptr == NULL) {
    fprintf(stderr, "Warning: Attempted to free null pointer in\n");
    fprintf(stderr, "File %s, Line: %d, Function: %s()\n", file, line, func);

    return; 
    }
    if (*ptr != NULL ) {
	free(*ptr);
	*ptr=NULL;
    }

}
