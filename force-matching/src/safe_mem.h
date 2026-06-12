/**
@file safe_mem.h 
@author Nicholas Dunn, Maria Lesniewski
@brief Functions related to safe memory management
*/

#ifndef SAFE_MEM
#define SAFE_MEM

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdlib.h>

/** This function replaces calloc in code, and causes the program to 
exit when memory cannot be allocated rather than permit undefined behavior.*/
void *ecalloc_impl(size_t number, size_t size, const char* file, int line, const char* func);

/*! This function replaces malloc in code, and causes the program to
exit when memory cannot be allocated rather than permit undefined behavior */
void *emalloc_impl(size_t size, const char* file, int line, const char* func);

/*! This function replaces realloc in code, and causes the program to
 * exit with an error when memory cannoy be allocated rather than 
 * permit undefined behavior. */
void *erealloc_impl(void *old_ptr, size_t size, const char* file, int line, const char* func);

/*! This function replaces free in code, and checks
 * if a pointer is NULL before trying to free it. */
void efree_impl(void **ptr, const char* file, int line, const char* func);

// MCL 02.25.26 - In the code we call ecalloc, emalloc, erealloc, efree. 
// The below macros link those calls to the above implementation functions 
// I did this so that the error messages can contain code location info
#define ecalloc(ptr, size) \
	ecalloc_impl((ptr), (size), __FILE__, __LINE__, __func__)

#define emalloc(size) \
	emalloc_impl((size), __FILE__, __LINE__, __func__)

#define erealloc(ptr, size) \
        erealloc_impl((ptr), (size), __FILE__, __LINE__, __func__)

#define efree(ptr)\
	efree_impl((void**)&(ptr), __FILE__, __LINE__, __func__)


#ifdef __cplusplus
}
#endif

#endif
