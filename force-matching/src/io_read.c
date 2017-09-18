/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file io_read.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Functions related to cgff and cgmap input 
*/

//c library includes
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//local includes
#include "cgff_types.h"
#include "gmx-interface.h"
#include "io_read.h"
#include "safe_mem.h"


/*****************************************************************************************
file_exists(): Checks whether the supplied filename can be opened for reading.  Returns
true if the file exists, returns false otherwise.
*****************************************************************************************/
int file_exists(const char *filename)
{
    FILE *file;
    if ((file = fopen(filename, "r"))) {
	fclose(file);
	return TRUE;
    }

    return FALSE;
}



/*****************************************************************************************
get_next_line(): Stores the next non-blank (excluding comments) line from fp_in into the
variable inp_line. Returns the number of characters in the line (including a newline
character if present). If the end of file is reached, l_inp is -1. The sscanf() call
tests if anything is on the line and the if statement ensures that blank lines are
ignored.
*****************************************************************************************/
int get_next_line(FILE * fp_in, tW_line inp_line)
{
    int l_inp;
    int test_line;
    tW_word test_word;

    do {
	l_inp = get_clean_line(fp_in, inp_line);

	test_line = sscanf(inp_line, "%s", test_word);
	if ((test_line == EOF) && (l_inp != -1)) {
	    l_inp = 0;
	}

    } while (l_inp == 0);

    return l_inp;
}


/*****************************************************************************************
get_clean_line(): Reads the next line of fp_in and removes commented sections from the
line. Copies the uncommented section into inp_line. Comments are define as all text after
the first element of COMMENTFLAG and are stripped using the function remove_comments().
Returns -1 if EOF has been reached and otherwise returns the length of the uncommented
string which may include the newline character.
*****************************************************************************************/
int get_clean_line(FILE * fp_in, tW_line inp_line)
{
    int l_inp = -1;

    char before[MAXLINELEN], after[MAXLINELEN];

    strcpy(inp_line, "");
    strcpy(before, "");
    strcpy(after, "");

    if (fgets(before, MAXLINELEN, fp_in) == NULL) {
	return -1;
    }

    l_inp = remove_comments(before, after);

    if (l_inp > 0) {
	strcpy(inp_line, after);
    }

    return l_inp;
}


/*****************************************************************************************
remove_comments(): Searches for characters that denote comments (COMMENTFLAG). Provides
a pointer to a string that has removed all sections after the first COMMENTFLAG. Returns
the length of the uncommented string which will include the new line character if no
comments are present.
*****************************************************************************************/
int remove_comments(char *before, char *after)
{
    int l_inp;

    l_inp = strcspn(before, COMMENTFLAG);

    strncpy(after, before, l_inp);

    after[l_inp] = '\0';

    return l_inp;
}


/****************************************************************************************
open_file(): Opens a file named "file_name" with option "mode." Returns the pointer to
the file. Prints error messages if file could not be opened and terminates execution.
****************************************************************************************/
FILE *open_file(tW_word file_name, char mode)
{
    FILE *fp = NULL;

    if (mode == 'r') {
	if ((fp = fopen(file_name, "r"))) {
	    return fp;
	}
    }

    else if (mode == 'w') {
	if ((fp = fopen(file_name, "w"))) {
	    return fp;
	}
    }

    printf("\nERROR: Cannont open %s for mode '%c'.\n", file_name, mode);
    exit(EXIT_FAILURE);

}


/*****************************************************************************************
match_word_list(): Determines if list_1 and list_2 match in either normal or reverse
order.
*****************************************************************************************/
int match_word_list(int n_words, tW_word list_1[], tW_word list_2[])
{
    int i;
    int test = 0;

    /* Do they match exactly? */
    for (i = 0; i < n_words; i++) {
	if (strcmp(list_1[i], list_2[i]) != 0) {
	    break;
	}
	test++;
    }
    if (test == n_words) {
	return TRUE;
    }

    /* Check reverse order. */
    for (i = 0; i < n_words; i++) {
	if (strcmp(list_1[i], list_2[n_words - i - 1]) != 0) {
	    return FALSE;
	}
    }
    return TRUE;
}



/*****************************************************************************************
delete_directive_from_line(): For lines that contain a directive marked by [...], the
ending brace is found and the line is truncated to the character to the right of the ']'.
It is assumed that the first occurance of ']' marks the end of the directive name.
*****************************************************************************************/
void delete_directive_from_line(tW_line inp_line)
{
    int flag;
    char *string;

    flag = 0;
    do {
	if (inp_line[flag] == ']') {
	    break;
	}
	flag++;
    } while (1);

    string = &(inp_line[flag + 1]);

//    strcpy(inp_line, string);
//  this is OK here because we're removing characters from inp_line, so it 
//  will never be longer than a tW_line should be as a result of our actions
//  here
    memmove(inp_line, string, strlen(string) + 1);

}



/*****************************************************************************************
get_word_list(): Reads a list of N_words words, each on separate lines, following a line
that contains [ keyword ]. Makes sure that the next non-blank line after the word list
contains [ End keyword ] to mark the end of the section. Returns the word list.
*****************************************************************************************/
tW_word *get_word_list(FILE * fp_in, tW_word keyword, int N_words)
{
    int i = 0;
    tW_line inp_line;
    tW_word *word_list;

    word_list = (tW_word *) emalloc(N_words * sizeof(tW_word));

    for (i = 0; i < N_words; i++) {
	get_next_line(fp_in, inp_line);

	check_inp_line("par.txt", keyword, inp_line);

	sscanf(inp_line, "%s", word_list[i]);
    }

    get_next_line(fp_in, inp_line);

    test_end_directive(inp_line, keyword, "par.txt");

    return word_list;
}


/*****************************************************************************************
test_end_directive(): Test to make sure the ending directive is found. Removes ']' if
no space between keyword and ending ']'.
*****************************************************************************************/
void test_end_directive(tW_line inp_line, tW_word keyword,
			tW_word input_file)
{
    tW_word test_keyword;

    strcpy(test_keyword, "\0");

    sscanf(inp_line, " [ End %s ]", test_keyword);

    if (test_keyword[strlen(test_keyword) - 1] == ']') {
	test_keyword[strlen(test_keyword) - 1] = '\0';
    }

    if ((strcmp(keyword, test_keyword) != 0)) {
	printf("\nERROR: Did not find End %s directive where expected.\n",
	       keyword);
	printf("  Found: \"%s\".\n", inp_line);
	printf("  Probably too many types listed under [ %s ] in %s.\n",
	       keyword, input_file);
	exit(EXIT_FAILURE);
    }
}


/*****************************************************************************************
check_inp_line(): Makes sure that the ending directive is not found were it is not 
supposed to be.
*****************************************************************************************/
void check_inp_line(tW_word input_file, tW_word keyword, tW_line inp_line)
{

    if (strstr(inp_line, "End") != NULL) {
	printf("\nERROR: Too few %s are listed in %s.\n", keyword,
	       input_file);
	exit(EXIT_FAILURE);
    }

}


/*****************************************************************************************
print_line_stars(): Prints a line of '*'. Useful for making output easier to read.
*****************************************************************************************/
void print_line_stars(FILE * fp)
{
    int i;

    for (i = 0; i < 91; i++) {
	fprintf(fp, "*");
    }
    fprintf(fp, "\n\n");

}


/*****************************************************************************************
check_word_list(): Checks to see if the n_words words in list1[] occur in list_2[] in 
any order. The flag bX_2[] ensures that there is a 1:1 matching between elements in 
list_1 and elements in list_2. (I think that this function is not realy useful anymore.)
*****************************************************************************************/
int check_word_list(int n_words, tW_word list_1[], tW_word list_2[])
{
    int i, j;
    bool bX_2[n_words];

    for (i = 0; i < n_words; i++) {
	bX_2[i] = TRUE;
    }

    for (i = 0; i < n_words; i++) {
	for (j = 0; j < n_words; j++) {
	    if ((strcmp(list_1[i], list_2[j]) == 0) && (bX_2[j])) {
		bX_2[j] = FALSE;
		break;
	    }
	}

	if (j == n_words) {
	    return FALSE;
	}
    }

    return TRUE;
}


/*****************************************************************************************
match_word(): Searches through N_words words in word_list to see if the string stored
in word is included in the list. Returns the index to word_list if there's a match
or returns -1 if there is no match.
*****************************************************************************************/
int match_word(int N_words, tW_word word, tW_word * word_list)
{
    int i, b_COMP;

    for (i = 0; i < N_words; i++) {
	b_COMP = strcmp(word, word_list[i]);

	if (b_COMP == 0) {
	    return i;
	}
    }

    return -1;
}
