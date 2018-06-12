/**
@file io_read.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef IO_READ
#define IO_READ

#include "cgff_types.h"

/****************************************************************************************/

int file_exists(const char *filename);

int get_next_line(FILE * fp_in, tW_line inp_line);

int get_clean_line(FILE * fp_in, tW_line inp_line);

int remove_comments(char *before, char *after);

FILE *open_file(const tW_word file_name, char mode);

int match_word_list(int n_words, tW_word list_1[], tW_word list_2[]);

void delete_directive_from_line(tW_line inp_line);

tW_word *get_word_list(FILE * fp_in, tW_word keyword, int N_words);

void test_end_directive(tW_line inp_line, tW_word keyword,
			tW_word input_file);

void check_inp_line(tW_word input_file, tW_word keyword, tW_line inp_line);

void print_line_stars(FILE * fp);

int check_word_list(int n_words, tW_word list_1[], tW_word list_2[]);

int match_word(int N_words, tW_word word, tW_word * word_list);

/****************************************************************************************/


#endif
