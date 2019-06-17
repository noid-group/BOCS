/**
@file io_read.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
*/

#ifndef IO_READ
#define IO_READ

#include "cgff_types.h"

/* MRD mimicking Gromacs' command line input */
enum {etINT, etBOOL, etSTRING, etREAL};

typedef struct t_pargs
{
  char *flag;
  bool bSet;
  int type;
  tW_word value;
  char *desc;
  bool bMandatory;
} t_pargs;

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

/* Gromacs-esque command line args */

const char * get_arg_type(t_pargs par);

void init_arg(t_pargs *param, const char * flag, int type, const char * desc, bool bMan);

void init_arg_def(t_pargs *param, const char * flag, int type, const char *desc, const char * init_val);

int get_command_line_args(int argc, char * argv[], int n_args, t_pargs *params);

bool print_arg_table(int n_args, t_pargs *params, bool HELP);

void check_mand_args(t_pargs *params, int nPars);

void build_filename(tW_word fnm, tW_word ipt_name, int ft, const char *def_end);

/****************************************************************************************/


#endif
