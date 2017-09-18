/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file gmx-interface.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef gromacs_interface
#define gromacs_interface


#if GMX == 45

// Gromacs libraries we want to keep
#include "typedefs.h"
#include "tpxio.h"		// for read_tps_conf
#include "statutil.h"		// for TRX_READ_XVF and trx_frames
#include "copyrite.h"
#include "tpxio.h"
#include "filenm.h"

// we want to remove these gmx dependencies:
#include "vec.h"		// for copy_dvec and math
#include "smalloc.h"		// for snew, sfree

#endif

// MRD
#if GMX == 46

#include "typedefs.h"
#include "tpxio.h"
#include "statutil.h"
#include "copyrite.h"
#include "tpxio.h"
#include "filenm.h"
#include "vec.h"
#include "smalloc.h"

#endif

#if GMX == 50

#include "legacyheaders/typedefs.h"
#include "fileio/tpxio.h"
#include <stdio.h>
#include "legacyheaders/readinp.h"
#include "fileio/pdbio.h"
#include "legacyheaders/types/oenv.h"
#include "legacyheaders/oenv.h"
#include "fileio/gmxfio.h"
#include "legacyheaders/copyrite.h"
#include "fileio/filenm.h"
#include "legacyheaders/vec.h"
#include "utility/smalloc.h"
#include "fileio/trxio.h"
#include "commandline/pargs.h"
#include "legacyheaders/macros.h"

#endif

#if GMX == 51

#include "legacyheaders/typedefs.h"
#include "fileio/tpxio.h"
#include <stdio.h>
#include "legacyheaders/readinp.h"
#include "fileio/pdbio.h"
#include "legacyheaders/types/oenv.h"
#include "legacyheaders/oenv.h"
#include "fileio/gmxfio.h"
#include "legacyheaders/copyrite.h"
#include "fileio/filenm.h"
#include "math/vec.h"
#include "utility/smalloc.h"
#include "fileio/trxio.h"
#include "commandline/pargs.h"


#include "legacyheaders/macros.h"

#endif

// MRD


#include "cgff_types.h"

typedef struct {
    bool b_Gromacs;		// struct. stored in gmx files ?
    bool b_Forces;		// forces present ?
    bool b_PBC;			// pbc ?
    double box[DIM][DIM];	// box for pbc
    bool b_Forces_1;		// forces for any struct?
    bool b_Forces_N;		// forces for all struct?
} tW_gmx_info;


typedef struct {
    int *bond_type_i;
    int *angle_type_i;
    int *dihedral_type_i;
    int *pair_type_i;
} tW_bond_type_indices;


/****************************************************************************************/
/* Structures and functions for gmx command line output *********************************/
/****************************************************************************************/

typedef struct tW_gmx_output {

	FILE *outstream;
	const char *prog;

	void (*copyright)(struct tW_gmx_output*);
	void (*thanks)(struct tW_gmx_output*);

} tW_gmx_output;

void copyright(tW_gmx_output *self);
void thanks(tW_gmx_output *self);

tW_gmx_output* init_tW_gmx_output(FILE *outstream, const char *prog);


/****************************************************************************************/
/* Structures and functions for t_filenames interface ***********************************/
/****************************************************************************************/

typedef struct tW_gmx_cgmap_input {

	#define N_ARGS_SOUGHT 5

	int N_files;

	#if GMX == 45
	output_env_t *oenv;
	t_filenm fnm[N_ARGS_SOUGHT];
	#elif GMX == 46
	output_env_t *oenv;
	t_filenm fnm[N_ARGS_SOUGHT];
        #elif GMX == 50
	output_env_t *oenv;
	t_filenm fnm[N_ARGS_SOUGHT];
        #elif GMX == 51 // MRD
	output_env_t *oenv; // MRD
	t_filenm fnm[N_ARGS_SOUGHT]; // MRD
	#endif


	const char *(*get_filename)(struct tW_gmx_cgmap_input*, const char*);
	bool (*arg_is_set)(struct tW_gmx_cgmap_input*, const char*);

} tW_gmx_cgmap_input;

const char *get_filename(tW_gmx_cgmap_input *self, const char *opt);
bool arg_is_set(tW_gmx_cgmap_input *self, const char *opt);

tW_gmx_cgmap_input* init_tW_gmx_cgmap_input(int argc, char *argv[]);


/****************************************************************************************/
/* Structures and functions for dummy t_filenames interface *****************************/
/****************************************************************************************/

//typedef struct tW_gmx_cgff_input {
//
//	output_env_t *oenv;
//
//} tW_gmx_cgff_input;
//
//tW_gmx_cgff_input* init_tW_gmx_cgff_input();

/****************************************************************************************/
/* Structures and functions for t_topology interface ************************************/
/****************************************************************************************/

typedef struct tW_gmx_topology {

	#if GMX == 45
	t_topology *contents;
	#elif GMX == 46
	t_topology *contents;
	#elif GMX == 50
	t_topology *contents;
	#elif GMX == 51 // MRD
	t_topology *contents; // MRD
	#endif

	bool b_tpr;

	int (*get_natoms)(struct tW_gmx_topology*);
	char *(*get_name)(struct tW_gmx_topology*);

	int (*get_nbonds)(struct tW_gmx_topology*);
	int (*get_nangles)(struct tW_gmx_topology*);
	int (*get_ndihs)(struct tW_gmx_topology*);
	int (*get_npairs)(struct tW_gmx_topology*);

	int *(*get_bond_list)(struct tW_gmx_topology*);
	int *(*get_angle_list)(struct tW_gmx_topology*);
	int *(*get_dih_list)(struct tW_gmx_topology*);
	int *(*get_pair_list)(struct tW_gmx_topology*);

	int (*get_resind)(struct tW_gmx_topology*, int);
	char *(*get_atomtype)(struct tW_gmx_topology*, int);
	int (*get_nexcl)(struct tW_gmx_topology*, int);
	int *(*get_excl_list)(struct tW_gmx_topology*, int);


} tW_gmx_topology;

int get_natoms(tW_gmx_topology *self);
char *get_name(tW_gmx_topology *self);

int get_nbonds(tW_gmx_topology *self);
int get_nangles(tW_gmx_topology *self);
int get_ndihs(tW_gmx_topology *self);
int get_npairs(tW_gmx_topology *self);

int *get_bond_list(tW_gmx_topology *self);
int *get_angle_list(tW_gmx_topology *self);
int *get_dih_list(tW_gmx_topology *self);
int *get_pair_list(tW_gmx_topology *self);

int get_resind(tW_gmx_topology *self, int i_atom);
char *get_atomtype(tW_gmx_topology *self, int i_atom);
int get_nexcl(tW_gmx_topology *self, int i_atom);
int *get_excl_list(tW_gmx_topology *self, int i_atom);

tW_gmx_topology* init_tW_gmx_topology();

void read_tpr(tW_gmx_topology *top, const char *tpr_fnm);
bool get_top(tW_word coord_fnm, tW_gmx_topology *top, bool TPR_flag, tW_word tpr_filename);
int print_tpr_file(FILE * fp_log, tW_gmx_topology *top);
void get_site_info(tW_CG_site *CG_struct, tW_system *sys,  tW_gmx_topology *top);

/****************************************************************************************/
/* Structures and functions for t_topology interface ************************************/
/****************************************************************************************/

typedef struct tW_gmx_trxframe {

	#if GMX == 45
	t_trxframe *contents;
	t_trxstatus *status;
	output_env_t oenv;
	#elif GMX == 46
	t_trxframe *contents;
	t_trxstatus *status;
	output_env_t oenv;
	#elif GMX == 50
	t_trxframe *contents;
	t_trxstatus *status;
	output_env_t *oenv;
	#elif GMX == 51 // MRD
	t_trxframe *contents; // MRD
	t_trxstatus *status; // MRD
	output_env_t *oenv; // MRD
	#endif


	void (*set_natoms)(struct tW_gmx_trxframe*, int);
	void (*set_atom_labels)(struct tW_gmx_trxframe*, struct tW_gmx_topology*);
	void (*print_frame_dump)(struct tW_gmx_trxframe*);

	void (*get_pos_of)(struct tW_gmx_trxframe*, int, dvec);
	void (*get_vel_of)(struct tW_gmx_trxframe*, int, dvec);
	void (*get_box)(struct tW_gmx_trxframe*, matrix);

} tW_gmx_trxframe;

void set_natoms(tW_gmx_trxframe *self, int natoms);
void set_atom_labels(tW_gmx_trxframe *self, tW_gmx_topology *top);
void print_frame_dump(tW_gmx_trxframe *self);

void get_pos_of(tW_gmx_trxframe *self, int i, dvec x);
void get_vel_of(tW_gmx_trxframe *self, int i, dvec v);
void get_box(tW_gmx_trxframe *self, matrix box);

tW_gmx_trxframe* init_tW_gmx_trxframe();

int read_first_trxframe(tW_gmx_trxframe *frame, const char *trx_fnm);
bool read_next_trxframe(tW_gmx_trxframe *frame);
void open_new_trxfile(tW_gmx_trxframe *frame, const char *trx_fnm);
void write_trxframe_to_file(tW_gmx_trxframe *frame);
void copy_trxframe_info(tW_gmx_trxframe *source, tW_gmx_trxframe *dest);
void map_trxframe(tW_gmx_trxframe *source, tW_gmx_trxframe *dest, tW_site_map *map);

int update_info_trr(tW_gmx_info * info, tW_gmx_trxframe *fr);


bool copy_trr_2_CGstruct(tW_gmx_trxframe *fr, tW_CG_site CG_struct[]);
bool copy_trr_2_CGstruct_ref(tW_gmx_trxframe *fr, tW_CG_site CG_struct[]);


bool read_trr_2_CGstruct(tW_gmx_info * info, tW_gmx_trxframe *fr, tW_CG_site CG_struct[]);	
bool read_trr_2_CGstruct_ref(tW_gmx_info * info, tW_gmx_trxframe *fr, tW_CG_site CG_struct[]);


/***************************************************************************************
 Functions that will live somewhere else in this file
****************************************************************************************/


// tW_trxframe
//int open_trr_file(char *trr_fnm, output_env_t oenv, t_trxstatus ** status, tW_gmx_info * info, t_trxframe * fr);	// 4.5.3



/* JFR - 07.16.12: for reading reference forces from a trr file */

/****************************************************************************************/
/*Functions called outside of this file**************************************************/
/****************************************************************************************/



/****************************************************************************************/
/*Functions associated with get_top()****************************************************/
/****************************************************************************************/

//bool read_tpr_file(char *tpr_fnm, t_topology * top, matrix box);

/****************************************************************************************/
/*Functions associated with setup_CG_struct()********************************************/
/****************************************************************************************/












/****************************************************************************************/
/*Functions associated with open_trr_file()**********************************************/
/****************************************************************************************/



/****************************************************************************************/

#endif
