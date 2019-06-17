
#ifndef GROMACS_TOPOLOGY_H
#define GROMACS_TOPOLOGY_H

#include <stdlib.h>
#include <stdio.h>
#include "rpc/xdr.h"

#include "cgff_types.h"
#include "wnoid_math.h"

/* ~ GROMACS #definitions ~ */

//need these from /gromacs/topology/idef.h
#define MAXATOMLIST 6
#define MAXFORCEPARAM   12
#define NR_RBDIHS   6
#define NR_CBTDIHS   6
#define NR_FOURDIHS     4

//need this from /gromacs/math/utilities.h
#define M_PI	3.14159265358979323846

//need this from /gromacs/math/vectypes.h
#define XX 0
#define YY 1
#define ZZ 2

// gromacs/fileio/trxio.h

 /* For trxframe.flags, used in trxframe read routines.
  * When a READ flag is set, the field will be read when present,
  * but a frame might be returned which does not contain the field.
  * When a NEED flag is set, frames not containing the field will be skipped.
  */

/* MRD while we don't actually use any of these, I copied them over
 * in case someone wants to in the future
 */

#define TRX_READ_X    (1<<0)
#define TRX_NEED_X    (1<<1)
#define TRX_READ_V    (1<<2)
#define TRX_NEED_V    (1<<3)
#define TRX_READ_F    (1<<4)
#define TRX_NEED_F    (1<<5)
/* Useful for reading natoms from a trajectory without skipping */
#define TRX_DONT_SKIP (1<<6)

/* For trxframe.not_ok */
#define HEADER_NOT_OK (1<<0)
#define DATA_NOT_OK   (1<<1)
#define FRAME_NOT_OK  (HEADER_NOT_OK | DATA_NOT_OK)
     
/* ~ End of GROMACS #definitions ~ */

enum { epbcXYZ, epbcNONE, epbcXY, epbcSCREW, epbcNR };

#define BOND_DIV 3
#define ANGLE_DIV 4
#define PDIH_DIV 5
#define RBDIH_DIV 5
#define TABDIH_DIV 5
#define LJ14_DIV 3

// structs relocated from gmx-interface.h

typedef struct {
    bool b_Gromacs;             // struct. stored in gmx files ?
    bool b_Forces;              // forces present ?
    bool b_PBC;                 // pbc ?
    double box[DIM][DIM];       // box for pbc
    bool b_Forces_1;            // forces for any struct?
    bool b_Forces_N;            // forces for all struct?
} tW_gmx_info;


typedef struct {
    int *bond_type_i;
    int *angle_type_i;
    int *dihedral_type_i;
    int *pair_type_i;
} tW_bond_type_indices;


/* Redone stuff for cgmap */

// GROMACS seems to have a spiffy method for parsing the command line flags/arguments
// Instead of having to look for each flag by hand and doing the appropriate things,
// I grouped the arguments you can pass into cgmap into three groups:
// files, time, and file types.
// I then wrote stupid little data structures and enums for each.
// In a loop over the command line args, I was then able to loop over the three different
// types of arguments, check for the flag, and do the appropriate stuff with the corresponding parameter

#define N_MAP_FILES 5
#define N_TIME_ARGS 3
#define N_FTYPE_ARGS 3

#define DUMP "dump"
#define BOCS "bocs"
#define GRO "gro"
#define TRJ "trj"
#define TRR "trr"
#define LMP "lmp"
#define LMPDATA "data"
#define LMPTRJ "lmp"

#define TOP "top"
#define TRAJ "traj"

enum {eDUMP, eBOCS, eGRO, eTRJ, eTRR, eLMP, eLMPDATA, eLMPTRJ};
enum {eTOP, eTRAJ};

enum {eAA, eCG, eCGTPR, eMAP, eGRO1};
enum {eBTIME, eETIME, eDELTAT};

/*
the flag is the thing you put on the command line to specify that argument.
	e.g. '-f' for a file

the shift is what I use to keep track of whether or not a flag has been specified
they work with the 3 unsigned shorts in the tW_gmx_cgmap_input struct
each struct has its own identifier, most conveniently represented in BINARY
they're all gonna be something like (1 << n) where n is 0-4 (for files, since there are 5 files)
which correspond to (in BINARY) 00000001 00000010 00000100 00001000 00010000
when they pass in, e.g. '-f',
i first check to see if that flag has already been seen:
(use_files & files[idx].file_shift) will evaluate to 0 if i have not seen that flag yet
& is the bitwise AND. it takes two BINARY numbers and compares each digit.
if a given digit in both numbers is 1, then the result has a 1 in that digit
if a given digit is 0 in either number, then the result has a 0 in that digit
once i have verified this is the first time I see that flag, i do 
use_files |= files[idx].file_shift
|= is the bitwise OR. it takes two BINARY numbers and compares each digit.
if a given digit in either number is 1, then the result has a 1 in that digit
if a given digit is 0 in both numbers, then the result has a 0 in that digit
I find this method of keeping track of what flags we've seen more elegant than
simply declaring a boolean variable for each flag. 

i dont think I actually use the idx... I should get rid of it

and then last is the corresponding parameter that follows the flag
*/
typedef struct tW_file_struct
{
  char *flag; 
  unsigned short file_shift; 
  int idx;
  tW_word fnm;
} tW_file_struct;

typedef struct tW_time_struct
{
  char *flag;
  unsigned short time_shift;
  int idx;
  float value;
} tW_time_struct;

typedef struct tW_fType_struct
{
  char *flag;
  unsigned short shift;
  int value;
} tW_fType_struct;

typedef struct tW_gmx_cgmap_input
{
  char **map_file_flags;
  char **time_flags;
  char **file_type_flags;

  tW_file_struct files[N_MAP_FILES];
  tW_time_struct times[N_TIME_ARGS];
  tW_fType_struct file_types[N_FTYPE_ARGS];

  unsigned short use_files;
  unsigned short file_times;  
  unsigned short ftflags;

  char *(*get_filename)(struct tW_gmx_cgmap_input*, const char*);
  bool (*arg_is_set)(struct tW_gmx_cgmap_input*, const char*);

} tW_gmx_cgmap_input;

char *get_filename(tW_gmx_cgmap_input *self, const char *opt);
bool arg_is_set(tW_gmx_cgmap_input *self, const char *opt);
tW_gmx_cgmap_input* init_tW_gmx_cgmap_input(int argc, char *argv[]);

/* end stuff for cgmap */


// A few of these are used as the index for getting counts of stuff, 
// e.g. top->contents->idef.il[F_BONDS].nr 
// from /gromacs/topology/idef.h
enum {
    F_BONDS,
    F_G96BONDS,
    F_MORSE,
    F_CUBICBONDS,
    F_CONNBONDS,
    F_HARMONIC,
    F_FENEBONDS,
    F_TABBONDS,
    F_TABBONDSNC,
    F_RESTRBONDS,
    F_ANGLES,
    F_G96ANGLES,
    F_RESTRANGLES,
    F_LINEAR_ANGLES,
    F_CROSS_BOND_BONDS,
    F_CROSS_BOND_ANGLES,
    F_UREY_BRADLEY,
    F_QUARTIC_ANGLES,
    F_TABANGLES,
    F_PDIHS,
    F_RBDIHS,
    F_RESTRDIHS,
    F_CBTDIHS,
    F_FOURDIHS,
    F_IDIHS,
    F_PIDIHS,
    F_TABDIHS,
    F_CMAP,
    F_GB12,
    F_GB13,
    F_GB14,
    F_GBPOL,
    F_NPSOLVATION,
    F_LJ14,
    F_COUL14,
    F_LJC14_Q,
    F_LJC_PAIRS_NB,
    F_LJ,
    F_BHAM,
    F_LJ_LR,
    F_BHAM_LR,
    F_DISPCORR,
    F_COUL_SR,
    F_COUL_LR,
    F_RF_EXCL,
    F_COUL_RECIP,
    F_LJ_RECIP,
    F_DPD,
    F_POLARIZATION,
    F_WATER_POL,
    F_THOLE_POL,
    F_ANHARM_POL,
    F_POSRES,
    F_FBPOSRES,
    F_DISRES,
    F_DISRESVIOL,
    F_ORIRES,
    F_ORIRESDEV,
    F_ANGRES,
    F_ANGRESZ,
    F_DIHRES,
    F_DIHRESVIOL,
    F_CONSTR,
    F_CONSTRNC,
    F_SETTLE,
    F_VSITE2,
    F_VSITE3,
    F_VSITE3FD,
    F_VSITE3FAD,
    F_VSITE3OUT,
    F_VSITE4FD,
    F_VSITE4FDN,
    F_VSITEN,
    F_COM_PULL,
    F_EQM,
    F_EPOT,
    F_EKIN,
    F_ETOT,
    F_ECONSERVED,
    F_TEMP,
    F_VTEMP_NOLONGERUSED,
    F_PDISPCORR,
    F_PRES,
    F_DVDL_CONSTR,
    F_DVDL,
    F_DKDL,
    F_DVDL_COUL,
    F_DVDL_VDW,
    F_DVDL_BONDED,
    F_DVDL_RESTRAINT,
    F_DVDL_TEMPERATURE, /* not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE) */
    F_NRE               /* This number is for the total number of energies      */
};

/* Simple GROMACS typedefs were redefined with a tW_ prefix. 
 * This was done so that the compiler didn't throw a fit due  
 * to numerous redefinitions while I was in the process of
 * removing our dependency on the GROMACS library.
 */

typedef int tW_gmx_bool;
typedef double tW_real;
typedef tW_real tW_rvec[DIM];
typedef double tW_dvec[DIM];
typedef tW_real tW_matrix[DIM][DIM];
typedef tW_real tW_tensor[DIM][DIM];
typedef int tW_ivec[DIM];
typedef int tW_imatrix[DIM][DIM];


/*
  This structure has a bunch of topology information for each molecule
  I use this as an intermediate when reading dumped topology files or .btp (bocstop) files.
  I read the file in, populate the stuff in this molecule struct, and once the file is done,
  I populate the stuff in tW_t_topology.
  I don't directly populate tW_t_topology because it has information for EACH molecule, whereas
  the topology files only store information for each molecule TYPE.
  Once i have that information for each molecule TYPE, I loop over the number of those molecules present
  to populate EACH molecule in tW_t_topology
*/
typedef struct tW_molecule{
  tW_word molname; 	// name of this molecule

  int n_apm; 		// number of atoms per molecule
  int *type_ids; 	// array of atom types (dim n_apm)
  double *m;		// array of atomic masses (dim n_apm)
  double *q;		// array of atomic charges (dim n_apm)
  tW_word *atom_names;	// array of atomic names (dim n_apm)
  tW_word *atom_types;  // array of atomic types (dim n_apm)
  tW_word *atom_Btypes; // array of atomic typesB (dim n_apm)
  int *residx; 		// array of residue index to which atom belongs (dim n_apm)

  int n_mols; 		// number of molecules

  int n_res; 		// number of "residues" in this molecule
  tW_word *resname;	// array of residue names (dim nres)

  int n_cg; 		// number of charge groups
  int *cg_start; 	// array of charge group starting indices (dim n_cg)
  int *cg_end;		// array of charge group ending indices (dim n_cg)

  int n_excls; 		// number of excl's
  int n_exclsa; 	// number of excla's
  int *n_epa;		// array of number of excls for each atom (dim n_apm)
  int **excls;		// array of lists of excls for each atom (dim n_apm x n_epa) - NOT NEC. RECTANGULAR

  int bond_nr;		// the value given in Bond: nr: X
  int n_bonds;		// number of bonds
  int *bond_types;	// types of bonds (dim n_bonds)
  int **bond_ij;	// array of atoms i and j for bonds (dim n_bonds x 2)

  int angle_nr;		// the value given in Angle: nr: X
  int n_angles;		// number of angles
  int *angle_types;	// types of angles (dim n_angles)
  int **angle_ijk;	// array of atoms i j and k for angles (dim n_angles x 3)
  
  int pdih_nr;		// the value given in Proper Dih.: nr: X
  int n_pdihs;		// number of proper dihedrals
  int *pdih_types;	// types of proper dihedrals (dim n_pdihs)
  int **pdih_ijkl;	// array of atoms i j k and l for pdihs (dim n_pdihs x 4)

  int rbdih_nr;		// the value given in RBdih: nr: X
  int n_rbdihs;		// number of rbdihs
  int *rbdih_types;	// types of rbdihs (dim n_rbdihs)
  int **rbdih_ijkl;	// array of atoms i j k and l for rbdihs (dim n_rbdihs x 4)

  int tabdih_nr;	// the value given in Tab. Dih.: nr: X
  int n_tabdihs;	// number of tabulated dihedrals
  int *tabdih_types;	// types of tabulated dihedrals (dim n_tbdihs)
  int **tabdih_ijkl;	// array of atoms i j k and l for tab dihs (dim n_tbdihs x 4)

  int lj14_nr;		// the value given in Lj14: nr: X
  int n_lj14s;		// number of lj14s
  int *lj14_types;	// types of lj14s (dim n_lj14s)
  int **lj14_ij;	// array of atoms i and j for lj14s (dim n_lj14s x 2)

} tW_molecule;

/* When our software linked against the GROMACS libraries, we (obviously) heavily
 * used the gromacs data structure t_topology. Our software was written with
 * "wrappers" to their data structures and functions. Although we no longer use 
 * their functions (save clear_dvec...), I copied the necessary type and structure
 * definitions from GROMACS and put them here.
 *
 * First, we have all of the types and structs required for t_topology.
 * t_topology contains elements of types: t_idef, t_atoms, t_atomtypes, t_block,
 * gmx_bool, t_blocka, and t_symtab. So, those struct definitons come before 
 * t_topology's definition, and any definitions required for those types preceed them.
 *
 * As with the simple gromacs typedefs above, a prefix of tW_ was appended to all
 * of the names, for the same reason as above.
 * 
 * The data structures are relatively untouched, other than appending a tW prefix to
 * their names and the types of their elements, even though we don't use all parts of them
 */


typedef int tW_t_functype;

typedef union tW_t_iparams
{
    /* Some parameters have A and B values for free energy calculations.
     * The B values are not used for regular simulations of course.
     * Free Energy for nonbondeds can be computed by changing the atom type.
     * The harmonic type is used for all harmonic potentials:
     * bonds, angles and improper dihedrals
     */
    struct { 
        tW_real a, b, c;
    } bham;
    struct { 
        tW_real rA, krA, rB, krB;
    } harmonic;
    struct { 
        tW_real klinA, aA, klinB, aB;
    } linangle;
    struct { 
        tW_real lowA, up1A, up2A, kA, lowB, up1B, up2B, kB;
    } restraint;
    /* No free energy supported for cubic bonds, FENE, WPOL or cross terms */
    struct { 
        tW_real b0, kb, kcub;
    } cubic;
    struct { 
        tW_real bm, kb;
    } fene;
    struct { 
        tW_real r1e, r2e, krr;
    } cross_bb;
    struct { 
        tW_real r1e, r2e, r3e, krt;
    } cross_ba;
    struct { 
        tW_real thetaA, kthetaA, r13A, kUBA, thetaB, kthetaB, r13B, kUBB;
    } u_b;
    struct { 
        tW_real theta, c[5];
    } qangle;
    struct { 
        tW_real alpha;
    } polarize;
    struct { 
        tW_real alpha, drcut, khyp;
    } anharm_polarize;
    struct { 
        tW_real al_x, al_y, al_z, rOH, rHH, rOD;
    } wpol;
    struct { 
        tW_real a, alpha1, alpha2, rfac;
    } thole;
    struct { 
        tW_real c6, c12;
    } lj;
    struct { 
        tW_real c6A, c12A, c6B, c12B;
    } lj14;
    struct { 
        tW_real fqq, qi, qj, c6, c12;
    } ljc14;
    struct { 
        tW_real qi, qj, c6, c12;
    } ljcnb;
    /* Proper dihedrals can not have different multiplicity when
     * doing free energy calculations, because the potential would not
     * be periodic anymore.
     */
    struct { 
        tW_real phiA, cpA; int mult; tW_real phiB, cpB;
    } pdihs;
    struct { 
        tW_real dA, dB;
    } constr;
    /* Settle can not be used for Free energy calculations of water bond geometry.
     * Use shake (or lincs) instead if you have to change the water bonds.
     */
    struct {
        tW_real doh, dhh;
    } settle;
    struct {
        tW_real b0A, cbA, betaA, b0B, cbB, betaB;
    } morse;
    struct {
        tW_real pos0A[DIM], fcA[DIM], pos0B[DIM], fcB[DIM];
    } posres;
    struct {
        tW_real pos0[DIM], r, k; int geom;
    } fbposres;
    struct {
        tW_real rbcA[NR_RBDIHS], rbcB[NR_RBDIHS];
    } rbdihs;
    struct {
        tW_real cbtcA[NR_CBTDIHS], cbtcB[NR_CBTDIHS];
    } cbtdihs;
    struct {
        tW_real a, b, c, d, e, f;
    } vsite;
    struct {
        int  n; tW_real a;
    } vsiten;
    /* NOTE: npair is only set after reading the tpx file */
    struct {
        tW_real low, up1, up2, kfac; int type, label, npair;
    } disres;
    struct {
        tW_real phiA, dphiA, kfacA, phiB, dphiB, kfacB;
    } dihres;
    struct {
        int  ex, power, label; tW_real c, obs, kfac;
    } orires;
    struct {
        int  table; tW_real kA; tW_real kB;
    } tab;
    struct {
        tW_real sar, st, pi, gbr, bmlt;
    } gb;
    struct {
        int cmapA, cmapB;
    } cmap;
    struct {
        tW_real buf[MAXFORCEPARAM];
    } generic;                                               /* Conversion */
} tW_t_iparams;

typedef struct
{
    tW_real *cmap; /* Has length 4*grid_spacing*grid_spacing, */
    /* there are 4 entries for each cmap type (V,dVdx,dVdy,d2dVdxdy) */
} tW_gmx_cmapdata_t;

typedef struct tW_gmx_cmap_t
{
    int             ngrid;        /* Number of allocated cmap (cmapdata_t ) grids */
    int             grid_spacing; /* Grid spacing */
    tW_gmx_cmapdata_t *cmapdata;     /* Pointer to grid with actual, pre-interpolated data */
} tW_gmx_cmap_t;

typedef int tW_atom_id;
typedef tW_atom_id tW_t_iatom;

typedef struct tW_t_ilist
{
    int      nr;
    int      nr_nonperturbed;
    tW_t_iatom *iatoms;
    int      nalloc;
} tW_t_ilist;


typedef struct tW_t_idef
{
    int         ntypes;			// how many elements are in functype and iparams
    int         atnr;			// number of atom types
    tW_t_functype *functype;		// dim ntypes. defines type of function to use for every force type. 
					// Note two bonds with different parameters are different force types
    tW_t_iparams  *iparams;
    tW_real        fudgeQQ;
    tW_gmx_cmap_t  cmap_grid;
    tW_t_iparams  *iparams_posres, *iparams_fbposres;
    int         iparams_posres_nalloc, iparams_fbposres_nalloc;

    tW_t_ilist     il[F_NRE];
    int         ilsort;
    int         nthreads;
    int        *il_thread_division;
    int         il_thread_division_nalloc;
} tW_t_idef;

/*
 * The struct t_idef defines all the interactions for the complete
 * simulation. The structure is setup in such a way that the multinode
 * version of the program  can use it as easy as the single node version.
 * General field description:
 *   int ntypes
 *      defines the number of elements in functype[] and param[].
 *   int nodeid
 *      the node id (if parallel machines)
 *   int atnr
 *      the number of atomtypes
 *   t_functype *functype
 *      array of length ntypes, defines for every force type what type of
 *      function to use. Every "bond" with the same function but different
 *      force parameters is a different force type. The type identifier in the
 *      forceatoms[] array is an index in this array.
 *   t_iparams *iparams
 *      array of length ntypes, defines the parameters for every interaction
 *      type. The type identifier in the actual interaction list
 *      (ilist[ftype].iatoms[]) is an index in this array.
 *   gmx_cmap_t cmap_grid
 *      the grid for the dihedral pair correction maps.
 *   t_iparams *iparams_posres, *iparams_fbposres
 *      defines the parameters for position restraints only.
 *      Position restraints are the only interactions that have different
 *      parameters (reference positions) for different molecules
 *      of the same type. ilist[F_POSRES].iatoms[] is an index in this array.
 *   t_ilist il[F_NRE]
 *      The list of interactions for each type. Note that some,
 *      such as LJ and COUL will have 0 entries.
 *   int ilsort
 *      The state of the sorting of il, values are provided above.
 *   int nthreads
 *      The number of threads used to set il_thread_division.
 *   int *il_thread_division
 *      The division of the normal bonded interactions of threads.
 *      il_thread_division[ftype*(nthreads+1)+t] contains an index
 *      into il[ftype].iatoms; thread th operates on t=th to t=th+1.
 *   int il_thread_division_nalloc
 *      The allocated size of il_thread_division,
 *      should be at least F_NRE*(nthreads+1).
 */


typedef struct tW_t_atom
{
    tW_real           m, q;        /* Mass and charge                      */
    tW_real           mB, qB;      /* Mass and charge for Free Energy calc */
    unsigned short type;        /* Atom type                            */
    unsigned short typeB;       /* Atom type for Free Energy calc       */
    int            ptype;       /* Particle type                        */
    int            resind;      /* Index into resinfo (in t_atoms)      */
    int            atomnumber;  /* Atomic Number or NOTSET              */
    char           elem[4];     /* Element name                         */
} tW_t_atom;

typedef struct tW_t_resinfo
{
    char          **name;       /* Pointer to the residue name          */
    int             nr;         /* Residue number                       */
    unsigned char   ic;         /* Code for insertion of residues       */
    int             chainnum;   /* Iincremented at TER or new chain id  */
    char            chainid;    /* Chain identifier written/read to pdb */
    char          **rtp;        /* rtp building block name (optional)   */
} tW_t_resinfo;

typedef struct tW_t_pdbinfo
{
    int      type;              /* PDB record name                      */
    int      atomnr;            /* PDB atom number                      */
    char     altloc;            /* Alternate location indicator         */
    char     atomnm[6];         /* True atom name including leading spaces */
    tW_real     occup;             /* Occupancy                            */
    tW_real     bfac;              /* B-factor                             */
    tW_gmx_bool bAnisotropic;      /* (an)isotropic switch                 */
    int      uij[6];            /* Anisotropic B-factor                 */
} tW_t_pdbinfo;


typedef struct tW_t_atoms
{
    int            nr;          /* Nr of atoms                          */
    tW_t_atom        *atom;        /* Array of atoms (dim: nr)             */
                                /* The following entries will not       */
                                /* always be used (nres==0)             */
    char          ***atomname;  /* Array of pointers to atom name       */
                                /* use: (*(atomname[i]))                */
    char          ***atomtype;  /* Array of pointers to atom types      */
                                /* use: (*(atomtype[i]))                */
    char          ***atomtypeB; /* Array of pointers to B atom types    */
                                /* use: (*(atomtypeB[i]))               */
    int              nres;      /* The number of resinfo entries        */
    tW_t_resinfo       *resinfo;   /* Array of residue names and numbers   */
    tW_t_pdbinfo       *pdbinfo;   /* PDB Information, such as aniso. Bfac */
} tW_t_atoms;

typedef struct tW_t_atomtypes
{
    int           nr;           /* number of atomtypes                          */
    tW_real         *radius;       /* GBSA radius for each atomtype                */
    tW_real         *vol;          /* GBSA efective volume for each atomtype       */
    tW_real         *surftens;     /* implicit solvent surftens for each atomtype  */
    tW_real         *gb_radius;    /* GB radius for each atom type                 */
    tW_real         *S_hct;        /* Overlap factors for HCT/OBC GB models        */
    int          *atomnumber;   /* Atomic number, used for QM/MM                */
} tW_t_atomtypes;

typedef struct tW_t_block
{
    int      nr;           /* The number of blocks          */
    tW_atom_id *index;        /* Array of indices (dim: nr+1)  */
    int      nalloc_index; /* The allocation size for index */
} tW_t_block;

typedef struct tW_t_blocka
{
    int      nr;    /* The number of blocks              */
    tW_atom_id *index; /* Array of indices in a (dim: nr+1) */
    int      nra;   /* The number of atoms               */
    tW_atom_id *a;     /* Array of atom numbers in each group  */
    /* (dim: nra)                           */
    /* Block i (0<=i<nr) runs from          */
    /* index[i] to index[i+1]-1. There will */
    /* allways be an extra entry in index   */
    /* to terminate the table               */
    int nalloc_index;           /* The allocation size for index        */
    int nalloc_a;               /* The allocation size for a            */
} tW_t_blocka;

typedef struct tW_t_symbuf
{
    int               bufsize;
    char            **buf;
    struct tW_t_symbuf  *next;
} tW_t_symbuf;

typedef struct tW_t_symtab
{
    int       nr;
    tW_t_symbuf *symbuf;
} tW_t_symtab;

typedef struct tW_t_topology
{
    tW_word            name;                        /* Name of the topology                 */
    tW_t_idef          idef;                        /* The interaction function definition  */
    tW_t_atoms         atoms;                       /* The atoms                            */
    tW_t_atomtypes     atomtypes;                   /* Atomtype properties                  */
    tW_t_block         cgs;                         /* The charge groups                    */
    tW_t_block         mols;                        /* The molecules                        */
    tW_gmx_bool        bIntermolecularInteractions; /* Inter.mol. int. ?   */
    tW_t_blocka        excls;                       /* The exclusions                       */
    tW_t_symtab        symtab;                      /* The symbol table                     */

//    tW_molecule       *molecules;
//    tW_word	      *force_names;
} tW_t_topology;


/* In additon to the t_topology data structure, we heavily use the t_trxframe
 * data structure. Below I have the definition for t_trxframe, similar to what 
 * I did above for t_topology, including tW prefixes
 *
 * Originally, we used t_trxstatus and gmx_output_env_t as well.
 * I eliminated use of those to make my life easier.
 * 
 * The way GROMACS does things, the FILE variable is in t_trxstatus.
 * Since I eliminated that, I added a FILE variable to tW_gmx_trxframe,
 * which is the container for tW_t_trxframe
 *
 */

typedef struct tW_t_trxframe
{
    int      flags;            /* flags for read_first/next_frame  */
    int      not_ok;           /* integrity flags                  */
    tW_gmx_bool bDouble;          /* Double precision?                */
    int      natoms;           /* number of atoms (atoms, x, v, f) */
    tW_real     t0;               /* time of the first frame, needed  *
                                * for skipping frames with -dt     */
    tW_real     tf;               /* internal frame time - DO NOT CHANGE */
    tW_real     tpf;              /* time of the previous frame, not  */
                               /* the read, but real file frames   */
    tW_real     tppf;             /* time of two frames ago           */
                               /* tpf and tppf are needed to       */
                               /* correct rounding errors for -e   */
    tW_gmx_bool        bTitle;
    char     *title;     /* title of the frame            */
    tW_gmx_bool        bStep;
    int             step;      /* MD step number                   */
    tW_gmx_bool        bTime;
    tW_real            time;      /* time of the frame                */
    tW_gmx_bool        bLambda;
    tW_gmx_bool        bFepState; /* does it contain fep_state?       */
    tW_real            lambda;    /* free energy perturbation lambda  */
    int             fep_state; /* which fep state are we in? */
    tW_gmx_bool        bAtoms;
    tW_t_atoms *atoms;     /* atoms struct (natoms)            */
    tW_gmx_bool        bPrec;
    tW_real            prec;      /* precision of x, fraction of 1 nm */
    tW_gmx_bool        bX;
    tW_rvec           *x;         /* coordinates (natoms)             */
    tW_gmx_bool        bV;
    tW_rvec           *v;         /* velocities (natoms)              */
    tW_gmx_bool        bF;
    tW_rvec           *f;         /* forces (natoms)                  */
    tW_gmx_bool        bBox;
    tW_matrix          box;       /* the 3 box vectors                */
    tW_gmx_bool        bPBC;
    int             ePBC;      /* the type of pbc                  */
//    tW_t_gmxvmdplugin* vmdplugin; // This is causing errors. We don't use it. I got rid of it.
} tW_t_trxframe;

//more from gmx-interface.h

// This is the modified wrapper for dealing with the t_topology data structure.
// All these function pointers likely aren't necesary now that we don't have to 
// deal with multiple versions of GROMACS. However, I kept using them anyways

/****************************************************************************************/
/*                   Structures and functions for t_topology interface                  */
/****************************************************************************************/

typedef struct tW_gmx_topology {

        tW_t_topology *contents; 
	int eFileType;

        tW_molecule       *molecules;
        tW_word           *force_names;
/*

When you dump a tpr file, the start of the topology section (the section we use) is:

topology:
   name="MIX"
   #atoms               = 804
   molblock (0):
      moltype              = 0 "BUT"
      #molecules           = 134
      #atoms_mol           = 2
      #posres_xA           = 0
      #posres_xB           = 0
   molblock (1):
      moltype              = 1 "DEC"
      #molecules           = 134
      #atoms_mol           = 4
      #posres_xA           = 0
      #posres_xB           = 0
   ffparams:
      atnr=2

This atnr value is deceptive. A topology file for a mixed system of CG butane and CG decane begins:

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
  1             1               no               1.0     1.0

[ atomtypes ]
;name       mass        charge   ptype       c6           c12
 CT         29.06200    0.000     A           0            1
 CM         42.08100    0.000     A           0            1

[ bondtypes ]
;  i     j    func    b0   kb
   CT    CT    1       0    1
   CT    CM    1       0    1
   CM    CM    1       0    1

[ angletypes ]
;  i    j    k     func   a0    cth
  CT   CM   CM     1      0     1

[ nonbond_params ]
  ; i      j     func          c6           c12
    CT     CT     1             0           1
    CM     CM     1             0           1
    CT     CM     1             0           1

And then goes on to define the molecule types, where both CT and CM are used.

One would EXPECT atnr to be 2. One would be WRONG. atnr = 1 in the dumped tpr file.

The only way I have found to fix this is by including an extra section at the end of the beginning of the top file.

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
  1             1               no               1.0     1.0

[ atomtypes ]
;name       mass        charge   ptype       c6           c12
 CT         29.06200    0.000     A           0            1
 CM         42.08100    0.000     A           0            1

[ bondtypes ]
;  i     j    func    b0   kb
   CT    CT    1       0    1
   CT    CM    1       0    1
   CM    CM    1       0    1

[ angletypes ]
;  i    j    k     func   a0    cth
  CT   CM   CM     1      0     1

[ nonbond_params ]
  ; i      j     func          c6           c12
    CT     CT     1             0           1
    CM     CM     1             0           1
    CT     CM     1             0           1

[ implicit_genborn_params ]
; atype    radius volume surftens gb_radius S_hct
  CT       0.3    1      1        0.29      0.6
  CM       0.4    1      1        0.38      0.5

including the implicit_genborn_params, atnr becomes 2 in the dumped tpr file!

I have absolutely NO idea what the genborn params are actually used for. But setting them is necessary
if you want to use the atnr value from the dumped tpr file as the number of atom types. otherwise, you'll have
to loop over the atoms in each molecule and see if you've found a new type name.

SO, that is what I do to populate the following two variables.

*/
        tW_word		  *atom_type_names;
	int		  n_atomtypes;
        bool b_tpr;

        int *int_map; // For writing LAMMPS files.

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

bool get_top(tW_word coord_fnm, tW_gmx_topology *top, bool TPR_flag, tW_word tpr_filename);
int print_tpr_file(FILE * fp_log, tW_gmx_topology *top);
void get_site_info(tW_CG_site *CG_struct, tW_system *sys,  tW_gmx_topology *top);


// Analogously, these are for dealing with the t_trxframe data structure

/****************************************************************************************/
/*                  Structures and functions for t_trxframe interface                   */
/****************************************************************************************/


typedef struct tW_gmx_trxframe {

        tW_t_trxframe *contents; 

/* These are the things I had to add once we got rid of GROMACS dependency.
 * the FILE variable is necessary because we don't use t_trxstatus, and that's
 * where the FILE variable is kept in the GROMACS data structures.
 * bInput tells whether you are reading from a file to this frame structure,
 * or writing from this frame structure to a file.
 * eFileType is an eNum for whether the frame in question is being read from/written to
 * a file of type dump, gro, bocs, trj, or trr.
 * counter is to count the number of frames read/written.
 * xdrmode and *xdr are necessary for reading/writing trr files.
 *
 */
        tW_word filename;
        FILE *fp;
        bool bInput;
	int eFileType;
	int counter;
        enum xdr_op xdrmode;
        XDR *xdr;
        bool bDouble;

/* LAMMPS stuff... */
        int *xids, *vids, *fids;
        int atid_id, attype_id, molid_id;

        void (*setup_xdr)(struct tW_gmx_trxframe*, const char *, bool);

        void (*set_natoms)(struct tW_gmx_trxframe*, int);
        void (*set_atom_labels)(struct tW_gmx_trxframe*, struct tW_gmx_topology*);
        void (*print_frame_dump)(struct tW_gmx_trxframe*); 

        void (*get_pos_of)(struct tW_gmx_trxframe*, int, dvec);
        void (*get_vel_of)(struct tW_gmx_trxframe*, int, dvec);
        void (*get_box)(struct tW_gmx_trxframe*, matrix);

} tW_gmx_trxframe;

void setup_xdr(tW_gmx_trxframe *self, const char *fnm, bool bRead);

void set_natoms(tW_gmx_trxframe *self, int natoms);
void set_atom_labels(tW_gmx_trxframe *self, tW_gmx_topology *top);
void print_frame_dump(tW_gmx_trxframe *self);

void get_pos_of(tW_gmx_trxframe *self, int i, dvec x);
void get_vel_of(tW_gmx_trxframe *self, int i, dvec v);
void get_box(tW_gmx_trxframe *self, matrix box);

tW_gmx_trxframe* init_tW_gmx_trxframe();


void copy_trxframe_info(tW_gmx_trxframe *source, tW_gmx_trxframe *dest);
void map_trxframe(tW_gmx_trxframe *source, tW_gmx_trxframe *dest, tW_site_map *map);

        
int update_info_trr(tW_gmx_info * info, tW_gmx_trxframe *fr);

bool copy_trr_2_CGstruct(tW_gmx_trxframe *fr, tW_CG_site CG_struct[]);
bool copy_trr_2_CGstruct_ref(tW_gmx_trxframe *fr, tW_CG_site CG_struct[]);

bool read_trr_2_CGstruct(tW_gmx_info * info, tW_gmx_trxframe *fr, tW_CG_site CG_struct[]);
bool read_trr_2_CGstruct_ref(tW_gmx_info * info, tW_gmx_trxframe *fr, tW_CG_site CG_struct[]);

//////////

void copyright();

// These two functions are used while reading files in
void test_line(tW_line inp_line, const char * term, tW_gmx_bool want, const char * errmsg);
void elim_char(tW_word word, char elim);

// These semi break up the reading in of topology files
void get_moltype_info(FILE *fp, tW_molecule * mol, tW_line * ret_inp_line);
void dump_molecule_info(tW_gmx_topology *top);

tW_gmx_bool read_box_dump_frame(FILE *fp, tW_matrix box);
void copy_atom_props(FILE *fp, tW_rvec *p, int natoms);
void get_prop(FILE *fp, tW_gmx_bool present, tW_rvec *p, char test_prop, int natoms, int step);


void check_count(int iarg, int argc, char * opt);

/* Functions for reading and writing various file types */

bool read_tpr_dump(tW_word fnm, tW_gmx_topology *top);

bool read_bocs_top(tW_word fnm, tW_gmx_topology * top);
void write_bocs_top(tW_word fnm, tW_gmx_topology * top);

/* Ever since I got the read/write trr functionality working, 
 * using dumped trajectories is very discouraged. */
int read_first_dump_frame(tW_gmx_trxframe *frame, const char *trx_fnm);
bool read_next_dump_frame(tW_gmx_trxframe *frame);
void write_dump_frame(tW_gmx_trxframe *fr);

int read_first_bocs_frame(tW_gmx_trxframe *frame, const char *trx_fnm);
bool read_next_bocs_frame(tW_gmx_trxframe *fr);
bool new_read_next_bocs_frame(tW_gmx_trxframe *fr);
void write_bocs_frame(tW_gmx_trxframe * fr);

int read_first_gro_frame(tW_gmx_trxframe *frame, const char *trx_fnm);
bool read_next_gro_frame(tW_gmx_trxframe *fr);
void write_gro_frame(tW_gmx_trxframe *fr);

/* 2 functions for comparing trajectories */

void comp_2_frames(tW_gmx_trxframe *fr1, tW_gmx_trxframe *fr2, tW_gmx_trxframe *fro);
void write_delta_frame(tW_gmx_trxframe *fr, FILE *fpbox, FILE *fpx, FILE *fpv, FILE *fpf);

/* Read/Write trj files */

// in gromacs, this was GROMACS_MAGIC. I have absolutely no idea what the point of
// it is, but as it's the first thing written in EVERY trr frame, we have to use it

#define BOCS_MAGIC 1993

enum {eioINT, eioFLOAT, eioDOUBLE, eioREAL, eioRVEC, eioNRVEC, eioSTRING};

static bool tW_do_binwrite(FILE *fp, const void *item, int nitem, int eio);
static bool tW_do_write_trjheader(tW_gmx_trxframe *fr);
static bool tW_do_write_trjstuff(tW_gmx_trxframe *fr);
bool tW_write_trj_frame(tW_gmx_trxframe *fr);
void wrap_tW_write_trj_frame(tW_gmx_trxframe *fr);

static bool tW_do_binread(FILE *fp, void *item, int nitem, int eio);
static bool tW_do_read_trjheader(tW_gmx_trxframe *fr);
static bool tW_do_read_trjstuff(tW_gmx_trxframe *fr);
bool tW_read_first_trj_frame(tW_gmx_trxframe *fr, const char *trx_fnm);
bool tW_read_next_trj_frame(tW_gmx_trxframe *fr);


/* and of course, after i get trj stuff working, I notice that
 * trj files are not supported in gromacs 5.1.4 */

/* Read/Write trr files */

static bool tW_do_write_xdr(tW_gmx_trxframe *fr, const void *item, int nitem, int eio);
static bool tW_do_write_trrheader(tW_gmx_trxframe *fr);
static bool tW_do_write_trrstuff(tW_gmx_trxframe *fr);
bool tW_write_trr_frame(tW_gmx_trxframe *fr);
void wrap_tW_write_trr_frame(tW_gmx_trxframe *fr);

static bool tW_do_read_xdr(tW_gmx_trxframe *fr, void *item, int nitem, int eio);
static bool get_trr_precision(tW_gmx_trxframe *fr);
static bool tW_do_read_trrheader(tW_gmx_trxframe *fr);
static bool tW_do_read_trrstuff(tW_gmx_trxframe *fr);
bool tW_read_first_trr_frame(tW_gmx_trxframe *fr, const char *trx_fnm);
bool tW_read_next_trr_frame(tW_gmx_trxframe *fr);

/* Read LAMMPS dump files */

int get_word_count(tW_line inp_line);
tW_word * get_words(tW_line inp_line);

int get_word_count_delim(tW_line inp_line, const char * delims);
void get_words_delim(tW_line inp_line, const char * delims, tW_word * word_list);

/* Unit conversion (assuming units real in LAMMPS)
Quantity	LAMMPS			GROMACS
Time		fs			ps			0.001
Distance	Angstrom		nm			0.1
Energy		kcal/mol		kJ/mol 			4.184
Velocity	Angstrom/fs		nm/ps			100.0
Force		kcal/(mol Ang)		kJ/(mol nm)		41.84
*/
#define T_LMP2GRO 0.001
#define X_LMP2GRO 0.1
#define V_LMP2GRO 100.0
#define F_LMP2GRO 41.84

int read_first_lammps_frame(tW_gmx_trxframe *fr, const char *trx_fnm);
bool read_next_lammps_frame(tW_gmx_trxframe *fr);

void write_lammps_data(tW_gmx_trxframe *fr, tW_gmx_topology *top);
void write_lammps_frame(tW_gmx_trxframe *fr, tW_gmx_topology *top);

/********************************************************************************
 these next functions are basically wrappers that call the appropriate function
 to do the requested operation based on the file type stored in fr->eFileType
 or top->eFileType
********************************************************************************/

int read_first_frame(tW_gmx_trxframe *fr, const char *fnm);
bool read_next_frame(tW_gmx_trxframe *fr, bool printCount );
void open_write_trajectory(tW_gmx_trxframe *fr, char *fnm);
void write_frame(tW_gmx_trxframe *fr, tW_gmx_topology *top);

bool read_topology(tW_gmx_topology *top, const char *fnm);

int get_type(const char *a_type, tW_gmx_topology *top);
void do_PBC(tW_gmx_trxframe *fr);

#endif
