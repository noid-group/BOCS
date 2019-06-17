/**
@file cgff_types.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
@brief Defines properties and structures of the cgff and cgmap calculations
*/

#ifndef NOID_TYPES
#define NOID_TYPES

#include <stdio.h>
#include <stdint.h>

#ifndef bool
#define bool int
#endif

//#ifndef DIM
#define DIM 3
typedef double dvec[DIM];
typedef double matrix[DIM][DIM];
//#endif

#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif

#define KILO 		 (1e3)			/* Thousand	*/
//#if GMX != 51
#define BOLTZMANN	 (1.380658e-23)		/* (J/K)	*/
#define AVOGADRO	 (6.0221367e23)		/* ()		*/
//#endif
#define RGAS             (BOLTZMANN*AVOGADRO)   /* (J/(mol K))  */
#define BOLTZ            (RGAS/KILO)    

/* Modes */
#define GMX_MODE "GROMACS"
#define PDB_MODE "PDB"

/* Basis function names. */
#define DELTA_NAME              "delta"
#define HARMONIC_NAME           "harmonic"
#define RYCKAERT_BELLEMANS_NAME "RB"
#define TOY_DIHED_NAME          "TOY"
#define POWER_NAME              "power"
/* START JFR */
#define LINEAR_NAME             "linear"
/*END JFR */
#define BSPLINE_NAME            "Bspline"

/* Basis function indices. */
#define DELTA_BASIS_INDEX              0
#define HARMONIC_BASIS_INDEX           1
#define RYCKAERT_BELLEMANS_BASIS_INDEX 2
#define TOY_DIHED_INDEX                3
#define POWER_INDEX                    4
/* START JFR */
#define LINEAR_BASIS_INDEX             5
/*END JFR */
#define BSPLINE_BASIS_INDEX            6	/* JFR - 07.22.12: Bspline basis set */

/* Non-bonded interaction names. */
#define NB_PAIR         "Pair_Interaction"

/* Bond interaction names. */
#define B_BOND_STRETCH  "BondStretch"
#define B_ANGLE         "Angle"
#define B_DIHEDRAL      "Dihedral"
#define B_NB_PAIR_BOND  "IntraMolec_NB_Pair"

/* Reference Forces */
#define REF_FORCE        "RefForces"

/* Interpolation names. */
#define REF_NO_INTERPOL_N     "none"
#define REF_LINEAR_INTERPOL_N "linear"

/* Indices to how to calc. ref forces. */
#define REF_NO_INTERPOL_I     0
#define REF_LINEAR_INTERPOL_I 1

/*JFR - how to calculate b_struct when using linear basis or the Bspline basis*/
#define LINEAR_STRUCT 1		// 0 -> same as delta basis, 1 -> take derivative of basis function analytically
#define BSPLINE_STRUCT 1	// 0 -> same as delta basis, 1 -> take derivative of basis function analytically

#define NOT_SET -1

#define MAX_NUM_BOND_COEFF 100
#define Max_Num_Sites_Prot 10000	/* 1000000 <- This used to be the max, we currently don't get anywhere close to this */
#define MAX_POWER_COEFF    5
#define MAX_BSPLINE_COEFF 20	/* JFR - 07.22.12: This is kspline */

#define MAXWORDLEN 1000
#define MAXLINELEN 5000

#define COMMENTFLAG "*!;"

#define FLOAT_EPS 1.0e-8		//NJD-9.24.14 - originally 5e-8

//#define N_PT 10 //JFR - 01/20/11 - calculates the force field up to N_PT orders by matrix perturbation theory 

#define DIHED_PREC   FLOAT_EPS	/* How close to zero to allow 1/sin(phi). */

//#define N_EIGEN 18 //JFR - 04/28/11 - Number of eigenvectors to print out from each side (pos and neg) of the spectrum 

/* CalcMODEs */
#define FULL        "FULL"
#define FIRST_HALF  "FIRST_HALF"
#define SECOND_HALF "SECOND_HALF"
#define TEST_INP    "TEST_INP"
#define IFULL        0
#define IFIRST_HALF  1
#define ISECOND_HALF 2
#define ITEST_INP    3

/* SOLN_METHs */
#define LU       "LU"
#define UU       "UU"
#define CHOLESKY "CHOLESKY"
#define SVD_SOLV  "SVD" // default SVD solution behavior

/* PC STUFF */
#define PC_DEFAULT "RPC"

/* Frame_Weighting selection */
#define NPT "NPT"
#define NONE "NONE"

typedef char tW_word[MAXWORDLEN];
typedef char tW_line[MAXLINELEN];

/* MRD 02.14.2019 bimask */
typedef uint64_t bitMask;

#define GET_N_SPOTS(x) ((int)((x)/64)+1)
#define SETSPOT1(bm, x) (bm[(int)((x)/64)] |= (1ULL << ((x) % 64)))
#define GETSPOT1(bm, x) ((bm[(int)((x)/64)] & (1ULL << ((x) % 64))) > 0)
#define CLEARSPOT1(bm, x) (bm[(int)((x)/64)] &= ~(1ULL << ((x) % 64)))
#define SETSPOT2(bm, n, x, y) (bm[(int)(((y)*(n)+(x))/64)] |= (1ULL << (((y)*(n)+(x)) % 64)))
#define GETSPOT2(bm, n, x, y) ((bm[(int)(((y)*(n)+(x))/64)] & (1ULL << (((y)*(n)+(x)) % 64))) > 0)
#define CLEARSPOT2(bm, n, x, y) (bm[(int)(((y)*(n)+(x))/64)] &= ~(1ULL << (((y)*(n)+(x)) % 64)))


typedef struct {
    char *cg_name;		/* name of site  */
    char *mol_name;		/* name of molecule */
    int res_no;			/* no.  of molecule */
    int n_atms;			/* no. of atms mapped to cg site */
    char **atm_name;		/* name of atoms */
    int *i_atm;			/* array of atm indices for atoms mapped to cg site */
    tW_word map_type;	        /* type of mapping.  usr, com, cog are supported */
    double *c_Ii;		/* mapping coefficients for site */
} tW_site_map;
/* The explicit mapping operator for a single site */

typedef struct {
    int total_atoms;
    char *molname;
    char **site_names;
    int n_sites;
    int *n_atoms;
    int **indices;
    char ***atom_names;
    char **atom_sites; // MRD 1.12.2018
    double *site_masses;
    double **weights;
    char **map_type;

} mapping_type;
/* Defines the mapping operator for a single molecule type */

typedef struct {
    int n_types;
    mapping_type *molecule_types;
    int *n_molecules;

} mapping_op;
/* Defines the mapping operator for a system of fixed composition */


typedef struct {
    tW_word inter_name;
    tW_word basis;
    tW_word name1;		// name of site1
    tW_word name2;		// name of site2
    int N_inter;		// no. of instances
    int N_pts;			// no. of grid points in table
    int N_coeff;		// no. of grid points in table
    int i_basis;		// type of basis functions
    int i_0;			// starting index for inter in sys->M,g2
    double dr;			// grid spacing
    double R_0;			// min grid pt
    double R_max;		// max grid pt
    int n_smooth;		// number of points in running average (use odd number since includes central point)
    int N_powers;
    int kspline;		/* JFR - 07.22.12: for Bspline basis, k is the order of the spline */
    int *powers;
    double *ptr_x;		// ptr to x values for truncated arrays
    double *ptr_u2;		// ptr to pair potential
    double *ptr_u2_forces;
    double *ptr_u2_struct;
    double *ptr_g2;		// ptr to rdf
    double *ptr_g2_cnt;		//
    double *ptr_L;		// ptr to laplacian
    double *ptr_b;		// ptr to deriv of pmf
    double *ptr_b_ref;		// ptr to deriv of pmf
    double *ptr_b_forces;
    double *ptr_b_struct;
} tW_type_inter2;

typedef struct {
    tW_word name;		// name of int. type (e.g. BondStretch)
    tW_word inter_name;		/* Name used to specify a specific Inter_Type[]. */
    tW_word basis;
    int N_Int_Sites;		// no. sites per interaction
    tW_word *Site_Types;	// names of types of interacting sites
    int i_basis;		// harmonic (i_basis=0); delta(bond) (i_basis=1);
    // delta(angle) (i_basis=2); harmonic angle (i_basis=3);
    // delta(dihedral) (i_basis=4); RB(dihedral) (i_basis=5)
    int N_coeff;		// no. of parameters describing force function, for delta assumes same as grid points
    int i_0;			// starting index in sys->M,g2,b,phi
    double dr;			// grid spacing
    double R_0;			// min grid pt
    double R_max;		// max grid pt
    int N_pts;			// no. of grid points in table for delta basis (bond)
    int n_smooth;
    int N_powers;
    int kspline;		/* JFR - 07.23.12: for Bspline basis, k is the order of the spline */
    int n_bonds_intramolec_pair_inter;
    int *powers;
    double *ptr_x;		// ptr to x - only used after trim_Mb
    double *ptr_g;		// ptr to pair distribution
    double *ptr_g_cnt;		// ptr to pair distribution
    double *ptr_L;		// ptr to laplacian
    double *ptr_b;		// ptr to force correlation fcn
    double *ptr_b_ref;		// ptr to force correlation fcn
    double *ptr_b_forces;
    double *ptr_b_struct;
    double *ptr_phi;		// ptr to force parameters
    double *ptr_phi_forces;
    double *ptr_phi_struct;
    /* below particular for each topology/system */
    int N_instances;		/* Number of times bond interactions appears. */
    int **Inter_List;		/* List of site indices for each of the N_instances occurances. */
} tW_Bonded_Inter;


typedef struct {
    tW_word inter_name;
    tW_word inter_type;
    tW_word basis;
    int i_basis;
    double dr;
    double R_min;
    double R_max;
    int N_pts;
    int N_coeff;
    int n_smooth;
    int i_0;
    int N_powers;
    int kspline;		/* JFR - 07.22.12: for Bspline basis, k is the order of the spline */
    int *powers;
    double *ptr_x;
    double *ptr_g;
    double *ptr_g_cnt;
    double *ptr_L;
    double *ptr_b;
    double *ptr_b_ref;
    double *ptr_b_forces;
    double *ptr_b_struct;
    double *ptr_phi;
    double *ptr_phi_forces;
    double *ptr_phi_struct;
} tW_Inter_Types;

typedef struct {		/* JFR - 08.10.11 */
    int flag_PT;		/* True => do the PT calculation */
    int N_PT;			/* Number of perturbation steps */
    int dPT;			/* spacing of output */
    int flag_MMOTF_SEP;		/* True => calc the decomposition of b at each perturbation step */
    int flag_eigen;		/* True => calc the eigenspectrum at each perturbation step */
    int N_eigen;		/* Number of eigenvectors to print out from each side of the spectrum (pos and neg) */
} tW_PT;

typedef struct {		/* JFR - 08.11.11 */
    int flag_Eigen;		/* True => calc the eigenspectrum of the G matrix */
    int N_Eigen;		/* Number of eigenvectors to print out from each side of the spectrum (pos and neg) */
    int flag_printn;		/* True => print fn, bn, Mn for various sets of eigenvectors */
    int DM;			/* spacing for grouping eigenvectors from the top */
    int DL;			/* spacing for grouping eigenvectors from the bottom */
    int flag_Gbar;		/* True => calc the eigenspectrum of the Gbar matrix (i.e., take out two body terms) */
    int flag_norm;		/* True => normalize the nb parts of the matrix by r^2 before calculating the eigenspectrum */
} tW_Eigen;

typedef struct {		/* JFR - 04.06.12 */
    int flag_SVD;		/* True => use the option below */
    double rcond;		/* Throw away eigenvalues under this number.  -1=> machine precision, FLOAT_EPS=> 1e-6. */
    int flag_printSV;		/* Explicitly calculate and print out singular values */
    int flag_printevecs;	/* If flag_printSV == TRUE, then TRUE => Print out evecs also */
    int flag_solve;		/* If flag_solve == TRUE, then use SVD algorithm to solve the equations */
} tW_SVD;

typedef struct {		/* JFR - 04.06.12 */
    int flag_TPR;		/* True => replace default tpr files */
    int N_TPR;			/* Number of tpr files, this should be the same as the number of trr files in inp.txt */
    tW_word *TPR_files;		/* name of each tpr file */
} tW_TPR;

typedef struct {		/* JFR - 06.27.12 */
    int flag_TPR_excl;		/* True => use alternative TPR file for exclusions */
    tW_word TPR_excl;		/* name of alternative TPR file */
} tW_TPR_EXCL;

typedef struct {		/* JFR - 04.06.12 */
    int flag_MT;		/* True => use the MT options */
    int flag_print;		/* True => print the metric tensor related files */
    int flag_norm;		/* True => normalize the metric tensor by r^2 for nonbonded interactions */
    int flag_Mcnt;		/* True => calculate the unweighted version of the matrix, Mcnt */
} tW_MT;

typedef struct {		/* JFR - 04.06.12 */
    int flag_MFD;		/* True => calculate the decomposition of b from each interaction in the system */
} tW_MFD;

typedef struct {		/* JFR - 04.06.12 */
    int flag_CalcMODE;		/* True => specify the calculation mode */
    int CalcMODE;		/* FULL(0) => do entire calculation, FIRST_HALF(1) => do the traj loop, calc correlation functions, print a save state */
    /*  SECOND_HALF(2) => read in save state, do the matrix inversion */
} tW_CalcMODE;

typedef struct {		/* JFR - 04.06.12 */
    int flag_PC;		/* True => specify the preconditioning for LU decomposition (and SVD if applicable) */
    tW_word RPC;		/* no, dimless, colnorm, MTvar, gcnt */
    int flag_normb;		/* True => rescale so that the norm of b is 1 */
    tW_word LPC;		/* no, dimless, rowmax, bvar, gcnt */
    int flag_normphi;		/* True => rescale so that the norm of phi is 1 */
} tW_PC;

typedef struct {		/* JFR - 04.06.12 */
    int flag_MEM;		/* True => Use the user memory options */
    int flag_LOWMEM;		/* True => run the calculation in low memory mode */
    tW_word info;		/* info = "forces" or "structures", which will be the only one used if flag_LOWMEM == TRUE */
    int flag_mult_top;		/* True => solve the problem for each topolgy separately */
} tW_MEM;

typedef struct {		/* JFR - 04.16.12 */
    int flag_SOLN;		/* True => specify the solution methods for the matrix inversion */
    tW_word SOLN_METH;		/* LU, UU, Cholesky, SVD */
} tW_SOLN;

typedef struct {		/* NJD - 03.10.15 */
    int flag_FRAMEWEIGHT;		/* True => Specify the ensemble of the input trr file */
    tW_word FRAMEWEIGHT;		/* NONE, NPT */
} tW_FRAMEWEIGHT;

typedef struct {		/* JFR - 04.16.12 */
    int flag_ERR;		/* True => calculate errors associated with the matrix inversion */
    char FACT;
} tW_ERR;

typedef struct {		/* JFR - 07.16.12 */
    int flag_REF;		/* True => use ref options specified in variable below */
    int flag_calcbref;		/* True => skip the triple loop, just calculate and print out bref */
    int flag_readbref;		/* True => use bref supplied in the file bref.dat (Note that [Reference_Potential] need not be specified) */
    int flag_splitfiles;	/* True => The MPI processes will be split between files instead of frames (np should divide n_files) */
    int flag_reftrr;		/* True => read in a trr file with the reference forces */
    int N_fnm;
    tW_word *reftrr_fnm;	/* reference force filename */
} tW_REF;

typedef struct {		/* JFR - 01.29.13 */
    int flag_TRIM;		/* True => trim the metric tensor below FE instead of the defult FLOAT_EPS */
    double FE;
} tW_TRIM;

typedef struct {		/* JFR - 01.29.13 */
    int flag_CHISQD;		/* True => read in a force file to use with the Chi2 calculation */
    tW_word force_fnm;		/* force filename */
} tW_CHISQD;

typedef struct {		/* JFR - 12.03.13 */
    int flag_REG;		/* True => Use Bayesian Inference to regularize the calculation */
    tW_word type;		/* type of REG (BAYES or UNCERT for now) */
    int Nmax;			/* Maximum number of iterations to optimize the regularization parameters (BAYES) */
    double tau_alpha;		/* Tolerence for convergence of the regularization matrix parameters (BAYES) */
    double tau_beta;		/* Tolerence for convergence of the precision parameter (BAYES or UNCERT) */
    int Nframes;		/* total number of data frames (BAYES) */
} tW_REG;

typedef struct {		/* JFR - 01.31.13 */
    int flag_RESCALE;		/* True => Use the force rescaling right preconditioning scheme */
    int Nmax;			/* Maximum number of iterations to optimize scaling */
    double tau_phi;		/* Tolerence for convergence of the scaling parameter */
} tW_RESCALE;

typedef struct {		/* JFR - 01.31.13 */
    int flag_CONSTRAIN;		/* True => Constrain the dihedral force coefficients to sum to zero */
    int Nmax;			/* Maximum number of iterations to optimize contraint scaling parameter, lambda */
    double tau_dih;		/* Tolerence for convergence of the constraint scaling parameter, lambda */
    double lambda;		/* initial constraint scaling parameter */
    double dlambda;		/*  amount to increase scaling parameter by after each step */
} tW_CONSTRAIN;

typedef struct {		/* JFR - 12.03.13 */
    int flag_ITER;		/* True => Apply options for iter-gYBG procedure */
    int flag_AAM2;		/* Use the AA 2-body contributions for the bond and angle interactions */
    tW_word AAM2_fnm;		/* AAM2 filename */
    int flag_bsolnerr;		/* TRUE => calculate the difference in the projected forces from the current and previous iterations */
    /*         using the current metric tensor                                                           */
} tW_ITER;

typedef struct {
    int N_Site_Types;		// no. of site types
    tW_word *Site_Types;	// list of names for site types

    tW_word **Inter_Map;	// list of lists of site names involved in pair interactions
    int **Inter_iMap;		// list of lists of indices for the above named sites
    int *Inter_Map_Len;		// number of pair interaction partners for each site type

    int N_Inter_Types;		/* Total number of interaction types. */
    tW_Inter_Types *Inter_Types;	/* Array of interaction types. */
    int N_Inter2_Types;		// no. of pair interaction types
    tW_type_inter2 *Inter2_Type_List;	// list of types of pair inter.
    int N_Bond_Int_Types;	// no. of bond interaction types
    tW_Bonded_Inter *Bonded_Inter_Types;	// array of bonded interaction types
    int N_coeff;		// no. of data_pts in g2, phi
    int N_pack;			// JFR - added 04.11.12: This is the number of elements in the matrix in packed form
    double Temperature;
    int nrexcl;			// Number of exclusion for FF; not needed for gmx loop
    double *x;			// grid of order parameters
    double *b;			// deriv. of pmf
    double *b_ref;
    double *b_forces;
    double *b_struct;
    double *g;			// rdf
    double *g_cnt;
    double *L;			// laplacian
    double *phi;		// ff coefficients
    double *phi_forces;
    double *phi_struct;

    /* MRD 02.14.2019 Eliminating the triple loop */
    bool SKIP_TRIPLE_LOOP;
    dvec *linear_half_matrix;
    dvec **half_matrix; 
    bitMask *bm_half_mat; 
    bool M_M2_proc;

    // JFR - added 04.11.12: put the matrix in packed form
    double *M;			// correlation kernel
    double *M2;			// 2-body contribution to M
    double *M_cnt;		// unweighted distributions
    //
    double *d2b;		// JFR - added 01.31.12: variance of b_forces for preconditioning
    double *d2M;		// JFR - added 01.31.12: variance of M for preconditioning
    double *rescale;		// JFR - added 01.31.12: rescaling factors for preconditioning (currently with gcnt)
    //
    /* NJD - add the extra G_wt, b_wt, etc. here */
    double *b_wt;		// accumulation of the frame-weighted b's
    double *M_wt;		// accumulation of the frame-weighted M's
    double *M2_wt;		// accumulation of the frame-weighted 2-body M's
    double wt_norm;             // the normalization for the frame-weighting
    //
    int flag_ref_potential;
    tW_PT PT_var;		// JFR - added 08.10.11 for solving the matrix equation using matrix perturbation theory
    tW_Eigen Eigen_var;		// JFR - added 08.11.11 for calculating eigenspectrum of the G matrix
    tW_SVD SVD_var;		// JFR - added 04.06.12 for solving the matric eqn with SVD and calculating condition numbers
    tW_TPR TPR_var;		// JFR - added 04.06.12 for over-writing default tpr filenames
    tW_TPR_EXCL TPR_EXCL_var;	// JFR - added 06.27.12 for using alternative TPR file for exclusions
    tW_MT MT_var;		// JFR - added 04.06.12 for options for printing and normalizing the metric tensor
    tW_MFD MFD_var;		// JFR - added 04.06.12 for calculating the mean force decomposition
    tW_CalcMODE CalcMODE_var;	// JFR - added 04.06.12 for splitting the FF calculation into two parts (for efficiency purposes) 
    tW_PC PC_var;		// JFR - added 04.06.12 for specifying a particular preconditioning for the matrix inversion
    tW_MEM MEM_var;		// JFR - added 04.13.12 for memory options
    tW_SOLN SOLN_var;		// JFR - added 04.16.12 for solution method options
    tW_FRAMEWEIGHT FRAMEWEIGHT_var;	// NJD - added 03.10.15 for ensemble selection options
    tW_ERR ERR_var;		// JFR - added 04.16.12 for error options
    double Chi2;		// JFR - added 06.27.12 calculate and print out Chi2
    tW_REF REF_var;		// JFR - added 07.16.12 options for precomputing and reading in bref
    tW_TRIM TRIM_var;		// JFR - added 01.29.13 options for trimming out basis vectors
    tW_CHISQD CHISQD_var;	// JFR - added 01.29.13 options for Chi2 calculation
    tW_REG REG_var;		// JFR - added 12.03.13 options for Regularization
    tW_RESCALE RESCALE_var;	// JFR - added 01.31.13 option for force rescale preconditioning
    tW_CONSTRAIN CONSTRAIN_var;	// JFR - added 01.31.13 option for constraining dihedrals
    tW_ITER ITER_var;		// JFR - added 12.03.13 options for the iter-gYBG procedure
} tW_system;

typedef struct {
    tW_word name;		/* Name from par.txt and atomtype in *.tpr. */
    int i_type;			/* Index to find site in tW_system.Site_Types[]. */
    int i_res;			/* Residue number read from *.tpr file. */
    int nr_excl;		/* Number of exlusions. */
    int *excl_list;		/* List of exlusions. */
    int nr_bonds;		/* Number of bond interactions. */
    int nr_bond_coeffs;		/* Number of coeffs. to determine for bond interactions. *//* NEEDED? */
    int *bond_type;		/* List of indices to tW_system.Bonded_Inter_Types[]. */
    int **bond_site;		/* List of pointers to tW_system.Bonded_Inter_Types[tW_CG_site.bond_type].Inter_List[]. */
    dvec r;			// coordinates for site
    dvec f;			// force on site
    dvec ref_f;			// forces from ref. pot.
} tW_CG_site;


typedef struct {
    FILE *fp_par;		/* ptr to parameter file                    */
    FILE *fp_log;		/* ptr to global log file                   */
    int N_struct_file;		/* no. of files with lists of structures    */
    tW_word *struct_file;	/* names of files listing structure files   */
    int N_struct;		/* no. of structure files                   */
    tW_word *structures;	/* names of structure files                 */
    double *p_struct;		/* p_gamma for topology gamma               */
    tW_word mode;		/* designates input and loop.               */
} tW_files;


typedef struct {
    int atom_position;
    int *dihedral_quartet;
    dvec grad_B[4];
    double cos_angle, sin_angle, angle;
} tW_dihedral;


typedef struct {
    tW_word name;
    int N_inter_sites;
    tW_word *site_types;
    int N_pts;
    int bonds_IntraMolecPairs;
    double *x;
    double *f;
    double dx;
    double x_0;
    double x_max;
} tW_ref_inter;


typedef struct {
    int N_pts;
    double *x;
    double *f;
    double dx;
    double x_0;
    double x_max;
} tW_ref_force;


typedef struct {
    int N_files;
    tW_word *fname;
    tW_ref_force *forces;
    int N_nb_pair_inter;
    tW_ref_inter *nb_pair_inter;
    int N_bondstretch_inter;
    tW_ref_inter *bondstretch_inter;
    int N_angle_inter;
    tW_ref_inter *angle_inter;
    int N_dihedral_inter;
    tW_ref_inter *dihedral_inter;
    int N_IntraMolecPairs_inter;
    tW_ref_inter *IntraMolecPairs_inter;
    int interpolation_index;
} tW_ref_potential;

#endif
