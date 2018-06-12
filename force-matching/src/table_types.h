/**
@file table_types.h 
@author Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
*/

#ifndef table_types
#define table_types

#define NB_INDEX       2
#define BOND_INDEX     3
#define ANGLE_INDEX    4
#define DIHED_INDEX    5
#define INTRA_INDEX    6

#define DELTA_INDEX   7
#define LINEAR_INDEX  8
#define BSPLINE_INDEX 9

#define NB_NAME    "nb"
#define BOND_NAME  "bond"
#define ANGLE_NAME "angle"
#define DIHED_NAME "dihedral"
#define INTRA_NAME "intra"

#define DELTA_BASIS   "delta"
#define LINEAR_BASIS  "linear"
#define BSPLINE_BASIS "Bspline"

#define DELTA_SM     3
#define LINEAR_SM    3
#define BSPLINE_SM   5

#define DELTA_INTERP_FACT   2
#define LINEAR_INTERP_FACT  2
#define BSPLINE_INTERP_FACT 1

#define BOND_SLOPE -200000
#define ANGLE_SLOPE -0.5
#define NB_SLOPE -10000
#define INTRA_SLOPE -10000

#define SWITCH_TOL 10

#define MAX_PTS 10000

typedef struct{
  double *r;
  double *f;
  double *u;
} Arrays;

typedef struct{
  int n_pts;
  int smn_pts;
  int type_index;
  int type_basis;
  double r_min;
  double r_max;
  double sloper;
  double slopel;
  int interp_fact;
  int sm;
} Input;

# endif
