#ifndef LDD_H
#define LDD_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdlib.h>
#include <stdio.h>

#include "cgff_types.h"
#include "wnoid_math.h"
#include "safe_mem.h"
#include "gromacs_topology.h"
#include "solv_lin_eqns.h"

/* Forward declarations for stuff */
typedef struct tW_gmx_topology tW_gmx_topology;
typedef struct tW_gmx_trxframe tW_gmx_trxframe;


/* ~ Indicator Function Stuff ~ */

/* enum for every different type of indicator function w(r) */
#define WTYPE_SHELL_NAME "shell"
#define WTYPE_SPHERE_NAME "sphere"
#define WTYPE_SMOOTH_NAME "smooth"
#define WTYPE_LUCY_NAME "lucy"
#define WTYPE_DPD_NAME "dpd"

enum{eSHELL,eSPHERE,eSMOOTH,eLUCY,eDPD};

/* Number of coefficients for each indicator and an enum for referencing each */
#define N_SHELL_COEFFS 4
enum{eshC6,eshC4,eshC2,eshC0};
#define N_SPHERE_COEFFS 4
enum{espC3, espC1, espC0, espCn1};
#define N_SMOOTH_COEFFS 6
enum{esmC5, esmC4, esmC3, esmC2, esmC1, esmC0};
#define N_LUCY_COEFFS 4
enum{eluC4, eluC3, eluC2, eluC0};
#define N_DPD_COEFFS 3
enum{edpC2, edpC1, edpC0};

typedef struct tW_indicator_function
{
  int type;             // takes the value eSHELL, eSPHERE, eSMOOTH, eLUCY, or eDPD
  double norm;          // the spatial integral of w(r)
  double r0;            // r0 (must be 0 for eLUCY and eDPD)
  double rC;            // rc
  /* R and r are only used for eSPHERE */
  double R;             // (rc + r0) / 2
  double r;             // (rc - r0) / 2

  int n_coeffs;         // takes the value N_SHELL_COEFFS, N_SPHERE_COEFFS, etc...
  double *coeffs;       // gets array of length n_coeffs.

  double (*get_w) (struct tW_indicator_function*, double);
  double (*get_wp1) (struct tW_indicator_function*, double);
  double (*get_wp2) (struct tW_indicator_function*, double);

} tW_indicator_function;

/* Functions for dealing with the indicator function */

double get_w_shell(struct tW_indicator_function *ind, double r);
double get_w_1_shell(struct tW_indicator_function *ind, double r);
double get_w_2_shell(struct tW_indicator_function *ind, double r);

double get_w_sphere(struct tW_indicator_function *ind, double r);
double get_w_1_sphere(struct tW_indicator_function *ind, double r);
double get_w_2_sphere(struct tW_indicator_function *ind, double r);

double get_w_smooth(struct tW_indicator_function *ind, double r);
double get_w_1_smooth(struct tW_indicator_function *ind, double r);
double get_w_2_smooth(struct tW_indicator_function *ind, double r);

double get_w_lucy(struct tW_indicator_function *ind, double r);
double get_w_1_lucy(struct tW_indicator_function *ind, double r);
double get_w_2_lucy(struct tW_indicator_function *ind, double r);

double get_w_dpd(struct tW_indicator_function *ind, double r);
double get_w_1_dpd(struct tW_indicator_function *ind, double r);
double get_w_2_dpd(struct tW_indicator_function *ind, double r);

char *get_indicator_type(struct tW_indicator_function *ind);
int wtype_str_to_int(tW_word type);
void pop_indicator(struct tW_indicator_function *ind, tW_word wtype, double r0, double rc);
tW_indicator_function * init_indicator_function();

/* Functions for dealing with the local densities themselves */

void calculate_local_densities(tW_gmx_topology *top, tW_gmx_trxframe *fr, bool bSelf,
                                    double **hiLDs, double **loLDs);
void init_LD_info(tW_gmx_trxframe *fr, tW_gmx_topology *top, tW_word fnm);

/* Functions for LDPDF specifically */
typedef double (*tally_func_t)(tW_gmx_trxframe*, int atomidx, int typeidx); // MCL setting up pointer to func type for generic tally
double ld_tally_impl(tW_gmx_trxframe *fr, int atomidx, int typidx);
double sg_tally_impl(tW_gmx_trxframe *fr, int atomidx, int typidx);

#ifdef __cplusplus
}
#endif

#endif
