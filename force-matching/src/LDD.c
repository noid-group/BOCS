
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "mkl_lapacke.h"

#include "omp.h"

#include "cgff_types.h"
#include "safe_mem.h"
#include "wnoid_math.h"
#include "io_read.h"
#include "LDD.h"
#include "gromacs_topology.h"

/* Note: Do not call any pop_coeffs_XXXX function from anywhere but
 * the switch block of pop_indicator.
 *
 * Additionally, do not directly call get_w_XXXX, get_w_1_XXXX, or
 * get_w_2_XXXX. Instead, use ind->get_w(ind, r), ind->get_wp1(ind, r), 
 * and ind->get_wp2(ind,r). That way, you're sure to call the right
 * function for the type of indicator you're using. This eliminates
 * the need to have switch blocks every time you need to call one
 * of these functions.
 * */

/* Indicator functions for Shell's w(r) */

/*  Shell's w(r):
 *  Taken from Sanyal, Shell. JCP 2016, 145.
 *  w(r) = 1  for r < r0
 *  w(r) = 0  for r > rC
 *  w(r) = c0 + c2 * r^2 + c4 * r^4 + c6 * r^6   for r0 <= r <= rC
 *  define t = r0 * r0 / rC / rC
 *  define d = (1-t)^3
 *  c0 = ind->coeffs[eshC0] = (1 - 3t) / d
 *  c2 = ind->coeffs[eshC2] = (6t)/(d rC^2)
 *  c4 = ind->coeffs[eshC4] = -(3(1+t))/(d rC^4)
 *  c6 = ind->coeffs[eshC6] = 2/(d rC^6)
 */

double get_w_shell(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 1.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return (ind->coeffs[eshC6] * pow(r,6) + ind->coeffs[eshC4] * pow(r,4) +
                 ind->coeffs[eshC2] * pow(r,2) + ind->coeffs[eshC0]); }
}
double get_w_1_shell(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( 6.0 * ind->coeffs[eshC6] * pow(r,5) +
                  4.0 * ind->coeffs[eshC4] * pow(r,3) +
                  2.0 * ind->coeffs[eshC2] * r); }
}
double get_w_2_shell(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( 30.0 * ind->coeffs[eshC6] * pow(r,4) +
                  12.0 * ind->coeffs[eshC4] * pow(r,2) +
                   2.0 * ind->coeffs[eshC2]); }
}

void pop_coeffs_shell(struct tW_indicator_function *ind)
{
  ind->coeffs = (double *) ecalloc(N_SHELL_COEFFS,sizeof(double));
  double r0 = ind->r0;
  double rc = ind->rC;
  ind->n_coeffs = N_SHELL_COEFFS;

  double denom = pow((1.0 - r0*r0/rc/rc),3);
  ind->coeffs[eshC0] = (1.0 - 3.0 * r0 * r0 / rc / rc) / denom;
  ind->coeffs[eshC2] = (6.0 * r0 * r0 / pow(rc,4)) / denom;
  ind->coeffs[eshC4] = (-3.0 * (1.0 + r0 * r0 / rc / rc) / pow(rc,4)) / denom;
  ind->coeffs[eshC6] = (2.0 / pow(rc,6)) / denom;
  ind->norm = 4.0 * 3.141592653589 * (pow(r0,3)/3.0 +
                  ind->coeffs[eshC6] * (pow(rc,9) - pow(r0,9)) / 9.0 +
                  ind->coeffs[eshC4] * (pow(rc,7) - pow(r0,7)) / 7.0 +
                  ind->coeffs[eshC2] * (pow(rc,5) - pow(r0,5)) / 5.0 +
                  ind->coeffs[eshC0] * (pow(rc,3) - pow(r0,3)) / 3.0);
}


/* Indicator functions for Sphere w(r) */

/* Sphere w(r):
 * It's like Shell's in that
 * w(r) = 1 for r < r0
 * w(r) = 0 for r > rC
 * It differs in how it gets there
 * w(r) = c3 * r^3 + c1 * r + c0 + cn1 / r    for r0 <= r <= rC
 * where
 * R = (rC + r0)/2
 * x = (rC - r0)/2
 * c3 = 1 / (16 * x^3)
 * c1 = (-3.0 * (R^2 + x^2)) / (8.0 * r^3)
 * c0 = (R^3 + x^3) / (2.0 * x^3)
 * cn1 = (-3.0 * (R^4+x^4))/(16 * r^3) + 3.0 * R^2 / (8.0 * r)
 * 
 * This is the equation that gives you the fraction of a smaller sphere
 * of radius x that lies inside a larger sphere or radius R when their
 * centers are distance r apart.
 * In other words, if a CG site j is considered a solid sphere of radius x,
 * and is a distance r away from CG site i, then w(r) tells you what fraction
 * of site j is inside a sphere of size R centered on site i.
 *
 *
 */

double get_w_sphere(struct tW_indicator_function *ind, double r)
{
  if (r <= ind->r0) { return 1.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( ind->coeffs[espC3] * pow(r,3) + ind->coeffs[espC1] * r +
                  ind->coeffs[espC0] + ind->coeffs[espCn1] / r); }
}
double get_w_1_sphere(struct tW_indicator_function *ind, double r)
{
  if (r <= ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( 3.0 * ind->coeffs[espC3] * r * r +
                        ind->coeffs[espC1] +
                 -1.0 * ind->coeffs[espCn1] / r / r); }
}
double get_w_2_sphere(struct tW_indicator_function *ind, double r)
{
  if (r <= ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( 6.0 * ind->coeffs[espC3] * r +
                  2.0 * ind->coeffs[espCn1] / pow(r,3)); }
}

void pop_coeffs_sphere(struct tW_indicator_function *ind)
{
  ind->coeffs = (double *) ecalloc(N_SPHERE_COEFFS,sizeof(double));
  ind->R = (ind->r0 + ind->rC)/2.0;
  ind->r = (ind->rC - ind->r0)/2.0;
  ind->n_coeffs = N_SPHERE_COEFFS;

  double r = ind->r;
  double R = ind->R;
  double rc = ind->rC;
  double r0 = ind->r0;

  ind->coeffs[espC3] = 1.0 / (16.0 * pow(r,3));
  ind->coeffs[espC1] = (-3.0 * (R * R + r * r)) / (8.0 * pow(r,3));
  ind->coeffs[espC0] = (pow(R,3) + pow(r,3)) / (2.0 * pow(r,3));
  ind->coeffs[espCn1] = (-3.0 * (pow(R,4) + pow(r,4)))/(16.0 * pow(r,3)) + (3.0 * pow(R,2)) / (8.0 * r);
  ind->norm = 4.0 * 3.141592653589 * (pow(r0,3)/3.0 +
                    ind->coeffs[espC3] * (pow(rc,6) - pow(r0,6)) / 6.0 +
                    ind->coeffs[espC1] * (pow(rc,4) - pow(r0,4)) / 4.0 +
                    ind->coeffs[espC0] * (pow(rc,3) - pow(r0,3)) / 3.0 +
                    ind->coeffs[espCn1] * (pow(rc,2) - pow(r0,2)) / 2.0);
}


/* Indicator functions for Smooth w(r) */
/*
 * Smooth w(r) is a 5th order polynomial parameterized such that
 * the second derivative is continuouos at r0 and rC
 *
 */

double get_w_smooth(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 1.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( ind->coeffs[esmC5] * pow(r,5) + ind->coeffs[esmC4] * pow(r,4) +
                  ind->coeffs[esmC3] * pow(r,3) + ind->coeffs[esmC2] * pow(r,2) +
                  ind->coeffs[esmC1] * r + ind->coeffs[esmC0]); }
}
double get_w_1_smooth(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( 5.0 * ind->coeffs[esmC5] * pow(r,4) +
                  4.0 * ind->coeffs[esmC4] * pow(r,3) +
                  3.0 * ind->coeffs[esmC3] * pow(r,2) +
                  2.0 * ind->coeffs[esmC2] * r +
                        ind->coeffs[esmC1]); }
}
double get_w_2_smooth(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( 20.0 * ind->coeffs[esmC5] * pow(r,3) +
                  12.0 * ind->coeffs[esmC4] * pow(r,2) +
                   6.0 * ind->coeffs[esmC3] * r +
                   2.0 * ind->coeffs[esmC2]); }
}

void pop_coeffs_smooth(struct tW_indicator_function *ind)
{
  ind->coeffs = (double *) ecalloc(N_SMOOTH_COEFFS,sizeof(double));
  ind->n_coeffs = N_SMOOTH_COEFFS;

  double rc = ind->rC;
  double r0 = ind->r0;
  double denom = (pow(r0,5) - pow(rc,5)) / 120.0 -
                 (r0 * rc * (pow(r0,3) - pow(rc,3))) / 24.0 +
                 (r0 * r0 * rc * rc * (r0 - rc)) / 12.0;

  ind->coeffs[esmC5] = 1.0 / 20.0 / denom;
  ind->coeffs[esmC4] = -(r0 + rc) / 8.0 / denom;
  ind->coeffs[esmC3] = (r0 * r0 + 4.0 * r0 * rc + rc * rc) / 12.0 / denom;
  ind->coeffs[esmC2] = (-r0 * rc * (r0 + rc)) / 4.0 / denom;
  ind->coeffs[esmC1] = (r0 * r0 * rc * rc) / 4.0 / denom;
  ind->coeffs[esmC0] = ((-pow(rc,5)) / 120.0 + r0 * pow(rc,4) / 24.0 - r0 * r0 * pow(rc,3) / 12.0) / denom;
  ind->norm = 4.0 * 3.141592653589 * (pow(r0,3) / 3.0 +
                    ind->coeffs[esmC5] * (pow(rc,8) - pow(r0,8)) / 8.0 +
                    ind->coeffs[esmC4] * (pow(rc,7) - pow(r0,7)) / 7.0 +
                    ind->coeffs[esmC3] * (pow(rc,6) - pow(r0,6)) / 6.0 +
                    ind->coeffs[esmC2] * (pow(rc,5) - pow(r0,5)) / 5.0 +
                    ind->coeffs[esmC1] * (pow(rc,4) - pow(r0,4)) / 4.0 +
                    ind->coeffs[esmC0] * (pow(rc,3) - pow(r0,3)) / 3.0);
}


/* Indicator functions for Lucy w(r) */

double get_w_lucy(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 1.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( ind->coeffs[eluC4] * pow(r,4) + ind->coeffs[eluC3] * pow(r,3) +
                  ind->coeffs[eluC2] * pow(r,2) + ind->coeffs[eluC0]); }
}
double get_w_1_lucy(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( 4.0 * ind->coeffs[eluC4] * pow(r,3) +
                  3.0 * ind->coeffs[eluC3] * pow(r,2) +
                  2.0 * ind->coeffs[eluC2] * r); }
}
double get_w_2_lucy(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return ( 12.0 * ind->coeffs[eluC4] * pow(r,2) +
                   6.0 * ind->coeffs[eluC3] * r +
                   2.0 * ind->coeffs[eluC2]); }
}

void pop_coeffs_lucy(struct tW_indicator_function *ind)
{
  ind->coeffs = (double *) ecalloc(N_LUCY_COEFFS,sizeof(double));
  ind->n_coeffs = N_LUCY_COEFFS;

  ind->r0 = 0.0;
  double rc = ind->rC;

  ind->coeffs[eluC4] = -3.0 / pow(rc,4);
  ind->coeffs[eluC3] = 8.0 / pow(rc,3);
  ind->coeffs[eluC2] = -6.0 / pow(rc,2);
  ind->coeffs[eluC0] = 1.0;
  ind->norm = 16 * 3.141592653589 / 105.0 * pow(rc,3);
}

/* DPD indicator functions */

double get_w_dpd(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 1.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return (ind->coeffs[edpC2] * r * r + ind->coeffs[edpC1] * r + ind->coeffs[edpC0]); }
}

double get_w_1_dpd(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return (2.0 * ind->coeffs[edpC2] * r + ind->coeffs[edpC1]); }
}

double get_w_2_dpd(struct tW_indicator_function *ind, double r)
{
  if (r < ind->r0) { return 0.0; }
  else if (r > ind->rC) { return 0.0; }
  else { return (2.0 * ind->coeffs[edpC2]); }
}

void pop_coeffs_dpd(struct tW_indicator_function *ind)
{
  ind->coeffs = (double *) ecalloc(N_DPD_COEFFS,sizeof(double));
  ind->n_coeffs = N_DPD_COEFFS;
  
  ind->r0 = 0.0;
  double rc = ind->rC;
  
  ind->coeffs[edpC2] = 1.0 / (rc * rc);
  ind->coeffs[edpC1] = -2.0 / rc;
  ind->coeffs[edpC0] = 1.0;
  ind->norm = (2.0 * 3.141592653589 * pow(rc,3))/(15.0);
}


/* Function to return the name of indicator function based on ind->type */

char *get_indicator_type(struct tW_indicator_function *ind)
{
  switch (ind->type)
  {
    case (eSHELL):
      return WTYPE_SHELL_NAME;
      break;
    case (eSPHERE):
      return WTYPE_SPHERE_NAME;
      break;
    case (eSMOOTH):
      return WTYPE_SMOOTH_NAME;
      break;
    case (eLUCY):
      return WTYPE_LUCY_NAME;
      break;
    case (eDPD):
      return WTYPE_DPD_NAME;
      break;
    default:
      /* ERROR */
      return "ERROR";
      break;
  }
}

/* Get indicator type idx from the word */
int wtype_str_to_int(tW_word type)
{
  if (strcmp(type,WTYPE_SHELL_NAME) == 0) { return eSHELL; }
  if (strcmp(type,WTYPE_SPHERE_NAME) == 0) { return eSPHERE; }
  if (strcmp(type,WTYPE_SMOOTH_NAME) == 0) { return eSMOOTH; }
  if (strcmp(type,WTYPE_LUCY_NAME) == 0) { return eLUCY; }
  if (strcmp(type,WTYPE_DPD_NAME) == 0) { return eDPD; }

  fprintf(stderr,"ERROR: unknown wtype: %s \n",type);
  fprintf(stderr,"Recognized indicator function types: %s %s %s %s %s \n",
                  WTYPE_SHELL_NAME, WTYPE_SPHERE_NAME, WTYPE_SMOOTH_NAME,
                  WTYPE_LUCY_NAME, WTYPE_DPD_NAME);
  exit(EXIT_FAILURE);


}

/* Function to populate the contents of the indicator function struct */
void pop_indicator(struct tW_indicator_function *ind, tW_word strwtype, double r0, double rc)
{
  int wtype = wtype_str_to_int(strwtype); 
  ind->type = wtype;
  ind->r0 = r0;
  ind->rC = rc;

  switch (wtype)
  {
    case (eSHELL):
      pop_coeffs_shell(ind);
      ind->get_w = get_w_shell;
      ind->get_wp1 = get_w_1_shell;
      ind->get_wp2 = get_w_2_shell;
      break;
    case (eSPHERE):
      pop_coeffs_sphere(ind);
      ind->get_w = get_w_sphere;
      ind->get_wp1 = get_w_1_sphere;
      ind->get_wp2 = get_w_2_sphere;
      break;
    case (eSMOOTH):
      pop_coeffs_smooth(ind);
      ind->get_w = get_w_smooth;
      ind->get_wp1 = get_w_1_smooth;
      ind->get_wp2 = get_w_2_smooth;
      break;
    case (eLUCY):
      pop_coeffs_lucy(ind);
      ind->get_w = get_w_lucy;
      ind->get_wp1 = get_w_1_lucy;
      ind->get_wp2 = get_w_2_lucy;
      break;
    case (eDPD):
      pop_coeffs_dpd(ind);
      ind->get_w = get_w_dpd;
      ind->get_wp1 = get_w_1_dpd;
      ind->get_wp2 = get_w_2_dpd;
      break;
    default:
      fprintf(stderr,"ERROR: Unrecognized wtype: %d \n",wtype);
      fprintf(stderr,"Type_Index  Type_Name\n----------  ---------\n");
      fprintf(stderr,"%10d  %9s\n",eSHELL,WTYPE_SHELL_NAME);
      fprintf(stderr,"%10d  %9s\n",eSPHERE,WTYPE_SPHERE_NAME);
      fprintf(stderr,"%10d  %9s\n",eSMOOTH,WTYPE_SMOOTH_NAME);
      fprintf(stderr,"%10d  %9s\n",eLUCY,WTYPE_LUCY_NAME);
      fprintf(stderr,"%10d  %9s\n",eDPD,WTYPE_DPD_NAME);
      break;
  }
}

/* Allocates memory for an indicator function struct and returns the address in memory */
tW_indicator_function *init_indicator_function()
{
  tW_indicator_function *ind;
  ind = (tW_indicator_function *) ecalloc(1,sizeof(tW_indicator_function));
  return ind;
}


/*
 * get_dx() does dx = xi - xj assuming the box has all 90 degree angles.
 * */

void get_dx(dvec xi, dvec xj, tW_matrix box, dvec dx)
{
  int i; 
  for (i = 0; i < DIM; ++i)
  { 
    if (fabs(xi[i] - xj[i]) <= box[i][i] / 2.0)
    {
      dx[i] = xi[i] - xj[i];
    }
    else
    {
      dx[i] = box[i][i] - fabs(xi[i] - xj[i]);
      if (xi[i] > xj[i]) { dx[i] *= -1.0; }
    }
  }
}

/* This function populates the local density array (matrix?) in fr */

/*
  Since we have to write each frame out after calculating local densities,
  parellelizing this by sending individual frames to individual processors
  would be annoying, because you then have to make sure they wait and 
  write their frame to the file at the proper time.

  Instead, use OpenMP to parallelize the local density calculation in each
  frame across processors.
*/
void calculate_local_densities(tW_gmx_topology *top, tW_gmx_trxframe *fr, bool bSelf, double **hiLDs, double **loLDs)
{
  int i;

  #pragma omp parallel for
  for (i = 0; i < fr->contents->natoms; ++i)
  {
    int jtype;
    /* Step 1: Initialize all of the local densities to 0 */
    for (jtype = 0; jtype < top->n_atomtypes; ++jtype)
    {
      fr->local_densities[i][jtype] = 0.0;
      fr->local_density_gradients[i][jtype][0] = 0.0;
      fr->local_density_gradients[i][jtype][1] = 0.0;
      fr->local_density_gradients[i][jtype][2] = 0.0;
    }

    int j, k; 
    double r;
    dvec xi, xj, dx, u_ij;
  
    int itype = get_type(*(top->contents->atoms.atomtype[i]),top);

    copy_vector(fr->contents->x[i],xi);
     
    // MCL 08.20.24  Get exclusion list for each atom before we compare to all other atoms
     int * excl_listi = top->get_excl_list(top,i);
     int nr_excl_for_i = top->get_nexcl(top,i);

    /* Step 2: Calculate the local densities */
    for (j = 0; j < fr->contents->natoms; ++j)
    {
      if ((bSelf == FALSE) && (i == j)) { continue; }
      
      if (((i != j) && skip_excl( nr_excl_for_i, excl_listi, j))) { continue; }

      jtype = get_type(*(top->contents->atoms.atomtype[j]),top);
      if (top->bLocalDens[itype][jtype])
      {
        copy_vector(fr->contents->x[j],xj);
        get_dx(xi,xj,fr->contents->box,dx);
        r = calc_norm(dx);
        if (r < top->indicators[itype][jtype]->rC)
        {
          fr->local_densities[i][jtype] += top->indicators[itype][jtype]->get_w(top->indicators[itype][jtype],r);
          if (r > 0.0)
          {
            scal_times_vect(1.0/r,dx,u_ij);
            fr->local_density_gradients[i][jtype][0] += top->indicators[itype][jtype]->get_wp1(top->indicators[itype][jtype],r) * u_ij[0];
            fr->local_density_gradients[i][jtype][1] += top->indicators[itype][jtype]->get_wp1(top->indicators[itype][jtype],r) * u_ij[1];
            fr->local_density_gradients[i][jtype][2] += top->indicators[itype][jtype]->get_wp1(top->indicators[itype][jtype],r) * u_ij[2];
          }
        }
      }
    }
    /* Step 2.5: normalize the local densities */
    for (jtype = 0; jtype < top->n_atomtypes; ++jtype) 
    { 
      fr->local_densities[i][jtype] /= top->indicators[itype][jtype]->norm;
      fr->local_density_gradients[i][jtype][0] /= top->indicators[itype][jtype]->norm;
      fr->local_density_gradients[i][jtype][1] /= top->indicators[itype][jtype]->norm;
      fr->local_density_gradients[i][jtype][2] /= top->indicators[itype][jtype]->norm;
    }

    efree(excl_listi); // MCL Added to avoid mem leaks while checking exclusions
    excl_listi = NULL;
  } // end for loop over i sites

  /* Step 3: update hi and lo LDs */
  #pragma omp critical
  {
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      int itype = get_type(*(top->contents->atoms.atomtype[i]),top);
      int jtype;
      for (jtype = 0; jtype < top->n_atomtypes; ++jtype)
      {
        if (fr->local_densities[i][jtype] > hiLDs[itype][jtype])
        {
          hiLDs[itype][jtype] = fr->local_densities[i][jtype];
        }
        if (fr->local_densities[i][jtype] < loLDs[itype][jtype])
        {
          loLDs[itype][jtype] = fr->local_densities[i][jtype];
        }
      }
    }
  }
}

void init_LD_info(tW_gmx_trxframe *fr, tW_gmx_topology *top, tW_word fnm)
{
  int i, j, n_types, test_sscanf, inp_int;
  float inp_float1, inp_float2;
  tW_line inp_line;
  tW_word *names, inp_word;

  FILE *fp = open_file(fnm,'r');

  get_next_line(fp,inp_line);
  test_sscanf = sscanf(inp_line," %d ",&n_types);
  if (n_types != top->n_atomtypes)
  {
    fprintf(stderr,"ERROR: number of atom types in file %s: %d\n",fnm,n_types);
    fprintf(stderr,"\tnumber of atom types in topology: %d\n",top->n_atomtypes);
    exit(1);
  }

  /* allocate memory */

  names = (tW_word *) ecalloc(n_types,sizeof(tW_word));
  fr->linear_lds = (double *) ecalloc(fr->contents->natoms * n_types, sizeof(double));
  fr->local_densities = (double **) ecalloc(fr->contents->natoms,sizeof(double *));
  fr->linear_ldgs = (dvec *) ecalloc(fr->contents->natoms * n_types, sizeof(dvec));
  fr->local_density_gradients = (dvec **) ecalloc(fr->contents->natoms,sizeof(dvec *));
  for (i = 0; i < fr->contents->natoms; ++i) 
  { 
    fr->local_densities[i] = &(fr->linear_lds[i*n_types]); 
    fr->local_density_gradients[i] = &(fr->linear_ldgs[i*n_types]);
  }


  top->int_map = (int *) ecalloc(n_types,sizeof(int));
  top->indicators = (tW_indicator_function ***) ecalloc(n_types, sizeof(tW_indicator_function **));
  top->bLocalDens = (bool **) ecalloc(n_types,sizeof(bool *));
  for (i = 0; i < n_types; ++i)
  {
    top->indicators[i] = (tW_indicator_function **) ecalloc(n_types, sizeof(tW_indicator_function *));
    top->bLocalDens[i] = (bool *) ecalloc(n_types,sizeof(bool));
    for (j = 0; j < n_types; ++j)
    {
      top->indicators[i][j] = init_indicator_function();
      top->bLocalDens[i][j] = FALSE;
    }
  }

  /* get map */
  for (i = 0; i < n_types; ++i)
  {
    get_next_line(fp,inp_line);
    test_sscanf = sscanf(inp_line, " %d %s ",&inp_int,&inp_word);
    if (test_sscanf != 2)
    {
      fprintf(stderr,"ERROR: expected to find site type index and site type name \n");
      fprintf(stderr,"\thowever, I found: %s",inp_line);
      exit(EXIT_FAILURE);
    }
    strcpy(names[inp_int],inp_word);
  }

  for (i = 0; i < n_types; ++i)
  {
    for (j = 0; j < n_types; ++j)
    {
      if (strcmp(top->atom_type_names[i],names[j]) == 0)
      {
        top->int_map[j] = i;
        j = n_types;
      }
    }
  }

  free(names);

  int counter = 0;
  while ((get_next_line(fp,inp_line) != -1) && (counter < n_types * n_types))
  {
    tW_word wtype;
    ++counter;
    test_sscanf = sscanf(inp_line," %d %d %s %f %f ",&i, &j, &wtype, &inp_float1, &inp_float2);
    if (test_sscanf != 5)
    {
      fprintf(stderr,"ERROR: unable to read 5 parameters: itype jtype wtype r0 rc from file: %s\n",fnm);
      fprintf(stderr,"\tread %d parameters from line: %s",test_sscanf,inp_line);
      exit(1);
    }
    int itype = top->int_map[i];
    int jtype = top->int_map[j];
    pop_indicator(top->indicators[itype][jtype],wtype,(double)inp_float1,(double)inp_float2);
    top->bLocalDens[top->int_map[i]][top->int_map[j]] = TRUE;
  }

  fclose(fp);

  top->hi_LDs = (double **) ecalloc(n_types, sizeof(double *));
  top->lo_LDs = (double **) ecalloc(n_types, sizeof(double *));
  for (i = 0; i < n_types; ++i)
  {
    top->hi_LDs[i] = (double *) ecalloc(n_types, sizeof(double));
    top->lo_LDs[i] = (double *) ecalloc(n_types, sizeof(double));
    for (j = 0; j < n_types; ++j)
    {
      top->hi_LDs[i][j] = 0.0;
      top->lo_LDs[i][j] = 9e30;
    }
  }

  fr->n_atomtypes = n_types;
  fr->bLDs = TRUE;
  fr->bLD_Grads = TRUE;
}

/* MCL LDPDF Extra Functions */ // Writing prototypes I'll define funcs with macros in LDPDF implementation 

/**
 * ld_tally_impl:
 * returns the local density of particles of type j (typeidx) around atom i (atomidx) 
 * @param fr The trajectory frame read  (Must contain LD Info)
 * @param atomidx the index of the atom to read
 * @param typeidx the index of the surrounding type to read
 * @return SG, the square of the gradient vector around i due to j
 * @see make_LDPDF.c
 * @note This is an impl function, calls to it are via macro definitions
 */
double ld_tally_impl(tW_gmx_trxframe *fr, int atomidx, int typidx)
{
return fr->local_densities[atomidx][typidx]; 
}
/**
 * sg_tally_impl:
 * returns the SG from particles of type j (typeidx) around atom i (atomidx)
 * @param fr The trajectory frame read  (Must contain Gradient Info)
 * @param atomidx the index of the atom to read
 * @param typeidx the index of the surrounding type to read
 * @return SG, the square of the gradient vector around i due to j
 * @see make_LDPDF.c
 * @note This is an impl function, calls to it are via macro definitions
 */
double sg_tally_impl(tW_gmx_trxframe *fr, int atomidx, int typeidx)
{
 double SG2; 
 SG2 = (fr->local_density_gradients[atomidx][typeidx][0])*(fr->local_density_gradients[atomidx][typeidx][0]) +
      fr->local_density_gradients[atomidx][typeidx][1]*(fr->local_density_gradients[atomidx][typeidx][1]) + 
      fr->local_density_gradients[atomidx][typeidx][2]*(fr->local_density_gradients[atomidx][typeidx][2]);

 return SG2;
}
