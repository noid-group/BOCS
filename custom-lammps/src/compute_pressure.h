/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(pressure,ComputePressure)

#else

#ifndef LMP_COMPUTE_PRESSURE_H
#define LMP_COMPUTE_PRESSURE_H

#include "compute.h"
#include <vector>




namespace LAMMPS_NS {

class ComputePressure : public Compute {
 public:
  ComputePressure(class LAMMPS *, int, char **);
  virtual ~ComputePressure();
  void init();
  double compute_scalar();
  double compute_cg_scalar(); //NJD 10/16/13
//  void compute_cg_vector(); //NJD 10/17/13
  double get_cg_p_corr(int, double *, int, double, double); // NJD 8/21/13
  double get_cg_fluct(double, double); //NJD 10/16/13
  void send_cg_info(int, int, double*, int, double); //NJD 10/17/13
  void send_cg_info(int, std::vector< std::vector<double> >); //NJD 12/17/13
  double get_cg_p_corr(std::vector< std::vector<double> >&, int, double ); //NJD 12/18/13
  double find_index(std::vector<double>& , double ); //NJD 12/18/13


  void compute_vector();
  void reset_extra_compute_fix(const char *);

 protected:
  double boltz,nktv2p,inv_volume;
  int nvirial,dimension;
  double **vptr;
  double *kspace_virial;
  Compute *temperature;
  char *id_temp;
  double virial[6];
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int fixflag,kspaceflag;


  int p_basis_type;		//NJD 12/17/13 - type of basis used for F(V)
  int p_match_flag;                // NJD 8/19/13 - flag for using pressure matching or not
  double vavg;			   // NJD 8/21/13 - the avg volume of the ref. aa system
  int N_mol;			   // NJD 8/21/13 - the number of particles in the system
  int N_basis;                   // NJD 8/19/13 - the number of basis functions to use in p_matching
  double *phi_coeff;          // NJD 8/19/13 - an array of the coefficients phi
//  double variance;                 // NJD 10/16/13 - the variance of the kinetic fluctuations we are replacing
  std::vector<std::vector<double> > splines; // NJD 12/17/13 - lookup table for interpolation


  void virial_compute(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute pressure must use group all

Virial contributions computed by potentials (pair, bond, etc) are
computed on all atoms.

E: Could not find compute pressure temperature ID

The compute ID for calculating temperature does not exist.

E: Compute pressure temperature ID does not compute temperature

The compute ID assigned to a pressure computation must compute
temperature.

E: Virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
