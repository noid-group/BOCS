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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "compute_pressure.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePressure::ComputePressure(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute pressure command");
  if (igroup) error->all(FLERR,"Compute pressure must use group all");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;
  pressflag = 1;
  timeflag = 1;
 
  p_match_flag = 0; // NJD

  // store temperature ID used by pressure computation
  // insure it is valid for temperature computation

  int n = strlen(arg[3]) + 1;
  id_temp = new char[n];
  strcpy(id_temp,arg[3]);

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Could not find compute pressure temperature ID");
  if (modify->compute[icompute]->tempflag == 0)
    error->all(FLERR,
               "Compute pressure temperature ID does not compute temperature");

  // process optional args

  if (narg == 4) {
    keflag = 1;
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = fixflag = 1;
  } else {
    keflag = 0;
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"ke") == 0) keflag = 1;
      else if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedralflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) improperflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg],"fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg],"virial") == 0) {
        pairflag = 1;
        bondflag = angleflag = dihedralflag = improperflag = 1;
        kspaceflag = fixflag = 1;
      } else error->all(FLERR,"Illegal compute pressure command");
      iarg++;
    }
  }

  vector = new double[6];
  nvirial = 0;
  vptr = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePressure::~ComputePressure()
{
  delete [] id_temp;
  delete [] vector;
  delete [] vptr;
}

/* ---------------------------------------------------------------------- */

void ComputePressure::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;

  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Could not find compute pressure temperature ID");
  temperature = modify->compute[icompute];

  // detect contributions to virial
  // vptr points to all virial[6] contributions

  delete [] vptr;
  nvirial = 0;
  vptr = NULL;

  if (pairflag && force->pair) nvirial++;
  if (bondflag && atom->molecular && force->bond) nvirial++;
  if (angleflag && atom->molecular && force->angle) nvirial++;
  if (dihedralflag && atom->molecular && force->dihedral) nvirial++;
  if (improperflag && atom->molecular && force->improper) nvirial++;
  if (fixflag)
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->virial_flag) nvirial++;

  if (nvirial) {
    vptr = new double*[nvirial];
    nvirial = 0;
    if (pairflag && force->pair) vptr[nvirial++] = force->pair->virial;
    if (bondflag && force->bond) vptr[nvirial++] = force->bond->virial;
    if (angleflag && force->angle) vptr[nvirial++] = force->angle->virial;
    if (dihedralflag && force->dihedral)
      vptr[nvirial++] = force->dihedral->virial;
    if (improperflag && force->improper)
      vptr[nvirial++] = force->improper->virial;
    if (fixflag)
      for (int i = 0; i < modify->nfix; i++)
        if (modify->fix[i]->virial_flag)
          vptr[nvirial++] = modify->fix[i]->virial;
  }

  // flag Kspace contribution separately, since not summed across procs

  if (kspaceflag && force->kspace) kspace_virial = force->kspace->virial;
  else kspace_virial = NULL;
}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputePressure::compute_scalar()
{
  if (p_match_flag) {scalar = compute_cg_scalar(); } // NJD
  else{
	  invoked_scalar = update->ntimestep;
	  if (update->vflag_global != invoked_scalar)
	    error->all(FLERR,"Virial was not tallied on needed timestep");


	  // invoke temperature it it hasn't been already
	  double t;
	  if (keflag) {
	    if (temperature->invoked_scalar != update->ntimestep)
	      t = temperature->compute_scalar();
	    else t = temperature->scalar;
	  }

	  if (dimension == 3) {
	    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
	    virial_compute(3,3);
	    if (keflag)
	      scalar = (temperature->dof * boltz * t +
		        virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
	    else
	      scalar = (virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
	  } else {
	    inv_volume = 1.0 / (domain->xprd * domain->yprd);
	    virial_compute(2,2);
	    if (keflag)
	      scalar = (temperature->dof * boltz * t +
		        virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
	    else
	      scalar = (virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
	  }
  }

  return scalar;
}

/* ----------------------------------------------------------------------
   compute required CG pressure correction for analytic basis //NJD 8/21/13
------------------------------------------------------------------------- */
 
double ComputePressure::get_cg_p_corr(int N_basis, double *phi_coeff, int N_mol, double vavg, double vCG)
{
	double correction = 0.0;

	for (int i=1; i<=N_basis; i++)
	{
		correction -= phi_coeff[i-1] * (N_mol * i / vavg) * pow( (1 / vavg) * (vCG - vavg), i-1);
	}

	return correction;
}

/* ----------------------------------------------------------------------
   finds the index of a point value within the given grid //NJD 12/18/13
------------------------------------------------------------------------- */

double ComputePressure::find_index(std::vector<double>& grid, double value)
{
	int i;
	double spacing;

	spacing = fabs(grid[1] - grid[0]);

	for(i = 0; i<(grid.size()-1); i++)
	{
		if (value>=grid[i] && value<=grid[i+1]) {return i;}
	}

	if (value>=grid[i] && value<=(grid[i]+spacing)) {return i;}

	fprintf(stderr, "ERROR: Value %f does not fall within the spline grid\n", value);
	for (int i = 0; i<grid.size(); i++)
	{
		fprintf(stderr, "grid %d: %f\n", i, grid[i]);
	}
	exit(1);
}


/* ----------------------------------------------------------------------
   compute required CG pressure correction for linear or cubic basis //NJD 12/18/13
------------------------------------------------------------------------- */

double ComputePressure::get_cg_p_corr(std::vector< std::vector<double> >& grid, int basis_type, double vCG )
{
	double correction, deltax;
	int i;

	i = find_index(grid[0], vCG);
	deltax = vCG - grid[0][i];

	if ( basis_type == 1 )
	{
		correction = grid[1][i] + ( grid[1][i+1] - grid[1][i] )/( grid[0][i+1] - grid[0][i] ) * ( deltax );
	} else if ( basis_type == 2 )
	{
		correction = grid[1][i] + (grid[2][i] * deltax) + (grid[3][i] * pow(deltax, 2)) + (grid[4][i] * pow(deltax, 3));
	} else
	{
		fprintf(stderr, "ERROR: bad spline type passed to get_cg_p_corr()\n");
		exit(1);
	}

	return correction;
}



/* ----------------------------------------------------------------------
   collect info for CG barostat from fix_NH //NJD 8/21/13
------------------------------------------------------------------------- */
void ComputePressure::send_cg_info(int basis_type, int sent_N_basis, double *sent_phi_coeff, int sent_N_mol, double sent_vavg)
{
//	fprintf(stderr, "cg pressure info sent\n");

	if(basis_type == 0)
	{
		p_basis_type = 0;
	}  else 
	{
		fprintf(stderr, "ERROR: Incorrect basis type passed to ComputePressure\n");
		exit(1);
 	}

	p_match_flag = 1;
	
	N_basis = sent_N_basis;
	phi_coeff = ((double *) calloc(N_basis, sizeof(double)) );
	for (int i=0; i<N_basis; i++)
	{
		phi_coeff[i] = sent_phi_coeff[i];
	}

	N_mol = sent_N_mol;
	vavg = sent_vavg;
}

/* ----------------------------------------------------------------------
   collect info for CG barostat from fix_NH for tabulated basis //NJD 12/17/13
------------------------------------------------------------------------- */
void ComputePressure::send_cg_info(int basis_type, std::vector< std::vector<double> > in_splines)
{
//	fprintf(stderr, "cg pressure info sent\n");

	if (basis_type == 1)
	{
		p_basis_type = 1;		
	}  else if(basis_type == 2)
	{
		p_basis_type = 2;		
	}   else
	{
		fprintf(stderr, "ERROR: Incorrect basis type passed to ComputePressure\n");
		exit(1);
 	}

	splines = in_splines;

	p_match_flag = 1;
	
}


/* ----------------------------------------------------------------------
   compute total corrected CG pressure, averaged over Pxx, Pyy, Pzz //NJD 8/21/13
------------------------------------------------------------------------- */
double ComputePressure::compute_cg_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  // invoke temperature it it hasn't been already

  double t;
  double volume, correction;
  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
      t = temperature->compute_scalar();
    else t = temperature->scalar;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    volume = (domain->xprd * domain->yprd * domain->zprd);

    if ( p_basis_type == 0 )   {correction = get_cg_p_corr(N_basis, phi_coeff, N_mol, vavg, volume);} // NJD 12/18/13
    else { correction = get_cg_p_corr(splines, p_basis_type, volume); } // NJD 12/18/13
    virial_compute(3,3);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
                virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p + (correction);
    else
      scalar = (virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p + (correction);
  } else {
	fprintf(stderr, "ERROR: Pressure matching not implemented in 2-d yet.\n");
	exit(1);
  }

  return scalar;
}


/* ----------------------------------------------------------------------
   compute pressure tensor
   assume KE tensor has already been computed
------------------------------------------------------------------------- */

void ComputePressure::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  double *ke_tensor;

  if (keflag) {
    if (temperature->invoked_vector != update->ntimestep)
      temperature->compute_vector();
    ke_tensor = temperature->vector;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6,3);
    if (keflag) {
      for (int i = 0; i < 6; i++)
	vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
    } else
      for (int i = 0; i < 6; i++)
	vector[i] = virial[i] * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4,2);
    if (keflag) {
      vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
      vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
      vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    } else {
      vector[0] = virial[0] * inv_volume * nktv2p;
      vector[1] = virial[1] * inv_volume * nktv2p;
      vector[3] = virial[3] * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    }
  }
}



/* ---------------------------------------------------------------------- */

void ComputePressure::virial_compute(int n, int ndiag)
{
  int i,j;
  double v[6],*vcomponent;

  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from forces and fixes

  for (j = 0; j < nvirial; j++) {
    vcomponent = vptr[j];
    for (i = 0; i < n; i++) v[i] += vcomponent[i];
  }

  // sum virial across procs

  MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,world);

  // KSpace virial contribution is already summed across procs

  if (kspace_virial)
    for (i = 0; i < n; i++) virial[i] += kspace_virial[i];

  // LJ long-range tail correction

  if (force->pair && force->pair->tail_flag)
    for (i = 0; i < ndiag; i++) virial[i] += force->pair->ptail * inv_volume;
}

/* ---------------------------------------------------------------------- */

void ComputePressure::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_temp;
  int n = strlen(id_new) + 1;
  id_temp = new char[n];
  strcpy(id_temp,id_new);
}
