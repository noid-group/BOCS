/*
@file gromacs_topology.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser, Maria Lesniewski
@brief Functions to interface with the copied Gromacs data structures defined in gromacs_topology.h
@	as well as read/write various file types
*/

//c library includes
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <tirpc/rpc/xdr.h>

//local includes
#include "cgff_types.h"
#include "safe_mem.h"
#include "wnoid_math.h"
#include "io_read.h"
#include "gromacs_topology.h"

// MCL generic error handler
#define THROW_TPR_ERROR(msg) do { \
	fprintf(stderr, "ERROR: Unexpected file format while translating dumped tpr\n"); \
	fprintf(stderr, "likely this is due to newer gromacs .tpr formats (post 2024.3)\n"); \
	fprintf(stderr, "please open a PR on github to alert us\n");\
	fprintf(stderr, "Be sure to copy/paste this error and let us know which version\n");\
	fprintf(stderr, "of gromacs you were using so that we can update our code\n");\
	fprintf(stderr, "%s::%d::%s()\nERROR:%s\n",__FILE__, __LINE__, __func__, msg);\
		exit(1);\
} while (0)


/*****************************************************************************************
copyright(): Prints the gromacs copyright and credits output
*****************************************************************************************/
void copyright()
{
  fprintf(stdout,"~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~\n");
  fprintf(stdout,"The \"BOCS\" software package was written by:\n");
  fprintf(stdout,"William Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, and Michael DeLyser\n");
  fprintf(stdout,"From the Noid lab at The Pennsylvania State University\n");
  fprintf(stdout,"Please cite our published work: \n\n");
  fprintf(stdout,"Dunn, NJH; Lebold, KM; DeLyser, MR; Rudzinski, JF; Noid, WG.\n");
  fprintf(stdout,"BOCS: Bottom-Up Open-Source Coarse-Graining Software. \n");
  fprintf(stdout,"J. Phys. Chem. B. 122, 13, 3363-3377.\n");
  fprintf(stdout,"\n\n");
  fprintf(stdout,"Our software uses multiple data structures and functions that have\n");
  fprintf(stdout,"been copied from Gromacs 4.5.3 or 5.1.4 and slightly edited.\n");
  fprintf(stdout,"The Gromacs development team was led by:\n");
  fprintf(stdout,"Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl\n");
  fprintf(stdout,"A full list of GROMACS contributors is available on their\n");
  fprintf(stdout,"website at http://www.gromacs.org/About_Gromacs/People \n");
  fprintf(stdout,"~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~\n");
    
}

/*****************************************************************************************
get_natoms(): Returns the number of atoms specified by topology 'self'
*****************************************************************************************/
int get_natoms(tW_gmx_topology *self)
{
  return self->contents->atoms.nr; 
}

/*****************************************************************************************
get_natoms(): Returns the name of topology 'self'
*****************************************************************************************/
char *get_name(tW_gmx_topology *self)
{
  return ((char *) self->contents->name); 
}


/*****************************************************************************************
get_nbonds(): Returns the number of bonds listed in self 
*****************************************************************************************/
int get_nbonds(tW_gmx_topology *self)
{
  return self->contents->idef.il[F_BONDS].nr; 
}

/*****************************************************************************************
get_bond_count(): Returns the number of bonds present in topology
*****************************************************************************************/
int get_bond_count(tW_gmx_topology *self)
{
  int i;
  int n_bonds = 0;
  for (i = 0; i < self->contents->mols.nr; ++i)
  {
    n_bonds += self->molecules[i].n_bonds * self->molecules[i].n_mols;
  }
  return n_bonds;
}


/*****************************************************************************************
get_n_bond_types(): Returns the number of different bond types
*****************************************************************************************/
int get_n_bond_types(tW_gmx_topology *self)
{
  int i;
  int n_bond_types = 0;
  for (i = 0; i < self->contents->idef.ntypes; ++i)
  {
    /* this should capture both BONDS and TABBONDS */
    if (strstr(self->force_names[i],"BONDS") != NULL) { ++n_bond_types; }
  }
  return n_bond_types;
}

/*****************************************************************************************
get_nangles(): Returns the number of angles listed in self 
*****************************************************************************************/
int get_nangles(tW_gmx_topology *self)
{
  return self->contents->idef.il[F_ANGLES].nr; 
}

/*****************************************************************************************
get_angle_count(): Returns the number of angles present in topology
*****************************************************************************************/
int get_angle_count(tW_gmx_topology *self)
{
  int i;
  int n_angles = 0;
  for (i = 0; i < self->contents->mols.nr; ++i)
  {
    n_angles += self->molecules[i].n_angles * self->molecules[i].n_mols;
  } 
  return n_angles;
}

/*****************************************************************************************
get_n_angle_types(): Returns the number of different angle interaction types
*****************************************************************************************/
int get_n_angle_types(tW_gmx_topology *self)
{
  int i;
  int n_angle_types = 0;
  for (i = 0; i < self->contents->idef.ntypes; ++i)
  {
    if (strstr(self->force_names[i],"ANGLES") != NULL) { ++n_angle_types; }
  }
  return n_angle_types;
}

/*****************************************************************************************
get_ndihs(): Returns the number of dihedrals listed in self 
*****************************************************************************************/
int get_ndihs(tW_gmx_topology *self)
{
 return self->contents->idef.il[F_PDIHS].nr; 
}


/*****************************************************************************************
get_dihedral_count(): Returns the number of dihedrals present in topology
*****************************************************************************************/
int get_dihedral_count(tW_gmx_topology *self)
{
  int i;
  int n_dihs = 0;
  for (i = 0; i < self->contents->mols.nr; ++i)
  {
    n_dihs += self->molecules[i].n_dihs * self->molecules[i].n_mols;
  } 
  return n_dihs;
}

/*****************************************************************************************
get_n_dihedral_types(): Returns the different number of dihedral interaction types
*****************************************************************************************/
int get_n_dihedral_types(tW_gmx_topology *self)
{
  int i;
  int n_dihedral_types = 0;
  for (i = 0; i < self->contents->idef.ntypes; ++i)
  {
    if (strstr(self->force_names[i],"DIHS") != NULL) { ++n_dihedral_types; }
  }
  return n_dihedral_types;
}

/*****************************************************************************************
get_npairs(): Returns the number of intramolecular nonbonded pairs listed in self 
*****************************************************************************************/
int get_npairs(tW_gmx_topology *self)
{
  return self->contents->idef.il[F_LJ14].nr; 
}

/*****************************************************************************************
get_bond_list() Returns a list of bonds in the t_ilist format:

+-----+--+--+-----+--+--+
|type1|i1|i2|type2|i1|i2|...
+-----+--+--+-----+--+--+
*****************************************************************************************/
int *get_bond_list(tW_gmx_topology *self)
{
  int i, nelements;
  int *bond_list;

  nelements = self->get_nbonds(self);

  bond_list = ecalloc(nelements, sizeof(int));

  for (i=0; i<nelements; i++)
  {
    bond_list[i] = self->contents->idef.il[F_BONDS].iatoms[i];
  }
   
  return bond_list; 
}

/*****************************************************************************************
get_angle_list() Returns a list of angles in the t_ilist format:

+-----+--+--+--+-----+--+--+--+
|type1|i1|i2|i3|type2|i1|i2|i3|...
+-----+--+--+--+-----+--+--+--+
*****************************************************************************************/
int *get_angle_list(tW_gmx_topology *self)
{
  int i, nelements;
  int *angle_list;

  nelements = self->get_nangles(self);

  angle_list = ecalloc(nelements, sizeof(int));

  for (i=0; i<nelements; i++)
  {
    angle_list[i] = self->contents->idef.il[F_ANGLES].iatoms[i];
  }
   
  return angle_list; 
}

/*****************************************************************************************
get_dih_list() Returns a list of dihedrals in the t_ilist format:

+-----+--+--+--+--+-----+--+--+--+--+
|type1|i1|i2|i3|i4|type2|i1|i2|i3|i4|...
+-----+--+--+--+--+-----+--+--+--+--+
*****************************************************************************************/
int *get_dih_list(tW_gmx_topology *self)
{
  int i, nelements;
  int *dih_list;
  nelements = self->get_ndihs(self);

  dih_list = ecalloc(nelements, sizeof(int));

  for (i=0; i<nelements; i++)
  {
    dih_list[i] = self->contents->idef.il[F_PDIHS].iatoms[i];
  }
  return dih_list; 
}

/*****************************************************************************************
get_pair_list() Returns a list of nonbonded intramolecular pairs in the t_ilist format:

+-----+--+--+-----+--+--+
|type1|i1|i2|type2|i1|i2|...
+-----+--+--+-----+--+--+
*****************************************************************************************/
int *get_pair_list(tW_gmx_topology *self)
{
  int i, nelements;
  int *pair_list;

  nelements = self->get_npairs(self);

  pair_list = ecalloc(nelements, sizeof(int));

  for (i=0; i<nelements; i++)
  {
    pair_list[i] = self->contents->idef.il[F_LJ14].iatoms[i];
  }
   
  return pair_list; 
}


/*****************************************************************************************
get_resind(): Returns the residue index of the i_atom'th atom in topology 'self'
*****************************************************************************************/
int get_resind(tW_gmx_topology *self, int i_atom)
{
  return self->contents->atoms.atom[i_atom].resind; 
}

/*****************************************************************************************
get_atomtype(): Returns the type of the i_atom'th atom in topology 'self'
*****************************************************************************************/
char *get_atomtype(tW_gmx_topology *self, int i_atom)
{
  return *(self->contents->atoms.atomtype[i_atom]); 
}

/*****************************************************************************************
get_nexcl(): Returns the number of exclusions for i_atom'th atom in topology 'self'
*****************************************************************************************/
int get_nexcl(tW_gmx_topology *self, int i_atom)
{
  int i_start, i_end, nr_excl;

  i_start = self->contents->excls.index[i_atom];
  i_end = self->contents->excls.index[i_atom + 1] - 1;
  nr_excl = i_end - i_start + 1;

  return nr_excl; 
}


/*****************************************************************************************
get_excl_list(): Returns the exclusion list for the i_atom'th atom in topolgy 'self'
*****************************************************************************************/
int *get_excl_list(tW_gmx_topology *self, int i_atom)
{
  int j;
  int i_start, i_end, nr_excl;
  int *excl_list = NULL;

  i_start = self->contents->excls.index[i_atom];
  i_end = self->contents->excls.index[i_atom + 1] - 1;
  nr_excl = i_end - i_start + 1;

  /* Generate exclusion list for cg site i_atom. */
  if (nr_excl > 0) 
  {
    /* Allocate memory for exclusion list for site i. */
    excl_list = (int *) ecalloc(nr_excl, sizeof(int));

    /* Copy exclusions (site indices) from top file for site i. */
    for (j = 0; j < nr_excl; j++) { excl_list[j] = self->contents->excls.a[i_start + j]; }

  }
  return excl_list; 
}


/*****************************************************************************************
init_tW_gmx_topology(): initialize the topology structure and align the function pointers
*****************************************************************************************/
tW_gmx_topology* init_tW_gmx_topology()
{
  tW_gmx_topology *top;

  top = emalloc(sizeof(tW_gmx_topology));
  top->contents = (tW_t_topology *) ecalloc(1,sizeof(tW_t_topology));
  top->eFileType = -1;
  top->b_tpr = FALSE;

  top->get_natoms = get_natoms; // Sets the pointer to the function get_natoms(self)
  top->get_name = get_name; // We use this pointer to function format so that our top structure is flexible to the input
  top->get_resind = get_resind; // e.g. we have access to the info in tW_t_top * contents through only the address at tW_gmx_top instead of having to pass input directly
  top->get_atomtype = get_atomtype; // so as long as we fill in tW_t_top the right way based on the input file type, we can use the tW_gmx_top generically for different kinds of topology input files
  top->get_nexcl = get_nexcl;
  top->get_excl_list = get_excl_list;

  top->get_nbonds = get_nbonds;
  top->get_bond_count = get_bond_count;
  top->get_n_bond_types = get_n_bond_types;
  top->get_nangles = get_nangles;
  top->get_angle_count = get_angle_count;
  top->get_n_angle_types = get_n_angle_types;
  top->get_ndihs = get_ndihs;
  top->get_dihedral_count = get_dihedral_count;
  top->get_n_dihedral_types = get_n_dihedral_types;
  top->get_npairs = get_npairs;

  top->get_bond_list = get_bond_list;
  top->get_angle_list = get_angle_list;
  top->get_dih_list = get_dih_list;
  top->get_pair_list = get_pair_list;

  // MCL 08.06.25 
  // Cleaning up pointers not yet set to avoid undefined behavior
  top->int_map = NULL;
  top->idx_file = NULL;
  top->bLocalDens = NULL;
  top->hi_LDs = NULL;
  top->lo_LDs = NULL;
  top->indicators = NULL;

  return top; 
}

/*****************************************************************************************
get_top(): Opens GROMACS topology (tpr) file and stores topology information in the 
data structure top.
*****************************************************************************************/
bool get_top(tW_word coord_fnm, tW_gmx_topology *top, bool TPR_flag, tW_word tpr_filename)
{
  bool bGromacs = FALSE;
  char *trr_ptr;		/* For time being we will only use trr files.  */
  char *tpr_fnm;		/* Name of tpr file determined from coord_fnm. */
  char *auto_tpr_fnm;
  matrix box;

  /* Get tpr file name. */
  if (TPR_flag == TRUE) 
  {
    if (file_exists(tpr_filename)) 
    {
      /* Get topology information. */
      //bGromacs = read_tpr_file(tpr_filename, top, box);
      //read_tpr(top, tpr_filename);
      top->b_tpr = read_topology(top, tpr_filename); // CLEANUP
    } 
    else 
    {
      fprintf(stderr, "ERROR: tpr file '%s' not found. Check par.txt.\n", tpr_filename);
      exit(EXIT_FAILURE);
    }

  } 
  else 
  {  /* do it the old way */
    auto_tpr_fnm = (char *) ecalloc(strlen(coord_fnm) + 1, sizeof(char));
    strcpy(auto_tpr_fnm, coord_fnm);
    trr_ptr = strstr(auto_tpr_fnm, ".trr\0");
    if (trr_ptr == NULL) 
    {
      auto_tpr_fnm = NULL;
      return bGromacs;
    }
    strcpy(trr_ptr, ".btp\0"); // MRD 9262017 changed tpr to btp "bocs top"

    if (file_exists(auto_tpr_fnm)) { top->b_tpr = read_topology(top,auto_tpr_fnm); }
    else 
    {
      fprintf(stderr, "ERROR: tpr file '%s' not found.  This file was "
                "generated automatically from the trr filename '%s'.\n",
                auto_tpr_fnm, coord_fnm);
      exit(EXIT_FAILURE);
    }

    free(auto_tpr_fnm);
  }

  return top->b_tpr; 
}

/*****************************************************************************************
print_tpr_file(): Prints information stored in top. Used for debugging.
*****************************************************************************************/
int print_tpr_file(FILE * fp_log, tW_gmx_topology *top)
{
  int i;
  int i_nr;
  int i_start, i_end, nr_excl;

  fprintf(fp_log, "\nIn print_tpr_file().\n");

  fprintf(fp_log, "  top->contents->name: %-10s\n\n", top->contents->name);
  fprintf(fp_log, "  No. of atomtypes: %-4d\n", top->contents->atomtypes.nr);
  fprintf(fp_log, "  No. of residues: %-4d\n", top->contents->atoms.nres);
  fprintf(fp_log, "  No. atoms: %-4d\n", top->contents->atoms.nr);
  fprintf(fp_log, "  Looping over atoms.\n");
  for (i = 0; i < top->contents->atoms.nr; i++) 
  {
    fprintf(fp_log, "    atm_no: %-4d   res_no: %-4d   ", i, top->contents->atoms.atom[i].resind);	// 4.5.3

    fprintf(fp_log, "atm_name: %-10s   atm_type: %-10s   type_id: %-4d  res_name: %-10s\n", 
                *(top->contents->atoms.atomname[i]), 
                *(top->contents->atoms.atomtype[i]), 
                top->contents->atoms.atom[i].type,
                *(top->contents->atoms.resinfo[top->contents->atoms.atom[i].resind].name));	// 4.5.3

  }
  fprintf(fp_log, "\n");

  fprintf(fp_log, "  Printing stuff from top->contents->atoms.excl.\n");
  fprintf(fp_log, "    No. of blocks; top->contents->excls.nr: %-4d \n",
        top->contents->excls.nr);
  fprintf(fp_log, "    No. of atoms; top->contents->excls.nra: %-4d \n",
        top->contents->excls.nra);
  fprintf(fp_log, "    Looping over top->contents->excls.nr (no. atoms+1 ).\n");
  for (i_nr = 0; i_nr <= top->contents->excls.nr; i_nr++) 
  {
    fprintf(fp_log, "      top->contents->excls.index[%-4d]: %-4d\n", i_nr,
        top->contents->excls.index[i_nr]);
  }
  fprintf(fp_log,
        "    Looping over top->contents->excls.nra (total no. of excl. indices).\n");
  for (i_nr = 0; i_nr < top->contents->excls.nra; i_nr++) 
  {
    fprintf(fp_log, "      top->contents->excls.a[%-4d]: %-4d\n", i_nr,
        top->contents->excls.a[i_nr]);
  }
  fprintf(fp_log, "    Printing excl. for each atom.\n");
  for (i_nr = 0; i_nr < top->contents->excls.nr; i_nr++) 
  {
    i_start = top->contents->excls.index[i_nr];
    i_end = top->contents->excls.index[i_nr + 1] - 1;
    nr_excl = i_end - i_start + 1;

    fprintf(fp_log, "      i_nr: %-4d   ", i_nr);
    fprintf(fp_log, "i_start: %-4d  i_end: %-4d  nr_excl: %-4d  ",
                     i_start, i_end, nr_excl);

    for (i = i_start; i < i_end + 1; i++) { fprintf(fp_log, "%-4d ", top->contents->excls.a[i]); }

    fprintf(fp_log, "\n");
  }
  fprintf(fp_log, "\n");

  fprintf(fp_log,
      "  Printing stuff about interactions from top->contents->idef. \n");
  fprintf(fp_log, "    top->contents->idef.ntypes: %-4d \n", top->contents->idef.ntypes);
  fprintf(fp_log, "    top->contents->idef.atnr  : %-4d \n\n", top->contents->idef.atnr);

  fprintf(fp_log, "  Printing F_BONDS interactions.\n");
  fprintf(fp_log, "    top->contents->idef.il[F_BONDS].nr: %-4d\n",
      top->contents->idef.il[F_BONDS].nr);
  for (i = 0; i < top->contents->idef.il[F_BONDS].nr; i += 3) 
  {
    fprintf(fp_log,
        "    top->contents->idef.il[F_BONDS].iatoms[%-4d]  inter_type: %-4d  atoms: %-4d %-4d\n",
        i, top->contents->idef.il[F_BONDS].iatoms[i],
        top->contents->idef.il[F_BONDS].iatoms[i + 1],
        top->contents->idef.il[F_BONDS].iatoms[i + 2]);
  }
  fprintf(fp_log, "\n");

  fprintf(fp_log, "  Printing F_ANGLES interactions.\n");
  fprintf(fp_log, "    top->contents->idef.il[F_ANGLES].nr: %-4d\n",
        top->contents->idef.il[F_ANGLES].nr);
  for (i = 0; i < top->contents->idef.il[F_ANGLES].nr; i += 4) 
  {
    fprintf(fp_log,
        "    top->contents->idef.il[F_ANGLES].iatoms[%-4d]  inter_type: %-4d  atoms: %-4d %-4d %-4d\n",
        i, top->contents->idef.il[F_ANGLES].iatoms[i],
        top->contents->idef.il[F_ANGLES].iatoms[i + 1],
        top->contents->idef.il[F_ANGLES].iatoms[i + 2],
        top->contents->idef.il[F_ANGLES].iatoms[i + 3]);
  }
  fprintf(fp_log, "\n");

  fprintf(fp_log, "  Printing F_PDIHS interactions.\n");
  fprintf(fp_log, "    top->contents->idef.il[F_PDIHS].nr: %-4d\n",
        top->contents->idef.il[F_PDIHS].nr);
  for (i = 0; i < top->contents->idef.il[F_PDIHS].nr; i += 5) 
  {
    fprintf(fp_log,
        "    top->contents->idef.il[F_PDIHS].iatoms[%-4d]  inter_type: %-4d  atoms: %-4d %-4d %-4d %-4d\n",
        i, top->contents->idef.il[F_PDIHS].iatoms[i],
        top->contents->idef.il[F_PDIHS].iatoms[i + 1],
        top->contents->idef.il[F_PDIHS].iatoms[i + 2],
        top->contents->idef.il[F_PDIHS].iatoms[i + 3],
        top->contents->idef.il[F_PDIHS].iatoms[i + 4]);
  }
  fprintf(fp_log, "\n");

  fprintf(fp_log, "  Printing F_LJ14 interactions.\n");
  fprintf(fp_log, "    top->contents->idef.il[F_LJ14].nr: %-4d\n",
        top->contents->idef.il[F_LJ14].nr);
  for (i = 0; i < top->contents->idef.il[F_LJ14].nr; i += 3) 
  {
    fprintf(fp_log,
        "    top->contents->idef.il[F_LJ14].iatoms[%-4d]  inter_type: %-4d  atoms: %-4d %-4d\n",
        i, top->contents->idef.il[F_LJ14].iatoms[i],
        top->contents->idef.il[F_LJ14].iatoms[i + 1],
        top->contents->idef.il[F_LJ14].iatoms[i + 2]);
  }
  fprintf(fp_log, "\n");

  fprintf(fp_log,
      "  Looping over interaction types: top->contents->idef.ntypes (%d)\n",
      top->contents->idef.ntypes);
  for (i = 0; i < top->contents->idef.ntypes; i++) 
  {
    fprintf(fp_log, "    i: %-4d  name: %-10s \n", i,
		//interaction_function[top->contents->idef.functype[i]].name);
    top->force_names[i]);
  }

  fprintf(fp_log, "End print_tpr_file().\n\n");

  return 0; 
}

/*****************************************************************************************
get_site_info(): Stores information from top into sys and CG_struct.
*****************************************************************************************/
void get_site_info(tW_CG_site *CG_struct, tW_system *sys, tW_gmx_topology *top)
{
  int i, j;
  int i_type;
  int i_start, i_end, nr_excl;
  int n_sites;

    /* JFR - 06.27.12: for bondref, add variables and read in new top file for exclusions */
  tW_gmx_topology *top_excl;
//    matrix box;
  top_excl = init_tW_gmx_topology();
  if (sys->TPR_EXCL_var.flag_TPR_excl == TRUE) 
  {
    //get_top(" ", &top_excl, box, sys->TPR_EXCL_var.flag_TPR_excl,sys->TPR_EXCL_var.TPR_excl);
    get_top(" ", top_excl, sys->TPR_EXCL_var.flag_TPR_excl, sys->TPR_EXCL_var.TPR_excl);
  }

  n_sites = top->get_natoms(top);

  for (i = 0; i < n_sites; i++) // Loop over every atom [instances of atoms, not types]
  {
    /* Get site name from GROMACS topology for site i. */
    strcpy(CG_struct[i].name, top->get_atomtype(top, i) );

    /* Does the site name from top match any site names read from par.txt? */
    i_type = match_word(sys->N_Site_Types, CG_struct[i].name, sys->Site_Types); // returns index of type that matched site_i
    if (i_type == -1) 
    {
      printf("\nERROR: Unknown site type found in GROMACS topology: %s.\n", top->get_name(top) );
      printf("i: %d  top_name: %s\n", i, CG_struct[i].name);
      exit(EXIT_FAILURE);
    }

    /* Stores index to find site type name in the list of site types: sys->Site_Types. */
    CG_struct[i].i_type = i_type;

    /* Determine residue number from top file. */
    //CG_struct[i].i_res = top.atoms.atom[i].resnr; // 4.0.7
    CG_struct[i].i_res = top->get_resind(top, i);	// 4.5.3

    /* Empty previous list of exculsions for site i. */
    if (CG_struct[i].excl_list != NULL) { efree(CG_struct[i].excl_list); }

    /* JFR - 06.27.12: for bondref, use top_excl for excl info */
    if (sys->TPR_EXCL_var.flag_TPR_excl == TRUE) 
    {
      /* Determine no. of exclusions for site i. */
      CG_struct[i].nr_excl = top_excl->get_nexcl(top_excl, i);
      CG_struct[i].excl_list = top_excl->get_excl_list(top_excl, i);
    } 
    else 
    {	/* do it the usual way */
      CG_struct[i].nr_excl = top->get_nexcl(top, i);
      CG_struct[i].excl_list = top->get_excl_list(top, i);
    }
  }  /* End loop over sites */

  /* MRD 06.18.20 */
  sys->n_atomtypes = top->n_atomtypes;
  sys->LD_type_names = (tW_word *) ecalloc(top->n_atomtypes, sizeof(tW_word));
  for (i = 0; i < n_sites; ++i)
  {
    CG_struct[i].LDs = (double *) ecalloc(top->n_atomtypes, sizeof(double));
    CG_struct[i].LD_type_idx = get_type(*(top->contents->atoms.atomtype[i]),top);
    strcpy(sys->LD_type_names[CG_struct[i].LD_type_idx],*(top->contents->atoms.atomtype[i])); // MRD 06.18.2020
    CG_struct[i].LDGrads = (dvec *) ecalloc(top->n_atomtypes, sizeof(dvec)); // MRD 06.18.2020
  }

}

/*****************************************************************************************
set_natoms(): sets the number of atoms in the frame self to natoms, and allocates memory
for those arrays that are flagged as present
*****************************************************************************************/
void set_natoms(tW_gmx_trxframe *self, int natoms)
{
  int i;
  self->contents->natoms = natoms;

  if (self->contents->bX) { self->contents->x = (tW_rvec *) ecalloc(natoms, sizeof(tW_rvec)); }
  if (self->contents->bV) { self->contents->v = (tW_rvec *) ecalloc(natoms, sizeof(tW_rvec)); }
  if (self->contents->bF) { self->contents->f = (tW_rvec *) ecalloc(natoms, sizeof(tW_rvec)); }

/* Local Density Stuff MRD 06.17.2020 */
  if (self->bLDs)
  {
    self->linear_lds = (double *) ecalloc(natoms * (self->n_atomtypes), sizeof(double));
    self->local_densities = (double **) ecalloc(natoms,sizeof(double *));
    for (i = 0; i < natoms; ++i)
    {
      self->local_densities[i] = &(self->linear_lds[i * (self->n_atomtypes)]);     
    }
    
    if (self->bLD_Grads)
    {
      self->linear_ldgs = (dvec *) ecalloc(natoms * (self->n_atomtypes), sizeof(dvec));
      self->local_density_gradients = (dvec **) ecalloc(natoms, sizeof(dvec *));   
      for (i = 0; i < natoms; ++i)
      {
        self->local_density_gradients[i] = &(self->linear_ldgs[i * (self->n_atomtypes)]);
      }
    } 
  }



}

/*****************************************************************************************
set_atom_labels(): sets the atom labels of the frame self if they are available in top
*****************************************************************************************/
void set_atom_labels(tW_gmx_trxframe *self, tW_gmx_topology *top)
{
  if (top->b_tpr) 
  {
    self->contents->bAtoms = TRUE;
    self->contents->atoms = &(top->contents->atoms);
  } 
  else 
  {
    self->contents->bAtoms = FALSE;
    self->contents->atoms = NULL;
  }

}


/*****************************************************************************************
print_frame_dump(): dumps trr frame details to the screen for reading
*****************************************************************************************/
void print_frame_dump(tW_gmx_trxframe *self)
{
  int i, j, i_atm, n_atms, res_no;

  tW_t_atom *atom_ptr;
  tW_t_atoms *atoms_ptr;

  printf("--- in print_trr_file ---\n");

  if (self->contents->bTitle) { printf("title: %s \n", self->contents->title); }

  n_atms = self->contents->natoms;
  printf("fr.natoms: %d \n", n_atms);

  printf("starting time.  fr.t0:   %f \n", self->contents->t0);
  printf("previous time.  fr.tpf:  %f \n", self->contents->tpf);
  printf("prevprev time.  fr.tppf:  %f \n", self->contents->tppf);

  printf("present  step.  fr.step: %d \n", self->contents->step);
  printf("present  time.  fr.time: %f \n", self->contents->time);

  printf("present  lambda.  fr.lambda: %f \n", self->contents->lambda);

  printf("printing atoms xvf.\n");
  for (i_atm = 0; i_atm < n_atms; i_atm++) 
  {
    printf("%i \t ", i_atm);

    if (self->contents->bAtoms) 
    {
      atoms_ptr = self->contents->atoms;
      atom_ptr = &(atoms_ptr->atom[i_atm]);
      res_no = atom_ptr->resind;
      printf("type_id: %d \t res_no: %d \t atm_name: %s \t atm_type: %s \t res_name: %s \n", 
        atom_ptr->type, res_no, *(atoms_ptr->atomname[i_atm]), 
        *(atoms_ptr->atomtype[i_atm]), *(atoms_ptr->resinfo[res_no].name));
    }

    if (self->contents->bX) 
    {
      printf("x: %f %f %f \t ", self->contents->x[i_atm][0], self->contents->x[i_atm][1],
       self->contents->x[i_atm][2]);
    }
    if (self->contents->bV) 
    {
      printf("v: %f %f %f \t ", self->contents->v[i_atm][0], self->contents->v[i_atm][1],
       self->contents->v[i_atm][2]);
    }
    if (self->contents->bF) 
    {
      printf("f: %f %f %f \t ", self->contents->f[i_atm][0], self->contents->f[i_atm][1],
       self->contents->f[i_atm][2]);
    }
    printf("\n");
  }

  if (self->contents->bBox) 
  {
    printf("printing box. \n");
    for (i = 0; i < DIM; i++) 
    {
      printf("i: %d \t", i);
      for (j = 0; j < DIM; j++) { printf(" %f ", self->contents->box[i][j]); }
      printf("\n");
    }
  }
}

/*****************************************************************************************
get_pos_of(): copies the dvec containing the position of the ith atom in the frame to x
*****************************************************************************************/
void get_pos_of(tW_gmx_trxframe *self, int i, dvec x)
{
  copy_vector(self->contents->x[i], x);
}


/*****************************************************************************************
get_vel_of(): copies the dvec containing the velocity of the ith atom in the frame to v
*****************************************************************************************/
void get_vel_of(tW_gmx_trxframe *self, int i, dvec v)
{
  copy_vector(self->contents->v[i], v);
}

/*****************************************************************************************
get_box(): copies the box matrix of this frame to box
*****************************************************************************************/
void get_box(tW_gmx_trxframe *self, matrix box)
{
  copy_matrix(self->contents->box, box);
}

/*****************************************************************************************
setup_xdr(): sets up the XDR data structure for reading/writing trr files MRD 92517
*****************************************************************************************/
void setup_xdr(tW_gmx_trxframe *self, const char *fnm, bool bRead)
{
  if (bRead) // MRD 02.11.2019 got rid of  if (!self->fp) { } around the fopen statements
  { // MCL 05.08.2026 - We still need to check if fp is valid before we try to operate on it.
    if ( fopen(fnm, "rb") )
    {
	    self->fp = fopen(fnm,"rb");
    	    self->xdrmode = XDR_DECODE; 
    }
    else { fprintf(stderr, "ERROR: Unable to open %s for mode 'rb'.\n", fnm);
   	   fprintf(stderr, "%s::%d::%s()\n",__FILE__, __LINE__, __func__); 
    	    exit("EXIT_FAILURE");}
  }
  else
  {
    self->fp = fopen(fnm,"wb");
    self->xdrmode = XDR_ENCODE;
  }

  self->xdr = (XDR *) ecalloc(1,sizeof(XDR));
  xdrstdio_create(self->xdr, self->fp, self->xdrmode); // This associates our xdr stream (self->xdr) with the file for I/O and the appt. translation mode (encoding/decoding)
						       //We can check if this is successful by checking the bool returned if we want
  if (bRead) // Grab the header and determine the precision here.
  {
    self->bDouble = get_trr_precision(self);    
    rewind(self->fp);
    if (self->bDouble) { fprintf(stderr,"trr file %s is double precision\n",fnm); }
    else { fprintf(stderr,"trr file %s is single precision\n",fnm); }
  }

}

/*****************************************************************************************
init_tW_gmx_trxframe(): initializes a tW_gmx_trxframe variable
*****************************************************************************************/
tW_gmx_trxframe* init_tW_gmx_trxframe()
{
  tW_gmx_trxframe *frame;

  frame = emalloc(sizeof(tW_gmx_trxframe));
  frame->counter = -1; // start at -1 so when we add 1, first frame index is 0.
  frame->eFileType = -1;
  frame->contents = emalloc(sizeof(tW_t_trxframe)); // MRD 82517

  frame->frame_size = (long) 0;
  frame->file_size = (long) 0;
  frame->n_frames = 0;
  frame->bInput = FALSE;
  frame->bDouble = TRUE;

/* Local Density Stuff MRD 06.17.2020 */
  frame->bLDs = FALSE; // assume no local densites, change later if needed
  frame->bLD_Grads = FALSE; // "

  frame->contents->t0 = (tW_real) 0.0;
  frame->contents->tf = (tW_real) 0.0;
  frame->contents->tpf = (tW_real) 0.0;
  frame->contents->tppf = (tW_real) 0.0;

  frame->setup_xdr = setup_xdr;

  frame->set_natoms = set_natoms;
  frame->set_atom_labels = set_atom_labels;
  frame->print_frame_dump = print_frame_dump;

  frame->get_pos_of = get_pos_of;
  frame->get_vel_of = get_vel_of;
  frame->get_box = get_box;

  //MCL 09.15.25 added to avoid undefined behavior when we check these variables (compiler dependent)
  frame->contents->title = NULL;
  frame->contents->ePBC = epbcXYZ; // Setting default to 3D periodic bounds

  //MCL 09.17.25 In the same spirit I am cleaning up the other non-initialized variables
  frame->fp= NULL; // Stuff directly in tW_gmx_trx
  frame->xids = NULL;
  frame->vids = NULL;
  frame->fids = NULL;

  frame->ld_id = NULL;
  frame->ldg_id = NULL;
  frame->linear_lds = NULL;
  frame->local_densities = NULL;
  frame->n_atomtypes = 0;
  frame->linear_ldgs = NULL;
  frame->local_density_gradients = NULL;

  //Stuff we use in contents (tW_t_trxframe)
  frame->contents->bDouble = (tW_gmx_bool) 0;
  frame->contents->bTitle = (tW_gmx_bool) 0;
  frame->contents->bStep = (tW_gmx_bool) 0;
  frame->contents->bX = (tW_gmx_bool) 0;
  frame->contents->bV = (tW_gmx_bool) 0;
  frame->contents->bF = (tW_gmx_bool) 0;
  frame->contents->bBox = (tW_gmx_bool) 0;
  frame->contents->bPrec = (tW_gmx_bool) 0;
  frame->contents->bAtoms = (tW_gmx_bool) 0;
  frame->contents->bBox = (tW_gmx_bool) 0;
  frame->contents->x = NULL;
  frame->contents->v = NULL;
  frame->contents->f = NULL;

  return frame;
}

/*****************************************************************************************
copy_trxframe_info(): copies the resolution-independent information from source to dest
*****************************************************************************************/
void copy_trxframe_info(tW_gmx_trxframe *source, tW_gmx_trxframe *dest)
{
  dest->counter = source->counter; // MRD 972017
  dest->contents->bTitle = source->contents->bTitle;
  dest->contents->title = source->contents->title;

  dest->contents->flags = source->contents->flags;

  dest->contents->bDouble = source->contents->bDouble;
  dest->contents->bTime = source->contents->bTime;
  if (dest->contents->bTime)
  {
    dest->contents->time = (tW_real) source->contents->time;
    dest->contents->t0 = (tW_real) source->contents->t0;
    dest->contents->tpf = (tW_real) source->contents->tpf;
    dest->contents->tppf = (tW_real) source->contents->tppf;
  }

  dest->contents->bStep = source->contents->bStep;
  if (dest->contents->bStep) { dest->contents->step = source->contents->step; }

  dest->contents->bLambda = source->contents->bLambda;
  if (source->contents->bLambda) { dest->contents->lambda = source->contents->lambda; }

  dest->contents->bPrec = source->contents->bPrec;
  if (source->contents->bPrec) { dest->contents->prec = source->contents->prec; }

  dest->contents->bBox = source->contents->bBox;
  if (source->contents->bBox) { copy_matrix(source->contents->box, dest->contents->box); }

  dest->contents->bX = source->contents->bX;
  dest->contents->bV = source->contents->bV;
  dest->contents->bF = source->contents->bF;
}

/*****************************************************************************************
clear_dvec(): This was copied from the GROMACS Library.
*****************************************************************************************/
void clear_dvec(dvec a)
{
  a[XX] = 0.0;
  a[YY] = 0.0;
  a[ZZ] = 0.0;
}


/*****************************************************************************************
map_trxframe(): maps source to dest according to the correspondence defined in map
*****************************************************************************************/
void map_trxframe(tW_gmx_trxframe *source, tW_gmx_trxframe *dest, tW_site_map *map)
{
  int I_site, i_atm, n_I, atm_no;

  dvec x, v, f, sc_x, sc_v;

  double c_Ii;

  for (I_site = 0; I_site < dest->contents->natoms; I_site++) 
  {
    n_I = map[I_site].n_atms;

    if (dest->contents->bX) { clear_dvec(dest->contents->x[I_site]); }
    if (dest->contents->bV) { clear_dvec(dest->contents->v[I_site]); }
    if (dest->contents->bF) { clear_dvec(dest->contents->f[I_site]); }


    for (i_atm = 0; i_atm < n_I; i_atm++) 
    {
      atm_no = map[I_site].i_atm[i_atm];

      c_Ii = map[I_site].c_Ii[i_atm];

      if (dest->contents->bX) 
      {
        copy_vector(source->contents->x[atm_no], x);
        scal_times_vect(c_Ii, x, sc_x);
        vect_inc(dest->contents->x[I_site], sc_x);
      }
      if (dest->contents->bV) 
      {
        copy_vector(source->contents->v[atm_no], v);
        scal_times_vect(c_Ii, v, sc_v);
        vect_inc(dest->contents->v[I_site], sc_v);
      }
      // NB at the moment we are assuming that no atom contributes 
      // to more than one cg sites, s.t. c_Ii = d_Ii and the 
      // \cal F_I = net atomistic force
      if (dest->contents->bF) 
      {
        copy_vector(source->contents->f[atm_no], f);
        vect_inc(dest->contents->f[I_site], f);
      }
    }
  }
}

/*****************************************************************************************
get_filename(): returns the filename corresponding to the argument opt. 
This replaces an old GROMACS function
*****************************************************************************************/
char *get_filename(tW_gmx_cgmap_input *self, const char *opt)
{
  int i;
  for (i = 0; i < N_MAP_FILES; ++i)
  {
    if (strcmp(self->files[i].flag,opt) == 0)
    {
      return self->files[i].fnm;
    }
  }  
  return NULL;
}

/*****************************************************************************************
arg_is_set(): returns TRUE if the argument corresponding to opt has been set, and FALSE otherwise
This replaces an old GROMACS function
*****************************************************************************************/
bool arg_is_set(tW_gmx_cgmap_input *self, const char *opt)
{
  int i;
  for (i = 0; i < N_MAP_FILES; ++i)
  {
    if (strcmp(opt,self->files[i].flag) == 0)
    {
      if ((self->use_files & self->files[i].file_shift) == 0) { return FALSE; }
      else { return TRUE; }
    }
  }
  for (i = 0; i < N_TIME_ARGS; ++i)
  {
    if (strcmp(opt,self->times[i].flag) == 0)
    {
      if ((self->file_times & self->times[i].time_shift) == 0) { return FALSE; }
      else { return TRUE; }
    }
  }
  for (i = 0; i < N_FTYPE_ARGS; ++i)
  {
    if (strcmp(opt,self->file_types[i].flag) == 0)
    {
      if ((self->ftflags & self->file_types[i].shift) == 0) { return FALSE; }
      else { return TRUE; }
    }
  }
  return FALSE;
}


/*****************************************************************************************
init_tW_gmx_cgmap_input(): Initializes a new instance of tW_gmx_cgmap_input
*****************************************************************************************/
tW_gmx_cgmap_input* init_tW_gmx_cgmap_input(int argc, char *argv[])
{
  int i, j, iarg;

  tW_gmx_cgmap_input *input;
//allocate memory
  input = ecalloc(1,sizeof(tW_gmx_cgmap_input));
//point structure's function pointers to functions
  input->get_filename = get_filename;
  input->arg_is_set = arg_is_set;

// Allocate memory for the data structures I made up. MRD
  input->map_file_flags = (char **) ecalloc(N_MAP_FILES,sizeof(char *));
  for (i = 0; i < N_MAP_FILES; ++i) { input->map_file_flags[i] = (char *) ecalloc(7,sizeof(char)); }

  input->time_flags = (char **) ecalloc(N_TIME_ARGS,sizeof(char *));
  for (i = 0; i < N_TIME_ARGS; ++i) { input->time_flags[i] = (char *) ecalloc(7,sizeof(char)); }

  input->file_type_flags = (char **) ecalloc(N_FTYPE_ARGS,sizeof(char *));
  for (i = 0; i < N_FTYPE_ARGS; ++i) { input->file_type_flags[i] = (char *) ecalloc(7,sizeof(char)); }

// Set the command line flags for the things we can read for cgmap
  strcpy(input->map_file_flags[eAA],"-f");
  strcpy(input->map_file_flags[eCG],"-o");
  strcpy(input->map_file_flags[eCGTPR],"-s");
  strcpy(input->map_file_flags[eMAP],"-p");
  strcpy(input->map_file_flags[eGRO1],"-c");

  strcpy(input->time_flags[eBTIME],"-b");
  strcpy(input->time_flags[eETIME],"-e");
  strcpy(input->time_flags[eDELTAT],"-dt");

  strcpy(input->file_type_flags[eAA],"-ipt");
  strcpy(input->file_type_flags[eCG],"-opt");
  strcpy(input->file_type_flags[eCGTPR],"-st");
// Initialize the file names
  strcpy(input->files[eAA].fnm,"unspecified");
  strcpy(input->files[eCG].fnm,"unspecified");
  strcpy(input->files[eCGTPR].fnm,"unspecified");
  strcpy(input->files[eMAP].fnm,"unspecified");
  strcpy(input->files[eGRO1].fnm,"unspecified");

//Copy flags over, set shifts and initial values for data structures' properties
  for ( i = 0; i < N_MAP_FILES; ++i)
  {
    input->files[i].flag = (char *) ecalloc(7,sizeof(char));
    strcpy(input->files[i].flag,input->map_file_flags[i]);
    input->files[i].file_shift = (1 << i);
    input->files[i].idx = i;
  }
  
  for ( i = 0; i < N_TIME_ARGS; ++i)
  {
    input->times[i].flag = (char *) ecalloc(7,sizeof(char));
    strcpy(input->times[i].flag,input->time_flags[i]);
    input->times[i].time_shift = (1 << i);
    input->times[i].idx = i;
    input->times[i].value = -1.0;
  }

  for ( i = 0; i < N_FTYPE_ARGS; ++i)
  {
    input->file_types[i].flag = (char *) ecalloc(7,sizeof(char));
    strcpy(input->file_types[i].flag,input->file_type_flags[i]);
    input->file_types[i].shift = (1 << i);
    input->file_types[i].value = -1;
  }

//These are the three elements that keep track of whether or not you used a specific cmd line arg
  input->use_files = 0;
  input->file_times = 0;
  input->ftflags = 0;

  iarg = 1;
  char **cfts;
  char *NA = "N/A";
  cfts = (char **) ecalloc(N_FTYPE_ARGS,sizeof(char *));
  for (i = 0; i < N_FTYPE_ARGS; ++i) { cfts[i] = NA; }
//Loop over command line arguments
  for (iarg = 1; iarg < argc; ++iarg)
  {
    for (i = 0; i < N_MAP_FILES; ++i) // Loop over possible files
    {
      if (strcmp(argv[iarg],input->files[i].flag) == 0) // check to see if flag provided matches this file's flag
      {
	if (iarg + 1 == argc) // Check to see if the flag was the last command line argument
	{
	  fprintf(stderr,"ERROR: command line flag %s listed last with no argument!\n",argv[iarg]);
	  exit(1);
	}
	if ((input->use_files & input->files[i].file_shift) != 0) // Check to see if you already specified this command line argument
	{
	  fprintf(stderr,"ERROR: specified command line flag %s twice!\n",argv[iarg]);
	  exit(1);
	}
	input->use_files |= input->files[i].file_shift; // Track that you have now seen this command line flag
	strcpy(input->files[i].fnm,argv[iarg+1]); // copy the filename that follows the flag
      }
    }
    for (i = 0; i < N_TIME_ARGS; ++i)
    {
      if (strcmp(argv[iarg],input->times[i].flag) == 0)
      {
	if (iarg + 1 == argc)
	{
	  fprintf(stderr,"ERROR: command line flag %s listed last with no argument!\n",argv[iarg]);
	  exit(1);
	}
	if ((input->file_times & input->times[i].time_shift) != 0)
	{
	  fprintf(stderr,"ERROR: specified command line flag %s twice!\n",argv[iarg]);
	  exit(1);
	}
	input->file_times |= input->times[i].time_shift;
	input->times[i].value = atof(argv[iarg+1]);
      }
    }
    for (i = 0; i < N_FTYPE_ARGS; ++i)
    {
      if (strcmp(argv[iarg],input->file_type_flags[i]) == 0)
      {
	if (iarg + 1 == argc)
	{
	  fprintf(stderr,"ERROR: command line flag %s listed last with no argument!\n",argv[iarg]);
	  exit(1);
	}
	if ((input->ftflags & input->file_types[i].shift) != 0)
	{
	  fprintf(stderr,"ERROR: specified command line flag %s twice!\n",argv[iarg]);
	  exit(1);
	}
	input->ftflags |= input->file_types[i].shift;
	cfts[i] = argv[iarg+1];
	if (strcmp(argv[iarg+1],DUMP) == 0) { input->file_types[i].value = eDUMP; }
	else if (strcmp(argv[iarg+1],BOCS) == 0) { input->file_types[i].value = eBOCS; }
	else if (strcmp(argv[iarg+1],GRO) == 0) { input->file_types[i].value = eGRO; }
	else if (strcmp(argv[iarg+1],TRJ) == 0) { input->file_types[i].value = eTRJ; }
	else if (strcmp(argv[iarg+1],TRR) == 0) { input->file_types[i].value = eTRR; }
	else { fprintf(stderr,"ERROR: unrecognized file type %s\n",argv[iarg+1]); exit(1); }
      }
    }
  }

  fprintf(stderr,"Option               Filename  Type         Description\n");
  fprintf(stderr,"----------------------------------------------------------------------\n");
  fprintf(stderr,"  -p     %20s  Input        Mapping Topology file\n",input->files[eMAP].fnm);
  fprintf(stderr,"  -f     %20s  Input        Trajectory: aa_trajectory.txt\n",input->files[eAA].fnm);
  fprintf(stderr,"  -o     %20s  Output, Opt  Trajectory: cg_mapped_traj.txt\n",input->files[eCG].fnm);
  fprintf(stderr,"  -s     %20s  Input, Opt.  CG Topology File\n",input->files[eCGTPR].fnm);
  fprintf(stderr,"  -c     %20s  Output, Opt. Coordinate file in Gromos-87 format\n",input->files[eGRO1].fnm);
  fprintf(stderr,"\n");
  fprintf(stderr,"Option       Type     Value     Description\n");
  fprintf(stderr,"-----------------------------------------------------------------\n");
  fprintf(stderr,"-[no]h       bool               Print help info and quit\n");
  fprintf(stderr,"-ipt       string   % 9s   Flag for what type of file -f is\n",cfts[eAA]);
  fprintf(stderr,"                                options: %s %s %s %s %s\n",DUMP,GRO,BOCS,TRJ,TRR);
  fprintf(stderr,"-opt       string   % 9s   Flag for what type of file -o is\n",cfts[eCG]);
  fprintf(stderr,"                                options: %s %s %s %s %s\n",DUMP,GRO,BOCS,TRJ,TRR);
  fprintf(stderr,"-st        string   % 9s   Flag for what type of file -s is\n",cfts[eCGTPR]);
  fprintf(stderr,"                                options: %s %s\n",DUMP,BOCS);
  fprintf(stderr,"-b           time   % 7f   First frame (ps) to read from trajectory\n",input->times[eBTIME].value);
  fprintf(stderr,"-e           time   % 7f   Last frame (ps) to read from trajectory\n",input->times[eETIME].value);
  fprintf(stderr,"-dt          time   % 7f   Only use frame when t MOD dt = first time (ps)\n",input->times[eDELTAT].value);

  // Check for errors
  // Check to see if you want to output a GRO file but forgot to specify a CG topology file
  if ((input->file_types[eCG].value == eGRO) && (! input->arg_is_set(input,"-s")))
  {
    fprintf(stderr,"ERROR: you requested the CG trajectory to be output in the .gro format, but did not provide a CG topology file\nPlease use the -s flag\n");
    exit(1);
  }
  // Check to see if you specified a file that needs a type specification, but forgot to provide said specification
  for (i = 0; i < N_FTYPE_ARGS; ++i)
  {
    if ((input->arg_is_set(input,input->files[i].flag)) && (!input->arg_is_set(input,input->file_types[i].flag)))
    {
      fprintf(stderr,"ERROR: you specified file %s using flag %s without specifying its type\nPlease use flag %s\n",input->files[i].fnm,input->files[i].flag,input->file_types[i].flag);
      exit(1);
    }
  }

  if ((input->use_files & input->files[eMAP].file_shift) == 0)
  {
    fprintf(stderr,"ERROR: Please provide a mapping topology file with flag -p\n");
    exit(1);
  }

  for (i=0;i<argc;i++) { if (strcmp(argv[i],"-h") == 0) { exit(0); } }


  return input;
}


/*****************************************************************************************
copy_trr_2_CGstruct(): copies position and force information from one data structure to another
*****************************************************************************************/
bool copy_trr_2_CGstruct(tW_gmx_trxframe *fr, tW_CG_site CG_struct[])
{
  int i, j;
  int n_atms = fr->contents->natoms;

  for (i = 0; i < n_atms; i++) 
  {
    if (fr->contents->bX) { copy_vector(fr->contents->x[i], CG_struct[i].r); }
    if (fr->contents->bF) { copy_vector(fr->contents->f[i], CG_struct[i].f); }
  }

  /* MRD 06.18.2020 */
  if (fr->bLDs)
  {
    for (i = 0; i < n_atms; ++i)
    {
      for (j = 0; j < fr->n_atomtypes; ++j)
      {
        CG_struct[i].LDs[j] = fr->local_densities[i][j];
      }
    }
    if (fr->bLD_Grads)
    {
      for (i = 0; i < n_atms; ++i)
      {
        for (j = 0; j < fr->n_atomtypes; ++j)
        {
          CG_struct[i].LDGrads[j][0] = fr->local_density_gradients[i][j][0]; //MCL i points the site index i, while j points to the type of particle surrounding site i, the third dimension is x y or z
          CG_struct[i].LDGrads[j][1] = fr->local_density_gradients[i][j][1];
          CG_struct[i].LDGrads[j][2] = fr->local_density_gradients[i][j][2];
        }
      }
    }
  }

  return fr->contents->bF;
}

/*****************************************************************************************
copy_trr_2_CGstruct_ref(): JFR - 07.16.12: for reading in reference forces from a trr file
*****************************************************************************************/
bool copy_trr_2_CGstruct_ref(tW_gmx_trxframe *fr, tW_CG_site CG_struct[])
{
  int i;
  int n_atms = fr->contents->natoms;

  for (i = 0; i < n_atms; i++) { if (fr->contents->bF) { copy_vector(fr->contents->f[i], CG_struct[i].ref_f); } }

  return fr->contents->bF;
}


/*****************************************************************************************
read_trr_2_CGstruct(): Tries to read another frame from trr file. If successful, copies
the info from frame data structure to CG_site data structure.
*****************************************************************************************/
bool read_trr_2_CGstruct(tW_gmx_info * info, tW_gmx_trxframe *fr, tW_CG_site CG_struct[])
{
  bool b_trr;

  b_trr = read_next_frame(fr, FALSE); // read_next_trxframe -> read_next_frame MRD 11.09.17

  if (b_trr) 
  {
    copy_trr_2_CGstruct(fr, CG_struct);
    update_info_trr(info, fr);
  }

  return b_trr;
}

/*****************************************************************************************
read_trr_2_CGstruct_ref(): JFR - 07.16.12: for reading in reference forces from a trr file
*****************************************************************************************/
bool read_trr_2_CGstruct_ref(tW_gmx_info * info, tW_gmx_trxframe *fr, tW_CG_site CG_struct[])
{
  bool b_trr;

  b_trr = read_next_frame(fr, FALSE); // read_next_trxframe -> read_next_frame MRD 11.09.17

  if (b_trr) 
  {
    copy_trr_2_CGstruct_ref(fr, CG_struct);
    update_info_trr(info, fr);
  }

  return b_trr;
}


/*****************************************************************************************
update_info_trr():
*****************************************************************************************/
int update_info_trr(tW_gmx_info *info, tW_gmx_trxframe *fr)
{
  int i, j;

  info->b_Gromacs = TRUE;
  info->b_Forces = fr->contents->bF;
  info->b_PBC = fr->contents->bBox;

  if (fr->contents->bF) { info->b_Forces_1 = TRUE; } 
  else { info->b_Forces_N = FALSE; }

  if (fr->contents->bBox) 
  {
    for (i = 0; i < DIM; i++) { for (j = 0; j < DIM; j++) { info->box[i][j] = fr->contents->box[i][j]; } }
  }

  return 0;
}

/* MRD FUNCTIONS */
/*
read_first_frame(): This function is a wrapper function. It calls the correct function to 
	read the header information from a simulation, based on the value of fr->eFileType.
	If for some catastrophic reason fr->eFileType has not been changed since calling
	init_tW_gmx_trxframe, it will still have it's initialized value of -1. If this is 
	the case, we check the filename for the default file extensions. If one is found,
	we set fr->eFileType. If one is not found, a long error message is printed. It returns
	the number of atoms/particles per frame of the trajectory file.
	Note: None of the read_first_XXX_frame functions actually store the box matrix
	nor the xvf vectors from the first frame. To get those, you must call read_next_frame
	after calling read_first_frame.
*/
int read_first_frame(tW_gmx_trxframe *fr, const char *fnm)
{
  strcpy(fr->filename,fnm);

  if (fr->eFileType == -1)
  {
    if (strstr(fnm,".dump") != NULL) { fr->eFileType = eDUMP; }
    else if (strstr(fnm,".gro") != NULL) { fr->eFileType = eGRO; }
    else if (strstr(fnm,".btj") != NULL) { fr->eFileType = eBOCS; }
    else if (strstr(fnm,".trj") != NULL) { fr->eFileType = eTRJ; }
    else if (strstr(fnm,".trr") != NULL) { fr->eFileType = eTRR; }
    else if (strstr(fnm,".lmp") != NULL) { fr->eFileType = eLMP; }
    else
    {
      fprintf(stderr,"ERROR: trajectory eFileType == -1\n");
      fprintf(stderr,"\tThis means the trajectory file type was never set\n");
      fprintf(stderr,"\tI checked filename \"%s\" for file extensions " 
		"\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\"\n"
		,fnm,".dump",".btj",".gro",".trj",".trr",".lmp");
      fprintf(stderr,"\twhich correspond to topology file types %s, %s, %s, %s, %s, and %s,"
		" respectively.\n",DUMP,BOCS,GRO,TRJ,TRR,LMP);
      fprintf(stderr,"\tHowever, none of the file extensions were found\n");
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
      exit(1);
    }
  } 

  switch (fr->eFileType)
  {
    case eDUMP:
      return read_first_dump_frame(fr,fnm);
      break;
    case eGRO:
      return read_first_gro_frame(fr,fnm);
      break;
    case eBOCS:
      return read_first_bocs_frame(fr,fnm);
      break;
    case eTRJ:
      return tW_read_first_trj_frame(fr,fnm);
      break;
    case eTRR:
      return tW_read_first_trr_frame(fr,fnm);
      break;
    case eLMP:
      return read_first_lammps_frame(fr,fnm);
      break;
    default:
      fprintf(stderr,"ERROR: unknown file type: %d\n",fr->eFileType);
      fprintf(stderr,"\tSupported types: %d %s, %d %s, %d %s, %d %s, %d %s, %d %s\n",
                        eDUMP,DUMP,eBOCS,BOCS,eGRO,GRO,eTRJ,TRJ,eTRR,TRR,eLMP,LMP);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
      exit(1);
      break;
  }
  return -1;
}


/*
read_next_frame(): This function is a wrapper, like read_first_frame. It calls the proper 
	function to read the next frame and store the information in fr->contents->XXXX,
	based (again) on the value of fr->eFileType. If that value is somehow -1, that means
	not only was it never changed after the frame was initialized with init_tW_gmx_trxframe,
	but read_first_frame was never called either, and accordingly no memory has been allocated
	for xvf vectors. Returns true if everything went ok.
*/
bool read_next_frame(tW_gmx_trxframe *fr, bool printCount)
{
  if (fr->eFileType == -1)
  {
    fprintf(stderr,"ERROR: trajectory frame eFileType = -1.\n");
    fprintf(stderr,"\tread_first_frame must have never been called\n");
    fprintf(stderr,"\tERROR: %s %d\n",__FILE__,__LINE__);
  }

  ++(fr->counter);
  if (printCount) { fprintf(stderr,"Reading frame %d\r",fr->counter); }

  switch (fr->eFileType)
  {
    case eDUMP:
      return read_next_dump_frame(fr);
      break;
    case eGRO:
      return read_next_gro_frame(fr);
      break;
    case eBOCS:
      return new_read_next_bocs_frame(fr);
      break;
    case eTRJ:
      return tW_read_next_trj_frame(fr);
      break;
    case eTRR:
      return tW_read_next_trr_frame(fr);
      break;
    case eLMP:
      return read_next_lammps_frame(fr);
      break;
    default:
      fprintf(stderr,"ERROR: unknown file type: %d\n",fr->eFileType);
      fprintf(stderr,"\tSupported types: %d %s, %d %s, %d %s, %d %s, %d %s, %d %s\n",
                        eDUMP,DUMP,eBOCS,BOCS,eGRO,GRO,eTRJ,TRJ,eTRR,TRR,eLMP,LMP);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
      exit(1);
      break;
  }
  return FALSE;
}

/*
open_write_trajectory(): This function opens the file fnm for writing the proper way
	depending on the value stored in fr->eFileType. If fr->eFileType == -1, it was
	never changed after being initialized in init_tW_gmx_trxframe, and we check for
	default file extensions to guess at what output type the user wanted. If none
	are found, we default to BOCS, and print a noisy warning message. 
*/

void open_write_trajectory(tW_gmx_trxframe *fr, char *fnm)
{
  strcpy(fr->filename,fnm);
  if (fr->eFileType == -1)
  {
    if (strstr(fnm,".trj") != NULL) { fr->eFileType = eTRJ; }
    else if (strstr(fnm,".trr") != NULL) { fr->eFileType = eTRR; }
    else if (strstr(fnm,".dump") != NULL) { fr->eFileType = eDUMP; }
    else if (strstr(fnm,".gro") != NULL) { fr->eFileType = eGRO; }
    else if (strstr(fnm,".btj") != NULL) { fr->eFileType = eBOCS; }
    else if (strstr(fnm,".data") != NULL) { fr->eFileType = eLMPDATA; }
    else if (strstr(fnm,".lmp") != NULL) { fr->eFileType = eLMPTRJ; }
    else
    {
      fr->eFileType = eBOCS;
      fprintf(stderr,"WARNING: eFileType for output file %s was never set!\n",fnm);
      fprintf(stderr,"\tI checked the file name for extensions \"%s\', \"%s\', \"%s\', \"%s\', and \"%s\', but none were found\n",".trj",".trr",".dump",".gro",".btj");
      fprintf(stderr,"\tAccordingly, I defaulted to the assignment of type %s\n",BOCS);
      fprintf(stderr,"\tIf you wanted a different trajectory type, either rerun the current program\n");
      fprintf(stderr,"\tand make sure your output filename has one of the default extensions, or\n");
      fprintf(stderr,"\ttranslate the %s trajectory output using the translator program\n",BOCS);
    }
  }
  switch (fr->eFileType)
  {
    case eTRR:
      fr->setup_xdr(fr,fnm,FALSE);
      break;
    case eTRJ:
      fr->fp = fopen(fnm,"wb");
      break;
    default:
      fr->fp = open_file(fnm, 'w');    
      break;
  }
}

/*
write_frame(): This function calls the proper function to write the current info in fr
	based on the file type in fr->eFileType. First, it checks if eFileType = -1,
	in which case open_write_trajectory has not been called, and it calls it after
	printing a noisy warning message.
*/
void write_frame(tW_gmx_trxframe *fr, tW_gmx_topology *top)
{
  if (fr->eFileType == -1)
  {
    fprintf(stderr,"WARNING: write_frame has been called before open_write_trajectory\n");
    fprintf(stderr,"\tThe desired output filename is unknown, so your output trajectory\n");
    fprintf(stderr,"\thas been given the default name \"output_trajectory.txt\". It is\n");
    fprintf(stderr,"\tin the %s trajectory format.\n",BOCS);
    open_write_trajectory(fr,"output_trajectory.txt");
  }

  switch (fr->eFileType)
  {
    case eDUMP:
      write_dump_frame(fr);
      break;
    case eGRO:
      write_gro_frame(fr);
      break;
    case eBOCS:
      write_bocs_frame(fr);
      break;
    case eTRJ:
      wrap_tW_write_trj_frame(fr);
      break;
    case eTRR:
      wrap_tW_write_trr_frame(fr);
      break;
    case eLMPTRJ:
      write_lammps_frame(fr, top);
      break;
    case eLMPDATA:
      write_lammps_data(fr, top);
      break;
    default:
      fprintf(stderr,"ERROR: unknown eFileType: %d\n",fr->eFileType);
      fprintf(stderr,"\tSupported types: %d %s, %d %s, %d %s, %d %s, %d %s, %d %s, %d %s\n",
			eDUMP,DUMP,eBOCS,BOCS,eGRO,GRO,eTRJ,TRJ,eTRR,TRR,eLMPTRJ,LMPTRJ,eLMPDATA,LMPDATA);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
      exit(1);
      break;
  }
}

/*
read_topology(): This function calls the proper function for reading the topology file fnm
	based on the value of top->eFileType. if eFileType = -1, it checks for the default
	file extensions, and assigns it. If none of the default file extensions are found, 
	it prints a noisy error message and terminates the program.
*/

bool read_topology(tW_gmx_topology *top, const char *fnm)
{
  
  // top->eFileType is initialized to -1. If it hasnt been changed, guess the type based on the filename.
  if (top->eFileType == -1)
  {
    if (strstr(fnm,".dump") != NULL) { top->eFileType = eDUMP; }
    else if (strstr(fnm,".btp") != NULL) { top->eFileType = eBOCS; }
    else
    {
      fprintf(stderr,"ERROR: topology eFileType == -1\n");
      fprintf(stderr,"\tThis means the topology file type was never set\n");
      fprintf(stderr,"\tI checked filename \"%s\" for file extensions \"%s\" and \"%s\"\n",fnm,".dump",".btp");
      fprintf(stderr,"\twhich correspond to topology file types %s and %s, respectively.\n",DUMP,BOCS);
      fprintf(stderr,"\tHowever, none of the file extensions were found\n");
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
      exit(1);
    }
  }

  switch (top->eFileType)
  {
    case eDUMP:
      return read_tpr_dump((char *)fnm,top);
      break;
    case eBOCS:
      return read_bocs_top((char *)fnm,top);
      break;
    default:
      fprintf(stderr,"ERROR: unrecognized topology eFileType: %d\n",top->eFileType);
      fprintf(stderr,"\tSupported types: %d %s, %d %s\n",eDUMP,DUMP,eBOCS,BOCS);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
      exit(1);
      break;
  }
  return FALSE;

}


/*****************************************************************************************
get_word_count_delim(): Gets the number of words separated by chars in delims in a line
*****************************************************************************************/

int get_word_count_delim(tW_line inp_line, const char * delims)
{
  int n_words = 0;
  int i, n_delims = 0;
  char **news, *new;
  char a_delim[2];
  char *p_delim;
  a_delim[0] = '!';
  a_delim[1] = '\0';
  p_delim = &(a_delim[0]);
  i = 0;
  while (delims[i] != '\0')
  {
    ++n_delims;
    ++i;
  }
  news = (char **) ecalloc(n_delims, sizeof(char *));
  new = &(inp_line[0]);
  do
  {
    ++n_words;
    bool flag;
    do
    {
      flag = FALSE;
      for (i = 0; i < n_delims; ++i) 
      { 
	if (*new == delims[i]) 
	{ 
	  flag = TRUE; 
	} 
      }
      if (flag) { ++new; }
    } while (flag);
    flag = TRUE;
    if ((*new == '\n') || (*new == '\0'))
    {
      flag = FALSE;
    }
    for (i = 0; i < n_delims; ++i) 
    { 
      a_delim[0] = delims[i]; 
      news[i] = strstr(new,p_delim); 
    }
    new = NULL;
    for (i = 0; i < n_delims; ++i) 
    { 
      if (news[i] != NULL) 
      { 
	new = news[i]; break; 
      } 
    }
    for (i = 0; i < n_delims; ++i) 
    { 
      if (news[i] != NULL) 
      { 
	if (news[i] < new) 
	{ 
	  new = news[i]; 
	} 
      } 
    }
    if (!flag) { new = NULL; }   
  } while (new != NULL);
  free(news);
  return n_words;
}

/*****************************************************************************************
get_words_delim(): supposed to be like python's string.split()
*****************************************************************************************/

void get_words_delim(tW_line inp_line, const char * delims, tW_word * word_list)
{
  int n_words = get_word_count_delim(inp_line, delims);
  int word_count = -1;
  int i;
  int n_delims = 0;
//  tW_word *word_list = (tW_word *) calloc(n_words,sizeof(tW_word));
  char **news;
  char *new;
  char *old;
  char *p_delim;
  char a_delim[2];
  a_delim[0] = '!';
  a_delim[1] = '\0';
  p_delim = &(a_delim[0]);
  new = &(inp_line[0]);

  i = 0;
  while (delims[i] != '\0')
  {
    ++n_delims;
    ++i;
  }
  
  news = (char **) ecalloc(n_delims,sizeof(char *));
  do
  {
    ++word_count;
    bool flag;
    do
    {
      flag = FALSE;
      for (i = 0; i < n_delims; ++i) { if (delims[i] == *(new)) { flag = TRUE; } }
      if (flag) { ++new; }
    } while (flag);
    flag = TRUE; 
    if ((*new == '\n') || (*new == '\0'))
    { 
      flag = FALSE;
    }
    old = new;
    for (i = 0; i < n_delims; ++i) 
    { 
      a_delim[0] = delims[i]; 
//      p_delim = &(a_delim[0]); 
      news[i] = strstr(new,p_delim); 
    }
    new = NULL;
    for (i = 0; i < n_delims; ++i) 
    { 
      if (news[i] != NULL) 
      { 
	new = news[i]; break; 
      } 
    }
    for (i = 0; i < n_delims; ++i) 
    { 
      if (news[i] != NULL) 
      { 
        if (news[i] < new) 
        { 
	  new = news[i]; 
	} 
      } 
    }
    
    if (new == NULL) 
    { 
      strcpy(word_list[word_count],old); 
    }
    else 
    {
      char hold = *new;
      *new = '\0';
      strcpy(word_list[word_count],old);
      *new = hold;
    }
    if (!flag) { new = NULL; }
  } while (new != NULL);
  free(news);
}

/*****************************************************************************************
check_count(): For command line argument parsing.
*****************************************************************************************/

void check_count(int iarg, int argc, char * opt)
{
  if (iarg >= argc)
  {
    fprintf(stderr,"ERROR: unable to find option following flag %s\n",opt);
    exit(1);
  }
}

/*****************************************************************************************
test_line(): Tests to see if a line contains "term" against whether or not it 
is expected to ("want"). Prints errmsg if the line contains the term and shouldn't,
or if the line doesn't contain the term and should.
*****************************************************************************************/
void test_line(tW_line inp_line, const char * term, tW_gmx_bool want, const char * errmsg)
{
  if (want)
  {
    if (strstr(inp_line,term) == NULL)
    {
      fprintf(stderr,"ERROR: %s\n",errmsg);
      fprintf(stderr,"\tline: %s\n",inp_line);
      exit(1);
    }
  }
  else
  {
    if (strstr(inp_line,term) != NULL)
    {
      fprintf(stderr,"ERROR: %s\n",errmsg);
      fprintf(stderr,"\tline: %s",inp_line);
      exit(1);
    }
  }
}

/*****************************************************************************************
elim_char(): Replaces occurrances of char "elim" with '\0' in word.
Useful because when reading some strings from some file formats, there's a comma
(or other character) at the end of the string and we want it gone.
*****************************************************************************************/
void elim_char(tW_word word, char elim)
{
  int let_idx = 0;
  while (word[let_idx] != '\0')
  {
    if (word[let_idx] == elim) { word[let_idx] = '\0'; }
    ++let_idx;
  }
}

/*****************************************************************************************
 * get_quoted_words: 
 * For an input line of text containing a known number of "", extract an expected number of words from between the quotes
 * e.g. If line = ' \"opls_136\",nameB=\"opls_135\" ' and expected_words = 2, the opls_13X will be checked
 * @param line a list of words
 * @param expected_words the number of quote pairs expected within lines
 * @return Error if words found != expected_words, OR pointer to list of found words [len = expected_words]
 ****************************************************************************************/
char **get_quoted_words(tW_line line, int expected_words)
{
  char *sst, *sse;
  char **the_words = (char **) ecalloc(expected_words,sizeof(char *));
  int n_quotes = 0, i;
  sse = &(line[0]);
  --sse;
  do // count how many times we find \" in string : [point to begining, find where \", search next part of string, stop when no more quotes]
  {
    ++sse;
    sse = strstr(sse,"\"");
    ++n_quotes;
  } while (sse != NULL);
  if (n_quotes/2 != expected_words)
  {
    fprintf(stderr,"ERROR: in string: %s\n",line);
    fprintf(stderr,"\tI found %d quotes, but I expected to find %d quotes (2 for each of %d words\n",n_quotes,expected_words*2,expected_words);
    exit(1);
  }
  sse = &(line[0]);
  for (i = 0; i < expected_words; ++i)
  {
    the_words[i] = (char *) ecalloc(50,sizeof(char));
    sse = strstr(sse,"\"");
    ++sse;
    sst = sse;
    sse = strstr(sse,"\"");
    *(sse) = '\0';
    strcpy(the_words[i],sst);
    *(sse) = '\"';
    ++sse;
  }
  return the_words; 
}

/*****************************************************************************************
dump_molecule_info(): Prints some basic info for each molecule type 
from the tW_molecule struct.
*****************************************************************************************/
void dump_molecule_info(tW_gmx_topology *top) 
{
  tW_molecule * my_mols = top->molecules;
  int n = top->contents->mols.nr;
  FILE *fp = fopen("dump_my_mols.txt","w");
  int i,j,k;
  for (i = 0; i < n; ++i)
  {
    fprintf(fp,"molecule index: %d\n",i);
    fprintf(fp,"molname: %s\n",my_mols[i].molname);
    fprintf(fp,"n_mols: %d\n",my_mols[i].n_mols);
    fprintf(fp,"n_apm: %d\n",my_mols[i].n_apm);

    fprintf(fp,"idx\tatom_names\tatom_types\ttype_ids\n");
    for (j = 0; j < my_mols[i].n_apm; ++j)
    {
       fprintf(fp,"%3d\t%10s\t%10s\t%8d\n",j,my_mols[i].atom_names[j],my_mols[i].atom_types[j],my_mols[i].type_ids[j]);
    }   
    
    fprintf(fp,"n_excls: %d\n",my_mols[i].n_excls);
    fprintf(fp,"n_exclsa: %d\n",my_mols[i].n_exclsa);
    fprintf(fp,"excls table:\n");
    for (j = 0; j < my_mols[i].n_apm; ++j)
    {
      for (k = 0; k < my_mols[i].n_epa[j]; ++k)
      {
	fprintf(fp,"%3d ",my_mols[i].excls[j][k]);
      }
      fprintf(fp,"\n");
    }    
    if (my_mols[i].n_bonds > 0)
    {
      fprintf(fp,"n_bonds: %d\n",my_mols[i].n_bonds);
      if (top->int_map)
      {
        fprintf(fp,"idx  type  at_i  at_j  LMP_ID\n");
        for (j = 0; j < my_mols[i].n_bonds; ++j)
        {
          fprintf(fp,"%3d  %4d  %4d  %4d  %6d\n",j,my_mols[i].bond_types[j],my_mols[i].bond_ij[j][0],my_mols[i].bond_ij[j][1],top->int_map[my_mols[i].bond_types[j]]);
        }
      }
      else
      {
        fprintf(fp,"idx  type  at_i  at_j\n");
        for (j = 0; j < my_mols[i].n_bonds; ++j)
        {
          fprintf(fp,"%3d  %4d  %4d  %4d\n",j,my_mols[i].bond_types[j],my_mols[i].bond_ij[j][0],my_mols[i].bond_ij[j][1]);
        }
      }
    }
    
    if (my_mols[i].n_angles > 0)
    {
      fprintf(fp,"n_angles: %d\n",my_mols[i].n_angles);
      if (top->int_map)
      {
        fprintf(fp,"idx  type  at_i  at_j  at_k  LMP_ID\n");
        for (j = 0; j < my_mols[i].n_angles; ++j)
        { 
          fprintf(fp,"%3d  %4d  %4d  %4d  %4d  %6d\n",j,my_mols[i].angle_types[j],my_mols[i].angle_ijk[j][0],my_mols[i].angle_ijk[j][1],my_mols[i].angle_ijk[j][2],top->int_map[my_mols[i].angle_types[j]]);
        }
      }
      else
      {
        fprintf(fp,"idx  type  at_i  at_j  at_k\n");
        for (j = 0; j < my_mols[i].n_angles; ++j)
        {
          fprintf(fp,"%3d  %4d  %4d  %4d  %4d\n",j,my_mols[i].angle_types[j],my_mols[i].angle_ijk[j][0],my_mols[i].angle_ijk[j][1],my_mols[i].angle_ijk[j][2]);
        }
      }
    }

    if (my_mols[i].n_dihs > 0)
    {
      fprintf(fp,"n_dihs: %d\n",my_mols[i].n_dihs);
      if (top->int_map)
      {
        fprintf(fp,"idx  type  at_i  at_j  at_k  at_l  LMP_ID\n");
        for (j = 0; j < my_mols[i].n_dihs; ++j)
        {
          fprintf(fp,"%3d  %4d  %4d  %4d  %4d  %4d  %6d\n",j,my_mols[i].dih_types[j],my_mols[i].dih_ijkl[j][0],my_mols[i].dih_ijkl[j][1],my_mols[i].dih_ijkl[j][2],my_mols[i].dih_ijkl[j][3],top->int_map[my_mols[i].dih_types[j]]);
        }
      }
      else
      {
        fprintf(fp,"idx  type  at_i  at_j  at_k  at_l\n");
        for (j = 0; j < my_mols[i].n_dihs; ++j)
        {
          fprintf(fp,"%3d  %4d  %4d  %4d  %4d  %4d\n",j,my_mols[i].dih_types[j],my_mols[i].dih_ijkl[j][0],my_mols[i].dih_ijkl[j][1],my_mols[i].dih_ijkl[j][2],my_mols[i].dih_ijkl[j][3]);
        }
      }
    }

    fprintf(fp,"~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  }
  fclose(fp);
}

/*****************************************************************************************
pop_contents(): pop here is short for populate. Maybe I should make that more explicit...
This function populates the top->contents->atoms data structure based on what was 
previously read and stored in the top->molecules data structure
*****************************************************************************************/
void pop_contents(tW_gmx_topology * top)
{
  int test_sscanf, inp_int, n_atoms;
  int i, j, k, l, let_idx;
  int n_moltypes = top->contents->mols.nr;

  n_atoms = top->contents->atoms.nr; 

  tW_molecule *my_mols;
  my_mols = (top->molecules);


  top->contents->atoms.atom = (tW_t_atom *) ecalloc(n_atoms, sizeof(tW_t_atom));

  top->contents->atoms.atomname = (char ***) ecalloc(n_atoms,sizeof(char **));
  top->contents->atoms.atomtype = (char ***) ecalloc(n_atoms,sizeof(char **));
  top->contents->atoms.atomtypeB = (char ***) ecalloc(n_atoms,sizeof(char **));
  for (i = 0; i < n_atoms; ++i)
  {
    (top->contents->atoms.atomname[i]) = (char **) ecalloc(1,sizeof(char *));
    (top->contents->atoms.atomtype[i]) = (char **) ecalloc(1,sizeof(char *));
    (top->contents->atoms.atomtypeB[i]) = (char **) ecalloc(1,sizeof(char *));
  }
  for (i = 0; i < n_atoms; ++i)
  {
    top->contents->atoms.atomname[i][0] = (char *) ecalloc(10,sizeof(char));
    top->contents->atoms.atomtype[i][0] = (char *) ecalloc(10,sizeof(char));
    top->contents->atoms.atomtypeB[i][0] = (char *) ecalloc(10,sizeof(char));
  }


// MRD 01.17.2020 This still works to get the total nres
  top->contents->atoms.nres = 0;
  for (i = 0; i < n_moltypes; ++i)
  {
    top->contents->atoms.nres += my_mols[i].n_res * my_mols[i].n_mols;
  }
  top->contents->atoms.resinfo = (tW_t_resinfo *) ecalloc(top->contents->atoms.nres,sizeof(tW_t_resinfo));
  
  for (i = 0; i < top->contents->atoms.nres; ++i)
  {
    top->contents->atoms.resinfo[i].name = (char **) ecalloc(1,sizeof(char *));    
  }
  for (i = 0; i < top->contents->atoms.nres; ++i)
  {
      *(top->contents->atoms.resinfo[i].name) = (char *) ecalloc(10,sizeof(char));
  }
  
  int resnr = 0;
  int it, im, ia, iattype;
  int molidx;
  int at_idx = 0;
  bool new_type;

  top->n_atomtypes = 0;
  char **at_type_names = (char **) calloc(1000,sizeof(char *));

  for (im = 0; im < top->n_molblocks; ++im)
  {
    molidx = top->molblocks[im].moltype;
    for (it = 0; it < top->molblocks[im].n_mols; ++it)
    {
      for (ia = 0; ia < my_mols[molidx].n_apm; ++ia)
      {
        top->contents->atoms.atom[at_idx].m = my_mols[molidx].m[ia];
        top->contents->atoms.atom[at_idx].q = my_mols[molidx].q[ia];
        top->contents->atoms.atom[at_idx].mB = my_mols[molidx].m[ia];
        top->contents->atoms.atom[at_idx].qB = my_mols[molidx].q[ia];
        top->contents->atoms.atom[at_idx].type = my_mols[molidx].type_ids[ia];
        top->contents->atoms.atom[at_idx].resind = resnr + my_mols[molidx].residx[ia];
        strcpy(*(top->contents->atoms.atomname[at_idx]),my_mols[molidx].atom_names[ia]);
        strcpy(*(top->contents->atoms.atomtype[at_idx]),my_mols[molidx].atom_types[ia]);
        strcpy(*(top->contents->atoms.atomtypeB[at_idx]),my_mols[molidx].atom_Btypes[ia]);
      
        new_type = TRUE;
        for (iattype = 0; iattype < top->n_atomtypes; ++iattype)
        {
          if (strcmp(my_mols[molidx].atom_types[ia],at_type_names[iattype]) == 0) { new_type = FALSE; }
        }
        if (new_type)
        {
          at_type_names[top->n_atomtypes] = my_mols[molidx].atom_types[ia];
          ++(top->n_atomtypes);
        }
        ++at_idx;
      }
      resnr += my_mols[molidx].n_res;
    }
  } 



////////////
  top->atom_type_names = (tW_word *) calloc(top->n_atomtypes,sizeof(tW_word));
  for (iattype = 0; iattype < top->n_atomtypes; ++iattype)
  {
    strcpy(top->atom_type_names[iattype],at_type_names[iattype]);
  }
  free(at_type_names);
////////////

  // Next we need excls stuff
  
// MRD 01.21.2020 This should still work for n_moltypes != n_molblocks
  top->contents->excls.nr = 0;
  top->contents->excls.nra = 0;
  for (i = 0; i < n_moltypes; ++i)
  {
    top->contents->excls.nr += my_mols[i].n_mols * my_mols[i].n_excls; 
    top->contents->excls.nra += my_mols[i].n_mols * my_mols[i].n_exclsa;
  }
  top->contents->excls.index = (int *) ecalloc(top->contents->excls.nr + 1,sizeof(int));
  top->contents->excls.index[0] = 0;
  at_idx = 1;

// MRD 01.21.2020 
  for (i = 0; i < top->n_molblocks; ++i)
  {
    molidx = top->molblocks[i].moltype;
    for (j = 0; j < top->molblocks[i].n_mols; ++j)
    {
      for (k = 0; k < my_mols[molidx].n_apm; ++k)
      {
        top->contents->excls.index[at_idx] = top->contents->excls.index[at_idx-1] + my_mols[molidx].n_epa[k];
        ++at_idx;
      }
    }
  } 
  top->contents->excls.a = (int *) ecalloc(top->contents->excls.nra,sizeof(int));
  int ea_idx = 0;
  int prev_molblocks = 0;
  for (i = 0; i < top->n_molblocks; ++i)
  {
    molidx = top->molblocks[i].moltype;
    for (j = 0; j < top->molblocks[i].n_mols; ++j)
    {
      for (k = 0; k < my_mols[molidx].n_apm; ++k)
      {
        for (l = 0; l < my_mols[molidx].n_epa[k]; ++l)
        {
          top->contents->excls.a[ea_idx] = my_mols[molidx].excls[k][l] + j * my_mols[molidx].n_apm + prev_molblocks;
          ++ea_idx;
        }
      }
    }
    prev_molblocks += top->molblocks[i].n_mols * my_mols[molidx].n_apm;
  }

// MRD 01.21.2020 This bit is okay still
  top->contents->idef.il[F_BONDS].nr = 0;
  top->contents->idef.il[F_ANGLES].nr = 0;
  top->contents->idef.il[F_PDIHS].nr = 0;
  top->contents->idef.il[F_LJ14].nr = 0;
  for (i = 0; i < n_moltypes; ++i)
  {
    top->contents->idef.il[F_BONDS].nr += my_mols[i].n_mols * my_mols[i].bond_nr;
    top->contents->idef.il[F_ANGLES].nr += my_mols[i].n_mols * my_mols[i].angle_nr;
    top->contents->idef.il[F_PDIHS].nr += my_mols[i].n_mols * my_mols[i].dih_nr;
    top->contents->idef.il[F_LJ14].nr += my_mols[i].n_mols * my_mols[i].imnb_nr;
  } 
  top->contents->idef.il[F_BONDS].iatoms = (tW_t_iatom *) ecalloc(top->contents->idef.il[F_BONDS].nr,sizeof(tW_t_iatom));
  top->contents->idef.il[F_ANGLES].iatoms = (tW_t_iatom *) ecalloc(top->contents->idef.il[F_ANGLES].nr,sizeof(tW_t_iatom));
  top->contents->idef.il[F_PDIHS].iatoms = (tW_t_iatom *) ecalloc(top->contents->idef.il[F_PDIHS].nr,sizeof(tW_t_iatom));
  top->contents->idef.il[F_LJ14].iatoms = (tW_t_iatom *) ecalloc(top->contents->idef.il[F_LJ14].nr,sizeof(tW_t_iatom));


// MRD 01.21.2020 
  prev_molblocks = 0;
  int b_idx = 0, a_idx = 0, pd_idx = 0, imnb_idx = 0;
//  MRD 01.22.2020
  for (i = 0; i < top->n_molblocks; ++i)
  {
    molidx = top->molblocks[i].moltype;
    for (j = 0; j < top->molblocks[i].n_mols; ++j)
    {
      for (k = 0; k < my_mols[molidx].n_bonds; ++k)
      {
	top->contents->idef.il[F_BONDS].iatoms[b_idx] = my_mols[molidx].bond_types[k];
	++b_idx;
	top->contents->idef.il[F_BONDS].iatoms[b_idx] = my_mols[molidx].bond_ij[k][0] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++b_idx;
        top->contents->idef.il[F_BONDS].iatoms[b_idx] = my_mols[molidx].bond_ij[k][1] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++b_idx;
      }
      for (k = 0; k < my_mols[molidx].n_angles; ++k)
      {
        top->contents->idef.il[F_ANGLES].iatoms[a_idx] = my_mols[molidx].angle_types[k];
        ++a_idx;
        top->contents->idef.il[F_ANGLES].iatoms[a_idx] = my_mols[molidx].angle_ijk[k][0] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++a_idx;
        top->contents->idef.il[F_ANGLES].iatoms[a_idx] = my_mols[molidx].angle_ijk[k][1] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++a_idx;
        top->contents->idef.il[F_ANGLES].iatoms[a_idx] = my_mols[molidx].angle_ijk[k][2] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++a_idx;
      }
      for (k = 0; k < my_mols[molidx].n_dihs; ++k)
      {
	top->contents->idef.il[F_PDIHS].iatoms[pd_idx] = my_mols[molidx].dih_types[k];
	++pd_idx;
	top->contents->idef.il[F_PDIHS].iatoms[pd_idx] = my_mols[molidx].dih_ijkl[k][0] + j * my_mols[molidx].n_apm + prev_molblocks;
	++pd_idx;
        top->contents->idef.il[F_PDIHS].iatoms[pd_idx] = my_mols[molidx].dih_ijkl[k][1] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++pd_idx;
        top->contents->idef.il[F_PDIHS].iatoms[pd_idx] = my_mols[molidx].dih_ijkl[k][2] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++pd_idx;
        top->contents->idef.il[F_PDIHS].iatoms[pd_idx] = my_mols[molidx].dih_ijkl[k][3] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++pd_idx;
      }
      for (k = 0; k < my_mols[molidx].n_imnbs; ++k)
      {
        top->contents->idef.il[F_LJ14].iatoms[imnb_idx] = my_mols[molidx].imnb_types[k];
        ++imnb_idx;
        top->contents->idef.il[F_LJ14].iatoms[imnb_idx] = my_mols[molidx].imnb_ij[k][0] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++imnb_idx;
        top->contents->idef.il[F_LJ14].iatoms[imnb_idx] = my_mols[molidx].imnb_ij[k][1] + j * my_mols[molidx].n_apm + prev_molblocks;
        ++imnb_idx;
      }
    }
    prev_molblocks = top->molblocks[i].n_mols * my_mols[molidx].n_apm;
  }
  top->contents->idef.functype = (tW_t_functype *) ecalloc(top->contents->idef.atnr,sizeof(tW_t_functype));
  for (i = 0; i < top->contents->idef.atnr; ++i)
  {
    top->contents->idef.functype[i] = i;
  }
}

/*********************************************************************************
find_line_in_file():
@arg fp pointer to file we are searching
@arg target string (char *) for pattern we are searching for
@arg ret_line string (char *) for storing the line that has the pattern
Reads a file looking for target pattern, if found, stores line in ret_line and 
stops file pointer after line containing pattern
Pass NULL if don't care. 
@return true if found
@return false if not
*********************************************************************************/
bool find_line_in_file(FILE *fp, const char * target, char * ret_line)
{
  tW_line inp_line;
  bool rew_flag = FALSE;
  int l_inp;

  l_inp = get_next_line(fp,inp_line);
  while (strstr(inp_line,target) == NULL)
  {
    l_inp = get_next_line(fp,inp_line);
    if (l_inp == -1) // EOF
    {
      if (rew_flag)
      {
	if (DEBUG_BUILD)
	{
	fprintf(stderr,"WARNING: went through entire file, rewound it, and went through entire file again.\n");
        fprintf(stderr,"\tNever found target: %s\n",target);
	}
	return FALSE; 
	//exit(1);
      }
      else
      {
//	fprintf(stderr,"WARNING: Reached end of file and was unable to find target: %s\n",target);
//	fprintf(stderr,"\tRewinding file and trying again incase target was passed previously\n");
        rewind(fp);
	rew_flag = TRUE;
      }
    } 
  }
  if (ret_line != NULL) {strcpy(ret_line,inp_line);}
  return TRUE;
}

/*****************************************************************************************
read_tpr_dump(): Function to read a text file that contains a dumped tpr file
*****************************************************************************************/
bool read_tpr_dump(tW_word fnm, tW_gmx_topology *top)
{

  if (DEBUG_BUILD) {fprintf(stderr, "DEBUG BUILD IS SET TO %d\n", DEBUG_BUILD);}
  top->contents = (tW_t_topology *) ecalloc(1,sizeof(tW_t_topology));

  FILE *fp = open_file(fnm,'r'); //fopen(fnm,"r");
  char * ret_line = (char *) ecalloc(MAXLINELEN,sizeof(char));
  tW_line inp_line;
  tW_line *molblock_lines;
  int test_sscanf, inp_int, n_atoms, n_molblocks = 0;
  int i, j, k, l, let_idx;
  bool bOk;
  tW_word inp_word;

  // top->contents->name is the first line under topology:
  bOk = dprocess_sysname_natoms(fp, top, &inp_line);
  if (!bOk) { THROW_TPR_ERROR("COULD NOT PROCESS TOPOLOGY: SECTION OF TPR\n");}
  // read ffparams section
  bOk = dget_top_ff_info(fp, top, &inp_line);
  if (!bOk) { THROW_TPR_ERROR("UNABLE TO PARSE TOP FF INFO\n"); }
  // Figure out total different molecule sections (blocks) in file
  bOk = dprocess_nmolblocks(fp, top, &inp_line);
  if (!bOk) { THROW_TPR_ERROR("COULD NOT PROCESS MOLTYPES-MOLBLOCKS\n");}

  // Have n_molblocks now. So now know the memory req. for storing block info
  top->molblocks = (tW_molblock *) ecalloc(top->n_molblocks, sizeof(tW_molblock));
  if (DEBUG_BUILD) {fprintf(stderr, "got %d molblocks \n", top->n_molblocks);}

  // Search molblocks for distinct type info 
  // e.g. my_mols[i].molname, my_mols[i].n_mols, n_molecule_types
  bOk = dget_moltypes_in_molblocks(fp, top, &inp_line);
  if (!bOk) { THROW_TPR_ERROR("COULD NOT PROCESS MOLBLOCK INFO FOR MOLTYPES\n"); }

  //Now we know how many distinct types to collect info for
  top->molecules = (tW_molecule *) ecalloc(top->n_molecule_types, sizeof(tW_molecule));
  tW_molecule *my_mols = top->molecules;
  top->contents->mols.nr = top->n_molecule_types;

  // So we'll collect that info now into our mol struct list
  bOk = dfind_all_molinfo(fp, my_mols, top->n_molecule_types, &inp_line);

  top->contents->atoms.nres = 0;
  n_atoms = 0;
 
  int n_moltypes = 0;
  // MRD 01.17.2020 I still want to check and make sure we have the correct number of atoms.
  n_atoms = 0;
  for (i = 0; i < top->n_molblocks; ++i)
  {
    n_atoms += top->molblocks[i].n_mols * my_mols[top->molblocks[i].moltype].n_apm;
    my_mols[top->molblocks[i].moltype].n_mols += top->molblocks[i].n_mols;
  }


  if (n_atoms != top->contents->atoms.nr)
  {
    fprintf(stderr,"WARNING: topology is supposed to contain %d atoms\n",top->contents->atoms.nr);
    for (i = 0; i < top->n_molblocks; ++i)
    {
fprintf(stderr,"\tmolblock: %d   moltype: %d  n_mols: %d   n_at_per_mol: %d   = %d atoms\n",i,
                      top->molblocks[i].moltype, top->molblocks[i].n_mols, my_mols[top->molblocks[i].moltype].n_mols,
                      top->molblocks[i].n_mols * my_mols[top->molblocks[i].moltype].n_mols);
    }
    fprintf(stderr,"Totals %d atoms from molblocks\n",n_atoms);
  }


  fclose(fp);

  pop_contents(top); // This is just a function that transfers the hard word from get_moletype_info to the topology data structure

  return TRUE;
}

/*****************************************************************************************
read_box_dump_frame(): This function reads the box information from a dumped trr file
*****************************************************************************************/
tW_gmx_bool read_box_dump_frame(FILE *fp, tW_matrix box)
{
  tW_line inp_line;
  int i, test_sscanf, alpha;
  float b1, b2, b3;
  for (i = 0; i < 3; ++i)
  {
    get_next_line(fp,inp_line); // box[   a]={ X.XXXXe+XX,  X.XXXXe+XX,   X.XXXXe+XX}
    test_sscanf = sscanf(inp_line," box[    %d]={ %f, %f, %f} ",&alpha,&b1,&b2,&b3);
    box[i][0] = b1;
    box[i][1] = b2;
    box[i][2] = b3;
  }
  return TRUE;
}

/*****************************************************************************************
read_first_dump_frame(): Reads the first frame of a dumped file. Gets "header" information
and finds out whether the positions/velocities/forces are present. 
*****************************************************************************************/
int read_first_dump_frame(tW_gmx_trxframe *frame, const char *trx_fnm)
{
  int n_atoms, step, test_sscanf;
  double time, lambda, d1, d2, d3;
  float f1, f2, f3;
  char prop;
  int i,j;
  tW_word inp_word;
  tW_line inp_line;

  frame->fp = open_file(trx_fnm,'r');
  
  get_next_line(frame->fp, inp_line); // cg.trr frame X:
  frame->contents->title = (char *) ecalloc(50,sizeof(char));
  test_sscanf = sscanf(inp_line," %s ",&inp_word);
  strcpy(frame->contents->title,inp_word);
  frame->contents->bTitle = TRUE;


  get_next_line(frame->fp, inp_line); // natoms=   %d  step=   %d  time=   %g  lambda= X
  test_sscanf = sscanf(inp_line," natoms= %d step= %d time= %g lambda= %g ",&n_atoms, &step, &time, &lambda);


  if (test_sscanf != 4) 
  { 
    fprintf(stderr,"ERROR: unable to read expected 4 arguments from line\n");
    fprintf(stderr,"\tnatoms, step, time, lambda\n");
    fprintf(stderr,"%s",inp_line);
    exit(1);
  }
  else
  {
    frame->contents->bStep = TRUE;
    frame->contents->bTime = TRUE;
    frame->contents->bLambda = TRUE;

    if (strstr(inp_line,"step") == NULL)
    { 
      fprintf(stderr,"ERROR: read 4 values from line %s",inp_line);
      fprintf(stderr,"\t but did not find \"step\"\n");
      frame->contents->bStep = FALSE;
    }
    if (strstr(inp_line,"time") == NULL)
    {
      fprintf(stderr,"ERROR: read 4 values from line %s",inp_line);
      fprintf(stderr,"\t but did not find \"time\"\n");
      frame->contents->bTime = FALSE;
    }
    if (strstr(inp_line,"lambda") == NULL)
    {
      fprintf(stderr,"ERROR: read 4 values from line %s",inp_line);
      fprintf(stderr,"\t but did not find \"lambda\"\n");
      frame->contents->bLambda = FALSE;
    }
    frame->contents->t0 = time;  
    frame->contents->tf = time;
    frame->contents->tpf = (tW_real) 0.0;
    frame->contents->tppf = (tW_real) 0.0;
  
    frame->contents->step = step;
    frame->contents->time = time;
    frame->contents->lambda = lambda;
  }
 
  get_next_line(frame->fp, inp_line); //box (3x3):

  frame->contents->bBox = read_box_dump_frame(frame->fp, frame->contents->box);

  get_next_line(frame->fp, inp_line); // x (n_atomsx3):
  test_sscanf = sscanf(inp_line," %c ", &prop);
  if (prop == 'x')
  {
    frame->contents->bX = TRUE;
  }
  else
  {
    fprintf(stderr,"ERROR: expected line \"     x (%dx3): \" in dumped trr file! \n",n_atoms);
    fprintf(stderr,"\tHOWEVER, we found line: %s\n",inp_line);
    fprintf(stderr,"\tWithout the positions, we can't do anything.\n");
    exit(1);
  }

  get_next_line(frame->fp, inp_line); // stealing a frame. My idea is to test if these numbers are doubles or floats

  test_sscanf = sscanf(inp_line," x[ %d]={ %g, %g, %g} ",&i,&d1,&d2,&d3);
  test_sscanf = sscanf(inp_line," x[ %d]={ %f, %f, %f} ",&i,&f1,&f2,&f3);

  if ((double) f1 == d1) { frame->contents->bDouble = TRUE; }
  else { frame->contents->bDouble = FALSE; }
  
  for (i = 1; i < n_atoms; ++i) { get_next_line(frame->fp, inp_line); }

  get_next_line(frame->fp, inp_line);
  test_sscanf = sscanf(inp_line," %c ", &prop);
  if (prop == 'v')
  {
    frame->contents->bV = TRUE;
    for (i = 0; i < n_atoms; ++i) { get_next_line(frame->fp, inp_line); }
    get_next_line(frame->fp, inp_line);
    test_sscanf = sscanf(inp_line," %c ", &prop);
  }

  if (prop == 'f') { frame->contents->bF = TRUE; }

  set_natoms(frame, n_atoms);

  rewind(frame->fp);

  return n_atoms;
}

/*****************************************************************************************
read_first_bocs_frame(): reads the first frame of a bocs.traj file. gets "header"
information
*****************************************************************************************/
int read_first_bocs_frame(tW_gmx_trxframe *frame, const char *trx_fnm)
{
  int natoms, test_sscanf, inp_int, i;
  int bDbl, bTtl, bStp, bTim, bAtm, bPrc, bBox, bX, bV, bF, bLD, bLDG;
  tW_line inp_line;
  tW_word inp_word;

  frame->fp = open_file(trx_fnm,'r');
  
  get_next_line(frame->fp,inp_line);
  get_next_line(frame->fp,inp_line);
  test_sscanf = sscanf(inp_line,"[natoms] %d ",&natoms);
  get_next_line(frame->fp,inp_line);
  test_line(inp_line,"[flags]",TRUE,"expected directive [flags] after [natoms]\n");
  get_next_line(frame->fp,inp_line);
  
  test_sscanf = sscanf(inp_line," %d %d %d %d %d %d %d %d %d %d %d %d ",&bDbl, &bTtl, &bStp, &bTim, &bAtm, &bPrc, &bBox, &bX, &bV, &bF, &bLD, &bLDG);
  if ((test_sscanf < 10) || (test_sscanf > 12)) 
  { 
    fprintf(stderr,"ERROR: read %d out of expected 10-12 arguments from line: %s\n",test_sscanf,inp_line);
    fprintf(stderr,"\tExpected: bDouble bTotal bStep bTime bAtom bPrecision bBox bX bV bF\n");
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
  }
  
  frame->contents->bDouble = (bDbl == 1 ? TRUE : FALSE);
  frame->contents->bTitle = (bTtl == 1 ? TRUE : FALSE);
  frame->contents->bStep = (bStp == 1 ? TRUE : FALSE);
  frame->contents->bTime = (bTim == 1 ? TRUE : FALSE);
  frame->contents->bAtoms = (bAtm == 1 ? TRUE : FALSE);
  frame->contents->bPrec = (bPrc == 1 ? TRUE : FALSE);
  frame->contents->bBox = (bBox == 1 ? TRUE : FALSE);
  frame->contents->bX = (bX == 1 ? TRUE : FALSE);
  frame->contents->bV = (bV == 1 ? TRUE : FALSE);
  frame->contents->bF = (bF == 1 ? TRUE : FALSE);

  if ((test_sscanf == 11) || (test_sscanf == 12)) { frame->bLDs = (bLD == 1 ? TRUE : FALSE); }
  if (test_sscanf == 12) { frame->bLD_Grads = (bLDG == 1 ? TRUE : FALSE); }

  if (frame->contents->bBox)
  {
    while (strstr(inp_line,"[box]") == NULL) { get_next_line(frame->fp,inp_line); }
    tW_word b0, b1, b2;
    for (i = 0; i < 3; ++i)
    {
      get_next_line(frame->fp, inp_line);
      test_sscanf = sscanf(inp_line," %s %s %s ",&b0,&b1,&b2);
      if (test_sscanf != 3)
      {
        fprintf(stderr,"ERROR reading box line %d\n",i+1);
        fprintf(stderr,"\tline: %s",inp_line);
        fprintf(stderr,"\ttest_sscanf: %d \n",test_sscanf);
        fprintf(stderr,"\t\' %s \'  \'  %s \'  \' %s \' \n",b0,b1,b2);
        exit(1);
      }
      frame->contents->box[i][0] = atof(b0);
      frame->contents->box[i][1] = atof(b1);
      frame->contents->box[i][2] = atof(b2);
    }
  }

  if (frame->bLDs)
  {
    while (strstr(inp_line,"[local_densities]") == NULL) { get_next_line(frame->fp,inp_line); }
    test_sscanf = sscanf(inp_line,"[local_densities] %d ",&inp_int);
    frame->n_atomtypes = inp_int;
  }


  set_natoms(frame,natoms);

  rewind(frame->fp);

  return natoms;
}

/*****************************************************************************************
copy_atom_props(): Reads "natoms" lines from a dump file and stores the numbers in 
the appropriate array
*****************************************************************************************/
void copy_atom_props(FILE *fp, tW_rvec *p, int natoms)
{
  int i, test_sscanf, idx;
  char prop;
  double p0, p1, p2;
  float f0, f1, f2;
  tW_line inp_line;
  
  for (i = 0; i < natoms; ++i)
  {
    get_next_line(fp,inp_line); //   p[  X]={ XXXXX,  XXXXXX,  XXXXXX}
    test_sscanf = sscanf(inp_line," %c[ %d]={ %f, %f, %f} ",&prop,&idx,&(f0),&(f1),&(f2));
    p[i][0] = (double) f0;
    p[i][1] = (double) f1;
    p[i][2] = (double) f2;
    if (test_sscanf != 5)
    {
      fprintf(stderr,"ERROR: unable to read 5 values from line: %s",inp_line);
      exit(1);
    }
  }

}

/*****************************************************************************************
get_prop(): makes sure the dumped trajectory file has the correct property next. If so,
calls a function to read it in
*****************************************************************************************/
void get_prop(FILE *fp, tW_gmx_bool present, tW_rvec *p, char test_prop, int natoms, int step)
{
  char prop;
  tW_line inp_line;
  int test_sscanf;
  if (present)
  {
    get_next_line(fp, inp_line);
    test_sscanf = sscanf(inp_line," %c ",&prop);
    if (prop != test_prop)
    {
      fprintf(stderr,"ERROR: missing property from frame %d!\n", step);
      fprintf(stderr,"\tshould have: %c (%dx3):\n",test_prop,natoms);
      fprintf(stderr,"\tinstead we have: %s",inp_line);
      exit(1);
    }
    copy_atom_props(fp, p, natoms);
  }

}

/*****************************************************************************************
read_next_dump_frame(): self explanatory?
*****************************************************************************************/
bool read_next_dump_frame(tW_gmx_trxframe *frame)
{
  int n_atoms, step, test_sscanf, linp;
  int lambda;
  float time, flambda;
  double dtime;
  tW_line inp_line;
  tW_word inp_word;
  int i;
  char prop;

  linp = get_next_line(frame->fp, inp_line); // cg.trr frame X:
  if (linp == -1)
  {
    fprintf(stderr,"Done after frame %d\n",frame->counter - 1);
    return FALSE;
  }
  get_next_line(frame->fp, inp_line); // natoms=   %d  step=   %d  time=  %g   lambda= X
  test_sscanf = sscanf(inp_line," natoms= %d step= %d %s lambda= %d ",&n_atoms, &step, &inp_word, &lambda);
  i = 0;
  while (inp_word[i+5] != '\0')
  {
    inp_word[i] = inp_word[i+5];
    ++i;
  }
  inp_word[i] = '\0';
  dtime = atof(inp_word);
  if (test_sscanf != 4) 
  { 
    fprintf(stderr,"test_sscanf: %d != 4\n",test_sscanf);
    fprintf(stderr,"WARNING: unable to read expected 4 arguments from line\n");
    fprintf(stderr,"\tnatoms, step, time, lambda\n");
    fprintf(stderr,"%s",inp_line);
    exit(1);
  }
  if (frame->contents->natoms != n_atoms)
  {
    fprintf(stderr,"ERROR: in frame at step: %d    time: %9.7e   we have %d atoms\n",step,dtime,n_atoms);
    fprintf(stderr,"\tHOWEVER, the first frame had %d atoms!\n",frame->contents->natoms);
    exit(1);
  }
  
  frame->contents->tppf = frame->contents->tpf;
  frame->contents->tpf = frame->contents->tf;
  frame->contents->tf = dtime;

  frame->contents->time = (tW_real) dtime;
  frame->contents->step = step;
  frame->contents->lambda = (tW_real)lambda;


  get_next_line(frame->fp, inp_line); //  box (3x3):
  frame->contents->bBox = read_box_dump_frame(frame->fp, frame->contents->box);

  get_prop(frame->fp, frame->contents->bX, frame->contents->x, 'x', frame->contents->natoms, step);
  get_prop(frame->fp, frame->contents->bV, frame->contents->v, 'v', frame->contents->natoms, step);
  get_prop(frame->fp, frame->contents->bF, frame->contents->f, 'f', frame->contents->natoms, step);
  
  return TRUE;
}

void write_dump_frame(tW_gmx_trxframe *fr)
{
  int i, j;
  fprintf(fr->fp,"%s frame %d: \n",fr->contents->title,fr->contents->step);
  fprintf(fr->fp,"   natoms= %9d  step= %9d  time=%9.7e  lambda= %9d\n",fr->contents->natoms, fr->contents->step, fr->contents->time, fr->contents->lambda);
  fprintf(fr->fp,"   box (3x3):\n");
  for (i = 0; i < 3; ++i)
  {
    fprintf(fr->fp,"      box[ %4d]={% 8.5e, % 9.5e, % 9.5e}\n",i,fr->contents->box[i][0],fr->contents->box[i][1],fr->contents->box[i][2]);
  }
  if (fr->contents->bX)
  {
    fprintf(fr->fp,"   x (%dx3): \n",fr->contents->natoms);
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      fprintf(fr->fp,"      x[%5d]={% 8.5e, % 9.5e, % 9.5e}\n",i,fr->contents->x[i][0],fr->contents->x[i][1],fr->contents->x[i][2]);
    }
  }
  if (fr->contents->bV)
  {
    fprintf(fr->fp,"   v (%dx3): \n",fr->contents->natoms);
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      fprintf(fr->fp,"      v[%5d]={% 8.5e, % 9.5e, % 9.5e}\n",i,fr->contents->v[i][0],fr->contents->v[i][1],fr->contents->v[i][2]);
    }
  }
  if (fr->contents->bF)
  {
    fprintf(fr->fp,"   f (%dx3): \n",fr->contents->natoms);
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      fprintf(fr->fp,"      f[%5d]={% 8.5e, % 9.5e, % 9.5e}\n",i,fr->contents->f[i][0],fr->contents->f[i][1],fr->contents->f[i][2]);
    }
  }

}

/*****************************************************************************************
read_first_gro_frame(): reads the first 3 lines of a .gro file, then rewinds it.
	does nothing with the first line. gets natoms from the second line.
	sets the "title" of the trajectory to the name of the gro file.
	sets bStep to true, and when reading the gro file, takes the value of
	fr->counter to be the step. from third line, determines if velocities are present	
*****************************************************************************************/
int read_first_gro_frame(tW_gmx_trxframe *frame, const char *fnm)
{
  int test_sscanf, natoms;
  tW_line inp_line;

  frame->fp = open_file(fnm,'r');
  
  get_next_line(frame->fp,inp_line);
  get_next_line(frame->fp,inp_line);
  test_sscanf = sscanf(inp_line," %d ",&natoms);

  frame->contents->bDouble = FALSE;
  frame->contents->bStep = TRUE;
  frame->contents->bTime = FALSE;
  frame->contents->bLambda = FALSE;
  frame->contents->bTitle = TRUE;
  frame->contents->title = (char *) ecalloc(50,sizeof(char));
  strcpy(frame->contents->title,fnm); 

  get_next_line(frame->fp,inp_line);
  int n_mol, n_atm;
  tW_word attype, moltype;
  float x0, x1, x2, v0, v1, v2;
  test_sscanf = sscanf(inp_line," %d%s %s %d %f %f %f %f %f %f ",&n_mol,&moltype,&attype,&n_atm,&x0,&x1,&x2,&v0,&v1,&v2);
 
  frame->contents->bBox = TRUE;
  frame->contents->bX = TRUE;
  if (test_sscanf == 10) { frame->contents->bV = TRUE; }
  else if (test_sscanf == 7) { frame->contents->bV = FALSE; }
  else
  {
    fprintf(stderr,"ERROR: in third line of gro file %s, I did not read 7 or 10 values\n",fnm);
    fprintf(stderr,"\tExpected #mol moltype attype #at x0 x1 x2 (v0 v1 v2)\n");
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    exit(1);
  }
  frame->contents->bF = FALSE;

  set_natoms(frame,natoms);

  rewind(frame->fp);
  return natoms;
}

/*
read_next_gro_frame(): reads the next frame from a gro file and puts x (and v) 
	vectors into fr->contents->x (v)
*/

bool read_next_gro_frame(tW_gmx_trxframe *fr)
{
  int step, test_sscanf, molnr, atnr;
  int i,j;
  float x0, x1, x2, v0, v1, v2;
  tW_word inp_word, moltype, attype;
  tW_line inp_line;

  if (get_next_line(fr->fp,inp_line) == -1) { fprintf(stderr,"Done after frame %d\n",fr->counter - 1); return FALSE; }
  fr->contents->step = fr->counter;


  get_next_line(fr->fp,inp_line);
  
  for (i = 0; i < fr->contents->natoms; ++i)
  {
    get_next_line(fr->fp,inp_line);
    if (i < 9999)
    {
      test_sscanf = sscanf(inp_line," %d%s %s %d %f %f %f %f %f %f ",&molnr,&moltype,&attype,&atnr,&x0,&x1,&x2,&v0,&v1,&v2);
      if (test_sscanf == 10 && fr->contents->bV)
      {
        fr->contents->x[i][0] = (double) x0;
        fr->contents->x[i][1] = (double) x1;
        fr->contents->x[i][2] = (double) x2;
        fr->contents->v[i][0] = (double) v0;
        fr->contents->v[i][1] = (double) v1;
        fr->contents->v[i][2] = (double) v2;
      }
      else if (test_sscanf == 7 && (!fr->contents->bV))
      {
        fr->contents->x[i][0] = (double) x0;
        fr->contents->x[i][1] = (double) x1;
        fr->contents->x[i][2] = (double) x2;
      }
      else if (test_sscanf == 10 && (!fr->contents->bV))
      {
	fprintf(stderr,"WARNING: read 10 items from line, but only expected 7.\n");
	fprintf(stderr,"\tline: %s",inp_line);
	fprintf(stderr,"\tMemory has not been allocated for velocities, so I will ignore them.\n");
	fprintf(stderr,"\tThe velocity flag is set based on the third line in the .gro file.\n");
        fr->contents->x[i][0] = (double) x0;
        fr->contents->x[i][1] = (double) x1;
        fr->contents->x[i][2] = (double) x2;
      }
      else if (test_sscanf == 7 && (fr->contents->bV))
      {
        fprintf(stderr,"WARNING: read 7 items from line, but expected 10.\n");
        fprintf(stderr,"\tline: %s",inp_line);
        fprintf(stderr,"\tThe velocity flag is set based on the third line in the .gro file.\n");
	fprintf(stderr,"\tTurning bV off.\n");
	fr->contents->bV = FALSE;
      }
      else
      {
        fprintf(stderr,"ERROR: unable to read expected 7 (10) arguments from line\n");
        fprintf(stderr,"\tresnr, resid, atid, atnr, x, y, z, (vx, vy, vz)\n");
        fprintf(stderr,"\t%s",inp_line);
        exit(1);
      }
    }
    else
    {
      char *ptr;
      ptr = &(inp_line[15]); // This SHOULD BE where the atom number starts.
      test_sscanf = sscanf(ptr,"%d %f %f %f %f %f %f ",&atnr, &x0, &x1, &x2, &v0, &v1, &v2);
      if (test_sscanf == 7 && fr->contents->bV)
      {
        fr->contents->x[i][0] = (double) x0;
        fr->contents->x[i][1] = (double) x1;
        fr->contents->x[i][2] = (double) x2;
        fr->contents->v[i][0] = (double) v0;
        fr->contents->v[i][1] = (double) v1;
        fr->contents->v[i][2] = (double) v2;
      }		
      else if (test_sscanf == 4 && (!fr->contents->bV))
      {
        fr->contents->x[i][0] = (double) x0;
        fr->contents->x[i][1] = (double) x1;
        fr->contents->x[i][2] = (double) x2;
      }
      else if (test_sscanf == 7 && (!fr->contents->bV))
      {
        fprintf(stderr,"WARNING: read positions and velocities from line, but only expected positions.\n");
        fprintf(stderr,"\tline: %s",inp_line);
        fprintf(stderr,"\tMemory has not been allocated for velocities, so I will ignore them.\n");
        fprintf(stderr,"\tThe velocity flag is set based on the third line in the .gro file.\n");
        fr->contents->x[i][0] = (double) x0;
        fr->contents->x[i][1] = (double) x1;
        fr->contents->x[i][2] = (double) x2;
      }
      else if (test_sscanf == 4 && fr->contents->bV)
      {
        fprintf(stderr,"WARNING: read positions from line, but expected velocities too.\n");
        fprintf(stderr,"\tline: %s",inp_line);
        fprintf(stderr,"\tThe velocity flag is set based on the third line in the .gro file.\n");
        fprintf(stderr,"\tTurning bV off.\n");
        fr->contents->bV = FALSE;
      }
      else
      {
	fprintf(stderr,"ERROR: unable to read expected 7 (4) arguments from end of line\n");
	fprintf(stderr,"\tatnr, x, y, z, (vx, vy, vz)\n");
	fprintf(stderr,"\t%s",ptr);
	exit(1);
      }
    }
  } 
  get_next_line(fr->fp,inp_line);
  test_sscanf = sscanf(inp_line," %f %f %f %f ",&x0, &x1, &x2, &v0);
  if (test_sscanf == 4) { fprintf(stderr,"ERROR: read 4 numbers for box. that means the box is not rectangular. \n"); }
 
  for (i = 0; i < DIM; ++i) { for (j = 0; j < DIM; ++j) { fr->contents->box[i][j] = (double) 0.0; }}
  fr->contents->box[0][0] = (double) x0;
  fr->contents->box[1][1] = (double) x1;
  fr->contents->box[2][2] = (double) x2;

  return TRUE;
}

/*
write_gro_frame(): This function prints the current frame in a .gro format. Plain and simple.
*/
void write_gro_frame(tW_gmx_trxframe *fr)
{
  int i, j;
  fprintf(fr->fp,"%s frame %d: \n", fr->contents->title, fr->contents->step);
  fprintf(fr->fp," %d\n",fr->contents->natoms);

  for (i = 0; i < fr->contents->natoms; ++i)
  {
    if (fr->contents->bX)
    {
      if (fr->contents->bV)
      {
	fprintf(fr->fp,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
		fr->contents->atoms->resinfo[fr->contents->atoms->atom[i].resind].nr,
		*(fr->contents->atoms->resinfo[fr->contents->atoms->atom[i].resind].name),
		*(fr->contents->atoms->atomname[i]),
		i+1,
		fr->contents->x[i][0],
		fr->contents->x[i][1],
		fr->contents->x[i][2],
		fr->contents->v[i][0],
		fr->contents->v[i][1],
		fr->contents->v[i][2]);
      }
      else
      {
        fprintf(fr->fp,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
                fr->contents->atoms->resinfo[fr->contents->atoms->atom[i].resind].nr,
                *(fr->contents->atoms->resinfo[fr->contents->atoms->atom[i].resind].name),
                *(fr->contents->atoms->atomname[i]),
                i+1,
                fr->contents->x[i][0],
                fr->contents->x[i][1],
                fr->contents->x[i][2]);
      }
    }
    else
    {
      fprintf(stderr,"ERROR: trying to print a frame in .gro format, but fr->contents->bX is false!\n");
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    }
  }
  fprintf(fr->fp,"%9.5f %9.5f %9.5f\n",fr->contents->box[0][0],fr->contents->box[1][1],fr->contents->box[2][2]);

}

/*
write_bocs_top(): The (not too awful) function to save the information that I've found is important
	in top in the bocstop format. If anything else ever becomes important in the topology,
	modify this function to write it as well. Likely, you'll add some [directive] at the end
	and then after that print the stuff. After that, be sure to modify read_bocs_top to check for
	the [directive] at the proper time in the file reading process.
*/

void write_bocs_top(tW_word fnm, tW_gmx_topology *top)
{
  FILE *fp = open_file(fnm,'w');
  int i, j, k;

//  fprintf(fp,"[name] %s\n\n",*(top->contents->name));
  fprintf(fp,"[name] %s\n\n",top->contents->name);
  fprintf(fp,"[natoms] %d\n\n",top->contents->atoms.nr);
  fprintf(fp,"[ffparams] %d\n",top->contents->idef.ntypes);
  fprintf(fp,"!ind \tname\n");
  for (i = 0; i < top->contents->idef.ntypes; ++i)
  {
    fprintf(fp,"%d \t%s\n",i,top->force_names[i]);
  }
  fprintf(fp,"!atnr  fudgeQQ\n");
  fprintf(fp,"%d %f\n\n",top->contents->idef.atnr,top->contents->idef.fudgeQQ);

// MRD 01.24.2020  
  fprintf(fp,"[molblocks] %d\n",top->n_molblocks);
  fprintf(fp,"!moltype  mol_name  n_mols\n");
  for (i = 0; i < top->n_molblocks; ++i)
  {
    fprintf(fp,"%d \t  %s \t  %d \n",top->molblocks[i].moltype,top->molblocks[i].molname,top->molblocks[i].n_mols);
  }
  fprintf(fp,"\n");

  fprintf(fp,"[moltypes] %d\n",top->contents->mols.nr);
  for (i = 0; i < top->contents->mols.nr; ++i)
  {
    fprintf(fp,"!ind \tname\tn_mols\tn_apm\tn_res\tn_cgs\texcls_nr \texcls_nra \tbond_nr \tangle_nr \tdih_nr   \tbigdih_nr  \timnb_nr\n");
    fprintf(fp,"%d \t%s\t%d\t%d\t%d\t%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\n",i,
		top->molecules[i].molname,
		top->molecules[i].n_mols,
		top->molecules[i].n_apm,
		top->molecules[i].n_res,
		top->molecules[i].n_cg,
		top->molecules[i].n_excls,
		top->molecules[i].n_exclsa,
		top->molecules[i].bond_nr,
		top->molecules[i].angle_nr,
		top->molecules[i].dih_nr,
		top->molecules[i].bigdih_nr,
                top->molecules[i].imnb_nr);
    if (top->molecules[i].n_res > 0)
    {
      fprintf(fp,"  [res]\n");
      for (j = 0; j < top->molecules[i].n_res; ++j)
      {
	fprintf(fp,"  %s\n",top->molecules[i].resname[j]);
      }
    }
    if (top->molecules[i].n_cg > 0)
    {
      fprintf(fp,"  [cg]\n  !%3s  %5s  %3s\n","idx","start","end");
      for (j = 0; j < top->molecules[i].n_cg; ++j)
      {
	fprintf(fp,"  %4d  %5d  %3d\n",j,top->molecules[i].cg_start[j],top->molecules[i].cg_end[j]);
      }
    }
    if (top->molecules[i].bond_nr > 0)
    {
      fprintf(fp,"  [bond]\n  !%3s  %4s  %3s  %3s\n","idx","type","i","j");
      for(j = 0; j < top->molecules[i].n_bonds; ++j)
      {
	fprintf(fp,"  %4d  %4d  %3d  %3d\n",j,top->molecules[i].bond_types[j],top->molecules[i].bond_ij[j][0],top->molecules[i].bond_ij[j][1]);
      }
    }
    if (top->molecules[i].angle_nr > 0)
    {
      fprintf(fp,"  [angle]\n  !%3s  %4s  %3s  %3s  %3s\n","idx","type","i","j","k");
      for (j = 0; j < top->molecules[i].n_angles; ++j)
      {
	fprintf(fp,"  %4d  %4d  %3d  %3d  %3d\n",j,top->molecules[i].angle_types[j],top->molecules[i].angle_ijk[j][0],top->molecules[i].angle_ijk[j][1],top->molecules[i].angle_ijk[j][2]);
      }
    }
    if (top->molecules[i].dih_nr > 0)
    {
      fprintf(fp,"  [dih]\n  !%3s  %4s  %3s  %3s  %3s  %3s\n","idx","type","i","j","k","l");
      for (j = 0; j < top->molecules[i].n_dihs; ++j)
      {
	fprintf(fp,"  %4d  %4d  %3d  %3d  %3d  %3d\n",j,top->molecules[i].dih_types[j],top->molecules[i].dih_ijkl[j][0],top->molecules[i].dih_ijkl[j][1],top->molecules[i].dih_ijkl[j][2],top->molecules[i].dih_ijkl[j][3]);
      }
    }
    if (top->molecules[i].bigdih_nr > 0)
    {
      fprintf(fp,"  [bigdih]\n  !%3s  %4s  %3s  %3s  %3s  %3s  %3s\n","idx","type","i","j","k","l","m");
      for (j = 0; j < top->molecules[i].n_bigdihs; ++j)
      {
        fprintf(fp,"  %4d  %4d  %3d  %3d  %3d  %3d  %3d\n",j,top->molecules[i].bigdih_types[j],top->molecules[i].bigdih_ijklm[j][0],top->molecules[i].bigdih_ijklm[j][1],top->molecules[i].bigdih_ijklm[j][2],top->molecules[i].bigdih_ijklm[j][3],top->molecules[i].bigdih_ijklm[j][4]);
      }
    }
    if (top->molecules[i].imnb_nr > 0)
    {
      fprintf(fp,"  [imnb]\n  !%3s  %4s  %3s  %3s\n","idx","type","i","j");
      for (j = 0; j < top->molecules[i].n_imnbs; ++j)
      {
        fprintf(fp,"  %4d  %4d  %3d  %3d\n",j,top->molecules[i].imnb_types[j],top->molecules[i].imnb_ij[j][0],top->molecules[i].imnb_ij[j][1]);
      }
    }
  }

  fprintf(fp,"\n[attypes]\n");
  for (i = 0; i < top->contents->mols.nr; ++i)
  {
    for (j = 0; j < top->molecules[i].n_apm; ++j)
    {
      fprintf(fp,"!moltype_ind  at_ind  type  m          q        resind  atomname  typename  n_epa\n");
      fprintf(fp,"%12d  %6d  %4d  %9.5f  %7.5f  %6d  %8s  %9s  %d\n",i,j,
			top->molecules[i].type_ids[j],
			top->molecules[i].m[j],
			top->molecules[i].q[j],
			top->molecules[i].residx[j],
			top->molecules[i].atom_names[j],
			top->molecules[i].atom_types[j],
			top->molecules[i].n_epa[j]);
      if (top->molecules[i].n_epa[j] > 0)
      {
	fprintf(fp,"  [epa]\n");
	for (k = 0; k < top->molecules[i].n_epa[j]; ++k)
	{
//	  fprintf(fp,"  %d\n",top->molecules[i].excls[j][k]); ONE LINE EPA
	  fprintf(fp,"  %d",top->molecules[i].excls[j][k]);
	}
        fprintf(fp,"\n");
      }
    }
  }
  
  fprintf(fp,"\n\n");
  fprintf(fp,"[Bocs_at_types] %d\n!idx  typename\n",top->n_atomtypes);
  for (i = 0; i < top->n_atomtypes; ++i)
  {
    fprintf(fp," %3d  %8s\n",i,top->atom_type_names[i]);
  }
 
  fclose(fp);
}

/*
read_bocs_top(): The function that reads a .btp file, and populates top accordingly.
	If write_bocs_top ever gets changed, make sure this function is changed properly
	to look for those changes.
*/

bool read_bocs_top(tW_word fnm, tW_gmx_topology * top)
{
  int test_sscanf, i, j, k, inp_int, i1, i2, i3, i4, i5, inp_int2;
  int n_mt, ind, n_mol, n_apm, n_res, n_cgs, excls_nr, excls_nra, bond_nr, angle_nr, dih_nr, bigdih_nr, imnb_nr, n_epa;
  float inp_flt1, inp_flt2;
  tW_word inp_word, name;
  tW_line inp_line;

  top->contents = (tW_t_topology *) ecalloc(1,sizeof(tW_t_topology));

  FILE *fp = open_file(fnm,'r');
//  top = init_tW_gmx_topology();  

/*
~~~~~~~~~~~~~~~~~~ [name] ~~~~~~~~~~~~~~~~~~~
*/
  get_next_line(fp,inp_line);
  test_sscanf = sscanf(inp_line,"[name] %s ",&inp_word);
  if (test_sscanf != 1) 
  { 
    fprintf(stderr,"ERROR: Expected to find \'[name]\' directive on first line of .btp file\n");
    fprintf(stderr,"\tAre you sure this file is a bocstop file?\n");
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
    exit(1); 
  }

  strcpy(top->contents->name,inp_word);  

/*
~~~~~~~~~~~~~~~~~~ [natoms] ~~~~~~~~~~~~~~~~~~
*/
  get_next_line(fp,inp_line);
  test_sscanf = sscanf(inp_line,"[natoms] %d ",&inp_int);
  if (test_sscanf != 1) 
  { 
    fprintf(stderr,"ERROR: Expected to find \'[natoms]\' directive after \'[name]\' directive in .btp file\n");
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
    exit(1); 
  }
  top->contents->atoms.nr = inp_int;

/*
~~~~~~~~~~~~~~~~~~ [ffparams] ~~~~~~~~~~~~~~~~~~
*/
  get_next_line(fp,inp_line);
  test_sscanf = sscanf(inp_line,"[ffparams] %d ",&inp_int);
  if (test_sscanf != 1) 
  { 
    fprintf(stderr,"ERROR: Expected to find \'[ffparams]\' directive after \'[natoms]\' directive in .btp file\n");
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
    exit(1); 
  }
  top->contents->idef.ntypes = inp_int;

  top->force_names = (tW_word *) ecalloc(inp_int,sizeof(tW_word));

  for (i = 0; i < top->contents->idef.ntypes; ++i)
  {
    get_next_line(fp,inp_line);
    inp_int = -1;
    strcpy(inp_word,"N/A");
    test_sscanf = sscanf(inp_line,"%d %s ",&inp_int,&inp_word);
    if (test_sscanf != 2) 
    { 
      fprintf(stderr,"ERROR: Expected index of %d th ffparam and name of that force\n",i+1);
      fprintf(stderr,"\tline: %s",inp_line);
      fprintf(stderr,"\t idx: %d\n",inp_int);
      fprintf(stderr,"\tname: %s\n",inp_word);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
      exit(1); 
    }
    strcpy(top->force_names[i],inp_word);
  }

  get_next_line(fp,inp_line);
  test_sscanf = sscanf(inp_line,"%d %f ",&i1,&inp_flt1);
  if (test_sscanf != 2) 
  { 
    fprintf(stderr,"ERROR: Expected to find atnr and fudgeQQ on\n");
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"\tatnr: %d\n",i1);
    fprintf(stderr,"\tfudgeQQ: %f\n",inp_flt1);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
    exit(1); 
  }
  top->contents->idef.atnr = i1;
  top->contents->atomtypes.nr = i1;
  top->contents->idef.fudgeQQ = (double) inp_flt1;

/* MRD 01.24.2020
~~~~~~~~~~~~~~~~~~ [molblocks] ~~~~~~~~~~~~~~~~~~
*/

  get_next_line(fp,inp_line);
  inp_int = -1;
  test_sscanf = sscanf(inp_line,"[molblocks] %d ",&(inp_int));
  if (test_sscanf != 1)
  {
    fprintf(stderr,"ERROR: expected to find directive \"[molblocks]\" followed by n_molblocks\n");
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"\tn_molblocks: %d\n",inp_int);
    fprintf(stderr,"ERROR: %s %d \n",__FILE__,__LINE__);
    exit(1);
  }
  top->n_molblocks = inp_int;
  top->molblocks = (tW_molblock *) ecalloc(top->n_molblocks, sizeof(tW_molblock));
  for (i = 0; i < top->n_molblocks; ++i)
  {
    get_next_line(fp,inp_line);
    test_sscanf = sscanf(inp_line,"%d %s %d ",&(inp_int),&(inp_word),&(inp_int2));
    if (test_sscanf != 3)
    {
      fprintf(stderr,"ERROR: expected to find molblock number %d (out of %d)\n",i+1,top->n_molblocks);
      fprintf(stderr,"\texpected moltype molname n_mols\n");
      fprintf(stderr,"\tline: %s",inp_line);
      fprintf(stderr,"\tmoltype: %d    molname: %s    n_mols: %d \n",inp_int,inp_word,inp_int2);
      fprintf(stderr,"ERROR: %s %d \n",__FILE__,__LINE__);
      exit(1);
    }
    top->molblocks[i].moltype = inp_int;
    strcpy(top->molblocks[i].molname,inp_word);
    top->molblocks[i].n_mols = inp_int2;
  }

/*
~~~~~~~~~~~~~~~~~~ [moltypes] ~~~~~~~~~~~~~~~~~~
*/
  get_next_line(fp,inp_line);
  inp_int = -1;
  test_sscanf = sscanf(inp_line,"[moltypes] %d ",&inp_int);
  if (test_sscanf != 1) 
  { 
    fprintf(stderr,"ERROR: Expected to find directive \"[moltypes]\" followed by n_moltypes\n");
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"\tn_moltypes: %d\n",inp_int);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
    exit(1); 
  }
  n_mt = inp_int;

  top->contents->mols.nr = n_mt; 
  top->n_molecule_types = n_mt;

  top->molecules = (tW_molecule *) ecalloc(n_mt, sizeof(tW_molecule));

/* MRD 11.7.2019
 * Due to the change around the 2019.4 new dump tpr format,
 * We now lump all dihedrals together. we used to have pdih, rbdih, and tabdih separate.
 * I'm going to check to see if we find 14 numbers, and if so, print a loud warning
 * that the .btp file is likely OLD and that they should regenerate it.
 */

  
  for (i = 0; i < n_mt; ++i)
  {
    tW_molecule *MM;
    MM = &(top->molecules[i]);
    get_next_line(fp,inp_line);
    strcpy(name,"N/A");
    ind = n_mol = n_apm = n_res = n_cgs = excls_nr = excls_nra = bond_nr = angle_nr = dih_nr = bigdih_nr = imnb_nr= -1;
    int testn1, testn2;
    test_sscanf = sscanf(inp_line,"%d %s %d %d %d %d %d %d %d %d %d %d %d %d ",&ind,&name,&n_mol,&n_apm,&n_res,&n_cgs,&excls_nr,&excls_nra,&bond_nr,&angle_nr,&dih_nr,&bigdih_nr,&imnb_nr,&testn1);
    if (test_sscanf == 14)
    {
      fprintf(stderr,"WARNING: We read 14 items from a moltypes line.\n");
      fprintf(stderr,"\tStarting on Nov 7, 2019, the .btp format changed. We now lump all dihedrals together in one,\n");
      fprintf(stderr,"\twhereas older .btp files will have separate numbers for pdih, rbdih, and tabdih.\n");
      fprintf(stderr,"\tPlease re-translate the dumped .tpr file into a new .btp file.\n");
      exit(EXIT_FAILURE);
    }
    if (test_sscanf != 13) 
    { 
      fprintf(stderr,"ERROR: Found %d out of expected 13 moltype properties\n",test_sscanf);
      fprintf(stderr,"\t%10s: %s\n","line",inp_line);
      fprintf(stderr,"\t%10s: %d\n","ind",ind);
      fprintf(stderr,"\t%10s: %s\n","name",name);
      fprintf(stderr,"\t%10s: %d\n","n_mol",n_mol);
      fprintf(stderr,"\t%10s: %d\n","n_apm",n_apm);
      fprintf(stderr,"\t%10s: %d\n","n_res",n_res);
      fprintf(stderr,"\t%10s: %d\n","n_cgs",n_cgs);
      fprintf(stderr,"\t%10s: %d\n","excls_nr",excls_nr);
      fprintf(stderr,"\t%10s: %d\n","excls_nra",excls_nra);
      fprintf(stderr,"\t%10s: %d\n","bond_nr",bond_nr);
      fprintf(stderr,"\t%10s: %d\n","angle_nr",angle_nr);
      fprintf(stderr,"\t%10s: %d\n","dih_nr",dih_nr);
      fprintf(stderr,"\t%10s: %d\n","bigdih_nr",bigdih_nr);
      fprintf(stderr,"\t%10s: %d\n","imnb_nr",imnb_nr);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
      exit(1); 
    }

/* Allocate memory or set values for all the elements of the tW_molecule data structure */
    strcpy(MM->molname,name);
    MM->n_mols = n_mol;
    MM->n_apm = n_apm;
    
    MM->type_ids = (int *) ecalloc(n_apm,sizeof(int));
    MM->m = (double *) ecalloc(n_apm,sizeof(double));
    MM->q = (double *) ecalloc(n_apm,sizeof(double));
    MM->atom_names = (tW_word *) ecalloc(n_apm,sizeof(tW_word));
    MM->atom_types = (tW_word *) ecalloc(n_apm,sizeof(tW_word));
    MM->atom_Btypes = (tW_word *) ecalloc(n_apm,sizeof(tW_word));
    MM->residx = (int *) ecalloc(n_apm,sizeof(int));
   
    MM->n_res = n_res;
    MM->resname = (tW_word *) ecalloc(n_res,sizeof(tW_word));

    MM->n_cg = n_cgs; 
    MM->cg_start = (int *) ecalloc(n_cgs,sizeof(int));
    MM->cg_end = (int *) ecalloc(n_cgs,sizeof(int));
    
    MM->n_excls = excls_nr;
    MM->n_exclsa = excls_nra;
    
    MM->bond_nr = bond_nr;
    MM->n_bonds = bond_nr / BOND_DIV;
    MM->bond_types = (int *) ecalloc(MM->n_bonds,sizeof(int));
    MM->bond_ij = (int **) ecalloc(MM->n_bonds,sizeof(int *));
    for (j = 0; j < MM->n_bonds; ++j) { MM->bond_ij[j] = (int *) ecalloc(2,sizeof(int)); }

    MM->angle_nr = angle_nr;
    MM->n_angles = angle_nr / ANGLE_DIV;
    MM->angle_types = (int *) ecalloc(MM->n_angles,sizeof(int));
    MM->angle_ijk = (int **) ecalloc(MM->n_angles,sizeof(int *));
    for (j = 0; j < MM->n_angles; ++j) { MM->angle_ijk[j] = (int *) ecalloc(3,sizeof(int)); }

    MM->dih_nr = dih_nr;
    MM->n_dihs = dih_nr / DIH_DIV; 
    MM->dih_types = (int *) ecalloc(MM->n_dihs,sizeof(int));
    MM->dih_ijkl = (int **) ecalloc(MM->n_dihs,sizeof(int *));
    for (j = 0; j < MM->n_dihs; ++j) { MM->dih_ijkl[j] = (int *) ecalloc(4,sizeof(int)); }

// MRD 01.24.2020
    MM->bigdih_nr = bigdih_nr;
    MM->n_bigdihs = bigdih_nr / BIGDIH_DIV; 
    MM->bigdih_types = (int *) ecalloc(MM->n_bigdihs,sizeof(int));
    MM->bigdih_ijklm = (int **) ecalloc(MM->n_bigdihs, sizeof(int *));
    for (j = 0; j < MM->n_bigdihs; ++j) { MM->bigdih_ijklm[j] = (int *) ecalloc(5,sizeof(int)); }

    MM->imnb_nr = imnb_nr;
    MM->n_imnbs = imnb_nr / IMNB_DIV;
    MM->imnb_types = (int *) ecalloc(MM->n_imnbs,sizeof(int));
    MM->imnb_ij = (int **) ecalloc(MM->n_imnbs,sizeof(int *));
    for (j = 0; j < MM->n_imnbs; ++j) { MM->imnb_ij[j] = (int *) ecalloc(2,sizeof(int)); }

/* If expected, read this molecule type's [res] directive */
    if (n_res > 0)
    {
      get_next_line(fp,inp_line);
      test_line(inp_line,"[res]",TRUE,"ERROR: unable to find [res]\n");
      for (j = 0; j < MM->n_res; ++j)
      {
	get_next_line(fp,inp_line);
 	test_sscanf = sscanf(inp_line," %s ",&inp_word);
	if (test_sscanf != 1) 
        { 
	  fprintf(stderr,"ERROR: expected resname on this line.\nline: %s",inp_line);	  
	  fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
	  exit(1); 
	}
	strcpy(MM->resname[j],inp_word);
      }     
    }

/* If expected, read this molecule type's [cg] directive */    
    if (n_cgs > 0)
    {
      get_next_line(fp,inp_line);
      test_line(inp_line,"[cg]",TRUE,"ERROR: unable to find [cg]\n");
      for (j = 0; j < MM->n_cg; ++j)
      {
	get_next_line(fp,inp_line);
	ind = i1 = i2 = -1;
	test_sscanf = sscanf(inp_line," %d %d %d ",&ind, &i1, &i2);
        if (test_sscanf != 3) 
	{ 
	  fprintf(stderr,"ERROR: expected to find idx, start, and end of cg group in line\n");
	  fprintf(stderr,"\t line: %s",inp_line);
	  fprintf(stderr,"\t  ind: %d\n\tstart: %d\n\t  end: %d\n",ind,i1,i2);
	  fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
	  exit(1); 
	}
        MM->cg_start[j] = i1;
        MM->cg_end[j] = i2;
      }
    }

/* If expected, read this molecule type's [bond] directive */
    if (bond_nr > 0)
    {
      get_next_line(fp,inp_line);
      test_line(inp_line,"[bond]",TRUE,"ERROR: unable to find [bond]\n");
      for (j = 0; j < MM->n_bonds; ++j)
      {
        get_next_line(fp,inp_line);
	ind = inp_int = i1 = i2 = -1;
        test_sscanf = sscanf(inp_line," %d %d %d %d ",&ind, &inp_int, &i1, &i2);
        if (test_sscanf != 4) 
	{ 
	  fprintf(stderr,"ERROR: Expected to find idx type i j for bond %d\n",j+1);
	  fprintf(stderr,"\tline: %s",inp_line);
	  fprintf(stderr,"\t idx: %d\n\ttype: %d\n\t   i: %d\n\t   j: %d\n",ind,inp_int,i1,i2);
	  fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
	  exit(1); 
	}
        MM->bond_types[j] = inp_int;
        MM->bond_ij[j][0] = i1;
        MM->bond_ij[j][1] = i2;
      }
    }

/* If expected, read this molecule type's [angle] directive */
    if (angle_nr > 0)
    {
      get_next_line(fp,inp_line);
      test_line(inp_line,"[angle]",TRUE,"ERROR: unable to find [angle]\n");
      for (j = 0; j < MM->n_angles; ++j)
      {
        get_next_line(fp,inp_line);
	ind = inp_int = i1 = i2 = i3 = -1;
        test_sscanf = sscanf(inp_line," %d %d %d %d %d ",&ind, &inp_int, &i1, &i2, &i3);
        if (test_sscanf != 5) 
	{ 
	  fprintf(stderr,"ERROR: expected to find idx type i j k for angle %d\n",j+1);
	  fprintf(stderr,"\tline: %s",inp_line);
	  fprintf(stderr,"\t idx: %d\n\ttype: %d\n\t   i: %d\n\t   j: %d\n\t   k: %d\n",ind,inp_int,i1,i2,i3);
	  fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
	  exit(1); 
	}
        MM->angle_types[j] = inp_int;
        MM->angle_ijk[j][0] = i1;
        MM->angle_ijk[j][1] = i2;
        MM->angle_ijk[j][2] = i3;
      }
    }

/* If expected, read this molecule type's [dih] directive */
    if (dih_nr > 0)
    {
      get_next_line(fp,inp_line);
      test_line(inp_line,"[dih]",TRUE,"ERROR: unable to find [dih]\n");
      for (j = 0; j < MM->n_dihs; ++j)
      {
        get_next_line(fp,inp_line);
	ind = inp_int = i1 = i2 = i3 = i4 = -1;
        test_sscanf = sscanf(inp_line," %d %d %d %d %d %d ",&ind, &inp_int, &i1, &i2, &i3, &i4);
        if (test_sscanf != 6) 
	{ 
	  fprintf(stderr,"ERROR: expected to find idx type i j k l for dih %d\n",j+1);
          fprintf(stderr,"\tline: %s",inp_line);
          fprintf(stderr,"\t idx: %d\n\ttype: %d\n\t   i: %d\n\t   j: %d\n\t   k: %d\n\t   l: %d\n",
									ind,inp_int,i1,i2,i3,i4);
	  fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
	  exit(1); 
	}
        MM->dih_types[j] = inp_int;
        MM->dih_ijkl[j][0] = i1;
        MM->dih_ijkl[j][1] = i2;
        MM->dih_ijkl[j][2] = i3;
        MM->dih_ijkl[j][3] = i4;
      }
    }

/* If expected, read this molecule type's [bigdih] directive MRD 01.24.2020 */
    if (bigdih_nr > 0)
    {
      get_next_line(fp,inp_line);
      test_line(inp_line,"[bigdih]",TRUE,"ERROR: unable to find [dih]\n");
      for (j = 0; j < MM->n_bigdihs; ++j)
      {
        get_next_line(fp,inp_line);
        ind = inp_int = i1 = i2 = i3 = i4 = i5 = -1;
        test_sscanf = sscanf(inp_line," %d %d %d %d %d %d %d ",&ind, &inp_int, &i1, &i2, &i3, &i4, &i5);
        if (test_sscanf != 7)
        {
          fprintf(stderr,"ERROR: expected to find idx type i j k l m for bigdih %d\n",j+1);
          fprintf(stderr,"\tline: %s",inp_line);
          fprintf(stderr,"\t idx: %d\n\ttype: %d\n\t   i: %d\n\t   j: %d\n\t   k: %d\n\t   l: %d\n\t   m: %d\n",
                                                                        ind,inp_int,i1,i2,i3,i4,i5);
          fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
          exit(1);
        }
        MM->bigdih_types[j] = inp_int;
        MM->bigdih_ijklm[j][0] = i1;
        MM->bigdih_ijklm[j][1] = i2;
        MM->bigdih_ijklm[j][2] = i3;
        MM->bigdih_ijklm[j][3] = i4;
        MM->bigdih_ijklm[j][4] = i5;
      }
    }

/* If expected, read this molecule type's [imnb] directive */
    if (imnb_nr > 0)
    {
      get_next_line(fp,inp_line);
      test_line(inp_line,"[imnb]",TRUE,"ERROR: unable to find [imnb]\n");
      for (j = 0; j < MM->n_imnbs; ++j)
      {
        get_next_line(fp,inp_line);
        ind = inp_int = i1 = i2 = -1;
        test_sscanf = sscanf(inp_line," %d %d %d %d ",&ind, &inp_int, &i1, &i2);
        if (test_sscanf != 4)
        {
          fprintf(stderr,"ERROR: expected to find idx type i j for imnb %d\n",j+1);
          fprintf(stderr,"\tline: %s",inp_line);
          fprintf(stderr,"\t idx: %d\n\ttype: %d\n\t   i: %d\n\t   j:%d\n",ind,inp_int,i1,i2);
          fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
          exit(1);
        }
        MM->imnb_types[j] = inp_int;
        MM->imnb_ij[j][0] = i1;
        MM->imnb_ij[j][1] = i2;
      }
    }
  }

/*
~~~~~~~~~~~~~~~~~~ [attypes] ~~~~~~~~~~~~~~~~~~
*/
  get_next_line(fp,inp_line);
  test_line(inp_line,"[attypes]",TRUE,"ERROR: unable to find [attypes]\n");
  for (i = 0; i < top->contents->mols.nr; ++i)
  {
    top->molecules[i].n_epa = (int *) ecalloc(top->molecules[i].n_apm,sizeof(int));
    top->molecules[i].excls = (int **) ecalloc(top->molecules[i].n_apm,sizeof(int *));
    for (j = 0; j < top->molecules[i].n_apm; ++j)
    {
      tW_molecule *MM;
      MM = &(top->molecules[i]);
      get_next_line(fp,inp_line);
      i1 = i2 = i3 = i4 = n_epa = -1;
      inp_flt1 = inp_flt2 = -1.0;
      strcpy(inp_word,"N/A");
      strcpy(name,"N/A");
      test_sscanf = sscanf(inp_line," %d %d %d %f %f %d %s %s %d ",&i1, &i2, &i3, &inp_flt1, &inp_flt2, &i4, &inp_word, &name, &n_epa);
      if (test_sscanf != 9) 
      { 
	fprintf(stderr,"ERROR: read %d out of 9 expected attype properties\n",test_sscanf);
	fprintf(stderr,"\tline: %s\n",inp_line);
	fprintf(stderr,"\t%12s: %d\n","moltype_ind",i1);
        fprintf(stderr,"\t%12s: %d\n","at_ind",i2);
        fprintf(stderr,"\t%12s: %d\n","type",i3);
        fprintf(stderr,"\t%12s: %f\n","mass",inp_flt1);
        fprintf(stderr,"\t%12s: %f\n","charge",inp_flt2);
        fprintf(stderr,"\t%12s: %d\n","resind",i4);
        fprintf(stderr,"\t%12s: %s\n","atomname",inp_word);
        fprintf(stderr,"\t%12s: %s\n","typename",name);
        fprintf(stderr,"\t%12s: %d\n","n_epa",n_epa);
	fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
	exit(1); 
      }
      (top->molecules[i].type_ids[j]) = i3;
      (top->molecules[i].m[j]) = (double) inp_flt1;
      (top->molecules[i].q[j]) = (double) inp_flt2;
      (top->molecules[i].residx[j]) = i4;
      strcpy((top->molecules[i].atom_names[j]),inp_word);
      strcpy((top->molecules[i].atom_types[j]),name);
      strcpy((top->molecules[i].atom_Btypes[j]),name);
      (top->molecules[i].n_epa[j]) = n_epa;      
      top->molecules[i].excls[j] = (int *) ecalloc(n_epa,sizeof(int));
      if (n_epa > 0)
      {
        get_next_line(fp,inp_line);
        test_line(inp_line,"[epa]",TRUE,"ERROR: unable to find [epa]\n");
	get_next_line(fp,inp_line);
	elim_char(inp_line,'\n');
        int n_words = get_word_count_delim(inp_line," ");
        tW_word * word_list = (tW_word *) ecalloc(n_words,sizeof(tW_word));
	get_words_delim(inp_line," ",word_list);
        for (k = 0; k < n_words; ++k) { (top->molecules[i].excls[j][k]) = atoi(word_list[k]); }
        efree(word_list); 
      }
    }
  }

/*
~~~~~~~~~~~~~~~~~~ [Bocs_at_types] ~~~~~~~~~~~~~~~~~~
*/
  get_next_line(fp,inp_line);
  test_line(inp_line,"[Bocs_at_types]",TRUE,"ERROR: unable to find [Bocs_at_types]\n");
  test_sscanf = sscanf(inp_line,"[Bocs_at_types] %d ",&inp_int);
  if (test_sscanf != 1)
  {
    fprintf(stderr,"ERROR: Unable to read number of bocs at types from directive line:\n%s",inp_line);
    exit(1);
  }
  
  int n_my_atom_types = inp_int;
  tW_word *Bocs_type_names = (tW_word *) calloc(n_my_atom_types,sizeof(tW_word));
  for (i = 0; i < n_my_atom_types; ++i)
  {
    get_next_line(fp,inp_line);
    test_sscanf = sscanf(inp_line," %d %s ",&inp_int,&(Bocs_type_names[i]));
    if (test_sscanf != 2)
    {
      fprintf(stderr,"ERROR: unable to read index and name of type idx %d under [Bocs_at_types] directive\n",i);
      fprintf(stderr,"\ttest_sscanf (should be 2): %d   idx (should be %d): %d  name: %s\n",test_sscanf,i,inp_int,Bocs_type_names[i]);
      fprintf(stderr,"\tline: %s",inp_line);
      exit(1);
    } 
  }

  fclose(fp);


  dump_molecule_info(top);
  pop_contents(top);

  if (top->n_atomtypes != n_my_atom_types)
  {
    fprintf(stderr,"ERROR: n_Bocs_at_types read from file is different from the value determined from the molecule information!\n");
    fprintf(stderr,"\tread: %d    determined: %d\n",n_my_atom_types,top->n_atomtypes);
  }

  for (i = 0; i < n_my_atom_types; ++i)
  {
    if (strcmp(Bocs_type_names[i],top->atom_type_names[i]) != 0)
    {
      fprintf(stderr,"ERROR: Read Bocs_at_type name %d as: %s but determined it to be: %s\n",i,Bocs_type_names[i],top->atom_type_names[i]);
      exit(1);
    }
  }
  free(Bocs_type_names);


  return TRUE;
}


void write_bocs_frame(tW_gmx_trxframe * fr)
{
  int i, j, k;
  fprintf(fr->fp,"[FRAME] %d\n\n",fr->counter);
  fprintf(fr->fp,"[natoms] %d\n\n",fr->contents->natoms);

  fprintf(fr->fp,"[flags]\n!double title step time atoms prec box x v f ld ldg\n");
  fprintf(fr->fp," %d      %d     %d    %d    %d     %d    %d   %d %d %d %d  %d\n\n",
		(fr->contents->bDouble ? 1 : 0),
		(fr->contents->bTitle ? 1 : 0),
                (fr->contents->bStep ? 1 : 0),
                (fr->contents->bTime ? 1 : 0),
                (fr->contents->bAtoms ? 1 : 0),
                (fr->contents->bPrec ? 1 : 0),
                (fr->contents->bBox ? 1 : 0),
                (fr->contents->bX ? 1 : 0),
                (fr->contents->bV ? 1 : 0),
                (fr->contents->bF ? 1 : 0),
                (fr->bLDs ? 1 : 0),
                (fr->bLD_Grads ? 1 : 0));

  if (fr->contents->bTitle) { fprintf(fr->fp,"[title] %s\n\n",fr->contents->title); }
  if (fr->contents->bStep) { fprintf(fr->fp,"[step] %d\n\n",fr->contents->step); }
  if (fr->contents->bTime) { fprintf(fr->fp,"[time] %9.5e\n\n",fr->contents->time); }
  if (fr->contents->bPrec) { fprintf(fr->fp,"[prec] %g\n\n",(double) fr->contents->prec); }

  if (fr->contents->bBox)
  {
    fprintf(fr->fp,"[box]\n");
    for (i = 0; i < 3; ++i)
    {
      for (j = 0; j < 3; ++j)
      {
        fprintf(fr->fp,"%9.7e  ",fr->contents->box[i][j]);
      }
      fprintf(fr->fp,"\n");
    }
    fprintf(fr->fp,"\n");
  }
  
  if (fr->contents->bX || fr->contents->bV || fr->contents->bF)
  {
    fprintf(fr->fp,"[%c%c%c]\n",(fr->contents->bX ? 'x' : '_'),(fr->contents->bV ? 'v' : '_'),(fr->contents->bF ? 'f' : '_'));
    for (i = 0; i < fr->contents->natoms; ++i)
    {
// TODO: 9.17.2020: find a way to trim a bunch of potential extra 0s off of the end of these numbers
      if (fr->contents->bX) { fprintf(fr->fp,"% 11.9e % 11.9e % 11.9e  ",fr->contents->x[i][0],fr->contents->x[i][1],fr->contents->x[i][2]); }
      if (fr->contents->bV) { fprintf(fr->fp,"% 7.4f % 7.4f % 7.4f  ",fr->contents->v[i][0],fr->contents->v[i][1],fr->contents->v[i][2]); }
      if (fr->contents->bF) { fprintf(fr->fp,"% 11.9e % 11.9e % 11.9e  ",fr->contents->f[i][0],fr->contents->f[i][1],fr->contents->f[i][2]); }
      fprintf(fr->fp,"\n");
    }
    fprintf(fr->fp,"\n");
  }

  if (fr->bLDs)
  {
    fprintf(fr->fp,"[local_densities] %d\n",fr->n_atomtypes);
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      for (j = 0; j < fr->n_atomtypes; ++j)
      {
// TODO: 9.17.2020: find a way to trim a bunch of potential extra 0s off of the end of these numbers
        fprintf(fr->fp," %14.10f ",fr->local_densities[i][j]);
      }
      fprintf(fr->fp,"\n");
    }
    fprintf(fr->fp,"\n");
    if (fr->bLD_Grads)
    {
      fprintf(fr->fp,"[local_density_gradients] \n");
      for (i = 0; i < fr->contents->natoms; ++i)
      {
        for (j = 0; j < fr->n_atomtypes; ++j)
        {
          for (k = 0; k < DIM; ++k)
          {
// TODO: 9.17.2020: find a way to trim a bunch of potential extra 0s off of the end of these numbers
            fprintf(fr->fp," %14.10f ",fr->local_density_gradients[i][j][k]);
          }
        }
        fprintf(fr->fp,"\n");
      }
      fprintf(fr->fp,"\n");
    } 
  }
 
  fprintf(fr->fp,"[END_FRAME]\n\n");
}

/*
check_dirs(): This function is called repeatedly when reading a bocs frame. It checks to make sure
	that the first three directives ([FRAME], [natoms], [flags]) were found, and is noisy if
	they haven't been found yet.
*/
void check_dirs(const char *current_dir, int ntest, bool bFRAME, bool bNATOMS, bool bFLAGS, int count)
{
  switch (ntest)
  {
    case 1:
      if (!bFRAME) { fprintf(stderr,"WARNING: found directive %s before [FRAME] in frame %d. [FRAME] should be first\n",current_dir,count); }
      break;
    case 2:
      if (!bFRAME) { fprintf(stderr,"WARNING: found directive %s before [FRAME] in frame %d. [FRAME] should be first\n",current_dir,count); }
      if (!bNATOMS) { fprintf(stderr,"WARNING: found directive %s before [natoms] in frame %d. [natoms] should be second\n",current_dir,count); }
      break;
    case 3:
      if (!bFRAME) { fprintf(stderr,"WARNING: found directive %s before [FRAME] in frame %d. [FRAME] should be first\n",current_dir,count); }
      if (!bNATOMS) { fprintf(stderr,"WARNING: found directive %s before [natoms] in frame %d. [natoms] should be second\n",current_dir,count); }
      if (!bFLAGS) { fprintf(stderr,"WARNING: found directive %s before [flags] in frame %d. [flags] should be third\n",current_dir,count); }
      break;
  }
}

/*
check_double_bocs_dir(): This function is called repeatedly when reading a bocs frame.
	It checks to see if you have already found a directive, and prints a warning
	message if you have.
*/
void check_double_bocs_dir(bool bTEST, const char *dir, int count)
{
  if (bTEST) { fprintf(stderr,"WARNING: found directive %s twice in frame %d\n",dir,count); }
}

/*
check_unexpected_dir(): This function is called repeatedly when reading a bocs frame.
	It checks to see if you found a directive that you weren't expecting based
	on the flags for the frame.
*/
void check_unexpected_dir(bool bTEST, const char *dir, int count)
{
  if (! bTEST) { fprintf(stderr,"WARNING: found directive %s in frame %d even though title flag is FALSE\n",dir,count); }
}

/*
read_bocs_frame_dir_XXXXX(): The next bunch of functions handle the reading 
	of the various directives and following data in a bocs trajectory file
*/
bool read_bocs_frame_dir_frame(tW_gmx_trxframe *fr, tW_line inp_line)
{
  int test_sscanf, frnr;
  test_sscanf = sscanf(inp_line,"[FRAME] %d ",&frnr);
  if (test_sscanf != 1)
  {
    fprintf(stderr,"ERROR: expected [FRAME] <frame_number> on line\n");
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"\tframe_number: %d\n",frnr);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    exit(1);
  }
  get_next_line(fr->fp,inp_line);
  return TRUE;
}

bool read_bocs_frame_dir_natoms(tW_gmx_trxframe *fr, tW_line inp_line)
{
  int test_sscanf, natoms;
  test_sscanf = sscanf(inp_line,"[natoms] %d ", &natoms);
  if (test_sscanf != 1)
  {
    fprintf(stderr,"ERROR: expected [natoms] <#atoms> on line\n");
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"\tn_atoms: %d\n",natoms);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    exit(1);
  }
  fr->contents->natoms = natoms;
  get_next_line(fr->fp, inp_line);
  return TRUE;
}

bool read_bocs_frame_dir_flags(tW_gmx_trxframe *fr, tW_line inp_line)
{
  int test_sscanf, dblFlag, ttlFlag, stpFlag, timFlag, atmFlag, prcFlag, boxFlag, xFlag, vFlag, fFlag;
 
  get_next_line(fr->fp,inp_line);

  test_sscanf = sscanf(inp_line," %d %d %d %d %d %d %d %d %d %d ",&dblFlag, &ttlFlag, &stpFlag, &timFlag, &atmFlag, &prcFlag, &boxFlag, &xFlag, &vFlag, &fFlag);
  if (test_sscanf != 10)
  {
    fprintf(stderr,"ERROR: expected to read 10 flags under directive [flags]\n");
    fprintf(stderr,"\tI was only able to read %d\n",test_sscanf);
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    exit(1);
  }

  if (fr->contents->bDouble && (dblFlag == 0)) { fprintf(stderr,"WARNING: bDouble flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bDouble = FALSE; }
  else if (! fr->contents->bDouble && (dblFlag == 1)) { fprintf(stderr,"WARNING: bDouble flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bDouble = TRUE; }

  if (fr->contents->bTitle && (ttlFlag == 0)) { fprintf(stderr,"WARNING: bTitle flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bTitle = FALSE; }
  else if (! fr->contents->bTitle && (ttlFlag == 1)) { fprintf(stderr,"WARNING: bTitle flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bTitle = TRUE; }

  if (fr->contents->bStep && (stpFlag == 0)) { fprintf(stderr,"WARNING: bStep flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bStep = FALSE; }
  else if (! fr->contents->bStep && (stpFlag == 1)) { fprintf(stderr,"WARNING: bStep flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bStep = TRUE; }

  if (fr->contents->bTime && (timFlag == 0)) { fprintf(stderr,"WARNING: bTime flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bTime = FALSE; }
  else if (! fr->contents->bTime && (timFlag == 1)) { fprintf(stderr,"WARNING: bTime flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bTime = TRUE; }

  if (fr->contents->bAtoms && (atmFlag == 0)) { fprintf(stderr,"WARNING: bAtoms flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bAtoms = FALSE; }
  else if (! fr->contents->bAtoms && (atmFlag == 1)) { fprintf(stderr,"WARNING: bAtoms flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bAtoms = TRUE; }

  if (fr->contents->bPrec && (prcFlag == 0)) { fprintf(stderr,"WARNING: bPrec flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bPrec = FALSE; }
  else if (! fr->contents->bPrec && (prcFlag == 1)) { fprintf(stderr,"WARNING: bPrec flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bPrec = TRUE; }

  if (fr->contents->bBox && (boxFlag == 0)) { fprintf(stderr,"WARNING: bBox flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bBox = FALSE; }
  else if (! fr->contents->bBox && (boxFlag == 1)) { fprintf(stderr,"WARNING: bBox flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bBox = TRUE; }

  if (fr->contents->bX && (xFlag == 0)) { fprintf(stderr,"WARNING: bX flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bX = FALSE; }
  else if (! fr->contents->bX && (xFlag == 1)) { fprintf(stderr,"WARNING: bX flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bX = TRUE; }

  if (fr->contents->bV && (vFlag == 0)) { fprintf(stderr,"WARNING: bV flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bV = FALSE; }
  else if (! fr->contents->bV && (vFlag == 1)) { fprintf(stderr,"WARNING: bV flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bV = TRUE; }

  if (fr->contents->bF && (fFlag == 0)) { fprintf(stderr,"WARNING: bF flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bF = FALSE; }
  else if (! fr->contents->bF && (fFlag == 1)) { fprintf(stderr,"WARNING: bF flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bF = TRUE; }

  get_next_line(fr->fp, inp_line);
  return TRUE;
}

bool read_bocs_frame_dir_title(tW_gmx_trxframe * fr, tW_line inp_line)
{
  int test_sscanf;
  tW_word inp_word;
  strcpy(inp_word,"N/A");
  test_sscanf = sscanf(inp_line,"[title] %s ",&inp_word);
  if (test_sscanf != 1)
  {
    fprintf(stderr,"ERROR: Expected [title] <TITLE> in line\n");
    fprintf(stderr,"\tline:%s\ttitle: %s\n",inp_line,inp_word);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    exit(1);
  }
  if (! fr->contents->title) { fr->contents->title = (char *) ecalloc(100,sizeof(char)); }
  strcpy(fr->contents->title,inp_word);

  get_next_line(fr->fp, inp_line);

  return TRUE;
}

bool read_bocs_frame_dir_step(tW_gmx_trxframe * fr, tW_line inp_line)
{
  int test_sscanf, step = -1;
  test_sscanf = sscanf(inp_line,"[step] %d ",&step);
  if (test_sscanf != 1)
  {
    fprintf(stderr,"ERROR: Expected [step] <STEP#> in line\n");
    fprintf(stderr,"\tline:%s\tstep#: %d\n",inp_line,step);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    exit(1);
  }
  fr->contents->step = step;
  get_next_line(fr->fp, inp_line);
  
  return TRUE;
}

bool read_bocs_frame_dir_time(tW_gmx_trxframe * fr, tW_line inp_line)
{
  int test_sscanf;
  float time;
  test_sscanf = sscanf(inp_line,"[time] %f ",&time);
  if (test_sscanf != 1)
  {
    fprintf(stderr,"ERROR: Expected [time] <TIME> in line\n");
    fprintf(stderr,"\tline:%s\ttime: %f\n",inp_line,time);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    exit(1);
  }
  fr->contents->time = (double) time;
  get_next_line(fr->fp, inp_line);

  return TRUE;
}

bool read_bocs_frame_dir_prec(tW_gmx_trxframe *fr, tW_line inp_line)
{
  int test_sscanf;
  float prec;
  test_sscanf = sscanf(inp_line,"[prec] %f ",&prec);
  if (test_sscanf != 1)
  {
    fprintf(stderr,"ERROR: expected [prec] <PREC> in line\n");
    fprintf(stderr,"\tline:%s\tprec: %f\n",inp_line,prec);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
    exit(1);
  }
  fr->contents->prec = (double) prec;
  get_next_line(fr->fp, inp_line);
  
  return TRUE;
}

bool read_bocs_frame_dir_box(tW_gmx_trxframe *fr, tW_line inp_line)
{
  int test_sscanf, i;
  float b0, b1, b2;
  for (i = 0; i < DIM; ++i)
  {
    get_next_line(fr->fp, inp_line);
    test_sscanf = sscanf(inp_line," %f %f %f ",&b0, &b1, &b2);
    if (test_sscanf != 3)
    {
      fprintf(stderr,"ERROR: unable to read line %d/3 of box dimensions!\n",i+1);
      fprintf(stderr,"\tline: %s\tb0: %f  b1: %f  b2: %f\n",inp_line,b0,b1,b2);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__);
      exit(1);
    }
    fr->contents->box[i][0] = (double) b0;
    fr->contents->box[i][1] = (double) b1;
    fr->contents->box[i][2] = (double) b2;
  }
  get_next_line(fr->fp, inp_line);
  return TRUE;
}

bool read_bocs_frame_dir_xvf(tW_gmx_trxframe *fr, tW_line inp_line)
{
  int test_sscanf, i;
  int nprops = 0;
  int pidx;
  tW_word props[9];
  tW_line orig_line;
  strcpy(orig_line,inp_line);

  if (strstr(inp_line,"x") != NULL) { nprops += 3; }
  if (strstr(inp_line,"v") != NULL) { nprops += 3; }
  if (strstr(inp_line,"f") != NULL) { nprops += 3; }

  for (i = 0; i < fr->contents->natoms; ++i)
  {
    get_next_line(fr->fp, inp_line);
    test_sscanf = sscanf(inp_line," %s %s %s %s %s %s %s %s %s ",
                                    &(props[0]), &(props[1]), &(props[2]),
                                    &(props[3]), &(props[4]), &(props[5]),
                                    &(props[6]), &(props[7]), &(props[8]));
    if (test_sscanf != nprops) 
    { 
      fprintf(stderr,"ERROR: expected %d properties per atom under directive %s\n",nprops,orig_line);
      fprintf(stderr,"       However, we found %s properties.\n",test_sscanf);
      exit(EXIT_FAILURE); 
    }
    pidx = 0;
    if (fr->contents->bX)
    {
      fr->contents->x[i][0] = atof(props[pidx]);
      fr->contents->x[i][1] = atof(props[pidx+1]);
      fr->contents->x[i][2] = atof(props[pidx+2]);
      pidx += 3;
    }
    if (fr->contents->bV)
    {
      fr->contents->v[i][0] = atof(props[pidx]);
      fr->contents->v[i][1] = atof(props[pidx+1]);
      fr->contents->v[i][2] = atof(props[pidx+2]);
      pidx += 3;
    }
    if (fr->contents->bF)
    {
      fr->contents->f[i][0] = atof(props[pidx]);
      fr->contents->f[i][1] = atof(props[pidx+1]);
      fr->contents->f[i][2] = atof(props[pidx+2]);
      pidx += 3;
    }
  }
  get_next_line(fr->fp,inp_line);

  return TRUE;

}

/* MRD 06.18.2020 */
bool read_bocs_frame_dir_local_densities(tW_gmx_trxframe *fr, tW_line inp_line)
{
  int i, test_sscanf;
  for (i = 0; i < fr->contents->natoms; ++i)
  {
    get_next_line(fr->fp, inp_line);
    pop_LD_vector(inp_line,fr->local_densities[i],fr->n_atomtypes);
  }
  get_next_line(fr->fp,inp_line);
}

bool read_bocs_frame_dir_local_density_gradients(tW_gmx_trxframe *fr, tW_line inp_line)
{
  int i, test_sscanf;
  for (i = 0; i < fr->contents->natoms; ++i)
  {
    get_next_line(fr->fp,inp_line);
    populate_LD_Gradients(inp_line, fr->local_density_gradients[i],fr->n_atomtypes);
  }
  get_next_line(fr->fp,inp_line);
}



/*
new_read_next_bocs_frame(): The newer function to read a bocs frame.
	it goes until it finds the directive [END_FRAME]. In the mean time,
	it tests to see if it found one of the various directives that are supported.
	If it has, it calls the corresponding function to read that directive (above).
*/
bool new_read_next_bocs_frame(tW_gmx_trxframe * fr)
{
  bool bEND, bFRAME, bNATOMS, bFLAGS, bTITLE, bSTEP, bTIME, bPREC, bBOX, bXVF, bLD, bLDG;
  bEND = bFRAME = bNATOMS = bFLAGS = bTITLE = bSTEP = bTIME = bPREC = bBOX = bXVF = bLD = bLDG = FALSE;
  char *xvf_dir = (char *) ecalloc(6,sizeof(char));
  tW_line inp_line;

//  fprintf(stderr,"Reading frame: %d\r",fr->counter);
  if (get_next_line(fr->fp,inp_line) == -1) { fprintf(stderr,"Done after frame %d\n",fr->counter - 1); return FALSE; }
  sprintf(xvf_dir,"[%c%c%c]",(fr->contents->bX ? 'x' : '_'),(fr->contents->bV ? 'v' : '_'),(fr->contents->bF ? 'f' : '_'));
  while (!bEND)
  {
    if (strstr(inp_line,"[FRAME]") != NULL) 
    { 
      check_double_bocs_dir(bFRAME,"[FRAME]",fr->counter);
      bFRAME = read_bocs_frame_dir_frame(fr,inp_line); 
    }
    else if (strstr(inp_line,"[natoms]") != NULL) 
    { 
      check_dirs("[natoms]",1,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bNATOMS,"[natoms]",fr->counter);
      bNATOMS = read_bocs_frame_dir_natoms(fr,inp_line); 
    }    
    else if (strstr(inp_line,"[flags]") != NULL)
    {
      check_dirs("[flags]",2,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bFLAGS,"[flags]",fr->counter);
      bFLAGS = read_bocs_frame_dir_flags(fr,inp_line);
      sprintf(xvf_dir,"[%c%c%c]",(fr->contents->bX ? 'x' : '_'),(fr->contents->bV ? 'v' : '_'),(fr->contents->bF ? 'f' : '_')); // Do this again incase any of them flipped
    }
    else if (strstr(inp_line,"[title]") != NULL)
    {
      check_dirs("[title]",3,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bTITLE,"[title]",fr->counter);
      check_unexpected_dir(fr->contents->bTitle,"[title]",fr->counter);
      bTITLE = read_bocs_frame_dir_title(fr,inp_line);
    }
    else if (strstr(inp_line,"[step]") != NULL)
    {
      check_dirs("[step]",3,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bSTEP,"[step]",fr->counter);
      check_unexpected_dir(fr->contents->bStep,"[step]",fr->counter);
      bSTEP = read_bocs_frame_dir_step(fr,inp_line);
    }
    else if (strstr(inp_line,"[time]") != NULL)
    {
      check_dirs("[time]",3,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bTIME,"[time]",fr->counter);
      check_unexpected_dir(fr->contents->bTime,"[time]",fr->counter);
      bTIME = read_bocs_frame_dir_time(fr,inp_line);
    }
    else if (strstr(inp_line,"[prec]") != NULL)
    {
      check_dirs("[prec]",3,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bPREC,"[prec]",fr->counter);
      check_unexpected_dir(fr->contents->bPrec,"[prec]",fr->counter);
      bPREC = read_bocs_frame_dir_prec(fr,inp_line);
    }
    else if (strstr(inp_line,"[box]") != NULL)
    {
      check_dirs("[box]",3,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bBOX,"[box]",fr->counter);
      check_unexpected_dir(fr->contents->bBox,"[box]",fr->counter);
      bBOX = read_bocs_frame_dir_box(fr,inp_line);
    } 
    else if (strstr(inp_line,xvf_dir) != NULL)
    {
      check_dirs(xvf_dir,3,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bXVF,xvf_dir,fr->counter);
      check_unexpected_dir((fr->contents->bX || fr->contents->bV || fr->contents->bF),xvf_dir,fr->counter);
      bXVF = read_bocs_frame_dir_xvf(fr,inp_line);
    }
    /* MRD 06.18.2020 */
    else if (strstr(inp_line,"[local_densities]") != NULL)
    {
      check_dirs("[local_densities]",3,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bLD,"[local_densities]",fr->counter);
      check_unexpected_dir(fr->bLDs,"[local_densities]",fr->counter);
      bLD = read_bocs_frame_dir_local_densities(fr,inp_line);
    }
    else if (strstr(inp_line,"[local_density_gradients]") != NULL)
    {
      check_dirs("[local_density_gradients]",3,bFRAME,bNATOMS,bFLAGS,fr->counter);
      check_double_bocs_dir(bLDG,"[local_density_gradients]",fr->counter);
      check_unexpected_dir(fr->bLD_Grads,"[local_density_gradients]",fr->counter);
      bLDG = read_bocs_frame_dir_local_density_gradients(fr,inp_line);
    }  
    else if (strstr(inp_line,"[END_FRAME]") != NULL) { bEND = TRUE; }
    else
    {
      fprintf(stderr,"WARNING: found the following line \"outside\" of a directive\n");
      fprintf(stderr,"\tline: %s",inp_line);
      fprintf(stderr,"\tI will simply keep reading lines until I find the next directive\n");
      get_next_line(fr->fp,inp_line);
    }
  }  
}

/*
read_next_bocs_frame(): This is the old function for reading the next bocs frame.
	I rewrote it so it can overcome more "errors" on its own, and the rewrite
	is above (new_read_next_bocs_frame). I don't think I actually call this
	function any more, but I didn't want to get rid of it just in case.
*/
bool read_next_bocs_frame(tW_gmx_trxframe * fr)
{
//  ++(fr->counter);
//  fprintf(stderr,"Reading frame: %d\r",fr->counter);
  tW_line inp_line;
  tW_word inp_word;
  int frnr, natoms, test_sscanf;
  int dblFlag, ttlFlag, stpFlag, timFlag, atmFlag, prcFlag, boxFlag, xFlag, vFlag, fFlag;
  int step;
  float time;
  float b[3];
  float v[3];
  tW_word xs[3], fs[3];
  double x[3], f[3]; 
  int i, j;

  if (get_next_line(fr->fp,inp_line) == -1) { fprintf(stderr,"Done after frame %d\n",fr->counter/* - 1*/); return FALSE; }
  test_line(inp_line, "[FRAME]", TRUE, "ERROR: Expected [FRAME] ####\n");
  frnr = -1;
  test_sscanf = sscanf(inp_line,"[FRAME] %d ",&frnr);
  if (test_sscanf != 1) 
  { 
    fprintf(stderr,"ERROR: expected [FRAME] <frame_number> on line\n");
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"\tframe_number: %d\n",frnr);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
    exit(1); 
  }
  
  get_next_line(fr->fp,inp_line);
  test_line(inp_line,"[natoms]", TRUE, "ERROR: Expected [natoms] ####\n");
  natoms = -1;
  test_sscanf = sscanf(inp_line,"[natoms] %d ",&natoms);
  if (test_sscanf != 1) 
  { 
    fprintf(stderr,"ERROR: Expected [natoms] <n_atoms> on line\n");
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"\tnatoms: %d\n",natoms);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
    exit(1); 
  }
  fr->contents->natoms = natoms;

  get_next_line(fr->fp,inp_line);
  test_line(inp_line,"[flags]",TRUE,"ERROR: Expected [flags]\n");
  get_next_line(fr->fp,inp_line);
  
  test_sscanf = sscanf(inp_line," %d %d %d %d %d %d %d %d %d %d ",&dblFlag, &ttlFlag, &stpFlag, &timFlag, &atmFlag, &prcFlag, &boxFlag, &xFlag, &vFlag, &fFlag);
  if (test_sscanf != 10) 
  { 
    fprintf(stderr,"ERROR: expected to read 10 flags under directive [flags]\n");
    fprintf(stderr,"\tI was only able to read %d\n",test_sscanf);
    fprintf(stderr,"\tline: %s",inp_line);
    fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
    exit(1); 
  }

  if (fr->contents->bDouble && (dblFlag == 0)) { fprintf(stderr,"WARNING: bDouble flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bDouble = FALSE; }
  else if (! fr->contents->bDouble && (dblFlag == 1)) { fprintf(stderr,"WARNING: bDouble flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bDouble = TRUE; }

  if (fr->contents->bTitle && (ttlFlag == 0)) { fprintf(stderr,"WARNING: bTitle flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bTitle = FALSE; }
  else if (! fr->contents->bTitle && (ttlFlag == 1)) { fprintf(stderr,"WARNING: bTitle flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bTitle = TRUE; }

  if (fr->contents->bStep && (stpFlag == 0)) { fprintf(stderr,"WARNING: bStep flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bStep = FALSE; }
  else if (! fr->contents->bStep && (stpFlag == 1)) { fprintf(stderr,"WARNING: bStep flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bStep = TRUE; }

  if (fr->contents->bTime && (timFlag == 0)) { fprintf(stderr,"WARNING: bTime flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bTime = FALSE; }
  else if (! fr->contents->bTime && (timFlag == 1)) { fprintf(stderr,"WARNING: bTime flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bTime = TRUE; }

  if (fr->contents->bAtoms && (atmFlag == 0)) { fprintf(stderr,"WARNING: bAtoms flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bAtoms = FALSE; }
  else if (! fr->contents->bAtoms && (atmFlag == 1)) { fprintf(stderr,"WARNING: bAtoms flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bAtoms = TRUE; }

  if (fr->contents->bPrec && (prcFlag == 0)) { fprintf(stderr,"WARNING: bPrec flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bPrec = FALSE; }
  else if (! fr->contents->bPrec && (prcFlag == 1)) { fprintf(stderr,"WARNING: bPrec flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bPrec = TRUE; }

  if (fr->contents->bBox && (boxFlag == 0)) { fprintf(stderr,"WARNING: bBox flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bBox = FALSE; }
  else if (! fr->contents->bBox && (boxFlag == 1)) { fprintf(stderr,"WARNING: bBox flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bBox = TRUE; }

  if (fr->contents->bX && (xFlag == 0)) { fprintf(stderr,"WARNING: bX flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bX = FALSE; }
  else if (! fr->contents->bX && (xFlag == 1)) { fprintf(stderr,"WARNING: bX flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bX = TRUE; }

  if (fr->contents->bV && (vFlag == 0)) { fprintf(stderr,"WARNING: bV flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bV = FALSE; }
  else if (! fr->contents->bV && (vFlag == 1)) { fprintf(stderr,"WARNING: bV flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bV = TRUE; }

  if (fr->contents->bF && (fFlag == 0)) { fprintf(stderr,"WARNING: bF flipped from 1 to 0 at counter %d\n",fr->counter); fr->contents->bF = FALSE; }
  else if (! fr->contents->bF && (fFlag == 1)) { fprintf(stderr,"WARNING: bF flipped from 0 to 1 at counter %d\n",fr->counter); fr->contents->bF = TRUE; }

/*  fr->contents->bDouble = ((dblFlag == 1) ? TRUE : FALSE);
  fr->contents->bTitle = ((ttlFlag == 1) ? TRUE : FALSE);
  fr->contents->bStep = ((stpFlag == 1) ? TRUE : FALSE);
  fr->contents->bTime = ((timFlag == 1) ? TRUE : FALSE);
  fr->contents->bAtoms = ((atmFlag == 1) ? TRUE : FALSE);
  fr->contents->bPrec = ((prcFlag == 1) ? TRUE : FALSE);
  fr->contents->bBox = ((boxFlag == 1) ? TRUE : FALSE);
  fr->contents->bX = ((xFlag == 1) ? TRUE : FALSE);
  fr->contents->bV = ((vFlag == 1) ? TRUE : FALSE);
  fr->contents->bF = ((fFlag == 1) ? TRUE : FALSE);*/


  if (fr->contents->bTitle) 
  { 
    get_next_line(fr->fp,inp_line);
    test_line(inp_line,"[title]",TRUE,"ERROR: title flag is 1 but found no [title]\n"); 
    strcpy(inp_word,"N/A");
    test_sscanf = sscanf(inp_line,"[title] %s ",&inp_word);
    if (test_sscanf != 1) 
    { 
      fprintf(stderr,"ERROR: Expected [title] <TITLE> in line\n");
      fprintf(stderr,"\tline:%s\ttitle: %s\n",inp_line,inp_word);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
      exit(1); 
    }
    if (! fr->contents->title) { fr->contents->title = (char *) ecalloc(100,sizeof(char)); }
    strcpy(fr->contents->title,inp_word);
  }

  if (fr->contents->bStep)
  {
    get_next_line(fr->fp,inp_line);
    test_line(inp_line,"[step]",TRUE,"ERROR: step flag is 1 but found no [step]\n");
    step = -1;
    test_sscanf = sscanf(inp_line,"[step] %d ",&step);
    if (test_sscanf != 1) 
    { 
      fprintf(stderr,"ERROR: Expected [step] <STEP_NUMBER> in line\n");
      fprintf(stderr,"\tline: %s\tstep: %d\n",inp_line,step);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
      exit(1); 
    }
    fr->contents->step = step;
  }
  
  if (fr->contents->bTime)
  {
    get_next_line(fr->fp,inp_line);
    test_line(inp_line,"[time]",TRUE,"ERROR: time flag is 1 but found no [time]\n");
    test_sscanf = sscanf(inp_line,"[time] %e ",&time);
    if (test_sscanf != 1) 
    { 
      fprintf(stderr,"ERROR: Expected [time] <TIME> in line\n");
      fprintf(stderr,"\tline: %s\ttime: %f\n",inp_line,time);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
      exit(1); 
    }
    fr->contents->time = time;
  }

  if (fr->contents->bPrec)
  {
    get_next_line(fr->fp,inp_line);
    test_line(inp_line,"[prec]",TRUE,"ERROR: prec flag is 1 but found no [prec]\n");
    test_sscanf = sscanf(inp_line,"[prec] %g ",&time);
    if (test_sscanf != 1) 
    { 
      fprintf(stderr,"ERROR: Expected [prec] <PREC> in line\n");
      fprintf(stderr,"\tline: %s\tprec: %f\n",inp_line,time);
      fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
      exit(1); 
    }
    fr->contents->prec = time;
  } 

  if (fr->contents->bBox)
  {
    get_next_line(fr->fp,inp_line);
    test_line(inp_line,"[box]",TRUE,"ERROR: box flag is 1 but found no [box]\n");
    for (i = 0; i < DIM; ++i)
    {
      get_next_line(fr->fp,inp_line);
      b[0] = b[1] = b[2] = -1.0;
      test_sscanf = sscanf(inp_line,"%f %f %f ",&b[0],&b[1],&b[2]);
      if (test_sscanf != 3) 
      { 
	fprintf(stderr,"ERROR: Expected 3 elements of box matrix in line\n");
	fprintf(stderr,"\tline: %s",inp_line);
	fprintf(stderr,"\tbox[%d][0]: %f\n\tbox[%d][1]: %f\n\tbox[%d][2]: %f\n",i,b[0],i,b[1],i,b[2]);
	fprintf(stderr,"ERROR: %s %d\n",__FILE__,__LINE__); 
	exit(1); 
      }
      fr->contents->box[i][0] = (double) b[0];
      fr->contents->box[i][1] = (double) b[1];
      fr->contents->box[i][2] = (double) b[2];
    }
  }

  if (fr->contents->bX || fr->contents->bV || fr->contents->bF)
  {
    char *dir = ecalloc(6,sizeof(char));
    sprintf(dir,"[%c%c%c]",(fr->contents->bX ? 'x' : '_'),(fr->contents->bV ? 'v' : '_'),(fr->contents->bF ? 'f' : '_'));
    get_next_line(fr->fp,inp_line);
    test_line(inp_line,dir,TRUE,"ERROR: wrong props [xvf]\n");

    for (i = 0; i < natoms; ++i)
    {
      get_next_line(fr->fp,inp_line);
      test_sscanf = sscanf(inp_line,"%s %s %s %f %f %f %s %s %s ",&xs[0],&xs[1],&xs[2],&v[0],&v[1],&v[2],&fs[0],&fs[1],&fs[2]);
      for (j = 0; j < DIM; ++j)
      {
        x[j] = atof(xs[j]);
        f[j] = atof(fs[j]);
      }
      if (test_sscanf != 9) 
      { 
	fprintf(stderr,"ERROR: Expected to find 9 values for atom %d (rx, ry, rz, vx, vy, vz, fx, fy, fz)\n",i);
        fprintf(stderr,"\tline: %s",inp_line);
	fprintf(stderr,"\tx[0]: %g\n\tx[1]: %g\n\tx[2]: %g\n",x[0],x[1],x[2]);
        fprintf(stderr,"\tv[0]: %f\n\tv[1]: %f\n\tv[2]: %f\n",v[0],v[1],v[2]);
        fprintf(stderr,"\tf[0]: %g\n\tf[1]: %g\n\tf[2]: %g\n",f[0],f[1],f[2]);
	fprintf(stderr,"ERROR: tssf: %d.  %s %d\n",test_sscanf,__FILE__,__LINE__); 
	exit(1); 
      }
      if (fr->contents->bX) { for (j = 0; j < DIM; ++j) { fr->contents->x[i][j] = (double) x[j]; } }
      if (fr->contents->bV) { for (j = 0; j < DIM; ++j) { fr->contents->v[i][j] = (double) v[j]; } }
      if (fr->contents->bF) { for (j = 0; j < DIM; ++j) { fr->contents->f[i][j] = (double) f[j]; } }
    } 
  }
  get_next_line(fr->fp,inp_line);
  test_line(inp_line,"[END_FRAME]",TRUE,"ERROR: Expected end of frame!\n");
  return TRUE;
}

/*
comp_2_frames(): I wrote this function to compare two frames (duh..). It was used
	early on when debugging. I'd read fr1 in from one file type, and fr2 in from
	a translated version of fr1 that was a different file type. I then check 
	to make sure they both have the same number of atoms. Then I take the difference
	between the first box matrix element, and the differences between all the xvf vectors.
	this all gets stored in fro.
*/
void comp_2_frames(tW_gmx_trxframe *fr1, tW_gmx_trxframe *fr2, tW_gmx_trxframe *fro)
{
  if (fr1->contents->natoms != fro->contents->natoms)
  {
    fprintf(stderr,"ERROR: in fr1 counter %d, natoms: %d, fro natoms: %d\n",fr1->counter,fr1->contents->natoms,fro->contents->natoms);
    exit(1);
  }
  if (fr2->contents->natoms != fro->contents->natoms)
  {
    fprintf(stderr,"ERROR: in fr2 counter %d, natoms: %d, fro natoms: %d\n",fr2->counter,fr2->contents->natoms,fro->contents->natoms);
    exit(1);
  }

  int i;

  if (fro->contents->bBox)
  {
    fro->contents->box[0][0] = fr2->contents->box[0][0] - fr1->contents->box[0][0];
  }


  if (fro->contents->bX)
  {
    for (i = 0; i < fro->contents->natoms; ++i)
    {
      vect_diff(fr2->contents->x[i], fr1->contents->x[i], fro->contents->x[i]);
    }
  }

  if (fro->contents->bV)
  {
    for (i = 0; i < fro->contents->natoms; ++i)
    {
      vect_diff(fr2->contents->v[i], fr1->contents->v[i], fro->contents->v[i]);
    }
  }
  
  if (fro->contents->bF)
  {
    for (i = 0; i < fro->contents->natoms; ++i)
    {
      vect_diff(fr2->contents->f[i], fr1->contents->f[i], fro->contents->f[i]);
    }
  }

}


/*
write_delta_frame(): This is the function called after comp_2_frames. It prints
	the differences in the box matrix element, as well as the difference of the 
	elements of all the xvf vectors to separate files. I then plotted these
	errors. Originally (before I figured out how to read/write trr files), the
	translated file would have a different precision then the original file. I needed
	to make sure the errors behaved as such. This, and the previous function,
	probably don't need to be used again. But if they do, they're here.
*/
void write_delta_frame(tW_gmx_trxframe *fr, FILE *fpbox, FILE *fpx, FILE *fpv, FILE *fpf)
{
  int i;
  if (fr->contents->bBox)
  {
    fprintf(fpbox,"%g\n",fr->contents->box[0][0]);
  }  

  if (fr->contents->bX)
  {
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      fprintf(fpx,"%g\n%g\n%g\n",fr->contents->x[i][0],fr->contents->x[i][1],fr->contents->x[i][2]);
    }
  }

  if (fr->contents->bV)
  {
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      fprintf(fpv,"%g\n%g\n%g\n",fr->contents->v[i][0],fr->contents->v[i][1],fr->contents->v[i][2]);
    }
  }

  if (fr->contents->bF)
  {
    for (i = 0; i < fr->contents->natoms; ++i)
    {
      fprintf(fpf,"%g\n%g\n%g\n",fr->contents->f[i][0],fr->contents->f[i][1],fr->contents->f[i][2]);
    }
  }
}

/*****************************************************************************************
The following are functions that have been copied and trimmed down from GROMACS.
They are used for reading/writing binary files (either trj files or trr (xdr) files).
*****************************************************************************************/

static bool tW_do_binwrite(FILE *fp, const void *item, int nitem, int eio)
{
  size_t size = 0, wsize;
  int ssize;

  if (!fp)
  {
    fprintf(stderr,"ERROR: fp is not open!\n");
    fp = fopen("error_open_trr_file.trr","wb");
    if (!fp) { fprintf(stderr,"ERROR: fp is still not open!\n"); exit(1); }
  }

  if (eio == eioSTRING)
  {
    size = ssize = strlen((char *) item) + 1;
    tW_do_binwrite(fp, &ssize, 1, eioINT);
  }
  else if (eio == eioINT) { size = sizeof(int); }
  else if (eio == eioFLOAT) { size = sizeof(float);} 
  else if (eio == eioDOUBLE) { size = sizeof(double);} 
  else if (eio == eioREAL) { size = sizeof(tW_real);} 
  else if (eio = eioRVEC) { size = sizeof(tW_rvec);}
  wsize = fwrite(item,size,nitem,fp);
  if (wsize != nitem)
  {
    fprintf(stderr,"ERROR: wsize == %d != %d == nitem\n",wsize,nitem);
  }
  return (wsize==nitem);
}


static bool tW_do_write_trjheader(tW_gmx_trxframe *fr)
{
  int magic = GROMACS_MAGIC;
  bool bOK = TRUE;

  bOK = tW_do_binwrite(fr->fp, &magic, 1, eioINT);
  if (!bOK) { fprintf(stderr,"ERROR: doesn't like binwrite magic\n"); exit(1); }

  char *buf = "GMX_trn_file";
  int zero = 0;
  tW_real rzero = 0;
  int box_size = (fr->contents->bBox ? sizeof(tW_matrix) : 0);
  int x_size = (fr->contents->bX ? fr->contents->natoms * sizeof(fr->contents->x[0]): 0);
  int v_size = (fr->contents->bV ? fr->contents->natoms * sizeof(fr->contents->v[0]): 0);
  int f_size = (fr->contents->bF ? fr->contents->natoms * sizeof(fr->contents->f[0]): 0);
  int natoms = fr->contents->natoms;
  int step = (fr->contents->bStep ? fr->contents->step : 0);
  tW_real lambda = (fr->contents->bLambda ? fr->contents->lambda : 0.0);
  tW_real time = (fr->contents->bTime ? fr->contents->time : 0.0);

  bOK = bOK && tW_do_binwrite(fr->fp, buf, 1, eioSTRING);

  bOK = bOK && tW_do_binwrite(fr->fp, &zero, 1, eioINT); // ir_size
  bOK = bOK && tW_do_binwrite(fr->fp, &zero, 1, eioINT); // e_size
  bOK = bOK && tW_do_binwrite(fr->fp, &box_size, 1, eioINT); // box_size
  bOK = bOK && tW_do_binwrite(fr->fp, &zero, 1, eioINT); // vir_size
  bOK = bOK && tW_do_binwrite(fr->fp, &zero, 1, eioINT); // pres_size
  bOK = bOK && tW_do_binwrite(fr->fp, &zero, 1, eioINT); // top_size
  bOK = bOK && tW_do_binwrite(fr->fp, &zero, 1, eioINT); // sym_size
  bOK = bOK && tW_do_binwrite(fr->fp, &x_size, 1, eioINT); // x_size
  bOK = bOK && tW_do_binwrite(fr->fp, &v_size, 1, eioINT); // v_size
  bOK = bOK && tW_do_binwrite(fr->fp, &f_size, 1, eioINT); // f_size
  bOK = bOK && tW_do_binwrite(fr->fp, &natoms, 1, eioINT); // natoms

  bOK = bOK && tW_do_binwrite(fr->fp, &step, 1, eioINT); // step
  bOK = bOK && tW_do_binwrite(fr->fp, &zero, 1, eioINT); // nre

  bOK = bOK && tW_do_binwrite(fr->fp, (fr->contents->bTime ? &(fr->contents->time) : &rzero), 1, eioREAL); // t
  bOK = bOK && tW_do_binwrite(fr->fp, (fr->contents->bLambda ? &(fr->contents->lambda) : &rzero), 1, eioREAL); // lambda


  return bOK;
}

static bool tW_do_write_trjstuff(tW_gmx_trxframe *fr)
{
  bool bOK = TRUE;
  
  if ( fr->contents->bBox ) bOK = bOK && tW_do_binwrite(fr->fp, fr->contents->box, DIM, eioRVEC);
  if ( fr->contents->bX   ) bOK = bOK && tW_do_binwrite(fr->fp, fr->contents->x, fr->contents->natoms, eioRVEC);
  if ( fr->contents->bV   ) bOK = bOK && tW_do_binwrite(fr->fp, fr->contents->v, fr->contents->natoms, eioRVEC);
  if ( fr->contents->bF   ) bOK = bOK && tW_do_binwrite(fr->fp, fr->contents->f, fr->contents->natoms, eioRVEC);

  return bOK;
}

bool tW_write_trj_frame(tW_gmx_trxframe *fr)
{
  bool bOK;

  bOK = tW_do_write_trjheader(fr);
  bOK = bOK && tW_do_write_trjstuff(fr);

  return bOK;
}

void wrap_tW_write_trj_frame(tW_gmx_trxframe *fr) // discards the bool return from tW_write_trr_frame
{
  tW_write_trj_frame(fr);
}

static bool tW_do_binread(FILE *fp, void *item, int nitem, int eio)
{
  size_t size=0, rsize;
  int ssize;
  if (eio == eioINT) { size = sizeof(int); }
  else if (eio == eioFLOAT) { size = sizeof(float); }
  else if (eio == eioDOUBLE) { size = sizeof(double); }
  else if (eio == eioREAL) { size = sizeof(tW_real); }
  else if (eio == eioRVEC) { size = DIM * sizeof(tW_real); }
  else if (eio == eioSTRING) { tW_do_binread(fp, &ssize, 1, eioINT); size = ssize; }
  
  if (item) { rsize = fread(item,size,nitem,fp); }
  else { fseek(fp, size*nitem, SEEK_CUR); rsize = nitem;}

  return (rsize == nitem);
}

static bool tW_do_read_trjheader(tW_gmx_trxframe *fr)
{
  int magic=GROMACS_MAGIC;
  int magic_test;
  char buf[256];
  bool bOK = TRUE;
  int ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size, x_size, v_size, f_size, natoms, step, nre;
  tW_real t, lambda;


  bOK = bOK && tW_do_binread(fr->fp, &magic_test, 1, eioINT); // magic

  bOK = bOK && tW_do_binread(fr->fp, buf, 1, eioSTRING); // buf


  bOK = bOK && tW_do_binread(fr->fp, &ir_size, 1, eioINT); // ir_size
  bOK = bOK && tW_do_binread(fr->fp, &e_size, 1, eioINT); // e_size
  bOK = bOK && tW_do_binread(fr->fp, &box_size, 1, eioINT); // box_size
  bOK = bOK && tW_do_binread(fr->fp, &vir_size, 1, eioINT); // vir_size
  bOK = bOK && tW_do_binread(fr->fp, &pres_size, 1, eioINT); // pres_size
  bOK = bOK && tW_do_binread(fr->fp, &top_size, 1, eioINT); // top_size
  bOK = bOK && tW_do_binread(fr->fp, &sym_size, 1, eioINT); // sym_size
  bOK = bOK && tW_do_binread(fr->fp, &x_size, 1, eioINT); // x_size
  bOK = bOK && tW_do_binread(fr->fp, &v_size, 1, eioINT); // v_size
  bOK = bOK && tW_do_binread(fr->fp, &f_size, 1, eioINT); // f_size
  bOK = bOK && tW_do_binread(fr->fp, &(fr->contents->natoms), 1, eioINT); // natoms

  fr->contents->bBox = ( box_size == 0 ? FALSE : TRUE ) ;
  fr->contents->bX = ( x_size == 0 ? FALSE : TRUE ) ;
  fr->contents->bV = ( v_size == 0 ? FALSE : TRUE ) ;
  fr->contents->bF = ( f_size == 0 ? FALSE : TRUE ) ; 

  if (!bOK) { return bOK; }

  bOK = bOK && tW_do_binread(fr->fp, &step, 1, eioINT); // step
  bOK = bOK && tW_do_binread(fr->fp, &nre, 1, eioINT); // nre
  bOK = bOK && tW_do_binread(fr->fp, &t, 1, eioREAL); // t
  bOK = bOK && tW_do_binread(fr->fp, &lambda, 1, eioREAL); // lambda

  fr->contents->bStep = ( step == 0 ? FALSE : TRUE ) ;
  if (fr->contents->bStep) { fr->contents->step = step; }
  fr->contents->bTime = ( t == 0.0 ? FALSE : TRUE ) ;
  if (fr->contents->bTime) { fr->contents->time = t; }
  fr->contents->bLambda = ( lambda == 0.0 ? FALSE : TRUE ) ;
  if (fr->contents->bLambda) { fr->contents->lambda = lambda; }

  return bOK;
}

static bool tW_do_read_trjstuff(tW_gmx_trxframe *fr)
{
  bool bOK = TRUE;

  if (fr->contents->bBox) { bOK = bOK && tW_do_binread(fr->fp,fr->contents->box,DIM,eioRVEC); }
  fprintf(stderr,"box[0][0]: %f %g %9.7e\n",fr->contents->box[0][0],fr->contents->box[0][0],fr->contents->box[0][0]);
  if (fr->contents->bX) { bOK = bOK && tW_do_binread(fr->fp,fr->contents->x,fr->contents->natoms, eioRVEC); }
  if (fr->contents->bV) { bOK = bOK && tW_do_binread(fr->fp,fr->contents->v,fr->contents->natoms, eioRVEC); }
  if (fr->contents->bF) { bOK = bOK && tW_do_binread(fr->fp,fr->contents->f,fr->contents->natoms, eioRVEC); }

  return bOK;
}

bool tW_read_first_trj_frame(tW_gmx_trxframe *fr, const char *trx_fnm)
{
  if (!fr->fp) { fr->fp = fopen(trx_fnm,"rb"); }

if (!fr->fp) { fprintf(stderr,"ERROR: unable to open file: %s\n",trx_fnm); exit(1); }

  bool bOK = tW_do_read_trjheader(fr);

  set_natoms(fr, fr->contents->natoms);

  rewind(fr->fp);

  return bOK;
}

bool tW_read_next_trj_frame(tW_gmx_trxframe *fr)
{
  bool bOK;
  bOK = tW_do_read_trjheader(fr);
  bOK = bOK && tW_do_read_trjstuff(fr);
  if (!bOK) { fprintf(stderr,"Done after frame %d\n",fr->counter - 1); }
  return bOK;  
}

/*****************************************************************************************
Functions analogous to the previous set, but this time for trr files, not trj files
*****************************************************************************************/

static bool tW_do_write_xdr(tW_gmx_trxframe *fr, const void *item, int nitem, int eio)
{
  unsigned char ucdum, *ucptr;
  bool		res = 0;
  float		fvec[DIM];
  double 	dvec[DIM];
  int		j, m, *iptr, idum;
  tW_real	*ptr;
  unsigned short us;
  double d = 0;
  float f = 0;
  char *cptr;
  int slen;
  switch (eio)
  {
    case eioREAL:
      if ( sizeof(tW_real) == sizeof(double) )
      {
	if (item) { d = *((tW_real *) item); }
        res = xdr_double(fr->xdr,&d);
        if (item) { *((tW_real *) item) = d; }
      }
      else
      {
	if (item) { f = *((tW_real *) item); }
	res = xdr_float(fr->xdr,&f);
	if (item) { *((tW_real *) item) = f; }
      }
      break;
    case eioFLOAT:
      if (item) { f = *((float *) item); }
      res = xdr_float(fr->xdr, &f);
      if (item) { *((float *) item) = f; }
      break;
    case eioDOUBLE:
      if (item) { d = *((double *) item); }
      res = xdr_double(fr->xdr, &d);
      if (item) { *((double *) item); }
      break;
    case eioINT:
      if (item) { idum = *(int *) item; }
      res = xdr_int(fr->xdr, &idum); 
      if (item) { *(int *) item = idum; }
      break;
    case eioRVEC:
      if ( sizeof(tW_real) == sizeof(double) )
      {
	if (item) { for (m = 0; m < DIM; ++m) { dvec[m] = ((tW_real *) item)[m]; } }
	res = xdr_vector(fr->xdr, (char *) dvec, DIM, (unsigned int) sizeof(double), (xdrproc_t) xdr_double);
	if (item) { for (m = 0; m < DIM; ++m) { ((tW_real *) item)[m] = dvec[m]; } }
      }
      else
      {
	if (item) { for (m = 0; m < DIM; ++m) { fvec[m] = ((tW_real *) item)[m]; } }
	res = xdr_vector(fr->xdr, (char *) fvec, DIM, (unsigned int) sizeof(float), (xdrproc_t) xdr_float);
	if (item) { for (m = 0; m < DIM; ++m) { ((tW_real *) item)[m] = fvec[m]; } }
      }
      break;
    case eioNRVEC:
      ptr = NULL;
      res = 1;
      for (j = 0; (j < nitem) && res; j++)
      {
        if (item) { ptr = ((tW_rvec *) item)[j];}
        res = tW_do_write_xdr(fr,ptr,1,eioRVEC);
      }
      break;
    case eioSTRING:
      if (item) { slen = strlen((char *) item) + 1; }
      else { slen = 0; }
      if (xdr_int(fr->xdr, &slen) <= 0) 
      { 
	fprintf(stderr,"ERROR: wrong string length %d for string %s\n",slen, (char *)item);
        fprintf(stderr,"FILE: %s   LINE: %d\n",__FILE__,__LINE__); 
      }
      cptr = (char *) item;
      res = xdr_string(fr->xdr, &cptr, slen);
      break;        
  }  
  return (res != 0);
}

static bool tW_do_write_trrheader(tW_gmx_trxframe *fr)
{
  int magic = GROMACS_MAGIC;
  bool bOK = TRUE;

  bOK = tW_do_write_xdr(fr, &magic, 1, eioINT);
  if (!bOK) { fprintf(stderr,"ERROR: doesn't like binwrite magic\n"); exit(1); }

  char *buf = "GMX_trn_file";
  int zero = 0; 
  tW_real rzero = 0;
  int box_size = (fr->contents->bBox ? sizeof(tW_matrix) : 0);
  int x_size = (fr->contents->bX ? fr->contents->natoms * sizeof(fr->contents->x[0]): 0);
  int v_size = (fr->contents->bV ? fr->contents->natoms * sizeof(fr->contents->v[0]): 0);
  int f_size = (fr->contents->bF ? fr->contents->natoms * sizeof(fr->contents->f[0]): 0);
  int natoms = fr->contents->natoms;
  int step = (fr->contents->bStep ? fr->contents->step : 0);
  tW_real lambda = (fr->contents->bLambda ? fr->contents->lambda : 0.0);
  tW_real time = (fr->contents->bTime ? fr->contents->time : 0.0);

  bOK = bOK && tW_do_write_xdr(fr, buf, 1, eioSTRING);

  bOK = bOK && tW_do_write_xdr(fr, &zero, 1, eioINT); // ir_size
  bOK = bOK && tW_do_write_xdr(fr, &zero, 1, eioINT); // e_size
  bOK = bOK && tW_do_write_xdr(fr, &box_size, 1, eioINT); // box_size
  bOK = bOK && tW_do_write_xdr(fr, &zero, 1, eioINT); // vir_size
  bOK = bOK && tW_do_write_xdr(fr, &zero, 1, eioINT); // pres_size
  bOK = bOK && tW_do_write_xdr(fr, &zero, 1, eioINT); // top_size
  bOK = bOK && tW_do_write_xdr(fr, &zero, 1, eioINT); // sym_size
  bOK = bOK && tW_do_write_xdr(fr, &x_size, 1, eioINT); // x_size
  bOK = bOK && tW_do_write_xdr(fr, &v_size, 1, eioINT); // v_size
  bOK = bOK && tW_do_write_xdr(fr, &f_size, 1, eioINT); // f_size
  bOK = bOK && tW_do_write_xdr(fr, &natoms, 1, eioINT); // natoms

  bOK = bOK && tW_do_write_xdr(fr, &step, 1, eioINT); // step
  bOK = bOK && tW_do_write_xdr(fr, &zero, 1, eioINT); // nre

  bOK = bOK && tW_do_write_xdr(fr, (fr->contents->bTime ? &(fr->contents->time) : &rzero), 1, eioREAL); // t
  bOK = bOK && tW_do_write_xdr(fr, (fr->contents->bLambda ? &(fr->contents->lambda) : &rzero), 1, eioREAL); // lambda


  return bOK;
  
}

static bool tW_do_write_trrstuff(tW_gmx_trxframe *fr)
{
  bool bOK = TRUE;

  if ( fr->contents->bBox ) bOK = bOK && tW_do_write_xdr(fr, fr->contents->box, DIM, eioNRVEC);
  if ( fr->contents->bX   ) bOK = bOK && tW_do_write_xdr(fr, fr->contents->x, fr->contents->natoms, eioNRVEC);
  if ( fr->contents->bV   ) bOK = bOK && tW_do_write_xdr(fr, fr->contents->v, fr->contents->natoms, eioNRVEC);
  if ( fr->contents->bF   ) bOK = bOK && tW_do_write_xdr(fr, fr->contents->f, fr->contents->natoms, eioNRVEC);

  return bOK;

}

bool tW_write_trr_frame(tW_gmx_trxframe *fr)
{
  bool bOK;
  
  bOK = tW_do_write_trrheader(fr);
  bOK = bOK && tW_do_write_trrstuff(fr);
  return bOK;
}
void wrap_tW_write_trr_frame(tW_gmx_trxframe *fr)
{
  tW_write_trr_frame(fr);
}
/*********************************************************************************
 * tW_do_read_xdr
 * Attempts to encode/decode various data types from the xdr data stream
 * [the data stream represents our trajectory], stores results in item
 * @param fr trajectory, not human readable
 * @param item stores decode results 
 * @return TRUE if encoding/decoding was successful FALSE if not
 *********************************************************************************/
static bool tW_do_read_xdr(tW_gmx_trxframe *fr, void *item, int nitem, int eio)
{
  unsigned char   ucdum, *ucptr;
  bool            res = 0; //res = result, represents success or failure
  float           fvec[DIM];
  double          dvec[DIM];
  int             j, m, *iptr, idum;
  tW_real         *ptr;
  unsigned short  us;
  double          d = 0;
  float           f = 0;
  char *cptr;
  int slen = 0;

  switch (eio)
  {
    case eioREAL:
/*      if ( sizeof(tW_real) == sizeof(double))
      {
        res = xdr_double(fr->xdr, &d);
        if (item) { *((tW_real *) item) = d; }
      }
      else
      {
        res = xdr_float(fr->xdr, &f);
        if (item) { *((tW_real *) item) = f; }
      }*/
      if (fr->bDouble)
      {
	res = xdr_double(fr->xdr, &d);
	if (item) { *((tW_real *) item) = d; }
      }
      else
      {
	res = xdr_float(fr->xdr, &f);
	if (item) { *((tW_real *) item) = (tW_real) f; }	
      }
      break;
    case eioFLOAT:
      res = xdr_float(fr->xdr, &f);
      if (item) { *((float *) item) = f; }
      break;
    case eioDOUBLE:
      res = xdr_double(fr->xdr, &d);
      if (item) { *((double *) item) = d; }
      break;
    case eioINT:
      res = xdr_int(fr->xdr, &idum);
      if (item) { *(int *) item = idum; }
      break;
    case eioRVEC:
      if ( sizeof(tW_real) == sizeof(double) )
      {
        if (fr->bDouble) // If detected double precision based on header
	{
	  res = xdr_vector(fr->xdr, (char *) dvec, DIM, (unsigned int) sizeof(double), (xdrproc_t) xdr_double);
	  if (item) { for (m = 0; m < DIM; ++m) { ((tW_real *) item)[m] = dvec[m]; } }
	}
	else // If detected single precision based on header
	{
	  res = xdr_vector(fr->xdr, (char *) fvec, DIM, (unsigned int) sizeof(float), (xdrproc_t) xdr_float);
	  if (item) { for (m = 0; m < DIM; ++m) { ((tW_real *) item)[m] = (tW_real) fvec[m]; } }
	}
      }
      else // This should never be called, tW_real is always double
      {
	res = xdr_vector(fr->xdr, (char *) fvec, DIM, (unsigned int) sizeof(float), (xdrproc_t) xdr_float);
	if (item) { for (m = 0; m < DIM; ++m) { ((tW_real *) item)[m] = fvec[m]; } }
      }
      break;
    case eioNRVEC:
      ptr = NULL;
      res = 1;
      for (j = 0; (j < nitem) && res; j++)
      {
	if (item) { ptr = ((tW_rvec *) item)[j];}
	res = tW_do_read_xdr(fr,ptr,1,eioRVEC);
      }
      break;
    case eioSTRING:
      if (xdr_int(fr->xdr, &slen) <= 0) 
      { 
	fprintf(stderr,"ERROR: wrong string length %d for string %s\n",slen,(char *)item);
        fprintf(stderr,"FILE: %s   LINE: %d\n",__FILE__,__LINE__); 
      }
      if (!item) { cptr = (char *) ecalloc(slen,sizeof(char)); }
      else { cptr = (char *) item; }
      if (cptr) { res = xdr_string(fr->xdr, &cptr, slen); }
      else { res = 1; }
      if (!item) { efree(cptr) ; }
      break;
  } 
 
  return (res != 0);
}

static bool get_trr_precision(tW_gmx_trxframe *fr)
{
  int magic=GROMACS_MAGIC;
  int magic_test;
  char buf[256];
  bool bOK = TRUE;
  bool bDouble;
  int nflsize_box = 0, nflsize_x = 0, nflsize_v = 0, nflsize_f = 0;
  int ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size, x_size, v_size, f_size, natoms, step, nre;
  tW_real t, lambda;


  bOK = bOK && tW_do_read_xdr(fr, &magic_test, 1, eioINT); // magic
  
  bOK = bOK && tW_do_read_xdr(fr, buf, 1, eioSTRING); // buf


  bOK = bOK && tW_do_read_xdr(fr, &ir_size, 1, eioINT); // ir_size
  bOK = bOK && tW_do_read_xdr(fr, &e_size, 1, eioINT); // e_size
  bOK = bOK && tW_do_read_xdr(fr, &box_size, 1, eioINT); // box_size
  bOK = bOK && tW_do_read_xdr(fr, &vir_size, 1, eioINT); // vir_size
  bOK = bOK && tW_do_read_xdr(fr, &pres_size, 1, eioINT); // pres_size
  bOK = bOK && tW_do_read_xdr(fr, &top_size, 1, eioINT); // top_size
  bOK = bOK && tW_do_read_xdr(fr, &sym_size, 1, eioINT); // sym_size
  bOK = bOK && tW_do_read_xdr(fr, &x_size, 1, eioINT); // x_size
  bOK = bOK && tW_do_read_xdr(fr, &v_size, 1, eioINT); // v_size
  bOK = bOK && tW_do_read_xdr(fr, &f_size, 1, eioINT); // f_size
  bOK = bOK && tW_do_read_xdr(fr, &(fr->contents->natoms), 1, eioINT); // natoms

//  if (!bOK) { return bOK; }

  bOK = bOK && tW_do_read_xdr(fr, &step, 1, eioINT); // step
  bOK = bOK && tW_do_read_xdr(fr, &nre, 1, eioINT); // nre
  bOK = bOK && tW_do_read_xdr(fr, &t, 1, eioREAL); // t 
  bOK = bOK && tW_do_read_xdr(fr, &lambda, 1, eioREAL); // lambda

  if (!bOK) // then at least one decoding got messed up
  {
    fprintf(stderr,"ERROR: unable to determine precision of trr file\n");
    exit(1);
  }

  nflsize_box = box_size / (DIM*DIM);
  nflsize_x = x_size / (fr->contents->natoms * DIM);
  nflsize_v = v_size / (fr->contents->natoms * DIM);
  nflsize_f = f_size / (fr->contents->natoms * DIM);
  
  if (nflsize_box != 0)
  {
    if (nflsize_box == sizeof(double)) { bDouble = TRUE; }
    else 
    { 
      bDouble = FALSE; 
      if (nflsize_box != sizeof(float))
      {
	fprintf(stderr,"ERROR: nflsize_box: %d   sizeof(double): %d   sizeof(float): %d\n",nflsize_box,sizeof(double),sizeof(float));
	fprintf(stderr,"\tUnable to determine precision of trr file: %s\n",fr->filename);
	exit(1);
      }
    }
  }
  else if (nflsize_x != 0)
  {
    if (nflsize_x == sizeof(double)) { bDouble = TRUE; }
    else
    {
      bDouble = FALSE; 
      if (nflsize_x != sizeof(float))
      { 
        fprintf(stderr,"ERROR: nflsize_x: %d   sizeof(double): %d   sizeof(float): %d\n",nflsize_x,sizeof(double),sizeof(float));
        fprintf(stderr,"\tUnable to determine precision of trr file: %s\n",fr->filename);
        exit(1);
      }
    }
  }
  else if (nflsize_v != 0)
  {
    if (nflsize_v == sizeof(double)) { bDouble = TRUE; }
    else
    {
      bDouble = FALSE; 
      if (nflsize_v != sizeof(float))
      { 
        fprintf(stderr,"ERROR: nflsize_v: %d   sizeof(double): %d   sizeof(float): %d\n",nflsize_v,sizeof(double),sizeof(float));
        fprintf(stderr,"\tUnable to determine precision of trr file: %s\n",fr->filename);
        exit(1);
      }
    }
  }
  else if (nflsize_f != 0)
  {
    if (nflsize_f == sizeof(double)) { bDouble = TRUE; }
    else
    {
      bDouble = FALSE; 
      if (nflsize_f != sizeof(float))
      { 
        fprintf(stderr,"ERROR: nflsize_f: %d   sizeof(double): %d   sizeof(float): %d\n",nflsize_f,sizeof(double),sizeof(float));
        fprintf(stderr,"\tUnable to determine precision of trr file: %s\n",fr->filename);
        exit(1);
      }
    }
  }

  return bDouble;
}

static bool tW_do_read_trrheader(tW_gmx_trxframe *fr)
{
  int magic=GROMACS_MAGIC;
  int magic_test;
  char buf[256];
  bool bOK = TRUE;
  int ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size, x_size, v_size, f_size, natoms, step, nre;
  tW_real t, lambda;


  bOK = bOK && tW_do_read_xdr(fr, &magic_test, 1, eioINT); // magic

  bOK = bOK && tW_do_read_xdr(fr, buf, 1, eioSTRING); // buf


  bOK = bOK && tW_do_read_xdr(fr, &ir_size, 1, eioINT); // ir_size
  bOK = bOK && tW_do_read_xdr(fr, &e_size, 1, eioINT); // e_size
  bOK = bOK && tW_do_read_xdr(fr, &box_size, 1, eioINT); // box_size
  bOK = bOK && tW_do_read_xdr(fr, &vir_size, 1, eioINT); // vir_size
  bOK = bOK && tW_do_read_xdr(fr, &pres_size, 1, eioINT); // pres_size
  bOK = bOK && tW_do_read_xdr(fr, &top_size, 1, eioINT); // top_size
  bOK = bOK && tW_do_read_xdr(fr, &sym_size, 1, eioINT); // sym_size
  bOK = bOK && tW_do_read_xdr(fr, &x_size, 1, eioINT); // x_size
  bOK = bOK && tW_do_read_xdr(fr, &v_size, 1, eioINT); // v_size
  bOK = bOK && tW_do_read_xdr(fr, &f_size, 1, eioINT); // f_size
  bOK = bOK && tW_do_read_xdr(fr, &(fr->contents->natoms), 1, eioINT); // natoms

  fr->contents->bBox = ( box_size == 0 ? FALSE : TRUE ) ;
  fr->contents->bX = ( x_size == 0 ? FALSE : TRUE ) ;
  fr->contents->bV = ( v_size == 0 ? FALSE : TRUE ) ;
  fr->contents->bF = ( f_size == 0 ? FALSE : TRUE ) ;
  if (!bOK) { return bOK; }

  bOK = bOK && tW_do_read_xdr(fr, &step, 1, eioINT); // step
  bOK = bOK && tW_do_read_xdr(fr, &nre, 1, eioINT); // nre
  bOK = bOK && tW_do_read_xdr(fr, &t, 1, eioREAL); // t 
  bOK = bOK && tW_do_read_xdr(fr, &lambda, 1, eioREAL); // lambda

  fr->contents->bStep = ( step == 0 ? FALSE : TRUE ) ;
  if (fr->contents->bStep) { fr->contents->step = step; }
  fr->contents->bTime = ( t == 0.0 ? FALSE : TRUE ) ;
  if (fr->contents->bTime) { fr->contents->time = t; }
  fr->contents->bLambda = ( lambda == 0.0 ? FALSE : TRUE ) ;
  if (fr->contents->bLambda) { fr->contents->lambda = lambda; }

  return bOK;
}
static bool tW_do_read_trrstuff(tW_gmx_trxframe *fr)
{
  bool bOK = TRUE;
  if (fr->bDouble)
  {
    if (fr->contents->bBox) { bOK = bOK && tW_do_read_xdr(fr,fr->contents->box,DIM,eioNRVEC); }
    if (fr->contents->bX) { bOK = bOK && tW_do_read_xdr(fr,fr->contents->x,fr->contents->natoms, eioNRVEC); }
    if (fr->contents->bV) { bOK = bOK && tW_do_read_xdr(fr,fr->contents->v,fr->contents->natoms, eioNRVEC); }
    if (fr->contents->bF) { bOK = bOK && tW_do_read_xdr(fr,fr->contents->f,fr->contents->natoms, eioNRVEC); }
  }
  else
  {
    if (fr->contents->bBox) { bOK = bOK && tW_do_read_xdr(fr,fr->contents->box,DIM,eioNRVEC); }
    if (fr->contents->bX) { bOK = bOK && tW_do_read_xdr(fr,fr->contents->x,fr->contents->natoms, eioNRVEC); }
    if (fr->contents->bV) { bOK = bOK && tW_do_read_xdr(fr,fr->contents->v,fr->contents->natoms, eioNRVEC); }
    if (fr->contents->bF) { bOK = bOK && tW_do_read_xdr(fr,fr->contents->f,fr->contents->natoms, eioNRVEC); }
  }
  return bOK;
}
bool tW_read_first_trr_frame(tW_gmx_trxframe *fr, const char *trx_fnm)
{
  bool bOK;
  fr->setup_xdr(fr,trx_fnm,TRUE); // link xdr object to file stream and determine trr precision
  bOK = tW_do_read_trrheader(fr);
  set_natoms(fr, fr->contents->natoms);
  rewind(fr->fp); 
  return bOK;
}
bool tW_read_next_trr_frame(tW_gmx_trxframe *fr)
{
  bool bOK;
  bOK = tW_do_read_trrheader(fr);
  bOK = bOK && tW_do_read_trrstuff(fr);
  if (!bOK) { fprintf(stderr,"Done after frame %d\n",fr->counter - 1); }
  return bOK;
}

/**********************************************************************
read_first_lammps_frame(): reads the first frame from a lammps custom dump file.
	that dump file MUST contain each atom's ID (since LAMMPS doesn't always
	print atoms in the same order from frame to frame, but we need to keep
	atoms straight for topology identification), and should contain
	x, y, z, (vx, vy, vz, fx, fy, fz) depending on what you're gonna do
	with the trajectory.
**********************************************************************/

int read_first_lammps_frame(tW_gmx_trxframe *fr, const char *fnm)
{
  int test_sscanf, n_atoms, n_words;
  tW_line inp_line;
  int i, n_LDs = 0, n_LDGs = 0;

  fr->fp = open_file(fnm,'r');

  get_next_line(fr->fp, inp_line); // ITEM: TIMESTEP
  get_next_line(fr->fp, inp_line); // #
  get_next_line(fr->fp, inp_line); // ITEM: NUMBER OF ATOMS
  get_next_line(fr->fp, inp_line); // N_ATOMS
  test_sscanf = sscanf(inp_line," %d ",&n_atoms);
  if (test_sscanf != 1) 
  { 
    fprintf(stderr,"ERROR: unable to read natoms from line 4\n");
    fprintf(stderr,"%s",inp_line);
    exit(1);
  } 
  
  get_next_line(fr->fp, inp_line); // ITEM: BOX BOUNDS pp pp pp
  get_next_line(fr->fp, inp_line); // xlo xhi
  get_next_line(fr->fp, inp_line); // ylo yhi
  get_next_line(fr->fp, inp_line); // zlo zhi

  get_next_line(fr->fp, inp_line); // ITEM: ATOMS xx xx xx xx xx xx xx xx

  elim_char(inp_line,'\n');
  n_words = get_word_count_delim(inp_line," ");
  tW_word *word_list = (tW_word *) ecalloc(n_words,sizeof(tW_word));
  get_words_delim(inp_line," ",word_list);  

  fr->xids = (int *) ecalloc(DIM,sizeof(int));
  fr->vids = (int *) ecalloc(DIM,sizeof(int));
  fr->fids = (int *) ecalloc(DIM,sizeof(int)); 

  unsigned short indices = 0;

  for (i = 0; i < n_words; ++i)
  {
    if (strcmp(word_list[i],"x") == 0) { fr->xids[0] = i-2; indices |= (1 << 0); }
    else if (strcmp(word_list[i],"y") == 0) { fr->xids[1] = i-2; indices |= (1 << 1); }
    else if (strcmp(word_list[i],"z") == 0) { fr->xids[2] = i-2; indices |= (1 << 2); }
    else if (strcmp(word_list[i],"vx") == 0) { fr->vids[0] = i-2; indices |= (1 << 3); }
    else if (strcmp(word_list[i],"vy") == 0) { fr->vids[1] = i-2; indices |= (1 << 4); }
    else if (strcmp(word_list[i],"vz") == 0) { fr->vids[2] = i-2; indices |= (1 << 5); }
    else if (strcmp(word_list[i],"fx") == 0) { fr->fids[0] = i-2; indices |= (1 << 6); }
    else if (strcmp(word_list[i],"fy") == 0) { fr->fids[1] = i-2; indices |= (1 << 7); }
    else if (strcmp(word_list[i],"fz") == 0) { fr->fids[2] = i-2; indices |= (1 << 8); }
    else if (strcmp(word_list[i],"id") == 0) { fr->atid_id = i-2; indices |= (1 << 9); }
    else if (strstr(word_list[i],"lddens") != NULL) { ++n_LDs; indices |= (1 << 10); }
    else if (strstr(word_list[i],"grad") != NULL) { ++n_LDGs; indices |= (1 << 11); }

//    else if (strcmp(word_list[i],"type") == 0) { fr->attype_id = i-2; indices |= (1 << 10); }
//    else if (strcmp(word_list[i],"mol") == 0) { fr->molid_id = i-2; indices |= (1 << 11); }
  }

  if (((indices & (1 << 0)) != 0) && ((indices & (1 << 1)) != 0) && ((indices & (1 << 2))) != 0)
  { fr->contents->bX = TRUE; }
  else { fr->contents->bX = FALSE; }
  if (((indices & (1 << 3)) != 0) && ((indices & (1 << 4)) != 0) && ((indices & (1 << 5))) != 0)
  { fr->contents->bV = TRUE; }
  else { fr->contents->bV = FALSE; }  
  if (((indices & (1 << 6)) != 0) && ((indices & (1 << 7)) != 0) && ((indices & (1 << 8))) != 0)
  { fr->contents->bF = TRUE; }
  else { fr->contents->bF = FALSE; } 

  if ((indices & ((1 << 9))) == 0) 
  {
    fprintf(stderr,"ERROR: Missing columns from lammps dump file.\n");
    if ((indices & (1 << 9)) == 0) { fprintf(stderr,"\tMissing id\n"); }
    exit(1);  
  } 

  /* Local Density 06.17.2020 MRD */
  fprintf(stderr,"n_LDs: %d \n",n_LDs);
  if (n_LDs > 0)
  {
    fr->bLDs = TRUE;
    fr->n_atomtypes = n_LDs;
    fr->ld_id = (int *) ecalloc(n_LDs+1,sizeof(int));
    int a_type;
    for (i = 0; i < n_words; ++i)
    {
      if (strstr(word_list[i],"lddens") != NULL)
      {
        test_sscanf = sscanf(word_list[i],"lddens%d",&(a_type));
        if (test_sscanf != 1)
        {
          fprintf(stderr,"LDDENS READ ERROR!\n");
        }
        fr->ld_id[a_type] = i-2;
        fprintf(stderr,"fr->ld_id[%d]: %d \n",a_type,fr->ld_id[a_type]);
      }
    }
    if (n_LDGs > 0)
    {
      char xyz;
      fr->bLD_Grads = TRUE;
      fr->ldg_id = (int *) ecalloc(n_LDGs+3,sizeof(int));
      for (i = 0; i < n_words; ++i)
      {
        if ((strstr(word_list[i],"grad") != NULL) && (strstr(word_list[i],"nrg") == NULL))
        {
          test_sscanf = sscanf(word_list[i],"grad%c%d",&(xyz),&(a_type));
          if (test_sscanf != 2)
          {
            fprintf(stderr,"LDGRAD READ ERROR! \n");
            fprintf(stderr,"line: %s",inp_line);
            fprintf(stderr,"word: %s\n",word_list[i]);
          }
          if (xyz == 'x') { fr->ldg_id[DIM*a_type] = i-2; }
          else if (xyz == 'y') { fr->ldg_id[DIM*a_type+1] = i-2; }
          else if (xyz == 'z') { fr->ldg_id[DIM*a_type+2] = i-2; }
          else
          {
            fprintf(stderr,"ERROR: did not understand word: %s \n",word_list[i]);
            fprintf(stderr,"Expected: gradx#, grady#, or gradz#\n");
            exit(1);
          }
        }
      }
    }
  }


// Generic LAMMPS stuff
  fr->contents->bBox = TRUE;
  
  fr->contents->bDouble = TRUE;
  fr->contents->bTime = TRUE;
  fr->contents->bStep = TRUE;
  fr->contents->bPrec = FALSE;
 
  rewind(fr->fp);

  set_natoms(fr,n_atoms);
  free(word_list);
  return n_atoms;
}

bool read_next_lammps_frame(tW_gmx_trxframe *fr)
{
  tW_line inp_line;
//  tW_word * word_list;
  int test_sscanf, time, n_atoms, n_words, at_idx, i, j;
  float blo, bhi;
  test_sscanf = get_next_line(fr->fp, inp_line); // ITEM: TIMESTEP
  if (test_sscanf == -1) { fprintf(stderr,"Done after frame %d\n",fr->counter - 1); return FALSE; }
  get_next_line(fr->fp, inp_line); // #
  test_sscanf = sscanf(inp_line," %d ",&time);
  if (test_sscanf != 1) { fprintf(stderr,"ERROR: unable to find int TIME listed on second line of frame\n%s",inp_line); exit(1);}
  fr->contents->time = (tW_real) (time * T_LMP2GRO); 
  fr->contents->step = (tW_real) (time);

  get_next_line(fr->fp, inp_line); // ITEM: NUMBER OF ATOMS
  get_next_line(fr->fp, inp_line); // N_ATOMS
  test_sscanf = sscanf(inp_line," %d ",&n_atoms);
  if (test_sscanf != 1) { fprintf(stderr,"ERROR: unable to find int N_ATOMS listed on fourth line of frame\n%s",inp_line); exit(1);}
  if (n_atoms != fr->contents->natoms)
  {
    fprintf(stderr,"ERROR: first frame has %d atoms, but frame %d has %d atoms!\n",fr->contents->natoms,fr->counter,n_atoms);
    exit(1);
  }

/* Get box info */
  get_next_line(fr->fp, inp_line); // ITEM: BOX BOUNDS pp pp pp
  for (i = 0; i < DIM; ++i) { for (j = 0; j < DIM; ++j) { fr->contents->box[i][j] = 0.0; } }
  get_next_line(fr->fp, inp_line); // xlo xhi
  test_sscanf = sscanf(inp_line," %f %f ",&blo, &bhi);
  if (test_sscanf != 2) { fprintf(stderr,"ERROR: unable to read xlo, xhi from sixth line of frame\n%s",inp_line); exit(1); }
  fr->contents->box[0][0] = (tW_real) ((bhi-blo) * X_LMP2GRO); 
  get_next_line(fr->fp, inp_line); // ylo yhi
  test_sscanf = sscanf(inp_line," %f %f ",&blo, &bhi);
  if (test_sscanf != 2) { fprintf(stderr,"ERROR: unable to read ylo, yhi from seventh line of frame\n%s",inp_line); exit(1); }
  fr->contents->box[1][1] = (tW_real) ((bhi-blo) * X_LMP2GRO); 
  get_next_line(fr->fp, inp_line); // zlo zhi
  test_sscanf = sscanf(inp_line," %f %f ",&blo, &bhi);
  if (test_sscanf != 2) { fprintf(stderr,"ERROR: unable to read zlo, zhi from eighth line of frame\n%s",inp_line); exit(1); }
  fr->contents->box[2][2] = (tW_real) ((bhi-blo) * X_LMP2GRO); 

/* Get atoms info */
  get_next_line(fr->fp, inp_line); // ITEM: ATOMS xx xx xx xx xx xx xx xx
  for (i=0; i < fr->contents->natoms; ++i)
  {
    get_next_line(fr->fp, inp_line); // xx xx xx xx xx xx xx xx
    n_words = get_word_count_delim(inp_line," ");
    tW_word *word_list = (tW_word *) ecalloc(n_words,sizeof(tW_word));
    get_words_delim(inp_line," ",word_list);
    at_idx = atoi(word_list[fr->atid_id]) - 1;
    if (fr->contents->bX)
    {
      for (j=0; j<DIM; ++j)
      {
        fr->contents->x[at_idx][j] = (tW_real) (atof(word_list[fr->xids[j]]) * X_LMP2GRO);
      }
    }
    if (fr->contents->bV)
    {
      for (j=0; j<DIM; ++j)
      {
        fr->contents->v[at_idx][j] = (tW_real) (atof(word_list[fr->vids[j]]) * V_LMP2GRO);
      }
    }
    if (fr->contents->bF)
    {
      for (j=0; j<DIM; ++j)
      {
        fr->contents->f[at_idx][j] = (tW_real) (atof(word_list[fr->fids[j]]) * F_LMP2GRO);
      }
    }
  /* Local Density */
    if (fr->bLDs)
    {
      for (j = 1; j <= fr->n_atomtypes; ++j)
      {
        fr->local_densities[at_idx][j-1] = (double) (atof(word_list[fr->ld_id[j]]) * LD_LMP2GRO);
      }
    }
  /* Local Density Gradients */
    if (fr->bLD_Grads)
    {
      for (j = 1; j <= fr->n_atomtypes; ++j)
      {
        fr->local_density_gradients[at_idx][j-1][0] = (double) (atof(word_list[fr->ldg_id[j*DIM]]) * LDG_LMP2GRO);
        fr->local_density_gradients[at_idx][j-1][1] = (double) (atof(word_list[fr->ldg_id[j*DIM+1]]) * LDG_LMP2GRO);
        fr->local_density_gradients[at_idx][j-1][2] = (double) (atof(word_list[fr->ldg_id[j*DIM+2]]) * LDG_LMP2GRO);
      }
    }
    free(word_list);
  }

  return TRUE;
}

/*
 write_lammps_frame(): this function writes a frame in a generic lammps dump text format
*/
void write_lammps_frame(tW_gmx_trxframe *fr, tW_gmx_topology *top)
{
  int i;
  fprintf(fr->fp,"ITEM: TIMESTEP\n%d\n",fr->counter);
  fprintf(fr->fp,"ITEM: NUMBER OF ATOMS\n%d\n",fr->contents->natoms);
  switch (fr->contents->ePBC)
  {
    case (epbcXYZ):
      fprintf(fr->fp,"ITEM: BOX BOUNDS pp pp pp\n");
      break;
    case (epbcNONE):
      fprintf(fr->fp,"ITEM: BOX BOUNDS ff ff ff\n");
      break;
    case (epbcXY):
      fprintf(fr->fp,"ITEM: BOX BOUNDS pp pp ff\n");
      break;
    default:
      fprintf(stderr,"ERROR: unsure of what ePBC: %d is\n",fr->contents->ePBC);
      fprintf(stderr,"\tepbcXYZ: %d   epbcNONE: %d   epbcXY: %d   epbcSCREW: %d   epbcNR: %d\n",epbcXYZ,epbcNONE,epbcXY,epbcSCREW,epbcNR);
      exit(1);
  }
  fprintf(fr->fp,"%g %g\n%g %g\n%g %g\n",0.0,(fr->contents->box[0][0] / X_LMP2GRO),
                                       0.0,(fr->contents->box[1][1] / X_LMP2GRO),
                                       0.0,(fr->contents->box[2][2] / X_LMP2GRO));
  fprintf(fr->fp,"ITEM: ATOMS id type ");
  if (fr->contents->bX) { fprintf(fr->fp,"x y z "); }
  if (fr->contents->bV) { fprintf(fr->fp,"vx vy vz "); }
  if (fr->contents->bF) { fprintf(fr->fp,"fx fy fz "); }
  fprintf(fr->fp,"\n");

  for (i = 0; i < fr->contents->natoms; ++i)
  {
    fprintf(fr->fp,"%d %d ",i+1,get_type(*(top->contents->atoms.atomtype[i]),top)+1);
    if (fr->contents->bX) { fprintf(fr->fp,"%16.12f %16.12f %16.12f ",fr->contents->x[i][0] / X_LMP2GRO,
                                   fr->contents->x[i][1] / X_LMP2GRO, fr->contents->x[i][2] / X_LMP2GRO); }
    if (fr->contents->bV) { fprintf(fr->fp,"%g %g %g ",fr->contents->v[i][0] / V_LMP2GRO,
                                    fr->contents->v[i][1] / V_LMP2GRO, fr->contents->v[i][2] / V_LMP2GRO); }
    if (fr->contents->bF) { fprintf(fr->fp,"%16.12f %16.12f %16.12f ",fr->contents->f[i][0] / F_LMP2GRO,
                                   fr->contents->f[i][1] / F_LMP2GRO, fr->contents->f[i][2] / F_LMP2GRO); }
    fprintf(fr->fp,"\n");
  }

}


/*
write_lammps_data(): this function can write a lammps data file for reading in to
	a lammps simulation to initialize the configuration.
*/
void write_lammps_data(tW_gmx_trxframe *fr, tW_gmx_topology *top)
{
  int molt_idx, mol_idx, at_idx, abs_at_idx, n_prev_mol = 0, n_prev_at = 0;
  int n_bonds = 0, n_angles = 0, n_dih = 0, n_imp = 0;
  int n_att = 0, n_bt = 0, n_at = 0, n_dt = 0;
  int i, j, k;
  FILE *fp = fr->fp;
  abs_at_idx = 0;

  fprintf(fp,"LAMMPS Description\n\n");
  fprintf(fp,"     %d  atoms\n",fr->contents->natoms);
  for (i = 0; i < top->contents->mols.nr; ++i)
  {
    n_bonds += top->molecules[i].n_bonds * top->molecules[i].n_mols;
    n_angles += top->molecules[i].n_angles * top->molecules[i].n_mols;
    n_dih += top->molecules[i].n_dihs * top->molecules[i].n_mols;
    for (j = 0; j < top->molecules[i].n_apm; ++j)
    {
      if (top->molecules[i].type_ids[j] + 1 > n_att ) { n_att = top->molecules[i].type_ids[j] + 1 ; }
    }
  }

  top->int_map = (int *) ecalloc(top->contents->idef.ntypes,sizeof(int));
  
  for (i = 0; i < top->contents->idef.ntypes; ++i)
  {
    top->int_map[i] = 0;
    if (strcmp(top->force_names[i],"BONDS") == 0) { ++n_bt; top->int_map[i] = n_bt; }
    else if (strcmp(top->force_names[i],"ANGLES") == 0) { ++n_at; top->int_map[i] = n_at; }
    else if (strcmp(top->force_names[i],"PDIHS") == 0) { ++n_dt; top->int_map[i] = n_dt; }
  }

  fprintf(fp,"     %d  bonds\n",n_bonds);
  fprintf(fp,"     %d  angles\n",n_angles);
  fprintf(fp,"     %d  dihedrals\n",n_dih);
  fprintf(fp,"     %d  impropers\n\n",n_imp);
  fprintf(fp,"     %d  atom types\n",n_att);
  fprintf(fp,"     %d  bond types\n",n_bt);
  fprintf(fp,"     %d  angle types\n",n_at);
  fprintf(fp,"     %d  dihedral types\n\n",n_dt);
  fprintf(fp,"     0  %f xlo xhi\n",fr->contents->box[0][0] / X_LMP2GRO);
  fprintf(fp,"     0  %f ylo yhi\n",fr->contents->box[1][1] / X_LMP2GRO);
  fprintf(fp,"     0  %f zlo zhi\n\n",fr->contents->box[2][2] / X_LMP2GRO);
  fprintf(fp,"Masses\n\n");

  double *masses = (double *) ecalloc(n_att,sizeof(double));
  double *charges = (double *) ecalloc(n_att,sizeof(double));
  for (i = 0; i < n_att; ++i)
  {
    for (j = 0; j < top->contents->mols.nr; ++j)
    {
      for (k = 0; k < top->molecules[j].n_apm; ++k)
      {
        if (top->molecules[j].type_ids[k] == i) 
	{ 
	  masses[i] = top->molecules[j].m[k]; 
	  charges[i] = top->molecules[j].q[k];
	  k = top->molecules[j].n_apm;
	  j = top->contents->mols.nr;
	}
      }
    }
    fprintf(fp,"    %d %g\n",i+1,masses[i]);
  }

  fprintf(fp,"\nAtoms\n\n");
  for (molt_idx = 0; molt_idx < top->contents->mols.nr; ++molt_idx)
  {
    for (mol_idx = 0; mol_idx < top->molecules[molt_idx].n_mols; ++mol_idx)
    {
      for (at_idx = 0; at_idx < top->molecules[molt_idx].n_apm; ++at_idx)
      {
	// AT_ID MOL_TYPE AT_TYPE CHARGE X Y Z
	// Note, LAMMPS starts indexing atom types, atoms, and molecules from 1
	// That is why there are +1s tacked on to abs_at_idx, mol_idx, and "atom".type
	fprintf(fp,"%d %d %d %g %f %f %f\n",abs_at_idx+1,
                                            mol_idx+1+n_prev_mol,
                                            top->contents->atoms.atom[abs_at_idx].type+1,
                                            top->contents->atoms.atom[abs_at_idx].q,
                                            fr->contents->x[abs_at_idx][0] / X_LMP2GRO,
                                            fr->contents->x[abs_at_idx][1] / X_LMP2GRO,
                                            fr->contents->x[abs_at_idx][2] / X_LMP2GRO);
        ++abs_at_idx;
      }
    }
    n_prev_mol += top->molecules[molt_idx].n_mols;
  }

  fprintf(fp,"\nBonds\n\n");
  abs_at_idx = 0; // reusing this instead of making abs_bond_idx
  n_prev_at = 0;
  for (molt_idx = 0; molt_idx < top->contents->mols.nr; ++molt_idx)
  {
    for (mol_idx = 0; mol_idx < top->molecules[molt_idx].n_mols; ++mol_idx)
    {
//reusing at_idx here instead of making bond_idx
      for (at_idx = 0; at_idx < top->molecules[molt_idx].n_bonds; ++at_idx)
      {
	fprintf(fp,"%d %d %d %d\n",abs_at_idx+1,
				   top->int_map[top->molecules[molt_idx].bond_types[at_idx]], 
				   top->molecules[molt_idx].bond_ij[at_idx][0] + n_prev_at + 1,
				   top->molecules[molt_idx].bond_ij[at_idx][1] + n_prev_at + 1);
	++abs_at_idx;
      }
      n_prev_at += top->molecules[molt_idx].n_apm;
    }
  }

  fprintf(fp,"\nAngles\n\n");
  abs_at_idx = 0; // reusing this instead of making abs_angle_idx
  n_prev_at = 0;
  for (molt_idx = 0; molt_idx < top->contents->mols.nr; ++molt_idx)
  {
    for (mol_idx = 0; mol_idx < top->molecules[molt_idx].n_mols; ++mol_idx)
    {
//reusing at_idx here instead of making angle_idx
      for (at_idx = 0; at_idx < top->molecules[molt_idx].n_angles; ++at_idx)
      {
        fprintf(fp,"%d %d %d %d %d\n",abs_at_idx+1,
                                      top->int_map[top->molecules[molt_idx].angle_types[at_idx]],
                                      top->molecules[molt_idx].angle_ijk[at_idx][0] + n_prev_at + 1,
                                      top->molecules[molt_idx].angle_ijk[at_idx][1] + n_prev_at + 1,
				      top->molecules[molt_idx].angle_ijk[at_idx][2] + n_prev_at + 1);
        ++abs_at_idx;
      }
      n_prev_at += top->molecules[molt_idx].n_apm;
    }
  }

  fprintf(fp,"\nDihedrals\n\n");
  abs_at_idx = 0; // reusing this instead of making abs_dih_idx
  n_prev_at = 0;
  for (molt_idx = 0; molt_idx < top->contents->mols.nr; ++molt_idx)
  {
    for (mol_idx = 0; mol_idx < top->molecules[molt_idx].n_mols; ++mol_idx)
    {
//reusing at_idx here instead of making dih_idx
      for (at_idx = 0; at_idx < top->molecules[molt_idx].n_dihs; ++at_idx)
      {
        fprintf(fp,"%d %d %d %d %d %d\n",abs_at_idx+1,
					 top->int_map[top->molecules[molt_idx].dih_types[at_idx]],
					 top->molecules[molt_idx].dih_ijkl[at_idx][0] + n_prev_at + 1,
					 top->molecules[molt_idx].dih_ijkl[at_idx][1] + n_prev_at + 1,
					 top->molecules[molt_idx].dih_ijkl[at_idx][2] + n_prev_at + 1,
					 top->molecules[molt_idx].dih_ijkl[at_idx][3] + n_prev_at + 1);
        ++abs_at_idx;
      }
      n_prev_at += top->molecules[molt_idx].n_apm;
    }
  }
  fclose(fp);
  dump_molecule_info(top);
}

/*
 get_type() finds the index for a_type as contained in the topology
 */
int get_type(const char *a_type, tW_gmx_topology *top)
{
  int i;
  for (i = 0; i < top->n_atomtypes; ++i)
  {
    if (strcmp(a_type,top->atom_type_names[i]) == 0) { return i; }
  }
  fprintf(stderr,"ERROR: Atom type: %s not found in topology\n",a_type);
  fprintf(stderr,"topology's atom types: \n");
  for (i = 0; i < top->n_atomtypes; ++i)
  {
    fprintf(stderr,"  %s  \n",top->atom_type_names[i]);
  }
  exit(1);
  return -1;
}

void do_PBC(tW_gmx_trxframe *fr)
{
  int i, j;
  for (i = 0; i < fr->contents->natoms; ++i)
  {
    for (j = 0; j < DIM; ++j)
    {
      if (fr->contents->x[i][j] >= fr->contents->box[j][j]) { fr->contents->x[i][j] -= fr->contents->box[j][j]; }
      else if (fr->contents->x[i][j] < 0.0) { fr->contents->x[i][j] += fr->contents->box[j][j]; }
    }
  }
}


/* Local Density Function 06.17.2020 MRD */

/* Populates the LD vector for an atom given an input line from a .btj file */
void pop_LD_vector(tW_line inp_line, double *lds, int ntypes)
{
  char *start;
//  float inp_float; MCL PREC FIX 02.28.25
  double inp_float;
  int i, test_sscanf;

  start = &(inp_line[0]);
  for (i = 0; i < ntypes; ++i)
  {
    test_sscanf = sscanf(start," %lf ",&inp_float); // MCL 02.05.25 - in sscanf this is a double, need to be careful reading in text data so that we don't round onto knot pointss for OG Bspline_deriv alg
    if (test_sscanf != 1)
    {
      fprintf(stderr,"ERROR: expected to read %dth float from line: %s",i+1,inp_line);
      fprintf(stderr,"\tCurrent location of line: %s",start);
      exit(1);
    }
    lds[i] = (double) inp_float;
    while ((*start) == ' ') { ++start; } // Finds the next number (non-space character, which is in inp_float)
    start = strstr(start," "); // Finds the space after that number
  }
}

void populate_LD_Gradients(tW_line inp_line, dvec *ldgs, int ntypes)
{
  char *start, *end;
  float inp_float;
  int i, j, test_sscanf;

  start = &(inp_line[0]);
  for (i = 0; i < ntypes; ++i)
  {
    for (j = 0; j < DIM; ++j)
    {
      test_sscanf = sscanf(start, " %f ", &inp_float);
      if (test_sscanf != 1)
      {
        fprintf(stderr,"ERROR: expected to read %dth float from line: %s",i*DIM+j+1,inp_line);
        fprintf(stderr,"\tCurrent location of line: %s",start);
        exit(1);
      }
      ldgs[i][j] = (double) inp_float;
      while ((*start) == ' ') { ++start; }
      start = strstr(start," ");
    }
  }
}

void copy_xvf(tW_gmx_trxframe *source, tW_gmx_trxframe *dest, bool bX, bool bV, bool bF)
{
  int i;

  dvec zero;
  zero[0] = zero[1] = zero[2] = 0.0;
  for (i = 0; i < source->contents->natoms; ++i)
  {
    if (bX && source->contents->bX && dest->contents->bX) { copy_vector(source->contents->x[i],dest->contents->x[i]); }
    else if (dest->contents->bX) { copy_vector(zero,dest->contents->x[i]); }
    if (bV && source->contents->bV && dest->contents->bV) { copy_vector(source->contents->v[i],dest->contents->v[i]); }
    else if (dest->contents->bV) { copy_vector(zero,dest->contents->v[i]); }
    if (bF && source->contents->bF && dest->contents->bF) { copy_vector(source->contents->f[i],dest->contents->f[i]); }
    else if (dest->contents->bF) { copy_vector(zero,dest->contents->f[i]); }
  }
}


void copy_lds(tW_gmx_trxframe *source, tW_gmx_trxframe *dest)
{
  int i = 0, j = 0;
  int n_atomtypes = source->n_atomtypes;
  if (source->bLDs)
  {
    dest->bLDs = source->bLDs;
    for (i = 0; i < source->contents->natoms; ++i)
    {
      for (j = 0; j < n_atomtypes; ++j)
      {
        dest->local_densities[i][j] = source->local_densities[i][j];
      }
    }
  }
}


void copy_ldgs(tW_gmx_trxframe *source, tW_gmx_trxframe *dest)
{
  int register i = 0;
  int register j = 0;
  int register n_atomtypes = source->n_atomtypes;
  int register natoms = source->contents->natoms;
  if (source->bLD_Grads)
  {
    dest->bLD_Grads = source->bLD_Grads;
    for (i = 0; i < natoms; ++i)
    {
      for (j = 0; j <= n_atomtypes; ++j)
      {
        dest->local_density_gradients[i][j][0] = source->local_density_gradients[i][j][0];
        dest->local_density_gradients[i][j][1] = source->local_density_gradients[i][j][1];
        dest->local_density_gradients[i][j][2] = source->local_density_gradients[i][j][2];
      }
    }
  }
}




/* MRD Index Files */

int get_n_words_from_line(tW_line inp_line)
{
  tW_line new_line;
  strcpy(new_line,inp_line);
  int n_words = 0;
  int n_chars = 0;
  int idx = 0;
  bool in_word;
  while ((new_line[idx] != '\n') && (new_line[idx] != '\0'))
  {
    if (new_line[idx] == '\t') { new_line[idx] = ' '; }
    ++idx;
  }
  n_chars = idx;
  idx = 0;
  in_word = FALSE;
  while (idx < n_chars)
  {
    if ((new_line[idx] != ' ') && (new_line[idx] != '\0') && (new_line[idx] != '\n'))
    {
      if (! in_word) { ++n_words; in_word = TRUE; }
    }
    else
    {
      in_word = FALSE;
    }
    ++idx;
  }
  return n_words;
}

void put_at_nums_in(tW_line inp_line, int *spots, int n_words_this_line)
{
  int at_idx = 0;
  int idx = 0;
  int n_chars;
  bool in_word = FALSE;
  char *p, *p2;
  tW_line new_line;
  strcpy(new_line,inp_line);
  while ((new_line[idx] != '\n') && (new_line[idx] != '\0'))
  {
    if (new_line[idx] == '\t') { new_line[idx] = ' '; }
    ++idx;
  }
  n_chars = idx;
  idx = 0;
  in_word = FALSE;
  while ((at_idx < n_words_this_line) && (idx < n_chars))
  {
    if ((new_line[idx] != ' ') && (new_line[idx] != '\0') && (new_line[idx] != '\n'))
    {
      if (! in_word)
      {
        p = &(new_line[idx]);
        p2 = strstr(p," ");
        if (p2 == NULL)
        {
          p2 = strstr(p,"\n");
        }
        *p2 = '\0';
        spots[at_idx] = atoi(p);
        ++at_idx;
        *p2 = ' ';
        in_word = TRUE;
      }
    }
    else
    {
      in_word = FALSE;
    }
    ++idx;
  }
}

void alloc_index_ng(tW_index_file *idx, int ng)
{
  idx->n_groups = ng;
  idx->group_names = (tW_word *) ecalloc(ng,sizeof(tW_word));
  idx->n_at_per_group = (int *) ecalloc(ng,sizeof(int));
  idx->atoms_in_group = (int **) ecalloc(ng,sizeof(int *));
}

void alloc_index_aig(tW_index_file *idx)
{
  int i;
  for (i = 0; i < idx->n_groups; ++i)
  {
    idx->atoms_in_group[i] = (int *) ecalloc(idx->n_at_per_group[i], sizeof(int));
  }
}

void read_index_file(tW_word fnm, tW_gmx_topology *top)
{
  FILE *fp = open_file(fnm,'r');
  tW_line inp_line;
  tW_word inp_word;
  int test_sscanf;
  int n_groups = 0;
  int i, j;
  while (get_next_line(fp,inp_line) != -1)
  {
    if ((strstr(inp_line,"[") != NULL) && (strstr(inp_line,"]") != NULL))
    {
      ++n_groups;
    }
  }
  rewind(fp);

  if (! top->idx_file) { top->idx_file = (tW_index_file *) ecalloc(1,sizeof(tW_index_file)); }

  alloc_index_ng(top->idx_file,n_groups);

  int group_idx = -1;
  int n_atoms_on_line;
  while (get_next_line(fp,inp_line) != -1)
  {
    if ((strstr(inp_line,"[") != NULL) && (strstr(inp_line,"]") != NULL))
    {
      ++group_idx;
      test_sscanf = sscanf(inp_line,"[ %s ]",&inp_word);
      top->idx_file->n_at_per_group[group_idx] = 0;
      if (test_sscanf != 1)
      {
        fprintf(stderr,"ERROR: unable to read name for index group idx %d\n",group_idx);
        fprintf(stderr,"\tline: %s",inp_line);
        exit(1);
      }
      strcpy(top->idx_file->group_names[group_idx],inp_word);
    }
    else
    {
      int n_words_this_line = get_n_words_from_line(inp_line);
      if (n_words_this_line > 0)
      {
        top->idx_file->n_at_per_group[group_idx] += n_words_this_line;
      }
    }
  }

  alloc_index_aig(top->idx_file);

  rewind(fp);
  group_idx = -1;
  int at_idx;
  while (get_next_line(fp,inp_line) != -1)
  {
    if ((strstr(inp_line,"[") != NULL) && (strstr(inp_line,"]") != NULL))
    {
      ++group_idx;
      at_idx = 0;
    }
    else
    {
      int n_words_this_line = get_n_words_from_line(inp_line);
      if (n_words_this_line > 0)
      {
        put_at_nums_in(inp_line,&(top->idx_file->atoms_in_group[group_idx][at_idx]),n_words_this_line);
        at_idx += n_words_this_line;
      }
    }
  }
  fclose(fp);
}

void copy_idx_file(tW_index_file *dest, tW_index_file *source)
{
  int i, j;
  alloc_index_ng(dest,source->n_groups);
  for (i = 0; i < dest->n_groups; ++i)
  {
    dest->n_at_per_group[i] = source->n_at_per_group[i];
    strcpy(dest->group_names[i],source->group_names[i]);
    dest->atoms_in_group[i] = (int *) ecalloc(dest->n_at_per_group[i],sizeof(int));
    for (j = 0; j < dest->n_at_per_group[i]; ++j)
    {
      dest->atoms_in_group[i][j] = source->atoms_in_group[i][j];
    }
  }
}

void generate_generic_index_file(tW_gmx_topology *top, bool bOld)
{
  int i, j, k, n_groups;
  tW_index_file old_idx_file;
  if (bOld)
  {
    copy_idx_file(&(old_idx_file),top->idx_file);
    efree(top->idx_file);
  }
  top->idx_file = (tW_index_file *) ecalloc(1,sizeof(tW_index_file));
  n_groups = 1 + top->contents->mols.nr +  top->n_atomtypes;
  bool *olds;
  if (bOld)
  {
    olds = (bool*) ecalloc(old_idx_file.n_groups, sizeof(bool));
    bool exists = FALSE;
    for (i = 0; i < old_idx_file.n_groups; ++i)
    {
      exists = FALSE;
      if (strcmp(old_idx_file.group_names[i],"System") == 0) { exists = TRUE; }
      for (j = 0; j < top->contents->mols.nr; ++j)
      {
        if (strcmp(old_idx_file.group_names[i],top->molecules[j].molname) == 0) { exists = TRUE; }
      }
      for (j = 0; j < top->n_atomtypes; ++j)
      {
        if (strcmp(old_idx_file.group_names[i],top->atom_type_names[j]) == 0) { exists = TRUE; }
      }

      if (! exists) { ++n_groups; }
      olds[i] = !exists;
    }
  }

  alloc_index_ng(top->idx_file,n_groups);
  int group_idx = 0;

  /* first group is always system and contains everything */
  strcpy(top->idx_file->group_names[group_idx],"System");
  top->idx_file->n_at_per_group[group_idx] = top->contents->atoms.nr;
  ++group_idx;
  /* next groups are the molecule groups */
  for (i = 0; i < top->contents->mols.nr; ++i)
  {
    strcpy(top->idx_file->group_names[group_idx],top->molecules[i].molname);
    top->idx_file->n_at_per_group[group_idx] = top->molecules[i].n_mols*top->molecules[i].n_apm;
    ++group_idx;
  }
  /* last groups are the atom groups */
  for (i = 0; i < top->n_atomtypes; ++i)
  {
    /* Make sure no atom type has same name as molecule type */
    bool mol_at_different = TRUE;
    for (j = 0; j < top->contents->mols.nr; ++j)
    {
      if (strcmp(top->molecules[j].molname,top->atom_type_names[i]) == 0) { mol_at_different = FALSE; }
    }
    if (mol_at_different)
    {
      strcpy(top->idx_file->group_names[group_idx],top->atom_type_names[i]);
    }
    else
    {
      sprintf(top->idx_file->group_names[group_idx],"%s_ATOMS",top->atom_type_names[i]);
    }
    int n_at = 0;
    for (j = 0; j < top->contents->mols.nr; ++j)
    {
      int n_in_mol = 0;
      for (k = 0; k < top->molecules[j].n_apm; ++k)
      {
        if (strcmp(top->atom_type_names[i],top->molecules[j].atom_types[k]) == 0)
        {
          ++n_in_mol;
        }
      }
      n_at += n_in_mol * top->molecules[j].n_mols;
    }
    top->idx_file->n_at_per_group[group_idx] = n_at;
    ++group_idx;
  }
  /* finally do old index file groups */
  if (bOld)
  {
    for (i = 0; i < old_idx_file.n_groups; ++i)
    {
      if (olds[i])
      {
        strcpy(top->idx_file->group_names[group_idx],old_idx_file.group_names[i]);
        top->idx_file->n_at_per_group[group_idx] = old_idx_file.n_at_per_group[i];
        ++group_idx;
      }
    }
  }

  alloc_index_aig(top->idx_file);
  /* Populate atoms_in_group s */
  /* System */
  group_idx = 0;
  for (i = 0; i < top->contents->atoms.nr; ++i)
  {
    top->idx_file->atoms_in_group[/*group_idx*/0][i] = i;
  }
  ++group_idx;
  /* Molecules */
  /* This relies on the assumption that all molecules of a given type are grouped together in the top/gro files */
  int at_idx = 0;
  for (i = 0; i < top->contents->mols.nr; ++i)
  {
    for (j = 0; j < top->molecules[i].n_mols * top->molecules[i].n_apm; ++j)
    {
      top->idx_file->atoms_in_group[group_idx][j] = at_idx;
      ++at_idx;
    }
    ++group_idx;
  }
  /* Atoms */
  int *idxes = (int *) ecalloc(top->n_atomtypes, sizeof(int));
  for (i = 0; i < top->n_atomtypes; ++i)
  {
    idxes[i] = 0;
  }
  for (i = 0; i < top->contents->atoms.nr; ++i)
  {
    int itype = get_type(*(top->contents->atoms.atomtype[i]),top);
    top->idx_file->atoms_in_group[group_idx+itype][idxes[itype]] = i;
    ++(idxes[itype]);
  }
  group_idx += top->n_atomtypes;
  /* Old ones */
  if (bOld)
  {
    for (i = 0; i < old_idx_file.n_groups; ++i)
    {
      if (olds[i])
      {
        for (j = 0; j < old_idx_file.n_at_per_group[i]; ++j)
        {
          top->idx_file->atoms_in_group[group_idx][j] = old_idx_file.atoms_in_group[i][j];
        }
        ++group_idx;
      }
    }
    efree(olds);
  }
  efree(idxes);
}

void print_index_names_and_atom_counts(tW_index_file *idx)
{
  int i;
  fprintf(stderr,"%3s    %20s    %10s\n","Idx","Index Group Name","Atom Count");
  fprintf(stderr,"---    --------------------    ----------\n");
  for (i = 0; i < idx->n_groups; ++i)
  {
    fprintf(stderr,"%3d    %20s    %6d\n",i,idx->group_names[i],idx->n_at_per_group[i]);
  }
}


/*****************************************************************************************
skip_excl():
*****************************************************************************************/
bool skip_excl(int nr_excl, int *excl_list, int j)
{
    int i;

    for (i = 0; i < nr_excl; i++) {
        if (j == excl_list[i]) {
            return TRUE;
        }
    }

    return FALSE;
}

/******************************************************************************************
dprocess_sysname_natoms():
@arg fp the pointer to the dumped tpr file
@arg top the pointer to the topology we are populating via tpr info
@arg lin_buf allocated space for working with lines (dummy space)
@return true (success) false (something bad happened)
This function reads dumped tpr files from the topology header and 2 lines below. 
In doing so it obtains the system name and number of atoms
it returns false if there was a problem and errors werent caught earlier
*******************************************************************************************/
bool dprocess_sysname_natoms(FILE* fp, tW_gmx_topology* top, tW_line* lin_buf)
{
 bool bOk;
 tW_word inp_word;
 int integer_holder;
// skip to topology:
 bOk = find_line_in_file(fp, "topology:", NULL);
 get_next_line(fp, *lin_buf); // expect name="name" undertop
 if (!bOk) { THROW_TPR_ERROR("Can't find `topology:`\n");}
 bOk =!(strstr(*lin_buf, "name=\"") == NULL); //NOT OKAY IF NO OVERLAP IS TRUE
 if (!bOk == TRUE) { THROW_TPR_ERROR("Can't find name=\"name\" under topology:\n"); }
 bOk = sscanf(*lin_buf, " name=\"%s ", &inp_word);
 if (!bOk) { THROW_TPR_ERROR("name=name fmt change\n"); }
 elim_char(inp_word, '"'); strcpy(top->contents->name, inp_word); //TOP INFO->SYSNAME
 get_next_line(fp,*lin_buf); // expect #atoms = XXXXX
 bOk = !((strstr(*lin_buf,"#atoms") == NULL));
 if (!bOk) { THROW_TPR_ERROR("#atoms = XXXX not found\n"); }
 bOk = sscanf(*lin_buf, " #atoms = %d ", &integer_holder);
 if (!bOk) { THROW_TPR_ERROR("atoms=XXX fmt change\n"); }
 top->contents->atoms.nr = integer_holder; //TOP INFO-> TOTAL ATOMS
 return TRUE;
}

/******************************************************************************************
dprocess_nmolblocks():
@arg fp the pointer to the dumped tpr file
@arg top the pointer to the topology we are populating via tpr info
@arg lin_buf allocated space for working with lines (dummy space)
@return true (success) false (something bad happened)
This function reads dumped tpr file to determine the number of molecule blocks
(It allows for there to be e.g. 3 molecule blocks where e.g. 2 different molecule types are listed)
(This feature is per MRD 01.17.2020 BOCS issue 11)
This sets us up for reading the molecule details later by figuring out how many are required.
it returns false if there was a problem and errors werent caught earlier
*******************************************************************************************/
bool dprocess_nmolblocks(FILE* fp, tW_gmx_topology* top, tW_line* lin_buf)
{
 bool bOk;
 tW_word inp_word;
 int EOF_sentinal = 0;
 int nmolblocks = 0;
 //Find molblock section in versions 5.1.4 and later
 bOk = find_line_in_file(fp, "#molblock", *lin_buf);
 if (!bOk)
 {
	  // We think this can happen in earlier versions of gromacs
	  // Scroll through file from here, count number times molblock comes up
	  find_line_in_file(fp, "topology:\n", NULL);
	  while (EOF_sentinal != -1) // EOF
	  {
		  EOF_sentinal = get_next_line(fp, *lin_buf);
		  if (strstr(*lin_buf, "molblock (") != NULL) { nmolblocks++; }
	  }
	  if (nmolblocks == 0) { THROW_TPR_ERROR("UNABLE TO COUNT \"molblock (X)\" SECTIONS\n");}
 }
 else 
 {
	  bOk = sscanf(*lin_buf, " #molblock = %d", &nmolblocks);
	  if (!bOk) { THROW_TPR_ERROR("UNABLE TO PARSE #molblock = # in tpr\n"); }
 }
 if (nmolblocks == 0) { THROW_TPR_ERROR("COULD NOT READ n in #molblock =n\n"); }
 top->n_molblocks = nmolblocks;
 
 return TRUE;
}

/**************************************************************************************
dget_moltypes_in_molblocks():
@arg fp the pointer to the dumped tpr file
@arg top the pointer to the topology we are populating via tpr info
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This function determines the number of unique molecule types listed in the tpr molblocks
It populates the list of unique moltypes, and coordinates them to each block listed
**************************************************************************************/
bool dget_moltypes_in_molblocks(FILE* fp, tW_gmx_topology* top, tW_line* lin_buf)
{
 bool bOk;
 tW_word inp_word; 
 int integer_holder;
 int nmoleculetypes = 0;
 tW_word molblock_header;
 for (int i = 0; i < top->n_molblocks ; i++)
 {
  sprintf(molblock_header, "molblock (%d):", i);
  // skip to moleblock (#):
  bOk = find_line_in_file(fp, molblock_header, *lin_buf);
  if (!bOk) { sprintf(molblock_header, "UNABLE TO FIND molblock (%d):",i);
              THROW_TPR_ERROR(molblock_header);}
  // We expect names and indices one line down from molblock header. e.g. moltype = 0 "BRANCH"
  get_next_line(fp, *lin_buf);
  bOk = sscanf(*lin_buf, " moltype = %d \"%s", &integer_holder, &inp_word);
  if (!bOk) { THROW_TPR_ERROR("COULD NOT PARSE MOLTYPE\n"); }
  elim_char(inp_word, '"');
  
  // Associate molecule with block
  top->molblocks[i].moltype = integer_holder;
  strcpy(top->molblocks[i].molname, inp_word);

  // If index read >= ntypes we know about, count total types as that index + 1 [eg i=0 is 1 type]
  // (This is an MRD innovation MCL ported here)
  if (integer_holder + 1 > nmoleculetypes) { nmoleculetypes = integer_holder + 1; }

  //figure out how many molecules in system are listed in this block. We expect it listed on the 
  //next line down e.g. #molecules = 
  get_next_line(fp, *lin_buf);
  bOk = sscanf(*lin_buf, " #molecules = %d ", &integer_holder);
  if (!bOk) { THROW_TPR_ERROR("UNABLE TO DETERMINE nmolecules from molblock\n");}
  top->molblocks[i].n_mols = integer_holder;
  //fprintf(stderr, "molblock[%d] n = %d\n", i, top->molblocks[i].n_mols);
 }
 top->n_molecule_types = nmoleculetypes;

 return TRUE;
}

/**************************************************************************************
dget_top_ff_info():
@arg fp the pointer to the dumped tpr file
@arg top the pointer to the topology we are populating via tpr info
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This function reads the list of different interactions that will be simulated if 
the tpr dumped is simulated in gromacs. [Does not account for distinct energy groups]
@note MCL 11.07.2025 To clarify, we don't do anything with this info for force-matching. The list of interactions read in in that case is determine by the par file and Bocs atom types.
@warning this function dynamically allocates memory for top structure members
**************************************************************************************/
bool dget_top_ff_info(FILE* fp, tW_gmx_topology* top, tW_line* lin_buf)
{
 bool bOk;
 tW_word interaction_name;
 tW_word functype_key;
 int integer_holder;
 float fudge_holder;

 bOk = find_line_in_file(fp, "ffparams", NULL);
 if (!bOk) { THROW_TPR_ERROR("Unable to read ffparams from dumped tpr\n");}
 get_next_line(fp, *lin_buf); //we expect atnr = X one line below
 bOk = sscanf(*lin_buf, " atnr=%d ", &integer_holder);
 if (!bOk) { THROW_TPR_ERROR("Unable to determine n_unique nonbonded atom params\n"); }
 top->contents->idef.atnr = integer_holder;
 top->contents->atomtypes.nr = top->contents->idef.atnr;
 get_next_line(fp, *lin_buf); // ntypes=X
 bOk = sscanf(*lin_buf, " ntypes=%d ", &integer_holder);
 if (!bOk) { THROW_TPR_ERROR("Unable to determine n_interactions defined by top\n");}
 top->contents->idef.ntypes = integer_holder;

 // Now we can make room for the list and search for names
 top->contents->idef.functype = (tW_t_functype *) ecalloc(top->contents->idef.ntypes,sizeof(tW_t_functype));
 top->contents->idef.iparams = (tW_t_iparams *) ecalloc(top->contents->idef.ntypes,sizeof(tW_t_iparams));
 top->force_names = (tW_word *) ecalloc(top->contents->idef.ntypes,sizeof(tW_word));
 
 // We expect lines like
 // functype[0]=LJ_SR, c6= 0.00000000e+00, c12= 1.00000000e+00, though some types take multiple lines
 for (int i = 0; i < top->contents->idef.ntypes; ++i)
  {
    sprintf(functype_key, "functype[%d]", i);
    bOk = find_line_in_file(fp, functype_key, *lin_buf);
    if (!bOk) { 
	       sprintf(functype_key, "Unable to find functype[%d] when nfunctypes is %d\n", i, top->contents->idef.ntypes);
	       THROW_TPR_ERROR(functype_key); }
    bOk = sscanf(*lin_buf," functype[%d]=%s ",&integer_holder,&interaction_name);
    if (!bOk) { THROW_TPR_ERROR("Unexpected functype format\n");}
    elim_char(interaction_name,',');
    strcpy(top->force_names[i],interaction_name);
  }
 bOk = find_line_in_file(fp, "fudgeQQ", *lin_buf);
 if (!bOk) { THROW_TPR_ERROR("Could not find mixing coeff\n"); }
 bOk = sscanf(*lin_buf, " fudgeQQ = %f ", &fudge_holder);
 if (!bOk) { THROW_TPR_ERROR("Unexpected fudgefmt \n"); }
 top->contents->idef.fudgeQQ = fudge_holder;
 return TRUE;
}
/**************************************************************************************
dfind_all_molinfo:
@arg fp the pointer to the dumped tpr file
@arg mol the list of molecule types in the top
@arg nmol the length of the list
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This function searches the tpr file for the various fields of the molecule structure.
It does this for every distinct molecule type in the system
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note This function replaces get_moltype_info
**************************************************************************************/
bool dfind_all_molinfo(FILE *fp, tW_molecule * mol, int nmol, tW_line * lin_buf)
{
 bool bOk; 
 tW_word word_holder;
 int integer_holder;
 tW_word moltype_header;

 for (int i = 0; i < nmol; i++) // Scan info for all dif mol types
 {
   sprintf(moltype_header, "moltype (%d):", i);
   bOk = find_line_in_file(fp, moltype_header, NULL);
   if (!bOk) {
	   sprintf(moltype_header, "Unable to find moltype (%d): when nmoltypes is %d\n", i, nmol);
	   THROW_TPR_ERROR(moltype_header); }
   get_next_line(fp, *lin_buf); //name = X
   bOk = sscanf(*lin_buf, " name=\"%s\" ", &word_holder);
   if (!bOk) { sprintf(moltype_header, "Could not parse name under moltype (%d):\n", i);
	   THROW_TPR_ERROR(moltype_header);}
   elim_char(word_holder,'"');
   strcpy(mol[i].molname, word_holder);

   get_next_line(fp, *lin_buf); // atoms:
   get_next_line(fp, *lin_buf); // atom (#):
   bOk = sscanf(*lin_buf, " atom (%d) ", &integer_holder);
   if (!bOk) { sprintf(moltype_header, "Unable to parse apm for moltype %d\n", i);
               THROW_TPR_ERROR(moltype_header); }
   // Now we have atoms per molecule
   mol[i].n_apm = integer_holder;
   bOk = dparse_atom_per_molecule_info(fp, &mol[i], lin_buf);
   if (!bOk) { sprintf(moltype_header, "Unable to parse atom info for moltype %d\n", i);
               THROW_TPR_ERROR(moltype_header);}
   bOk = dparse_residue_per_molecule_info(fp, &mol[i], lin_buf);
   if (!bOk) { sprintf(moltype_header, "Unable to parse residue info for moltype %d\n", i);
               THROW_TPR_ERROR(moltype_header);}
   bOk = dparse_charge_groups_per_molecule(fp, &mol[i], lin_buf, i);
   if (!bOk) { sprintf(moltype_header, "Unable to parse charge_group info for moltype %d\n", i);
               THROW_TPR_ERROR(moltype_header);}
   bOk = dparse_atomic_exclusions_per_molecule(fp, &mol[i], lin_buf, i);
   if (!bOk) { sprintf(moltype_header, "Unable to parse excl. per. atom info for moltype %d\n", i);
               THROW_TPR_ERROR(moltype_header);}
   bOk = dparse_bonded_int_lists_per_molecule(fp, &mol[i], lin_buf, i);

 }
 
}
/**************************************************************************************
dparse_atom_per_molecule_info:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This function is meant to parse the section of text following moltype(i)-> atom (apm)
It fills in the atomic info for the passed molecule structure
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note This function is a helper function @see dfind_all_molinfo
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
bool dparse_atom_per_molecule_info(FILE* fp, tW_molecule* moli, tW_line* lin_buf)
{
 bool bOk; 
 int at_per_mi = moli->n_apm;
 int atype_id, residx;
 tW_word atname;
 tW_word atype, atomBtype;
 float mass, charge; // doubles > float but sscanf works with float


 // MRD 11.05.2019 now we finally have n_apm, allocate memory for stuff
  moli->type_ids = (int *) ecalloc(at_per_mi,sizeof(int));
  moli->residx = (int *) ecalloc(at_per_mi,sizeof(int));
  moli->atom_names = (tW_word *) ecalloc(at_per_mi,sizeof(tW_word));
  moli->atom_types = (tW_word *) ecalloc(at_per_mi,sizeof(tW_word));
  moli->atom_Btypes = (tW_word *) ecalloc(at_per_mi,sizeof(tW_word));
  moli->m = (double *) ecalloc(at_per_mi,sizeof(double));
  moli->q = (double *) ecalloc(at_per_mi,sizeof(double));

  // off target variables (dummies for reading/checking)
  int atom_index, dumi, dumi2, dumi3;
  tW_word dumword, dumword2;
  float dumf, dumf1;

  for (int i = 0; i < at_per_mi; i++) // Scan for mass/charge/residx
  {
    get_next_line(fp,*lin_buf); 
    bOk = sscanf(*lin_buf," atom[ %d]={type= %d, typeB= %d, ptype= %s m= %f, q= %f, mB= %g, qB= %g, resind= %d, atomnumber= %d} ",&atom_index, &atype_id, &dumi, &dumword, &mass, &charge, &dumf, &dumf1, &residx, &dumi3);
    if (bOk != 10) { THROW_TPR_ERROR("COULDNT READ MASS CHARGE RESIDX PROPERLY\n"); }
    bOk=(atom_index == i);//if file alligned correctly, index should correspond to line we're processing
    if (!bOk) { THROW_TPR_ERROR("DID NOT PROPERLY READ ATOMS_PER_MOLECULE\n"); } 
    moli->type_ids[i] = atype_id;
    moli->m[i] = mass;
    moli->q[i] = charge;
    moli->residx[i] = residx;
  }
  get_next_line(fp, *lin_buf); // skip second atom (apm):
  
  for (int i = 0; i < at_per_mi; i++) // Scan for atom names
  {
    get_next_line(fp, *lin_buf);
    bOk = sscanf(*lin_buf," atom[%d]={name=\"%s\"} ",&atom_index, &atname);
    if (bOk != 2) { THROW_TPR_ERROR("DID NOT PARSE ATOM NAMES\n"); }
    if (atom_index != i) { THROW_TPR_ERROR("APM READER MISALIGNED IN ATOM NAMES\n"); }
    elim_char(atname, '"');
    strcpy(moli->atom_names[i], atname);
  }
  get_next_line(fp, *lin_buf); // skip type (apm):
  
  for (int i = 0; i < at_per_mi; i++) // Scan for atom types, they list 2, we take first
  {
    get_next_line(fp, *lin_buf);
   // bOk = sscanf(*lin_buf, " type[%d]=\"%s,nameB=\"%s}", &atom_index, &atype, &dumword);
    bOk = sscanf(*lin_buf, " type[%d]={name=\"%s", &atom_index, &atype);// &dumword);
   // fprintf(stderr, "%s %d\n", *lin_buf, bOk);
    if (bOk != 2) { THROW_TPR_ERROR("DID NOT PARSE ATOM TYPES\n"); }
    if (atom_index != i) { THROW_TPR_ERROR("APM READER MISALIGNED IN ATOM TYPES\n"); }
    elim_char(atype, '"');
    strcpy(moli->atom_types[i], atype);
  }
  return TRUE;
}

/**************************************************************************************
dparse_residue_per_molecule_info:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This function is meant to parse the section of text following moltype(i)-> residue(nresidue)
It fills in the residue info for the passed molecule structure
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note This function is a helper function @see dfind_all_molinfo
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
bool dparse_residue_per_molecule_info(FILE* fp, tW_molecule* moli, tW_line* lin_buf)
{
 bool bOk;
 int nresidues, residue_idx;
 tW_word resname;

 int dumi;
// tW_word dum_word;
 char dum_word[50];
 char atype[50];

 // Find nresidue and make space
 get_next_line(fp, *lin_buf); // should be residue (nres):, follow by the list
 bOk = sscanf(*lin_buf, " residue (%d): ", &nresidues);
 if (!bOk) { sprintf(resname, "COULD NOT READ NRESIDUE FROM LINE: %s", lin_buf);
	 THROW_TPR_ERROR(resname); }
 moli->n_res = nresidues;
 moli->resname = ( tW_word *) ecalloc(nresidues, sizeof(tW_word));
 for (int i = 0; i < nresidues; i++) // Scan for resnames
 {
  get_next_line(fp, *lin_buf);
  bOk = sscanf(*lin_buf," residue[%d]={name=\"%s\"", &residue_idx, &resname);
		  //, nr=%d, ic\'%s\'} ",&residue_idx, &resname, &dumi, dum_word);
  //fprintf(stderr, "bOk %d string: %s\n", bOk, lin_buf);
  if (bOk != 2) { THROW_TPR_ERROR("RESIDUE INFO FMT CHANGE\n");}
  if (residue_idx != i) { THROW_TPR_ERROR("RESIDUE READER MISALIGNED\n"); }
  elim_char(resname, '"');
  strcpy(moli->resname[i], resname);
 }
 return TRUE;
}
/**************************************************************************************
dparse_charge_groups_per_molecule:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This function is meant to parse the section of text following moltype(i)-> cgs:
It fills in the number of charge groups/members involved for the passed molecule structure
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note newer tprs do not have charge groups. Also charge groups are blocked by contiguous atoms with the same type in the .top file. (not the indices passed into the .top e.g. 1 2 1 2 is 4 dif charge groups)
@note This function is a helper function @see dfind_all_molinfo
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
bool dparse_charge_groups_per_molecule(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeid)
{
 bool bOk;
 int n_cg, cg_idx, idxatomstart, idxatomstop;
 tW_word moltype_header;
 sprintf(moltype_header, "moltype (%d):", moltypeid);
 find_line_in_file(fp, moltype_header, NULL); // go to correct section
 bOk = find_line_in_file(fp, "cgs:", NULL); // see if found
 if (bOk) // parse ncgs found
 {
   get_next_line(fp,*lin_buf); // nr=X
   bOk = sscanf(*lin_buf," nr=%d ",&n_cg);
   if (!bOk) { THROW_TPR_ERROR("Found charge group info but could not parse cgs\n");} 
 }
 if(!bOk) { n_cg = 0; } // then there are none to find
 moli->n_cg = n_cg;

 moli->cg_start = (int *) ecalloc(n_cg, sizeof(int));
 moli->cg_end = (int *) ecalloc(n_cg, sizeof(int));

 for (int i = 0; i < n_cg; i++) // Scan for cg members cgs[X]={start..stop}
 {
   get_next_line(fp, *lin_buf);
   bOk = sscanf(*lin_buf," cgs[%d]={%d..%d} ",&cg_idx, &idxatomstart, &idxatomstop);
   if (bOk !=3) { THROW_TPR_ERROR("charge group iatom fmt change\n"); }
   moli->cg_start[i] = idxatomstart;
   moli->cg_end[i] = idxatomstop;
 }
 return TRUE;
}
/**************************************************************************************
dparse_atomic_exclusions_per_molecule:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This function is meant to parse the section of text following moltype(i)-> excls:
It is *VERY* important for cgff/ldcalc that this function works.
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note older and newer tprs will format the exclusion list differently
@note This function is a helper function @see dfind_all_molinfo
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
bool dparse_atomic_exclusions_per_molecule(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeid)
{
 bool bOk;
 bool bnew_excl_fmt = FALSE;
 int nr_exclusion_lists, total_members_across_lists, lstid;
 int exc_for_ati_startidx, exc_for_ati_stopidx, lenlist_ati;
 char* excl_mem_parser;
 char* excl_mem_parser2;

 int members_found, member_id;
 int EOF_sentinal;
 
 tW_word moltype_header;
 sprintf(moltype_header, "moltype (%d):", moltypeid);
 find_line_in_file(fp, moltype_header, NULL); // go to correct section
 bOk = find_line_in_file(fp, "excls:", NULL); // see if found
 if (!bOk) { THROW_TPR_ERROR("TRIED TO PARSE EXCLUSIONS but info not found excls: SECTION\n"); }
 get_next_line(fp,*lin_buf); // nr=X or numLists=X
 //fprintf(stderr, "found %s\n", *lin_buf); 
 bOk = sscanf(*lin_buf," nr=%d ",&nr_exclusion_lists);
 if (!bOk) 
 { 
   bOk = sscanf(*lin_buf, " numLists=%d", &nr_exclusion_lists);
   if (!bOk) {THROW_TPR_ERROR("nr=X/numLists=X NOT FOUND UNDER exclusions\n"); }
 }

 moli->n_excls = nr_exclusion_lists;
 if (nr_exclusion_lists != moli->n_apm)
 { THROW_TPR_ERROR("NUMBER OF ATOMS THAT HAVE EXCLUSION LISTS != NUMBER ATOMS IN MOLECULE\n"); }
 get_next_line(fp, *lin_buf); // nra=X or numElements
 bOk = sscanf(*lin_buf, " nra=%d ", &total_members_across_lists);
 if (!bOk)
 {
   bOk = sscanf(*lin_buf, " numElements=%d ", &total_members_across_lists); 
   if (!bOk) {THROW_TPR_ERROR("nra=X/numElements=X NOT FOUND UNDER nr-X:\n");}
   bnew_excl_fmt= TRUE;
 }
 moli->n_exclsa = total_members_across_lists; // exclusions are a jagged array

 // Now we can allocate space for each exclusion list
 moli->excls = (int **) ecalloc(nr_exclusion_lists, sizeof(int *));
 moli->n_epa = (int *) ecalloc(nr_exclusion_lists, sizeof(int *)); // will hold each list len

 for (int i = 0; i < nr_exclusion_lists; i++) // Scan for list members & parse
 {
   //Determine number of entries in each list
   get_next_line(fp, *lin_buf);
   bOk = sscanf(*lin_buf, " excls[%d][%d..%d]", &lstid, &exc_for_ati_startidx, &exc_for_ati_stopidx);
   if (bOk == 3) // old fmt
   {
     lenlist_ati = exc_for_ati_stopidx - exc_for_ati_startidx + 1;
   }
   if (bOk !=3) // new fmt
   { 
     bOk = sscanf(*lin_buf, " excls[%d][num=%d]", &lstid, &lenlist_ati);
     if (bOk !=2) {THROW_TPR_ERROR("UNEXPECTED EXCLUSION LIST FMT\n"); }
   }
   if (lstid != i) { THROW_TPR_ERROR("EXCLUSION PARSER MISALIGNED\n"); }
   moli->n_epa[i] = lenlist_ati;

   // Make space & fill
   moli->excls[i] = (int *) ecalloc(lenlist_ati, sizeof(int));

   // Expect lines like excls[X][start..stop]{atomid1, atomid2}, gotta store and parse the brackets
   excl_mem_parser = strstr(*lin_buf, "=");
   if (bnew_excl_fmt) // then we gotta split by = twice
   { 
     excl_mem_parser2 = excl_mem_parser+1; 
     excl_mem_parser = strstr(excl_mem_parser2, "=");
   }
   if (DEBUG_BUILD) {fprintf(stderr, "EXCL PARSER %s", excl_mem_parser);}
   members_found = 0;
   while (members_found < lenlist_ati)
   {
     bOk = sscanf(excl_mem_parser, "%d", &member_id);
     if (!bOk) // Then keep searching string subsections until it is okay
     {
       excl_mem_parser = (excl_mem_parser + 1); // address of next char in
       if ( (strstr(excl_mem_parser, ",") == NULL) &&
	  (strstr(excl_mem_parser, "}") == NULL) &&
	  (members_found < lenlist_ati) ) // { } This list can take multiple lines
          {
            EOF_sentinal = get_next_line(fp, *lin_buf);
	    if (EOF_sentinal == -1 ) { THROW_TPR_ERROR("Unable to complete reading excl list\n"); }
	    excl_mem_parser = *lin_buf;
          }
     }
     else // we probably found an id
     {
       moli->excls[i][members_found] = member_id;
       members_found++;
       excl_mem_parser = strstr(excl_mem_parser, ","); // split remaining line by ,
     }
   }
 }
 return TRUE;
}

/**************************************************************************************
dparse_Bonded_int_lists_per_molecule:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
@arg moltypeidx the index of the molecule type we should fill this info for
@return TRUE if successful FALSE (something bad happened)
This is a wrapper function.
The goal of this function is to generate a list of types of intra (Bonded) interactions in the topology
and the indices of the molecule's atoms involved in each.
It is *VERY* important for cgff/ldcalc that this function works.
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note This function is a helper function @see dfind_all_molinfo
**************************************************************************************/
bool dparse_bonded_int_lists_per_molecule(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeidx)
{
  int n_found;
  
  n_found = dfill_bonds_from_gmx_bondtype_headers(fp, moli, lin_buf, moltypeidx);
  if (DEBUG_BUILD) { fprintf(stderr, "Found %d bonds in moleculetype[%d]\n", n_found, moltypeidx);}
  n_found = dfill_angles_from_gmx_angletype_headers(fp, moli, lin_buf, moltypeidx);
  if (DEBUG_BUILD) {fprintf(stderr, "Found %d angles in moleculetype[%d]\n", n_found, moltypeidx);}
  n_found = dfill_dihedrals_from_gmx_dihedraltype_headers(fp, moli, lin_buf, moltypeidx);
  if (DEBUG_BUILD) {fprintf(stderr, "Found %d dihedrals in moleculetype[%d]\n", n_found, moltypeidx);}
  n_found = dfill_bigdihedrals_from_gmx_bigdihedraltype_headers(fp, moli, lin_buf, moltypeidx);
  if (DEBUG_BUILD) {fprintf(stderr, "Found %d bigdihedrals in moleculetype[%d]\n", n_found, moltypeidx);}
  n_found = dfill_14nb_from_gmx_IMNB_headers(fp, moli, lin_buf, moltypeidx);
  if (DEBUG_BUILD) {fprintf(stderr, "Found %d 14nb interactions in moleculetype[%d]\n", n_found, moltypeidx);}
  n_found =  dproc_other_headers(fp, moli, lin_buf, moltypeidx);

  return TRUE;
}

/**************************************************************************************
dfill_bonds_from_gmx_bondtype_headers:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This is one of the functions we expect to update over time
It stores a list of gromacs headers that signal bond member lists and fills the interactions
in the topology
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note older and newer tprs will format the exclusion list differently
@note This function is a helper function to  @see dfind_all_molinfo -> @see dparse_bond_int_lists_per_molecule
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
int dfill_bonds_from_gmx_bondtype_headers(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeid)
{
  bool bOk;

  int nBondTypes = 10; // increase this and add headers for new gromacs versions
  const tW_word BondTypes[] = {"Bond:","G96Bond:","Morse:","Cubic Bonds:","Connect Bonds:",
	                       "Harmonic Pot.:","FENE Bonds:","Tab. Bonds:","Tab. Bonds NC:",
			       "Restraint Pot.:"};
  int integer_holder;
  int change_in_bond_count, prev_n_bonds, n_words;
  int BIDX = 0;

  moli->bond_nr = 0; // total entries for bond lists (memberids + interaction indices)
  moli->n_bonds = 0; // Instances of bonds
  tW_word moltype_header;
  sprintf(moltype_header, "moltype (%d):", moltypeid);

  for (int i = 0; i < nBondTypes; i++) // Scan all known bond headers
  {
   find_line_in_file(fp, moltype_header, NULL); //make sure we're always in the right section to start
   find_line_in_file(fp, BondTypes[i], NULL); // If exists find next header instance
   get_next_line(fp, *lin_buf); // nr: X should always be under the header
   bOk = sscanf(*lin_buf, " nr: %d ", &integer_holder); //space needed to hold n_members + int_idx
   if (!bOk) { THROW_TPR_ERROR("Could not find nr: X under Bondtype header\n"); }
   if (integer_holder == 0) { continue; }

   // Mark prior found info and figure out new space requirements
   prev_n_bonds = moli->n_bonds;
   moli->bond_nr += integer_holder;
   moli->n_bonds += integer_holder/BOND_DIV; // 3 entries, 2 members + 1 idx
   change_in_bond_count = integer_holder/BOND_DIV;

   // Make room
   moli->bond_types = (int *) erealloc(moli->bond_types, moli->n_bonds * sizeof(int)); //list of instances 
   moli->bond_type_top_cat = (int *) erealloc(moli->bond_type_top_cat, moli->n_bonds * sizeof(int));
   moli->bond_ij = (int **) erealloc(moli->bond_ij, moli->n_bonds * sizeof(int *)); // members per inst only

   for (int j = prev_n_bonds; j < moli->n_bonds; j++) {moli->bond_ij[j] = (int *) ecalloc(2,sizeof(int));}

   // skip header
   get_next_line(fp, *lin_buf); // if != 0, iatoms:

   for (int j = 0; j < change_in_bond_count; j++) 
   {
     get_next_line(fp,*lin_buf);
     elim_char(*lin_buf,'\n');
     n_words = get_word_count_delim(*lin_buf," ");
     tW_word * word_list = (tW_word *) ecalloc(n_words,sizeof(tW_word));
     get_words_delim(*lin_buf," ",word_list);

     bOk = sscanf(word_list[1],"type=%d",&integer_holder);
     moli->bond_types[BIDX] = integer_holder;
     moli->bond_type_top_cat[BIDX] = i;
     moli->bond_ij[BIDX][0] = atoi(word_list[3]);
     moli->bond_ij[BIDX][1] = atoi(word_list[4]);
     ++BIDX;
     efree(word_list);
   }
  }

  return moli->n_bonds;
}
/***************************************************************************************
dfill_angles_from_gmx_angletype_headers:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This is one of the functions we expect to update over time
It stores a list of gromacs headers that signal angle member lists and fills the interactions
in the topology
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note older and newer tprs will format the exclusion list differently
@note This function is a helper function to  @see dfind_all_molinfo -> @see dparse_bond_int_lists_per_molecule
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
int dfill_angles_from_gmx_angletype_headers(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeid)
{
  tW_word moltype_header;
  sprintf(moltype_header, "moltype (%d):", moltypeid);
  bool bOk;
  int AIDX = 0;
  int nAngleTypes = 10;
  const tW_word AngleTypes[] = {"Angle:","G96Angle:","Restricted Angles:","Lin. Angle:",
	                        "Bond-Cross:","BA-Cross:","U-B:","Quartic Angles:","Tab. Angles:","Restr. Angles:"};
  moli->angle_nr = moli->n_angles = 0;
  int integer_holder;
  int prev_n_angles, n_new_angles;

  for (int i = 0; i < nAngleTypes; i++)
  {
    find_line_in_file(fp, moltype_header, NULL); //make sure we're always in the right section to start
    if (!find_line_in_file(fp, AngleTypes[i], NULL)) { continue; } // If exists find next header instance, if not, continue to next one
    get_next_line(fp, *lin_buf); // nr: X should always be under the header
    bOk = sscanf(*lin_buf, " nr: %d ", &integer_holder); //space needed to hold n_members + int_idx
    if (!bOk) { THROW_TPR_ERROR("Could not find nr: X under Angletype header\n"); }
    if (integer_holder == 0) { continue; }
    get_next_line(fp,*lin_buf); // iatoms:
    moli->angle_nr += integer_holder;
    prev_n_angles = moli->n_angles;
    moli->n_angles += integer_holder/ANGLE_DIV;
    n_new_angles = integer_holder/ANGLE_DIV;
    moli->angle_types = (int *) erealloc(moli->angle_types, moli->n_angles * sizeof(int));
    moli->angle_type_top_cat = (int *) erealloc(moli->angle_type_top_cat, moli->n_angles * sizeof(int));
    moli->angle_ijk = (int **) erealloc(moli->angle_ijk, moli->n_angles * sizeof(int *));
    for (int j = prev_n_angles; j < moli->n_angles; ++j) {moli->angle_ijk[j] = (int *) ecalloc(3,sizeof(int));}
    for (int j = 0; j < n_new_angles; ++j)
    {
      get_next_line(fp,*lin_buf);
      
      elim_char(*lin_buf,'\n');
      int n_words = get_word_count_delim(*lin_buf," ");
      tW_word * word_list = (tW_word *) ecalloc(n_words,sizeof(tW_word));
      get_words_delim(*lin_buf," ",word_list);
      
      bOk = sscanf(word_list[1],"type=%d",&integer_holder);
      moli->angle_types[AIDX] = integer_holder;
      moli->angle_type_top_cat[AIDX] = i;
      moli->angle_ijk[AIDX][0] = atoi(word_list[3]);
      moli->angle_ijk[AIDX][1] = atoi(word_list[4]);
      moli->angle_ijk[AIDX][2] = atoi(word_list[5]);
      ++AIDX;
      efree(word_list);
    }
  }
  return moli->n_angles;
}

/***************************************************************************************
dfill_dihedrals_from_gmx_dihedraltype_headers:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
@return TRUE if successful FALSE (something bad happened)
This is one of the functions we expect to update over time
It stores a list of gromacs headers that signal dihedral member lists and fills the interactions
in the topology
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note older and newer tprs will format the exclusion list differently
@note This function is a helper function to  @see dfind_all_molinfo -> @see dparse_bond_int_lists_per_molecule
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
int dfill_dihedrals_from_gmx_dihedraltype_headers(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeid)
{
  tW_word moltype_header;
  sprintf(moltype_header, "moltype (%d):", moltypeid);
  bool bOk;
  bool bsecond_Improper_Dih_head = FALSE;
  bool checked_for_second = FALSE;
  int DIDX = 0;
  int  nDihedralTypes = 9;
  const tW_word DihedralTypes[] = {"Proper Dih.:","Ryckaert-Bell.:", "Restr. Dih.:","Restricted Dih.:",
	                           "CBT Dih.:","Fourier Dih.:","Improper Dih.:","Per. Imp. Dih.:","Tab. Dih.:"}; 

  moli->dih_nr = moli->n_dihs = 0;
  int prev_n_dihs, n_new_dihs, n_words, integer_holder;

  int improp_i = -1; // Default index of improper dihedral type in DihedralTypes array
  for (int i = 0; i < nDihedralTypes; i++) 
  { 
    if (strcmp(DihedralTypes[i], "Improper Dih.:") == 0)
    { improp_i = i; 
      if (DEBUG_BUILD) printf("Found Improper Dih. at index %d\n", improp_i);
    }
  }  

  for (int i = 0; i < nDihedralTypes; i++)
  {
    find_line_in_file(fp, moltype_header, NULL); //make sure we're always in the right section to start
    if(!find_line_in_file(fp, DihedralTypes[i], NULL)){continue;} // If exists find next header instance, if not, continue to next one
    if (bsecond_Improper_Dih_head && (i == improp_i)) // Check if processing Improper Dih. The integer check is more straightforward to implement than the string
    {find_line_in_file(fp, DihedralTypes[improp_i], NULL);} // Go to second header if needed
    get_next_line(fp, *lin_buf); // nr: X should always be under the header
    bOk = sscanf(*lin_buf, " nr: %d ", &integer_holder); //space needed to hold n_members + int_idx
    if (!bOk) { THROW_TPR_ERROR("Could not find nr: X under Dihedraltype header\n"); }
    if (integer_holder == 0) 
    { 
      if (i != improp_i) {continue;} // Check for second if first improp was 0
      if (bsecond_Improper_Dih_head) { continue; }
      else { get_next_line(fp, *lin_buf); 
	     bsecond_Improper_Dih_head=((strstr(*lin_buf, DihedralTypes[i]) != NULL));
	     if (!checked_for_second){ i=i-1;}
      	     checked_for_second=TRUE;
           } 
    }
    get_next_line(fp, *lin_buf); // iatoms:

    moli->dih_nr += integer_holder;
//    fprintf(stderr, "Got %d dih for %s of type %s\n", integer_holder, moltype_header, DihedralTypes[i]);
    prev_n_dihs = moli->n_dihs;
    moli->n_dihs += integer_holder/DIH_DIV;
    n_new_dihs = integer_holder/DIH_DIV;
    moli->dih_types = (int *) erealloc(moli->dih_types, moli->n_dihs * sizeof(int));
    moli->dih_type_top_cat = (int *) erealloc(moli->dih_type_top_cat, moli->n_dihs * sizeof(int));
    moli->dih_ijkl = (int **) erealloc(moli->dih_ijkl, moli->n_dihs * sizeof(int *));
    for (int j = prev_n_dihs; j < moli->n_dihs; ++j) { moli->dih_ijkl[j] = (int *) ecalloc(4,sizeof(int)); }
    for (int j = 0; j < n_new_dihs; ++j)
    {
      get_next_line(fp,*lin_buf);
     // fprintf(stderr, "%s\n", lin_buf);
      elim_char(*lin_buf,'\n');
      n_words = get_word_count_delim(*lin_buf," ");
      tW_word * word_list = (tW_word *) ecalloc(n_words,sizeof(tW_word));
      get_words_delim(*lin_buf," ",word_list);
     
      bOk = sscanf(word_list[1],"type=%d",&integer_holder);
      moli->dih_types[DIDX] = integer_holder;
      moli->dih_type_top_cat[DIDX] = i;
      moli->dih_ijkl[DIDX][0] = atoi(word_list[3]);
      moli->dih_ijkl[DIDX][1] = atoi(word_list[4]);
      moli->dih_ijkl[DIDX][2] = atoi(word_list[5]);
      moli->dih_ijkl[DIDX][3] = atoi(word_list[6]);
      ++DIDX;
      efree(word_list);
    }
    if (i == improp_i) // Then we need to check if we're in the version regime where Improp Dih is listed 2x per molecule.
    {
      get_next_line(fp, *lin_buf); // It will be one line down from the first one if its here.
      if ((strstr(*lin_buf, "Improper Dih.:") != NULL) && !bsecond_Improper_Dih_head)
      {
//        fprintf(stderr, "In the improper check\n");
	bsecond_Improper_Dih_head = TRUE;
	i=i-1; // We'll process this kind of header again if we need to (but not if we've done this already)
      }
    }
  }
  return moli->n_dihs;
}

/***************************************************************************************
dfill_bigdihedrals_from_gmx_bigdihedraltype_headers:
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
This is one of the functions we expect to update over time
It stores a list of gromacs headers that signal bigdihedral (Charm) member lists and fills the interactions
in the topology
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note older and newer tprs will format the exclusion list differently
@note This function is a helper function to  @see dfind_all_molinfo -> @see dparse_bond_int_lists_per_molecule
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
int dfill_bigdihedrals_from_gmx_bigdihedraltype_headers(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeid)
{
  
  tW_word moltype_header;
  sprintf(moltype_header, "moltype (%d):", moltypeid);
  bool bOk;

  const tW_word BigDihedralTypes[] = {"CMAP Dih.:"}; // MRD 01.24.2020
  int  nBigDihedralTypes = 1;
  int BDIDX, integer_holder;
  int prev_n_bigdihs, n_new_bigdihs, n_words;
  moli->n_bigdihs = 0;
  for (int i=0; i < nBigDihedralTypes; i++)
  {
    find_line_in_file(fp, moltype_header, NULL); //make sure we're always in the right section to start
    find_line_in_file(fp, BigDihedralTypes[i], NULL); // If exists find next header instance
    get_next_line(fp, *lin_buf); // nr: X should always be under the header
    bOk = sscanf(*lin_buf, " nr: %d ", &integer_holder); //space needed to hold n_members + int_idx
    if (!bOk) { THROW_TPR_ERROR("Could not find nr: X under BigDihedraltype header\n"); }
    if (integer_holder == 0) { continue; }
    get_next_line(fp,*lin_buf); // iatoms:

    moli->bigdih_nr += integer_holder;
    prev_n_bigdihs = moli->n_bigdihs;
    moli->n_bigdihs += integer_holder/BIGDIH_DIV;
    n_new_bigdihs = integer_holder/BIGDIH_DIV;
    moli->bigdih_types = (int *) erealloc(moli->bigdih_types, moli->n_bigdihs * sizeof(int));
    moli->bigdih_type_top_cat = (int *) erealloc(moli->bigdih_type_top_cat, moli->n_bigdihs * sizeof(int));
    moli->bigdih_ijklm = (int **) erealloc(moli->bigdih_ijklm, moli->n_bigdihs * sizeof(int *));
    for (int j = prev_n_bigdihs; j < moli->n_bigdihs; ++j) { moli->bigdih_ijklm[j] = (int *) ecalloc(5,sizeof(int)); }
    for (int j = 0; j < n_new_bigdihs; ++j)
    {
      get_next_line(fp,*lin_buf);
     
      elim_char(*lin_buf,'\n');
      n_words = get_word_count_delim(*lin_buf," ");
      tW_word * word_list = (tW_word *) ecalloc(n_words,sizeof(tW_word));
      get_words_delim(*lin_buf," ",word_list);
     
      bOk = sscanf(word_list[1],"type=%d",&integer_holder);
      moli->bigdih_types[BDIDX] = integer_holder;
      moli->bigdih_type_top_cat[BDIDX] = i;
      moli->bigdih_ijklm[BDIDX][0] = atoi(word_list[3]);
      moli->bigdih_ijklm[BDIDX][1] = atoi(word_list[4]);
      moli->bigdih_ijklm[BDIDX][2] = atoi(word_list[5]);
      moli->bigdih_ijklm[BDIDX][3] = atoi(word_list[6]);
      moli->bigdih_ijklm[BDIDX][4] = atoi(word_list[7]);
      ++BDIDX;
      efree(word_list);
    }
  
  }
  return moli->n_bigdihs;
}

/****************************************************************************
dfill_14nb_from_gmx_IMNB_headers
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
This is one of the functions we expect to update over time
It stores a list of gromacs headers that signal bigdihedral (Charm) member lists and fills the interactions
in the topology
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note older and newer tprs will format the exclusion list differently
@note This function is a helper function to  @see dfind_all_molinfo -> @see dparse_bond_int_lists_per_molecule
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
int dfill_14nb_from_gmx_IMNB_headers(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeid)
{
  tW_word moltype_header;
  sprintf(moltype_header, "moltype (%d):", moltypeid);
  const tW_word IMNBTypes[] = {"LJ-14:","Coulomb-14:","LJC-14 q:","LJC Pairs NB:"};
  int nIMNBTypes = 4;
  int IMNBIDX = 0;
  int integer_holder, prev_n_imnbs, n_words, n_new_imnbs;
  moli->imnb_nr = moli->n_imnbs = 0;
  bool bOk;
  
   for (int i=0; i < nIMNBTypes; i++)
  {
    find_line_in_file(fp, moltype_header, NULL); //make sure we're always in the right section to start
    find_line_in_file(fp, IMNBTypes[i], NULL); // If exists find next header instance
    get_next_line(fp, *lin_buf); // nr: X should always be under the header
    bOk = sscanf(*lin_buf, " nr: %d ", &integer_holder); //space needed to hold n_members + int_idx
    if (!bOk) { THROW_TPR_ERROR("Could not find nr: X under 14NB header\n"); }
    if (integer_holder == 0) { continue; }
    get_next_line(fp,*lin_buf); // iatoms:
    moli->imnb_nr += integer_holder;
    prev_n_imnbs = moli->n_imnbs;
    moli->n_imnbs += integer_holder / IMNB_DIV;
    n_new_imnbs = integer_holder / IMNB_DIV;
    moli->imnb_types = (int *) erealloc(moli->imnb_types, moli->n_imnbs * sizeof(int));
    moli->imnb_type_top_cat = (int *) erealloc(moli->imnb_type_top_cat, moli->n_imnbs * sizeof(int));
    moli->imnb_ij = (int **) erealloc(moli->imnb_ij, moli->n_imnbs * sizeof(int *));
    for (int j = prev_n_imnbs; j < moli->n_imnbs; ++j) { moli->imnb_ij[j] = (int *) ecalloc(2,sizeof(int)); }
    for (int j = 0; j < n_new_imnbs; ++j)
    {
      get_next_line(fp,*lin_buf);
      elim_char(*lin_buf,'\n');
      n_words = get_word_count_delim(*lin_buf," ");
      tW_word * word_list = (tW_word *) ecalloc(n_words, sizeof(tW_word));
      get_words_delim(*lin_buf," ",word_list);
    
      bOk = sscanf(word_list[1],"type=%d",&integer_holder);
      moli->imnb_types[IMNBIDX] = integer_holder;
      moli->imnb_type_top_cat[IMNBIDX] = i;
      moli->imnb_ij[IMNBIDX][0] = atoi(word_list[3]);
      moli->imnb_ij[IMNBIDX][1] = atoi(word_list[4]);
      ++IMNBIDX;
      efree(word_list);
    } 
  }
   return moli->n_imnbs;
}

/****************************************************************************
dproc_other_headers
@arg fp the pointer to the dumped tpr file
@arg moli molecule_typei
@arg lin_buf allocated space for working with lines (something like address of (char*))
This is one of the functions we expect to update over time
It stores a list of gromacs headers that contain interaction lists BOCS doesnt handle
@note in legacy versions of BOCS this functionality was handled one at a time in get_moltype_info
@note older and newer tprs will format the exclusion list differently
@note This function is a helper function to  @see dfind_all_molinfo -> @see dparse_bond_int_lists_per_molecule
@warning This function dynamically allocates memory in the molecule structure (in top structure)
**************************************************************************************/
int dproc_other_headers(FILE* fp, tW_molecule* moli, tW_line* lin_buf, int moltypeid)
{
  tW_word moltype_header;
  sprintf(moltype_header, "moltype (%d):", moltypeid);
  const tW_word OtherTypes[] = {"GB 1-2 Pol.:","GB 1-2 Pol. (unused):","GB 1-3 Pol.:",
              "GB 1-3 Pol. (unused):","GB 1-4 Pol. (unused):", "GB 1-4 Pol.:",
  	          "GB Polariz.:","GB Polarization (unused):","Nonpolar Sol.:","Nonpolar Sol. (unused):",
              "LJ (SR):","LJ:","B.ham:","Coul:", "Flat-b. P-R.:", "Dih. Rest. Vi.:",
			        "Buck.ham (SR):","LJ (unused):","B.ham (unused):","Disper. corr.:",
			        "Coulomb (SR):","Coul (unused):","RF excl.:","Coul. recip.:",
			        "LJ recip.:","DPD:","Polarization:","Water Pol.:","Thole Pol.:",
			        "Anharm. Pol.:","Position Rest.:","Flat-bottom posres:","Dis. Rest.:",
			        "D.R.Viol. (nm):","Orient. Rest.:","Ori. R. RMSD:","Angle Rest.:",
			        "Angle Rest. Z:","Dih. Rest.:","Dih. Rest. Viol.:","Constraint:",
              "Constr. No Co.:", "Virtual site 1:", "Virt. site 2fd:", "Virt. site 3fd:",
			        "Constr. No Conn.:","Settle:","Virtual site 2:","Virtual site 3:",
			        "Virtual site 3fd:","Virtual site 3fad:","Virtual site 3out:",
			        "Virtual site 4fd:","Virtual site 4fdn:","Virtual site N:",
              "Vir. site 3fad:","Vir. site 3out:","Virt. site 4fd:","Vir. site 4fdn:",
			        "COM Pull En.:","Quantum En.:","Potential:","Kinetic En.:",
			        "Total Energy:","Conserved En.:","Temperature:","Vir. Temp. (not used):",
			        "Pres. DC:","Pressure:","dH/dl constr.:","dVremain/dl:","dEkin/dl:",
			        "dVcoul/dl:","dVvdw/dl:","dVbonded/dl:","dVrestraint/dl:",
			        "dVtemperature/dl:","dVtemp/dl:","Vir. Temp.:", "Dens. fitting:"};
  int nOtherTypes = 59;
  int IMNBIDX = 0;
  bool bOk;
  int integer_holder, prev_n_imnbs, n_words;
  
  moli->n_other = (int *) ecalloc(nOtherTypes,sizeof(int));
  moli->other_nr = (int *) ecalloc(nOtherTypes,sizeof(int));
  moli->other_types = (int **) ecalloc(nOtherTypes,sizeof(int *));
  moli->other_ijklmn = (int ***) ecalloc(nOtherTypes,sizeof(int **));

  for (int i=0; i < nOtherTypes; i++)
  {
    find_line_in_file(fp, moltype_header, NULL); //make sure we're always in the right section to start
    find_line_in_file(fp, OtherTypes[i], NULL); // If exists find next header instance
    get_next_line(fp, *lin_buf); // nr: X should always be under the header
    bOk = sscanf(*lin_buf, " nr: %d ", &integer_holder); //space needed to hold n_members + int_idx
    if (!bOk) { THROW_TPR_ERROR("Could not find nr: X under 14NB header\n"); }
    if (integer_holder == 0) { continue; }
    get_next_line(fp, *lin_buf); // iatoms:
    fprintf(stderr,"WARNING: nonzero number of interactions of type %s found\n",OtherTypes[i]);
    fprintf(stderr,"These types of interactions do not correspond to BondStretch, Angle, Dihedral, nor IntraMolec_NB_Pair types\n");
    fprintf(stderr,"Accordingly, they will not be used nor transferred to the .btp file\n");
    get_next_line(fp,*lin_buf); // iatoms:
    get_next_line(fp,*lin_buf);
    // split this line up into words. is the number of words minus 2
    // X type=XX (XXXX) I J K ...
    elim_char(*lin_buf,'\n');
    n_words = get_word_count_delim(*lin_buf," ");
    tW_word * word_list = (tW_word *) ecalloc(n_words,sizeof(tW_word));
    get_words_delim(*lin_buf," ",word_list);

    moli->n_other[i] = moli->other_nr[i] / (n_words - 2);
    moli->other_types[i] = (int *) ecalloc(moli->n_other[i], sizeof(int));
    moli->other_ijklmn[i] = (int **) ecalloc(moli->n_other[i], sizeof(int *));
    for (int j = 0; j < moli->n_other[i]; ++j)
    { moli->other_ijklmn[i][j] = (int *) ecalloc(n_words-3,sizeof(int)); }

    bOk = sscanf(word_list[1],"type=%d",&integer_holder);
    moli->other_types[i][0] = integer_holder;
    for (int j = 3; j < n_words; ++j) { moli->other_ijklmn[i][0][j-3] = atoi(word_list[j]); }
    efree(word_list);
    for (int j = 1; j < moli->n_other[i]; ++j)
    {
      get_next_line(fp,*lin_buf);
      elim_char(*lin_buf,'\n');
      int n_words2 = get_word_count_delim(*lin_buf," ");
      tW_word * word_list2 = (tW_word *) ecalloc(n_words2,sizeof(tW_word));
      get_words_delim(*lin_buf," ",word_list2);
      moli->other_types[i][j] = integer_holder;
      for (int k = 3; k < n_words; ++k) { moli->other_ijklmn[i][j][k-3] = atoi(word_list2[k]); }
      efree(word_list2);
     }
  }
  // Need to consider how to scan for unknown headers
  return 0;
}
