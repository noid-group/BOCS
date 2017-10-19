/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file gmx-interface_46.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
@brief Functions to interface cgff and cgmap with the gmx library version 4.6
*/

//c library includes
#include <stdio.h>
#include <string.h>

//local includes
#include "gmx-interface.h"

#include "safe_mem.h"
#include "wnoid_math.h"
#include "io_read.h"


/*****************************************************************************************
copyright(): Prints the gromacs copyright and credits output
*****************************************************************************************/
void copyright(tW_gmx_output *self)
{
    CopyRight(self->outstream, self->prog);
}


/*****************************************************************************************
thanks(): Prints the gromacs funny end quote
*****************************************************************************************/
void thanks(tW_gmx_output *self)
{
    thanx(self->outstream);
}


/*****************************************************************************************
init_tW_gmx_output(): Initializes the interface to the gromacs output-to-screen functions
*****************************************************************************************/
tW_gmx_output* init_tW_gmx_output(FILE *outstream, const char *prog)
{
    tW_gmx_output *out;

    out = emalloc(sizeof(tW_gmx_output));

    out->outstream = outstream;
    out->prog = strdup(prog);

    out->copyright = copyright;
    out->thanks = thanks;

    return out;
}


/*****************************************************************************************
get_natoms(): Returns the number of atoms specified by topology 'self'
*****************************************************************************************/
int get_natoms(tW_gmx_topology *self)
{
    //top.atoms.nr

    return self->contents->atoms.nr;
}

/*****************************************************************************************
get_natoms(): Returns the name of topology 'self'
*****************************************************************************************/
char *get_name(tW_gmx_topology *self)
{
    //*(top.name)

    return *(self->contents->name);
}


/*****************************************************************************************
get_nbonds(): Returns the number of bonds listed in self 
*****************************************************************************************/
int get_nbonds(tW_gmx_topology *self){
   return self->contents->idef.il[F_BONDS].nr;
}

/*****************************************************************************************
get_nangles(): Returns the number of angles listed in self 
*****************************************************************************************/
int get_nangles(tW_gmx_topology *self){
   return self->contents->idef.il[F_ANGLES].nr;
}

/*****************************************************************************************
get_ndihs(): Returns the number of dihedrals listed in self 
*****************************************************************************************/
int get_ndihs(tW_gmx_topology *self){
   return self->contents->idef.il[F_RBDIHS].nr;
}

/*****************************************************************************************
get_npairs(): Returns the number of intramolecular nonbonded pairs listed in self 
*****************************************************************************************/
int get_npairs(tW_gmx_topology *self){
   return self->contents->idef.il[F_LJ14].nr;
}

/*****************************************************************************************
get_bond_list() Returns a list of bonds in the t_ilist format:

+-----+--+--+-----+--+--+
|type1|i1|i2|type2|i1|i2|...
+-----+--+--+-----+--+--+
*****************************************************************************************/
int *get_bond_list(tW_gmx_topology *self){
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
int *get_angle_list(tW_gmx_topology *self){
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
int *get_dih_list(tW_gmx_topology *self){
   int i, nelements;
   int *dih_list;

   nelements = self->get_ndihs(self);

   dih_list = ecalloc(nelements, sizeof(int));

   for (i=0; i<nelements; i++)
   {
       dih_list[i] = self->contents->idef.il[F_RBDIHS].iatoms[i];
   }
   
   return dih_list;
}

/*****************************************************************************************
get_pair_list() Returns a list of nonbonded intramolecular pairs in the t_ilist format:

+-----+--+--+-----+--+--+
|type1|i1|i2|type2|i1|i2|...
+-----+--+--+-----+--+--+
*****************************************************************************************/
int *get_pair_list(tW_gmx_topology *self){
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
    //top.atoms.atom[i].resind
    
    return self->contents->atoms.atom[i_atom].resind;
}

/*****************************************************************************************
get_atomtype(): Returns the type of the i_atom'th atom in topology 'self'
*****************************************************************************************/
char *get_atomtype(tW_gmx_topology *self, int i_atom)
{

    //*(top.atoms.atomtype[i])

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
    if (nr_excl > 0) {
	/* Allocate memory for exclusion list for site i. */
	excl_list = (int *) ecalloc(nr_excl, sizeof(int));

	/* Copy exclusions (site indices) from top file for site i. */
	for (j = 0; j < nr_excl; j++) {
	    excl_list[j] = self->contents->excls.a[i_start + j];
        }

    }
    return excl_list;
}


/*****************************************************************************************
read_tpr(): Wrapper function for the tpr reading functionality in the gromacs library
*****************************************************************************************/
tW_gmx_topology* init_tW_gmx_topology()
{
	tW_gmx_topology *top;

	top = emalloc(sizeof(tW_gmx_topology));

	top->b_tpr = FALSE;

	top->get_natoms = get_natoms;
	top->get_name = get_name;
	top->get_resind = get_resind;
	top->get_atomtype = get_atomtype;
	top->get_nexcl = get_nexcl;
	top->get_excl_list = get_excl_list;

	top->get_nbonds = get_nbonds;
	top->get_nangles = get_nangles;
	top->get_ndihs = get_ndihs;
	top->get_npairs = get_npairs;

	top->get_bond_list = get_bond_list;
	top->get_angle_list = get_angle_list;
	top->get_dih_list = get_dih_list;
	top->get_pair_list = get_pair_list;

	return top;
}







/*****************************************************************************************
read_tpr(): Wrapper function for the tpr reading functionality in the gromacs library
*****************************************************************************************/
void read_tpr(tW_gmx_topology *top, const char *tpr_fnm)
{
    int ePBC;
    char *title;
    dvec *xtop;
    matrix box;

    title = ecalloc(100, sizeof(char));

    top->b_tpr = FALSE;

    top->contents = emalloc(sizeof(t_topology));
    top->b_tpr = read_tps_conf(tpr_fnm, title, top->contents, &ePBC, NULL, NULL, box, FALSE);

    efree(title);

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
    if (TPR_flag == TRUE) {

	if (file_exists(tpr_filename)) {

	    /* Get topology information. */
	    //bGromacs = read_tpr_file(tpr_filename, top, box);
	    read_tpr(top, tpr_filename);
	} else {
	    fprintf(stderr, "ERROR: tpr file '%s' not found. Check par.txt.\n", tpr_filename);
	    exit(EXIT_FAILURE);
	}

    } else {			/* do it the old way */

	auto_tpr_fnm = (char *) ecalloc(strlen(coord_fnm) + 1, sizeof(char));
	strcpy(auto_tpr_fnm, coord_fnm);
	trr_ptr = strstr(auto_tpr_fnm, ".trr\0");
	if (trr_ptr == NULL) {
	    auto_tpr_fnm = NULL;
	    return bGromacs;
	}
	strcpy(trr_ptr, ".tpr\0");

	if (file_exists(auto_tpr_fnm)) {
	    /* Get topology information. */

	    //bGromacs = read_tpr_file(auto_tpr_fnm, top, box);
	    read_tpr(top, auto_tpr_fnm);
	} else {
	    fprintf(stderr,
		    "ERROR: tpr file '%s' not found.  This file was generated automatically from the trr filename '%s'.\n",
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

    fprintf(fp_log, "  top->contents->name: %-10s\n\n", *(top->contents->name));
    fprintf(fp_log, "  No. of atomtypes: %-4d\n", top->contents->atomtypes.nr);
    fprintf(fp_log, "  No. of residues: %-4d\n", top->contents->atoms.nres);
    fprintf(fp_log, "  No. atoms: %-4d\n", top->contents->atoms.nr);
    fprintf(fp_log, "  Looping over atoms.\n");
    for (i = 0; i < top->contents->atoms.nr; i++) {
	//fprintf( fp_log, "    atm_no: %-4d   res_no: %-4d   ", i, top->contents->atoms.atom[i].resnr ); // 4.0.7
	fprintf(fp_log, "    atm_no: %-4d   res_no: %-4d   ", i, top->contents->atoms.atom[i].resind);	// 4.5.3

	//fprintf( fp_log, "atm_name: %-10s   atm_type: %-10s   type_id: %-4d   res_name: %-10s\n",
	//        *(top->contents->atoms.atomname[i]), *(top->contents->atoms.atomtype[i]), top->contents->atoms.atom[i].type, 
	fprintf(fp_log, "atm_name: %-10s   atm_type: %-10s   type_id: %-4d   res_name: %-10s\n", 
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
    for (i_nr = 0; i_nr <= top->contents->excls.nr; i_nr++) {
	fprintf(fp_log, "      top->contents->excls.index[%-4d]: %-4d\n", i_nr,
		top->contents->excls.index[i_nr]);
    }
    fprintf(fp_log,
	    "    Looping over top->contents->excls.nra (total no. of excl. indices).\n");
    for (i_nr = 0; i_nr < top->contents->excls.nra; i_nr++) {
	fprintf(fp_log, "      top->contents->excls.a[%-4d]: %-4d\n", i_nr,
		top->contents->excls.a[i_nr]);
    }
    fprintf(fp_log, "    Printing excl. for each atom.\n");
    for (i_nr = 0; i_nr < top->contents->excls.nr; i_nr++) {
	i_start = top->contents->excls.index[i_nr];
	i_end = top->contents->excls.index[i_nr + 1] - 1;
	nr_excl = i_end - i_start + 1;

	fprintf(fp_log, "      i_nr: %-4d   ", i_nr);
	fprintf(fp_log, "i_start: %-4d  i_end: %-4d  nr_excl: %-4d  ",
		i_start, i_end, nr_excl);

	for (i = i_start; i < i_end + 1; i++) {
	    fprintf(fp_log, "%-4d ", top->contents->excls.a[i]);
	}

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
    for (i = 0; i < top->contents->idef.il[F_BONDS].nr; i += 3) {
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
    for (i = 0; i < top->contents->idef.il[F_ANGLES].nr; i += 4) {
	fprintf(fp_log,
		"    top->contents->idef.il[F_ANGLES].iatoms[%-4d]  inter_type: %-4d  atoms: %-4d %-4d %-4d\n",
		i, top->contents->idef.il[F_ANGLES].iatoms[i],
		top->contents->idef.il[F_ANGLES].iatoms[i + 1],
		top->contents->idef.il[F_ANGLES].iatoms[i + 2],
		top->contents->idef.il[F_ANGLES].iatoms[i + 3]);
    }
    fprintf(fp_log, "\n");

    fprintf(fp_log, "  Printing F_RBDIHS interactions.\n");
    fprintf(fp_log, "    top->contents->idef.il[F_RBDIHS].nr: %-4d\n",
	    top->contents->idef.il[F_RBDIHS].nr);
    for (i = 0; i < top->contents->idef.il[F_RBDIHS].nr; i += 5) {
	fprintf(fp_log,
		"    top->contents->idef.il[F_RBDIHS].iatoms[%-4d]  inter_type: %-4d  atoms: %-4d %-4d %-4d %-4d\n",
		i, top->contents->idef.il[F_RBDIHS].iatoms[i],
		top->contents->idef.il[F_RBDIHS].iatoms[i + 1],
		top->contents->idef.il[F_RBDIHS].iatoms[i + 2],
		top->contents->idef.il[F_RBDIHS].iatoms[i + 3],
		top->contents->idef.il[F_RBDIHS].iatoms[i + 4]);
    }
    fprintf(fp_log, "\n");

    fprintf(fp_log, "  Printing F_LJ14 interactions.\n");
    fprintf(fp_log, "    top->contents->idef.il[F_LJ14].nr: %-4d\n",
	    top->contents->idef.il[F_LJ14].nr);
    for (i = 0; i < top->contents->idef.il[F_LJ14].nr; i += 3) {
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
    for (i = 0; i < top->contents->idef.ntypes; i++) {
	fprintf(fp_log, "    i: %-4d  name: %-10s \n", i,
		interaction_function[top->contents->idef.functype[i]].name);
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
    if (sys->TPR_EXCL_var.flag_TPR_excl == TRUE) {
	//get_top(" ", &top_excl, box, sys->TPR_EXCL_var.flag_TPR_excl,	sys->TPR_EXCL_var.TPR_excl);
	get_top(" ", top_excl, sys->TPR_EXCL_var.flag_TPR_excl, sys->TPR_EXCL_var.TPR_excl);
    }

    n_sites = top->get_natoms(top);

    for (i = 0; i < n_sites; i++) {
	/* Get site name from GROMACS topology for site i. */
	strcpy(CG_struct[i].name, top->get_atomtype(top, i) );

	/* Does the site name from top match any site names read from par.txt? */
	i_type = match_word(sys->N_Site_Types, CG_struct[i].name, sys->Site_Types);
	if (i_type == -1) {
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
	if (CG_struct[i].excl_list != NULL) {
	    efree(CG_struct[i].excl_list);
	}

	/* JFR - 06.27.12: for bondref, use top_excl for excl info */
	if (sys->TPR_EXCL_var.flag_TPR_excl == TRUE) {
	    /* Determine no. of exclusions for site i. */
  	    CG_struct[i].nr_excl = top_excl->get_nexcl(top_excl, i);
	    CG_struct[i].excl_list = top_excl->get_excl_list(top_excl, i);

	} else {		/* do it the usual way */

	    CG_struct[i].nr_excl = top->get_nexcl(top, i);
	    CG_struct[i].excl_list = top->get_excl_list(top, i);
	
	}

    }				/* End loop over sites */

}




/*****************************************************************************************
set_natoms(): sets the number of atoms in the frame self to natoms, and allocates memory
for those arrays that are flagged as present
*****************************************************************************************/
void set_natoms(tW_gmx_trxframe *self, int natoms)
{
    self->contents->natoms = natoms;

    if (self->contents->bX) {
	self->contents->x = ecalloc(natoms, 3*sizeof(real)); 
    }

    if (self->contents->bV) {
	self->contents->v = ecalloc(natoms, 3*sizeof(real));
    }

    if (self->contents->bF) {
	self->contents->f = ecalloc(natoms, 3*sizeof(real));
    }
}

/*****************************************************************************************
set_atom_labels(): sets the atom labels of the frame self if they are available in top
*****************************************************************************************/
void set_atom_labels(tW_gmx_trxframe *self, tW_gmx_topology *top)
{
    if (top->b_tpr) {
	self->contents->bAtoms = TRUE;
	self->contents->atoms = &(top->contents->atoms);
    } else {
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

    t_atom *atom_ptr;
    t_atoms *atoms_ptr;

    printf("--- in print_trr_file ---\n");

    if (self->contents->bTitle) {
	printf("title: %s \n", self->contents->title);
    }

    n_atms = self->contents->natoms;
    printf("fr.natoms: %d \n", n_atms);

    printf("starting time.  fr.t0:   %f \n", self->contents->t0);
    printf("previous time.  fr.tpf:  %f \n", self->contents->tpf);
    printf("prevprev time.  fr.tpf:  %f \n", self->contents->tppf);

    printf("present  step.  fr.step: %d \n", self->contents->step);
    printf("present  time.  fr.time: %f \n", self->contents->time);

    printf("present  lambda.  fr.lambda: %f \n", self->contents->lambda);

    printf("printing atoms xvf.\n");
    for (i_atm = 0; i_atm < n_atms; i_atm++) {
	printf("%i \t ", i_atm);

	if (self->contents->bAtoms) {
	    atoms_ptr = self->contents->atoms;
	    atom_ptr = &(atoms_ptr->atom[i_atm]);
	    res_no = atom_ptr->resind;
	    printf("type_id: %d \t res_no: %d \t atm_name: %s \t atm_type: %s \t res_name: %s \n", 
		atom_ptr->type, res_no, *(atoms_ptr->atomname[i_atm]), 
		*(atoms_ptr->atomtype[i_atm]), *(atoms_ptr->resinfo[res_no].name));
	}

	if (self->contents->bX) {
	    printf("x: %f %f %f \t ", self->contents->x[i_atm][0], self->contents->x[i_atm][1],
		   self->contents->x[i_atm][2]);
	}
	if (self->contents->bV) {
	    printf("v: %f %f %f \t ", self->contents->v[i_atm][0], self->contents->v[i_atm][1],
		   self->contents->v[i_atm][2]);
	}
	if (self->contents->bF) {
	    printf("f: %f %f %f \t ", self->contents->f[i_atm][0], self->contents->f[i_atm][1],
		   self->contents->f[i_atm][2]);
	}
	printf("\n");
    }

    if (self->contents->bBox) {
	printf("printing box. \n");
	for (i = 0; i < DIM; i++) {
	    printf("i: %d \t", i);
	    for (j = 0; j < DIM; j++) {
		printf(" %f ", self->contents->box[i][j]);
	    }
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
init_tW_gmx_trxframe(): initializes a tW_gmx_trxframe variable
*****************************************************************************************/
tW_gmx_trxframe* init_tW_gmx_trxframe()
{
	tW_gmx_trxframe *frame;

	frame = emalloc(sizeof(tW_gmx_trxframe));

	frame->contents = emalloc(sizeof(t_trxframe));
	//frame->oenv = emalloc(sizeof(output_env_t));
	snew(frame->oenv,1);
        output_env_init_default(frame->oenv);

	frame->set_natoms = set_natoms;
	frame->set_atom_labels = set_atom_labels;
	frame->print_frame_dump = print_frame_dump;

	frame->get_pos_of = get_pos_of;
	frame->get_vel_of = get_vel_of;
	frame->get_box = get_box;

	return frame;
}



/*****************************************************************************************
read_first_trxframe(): Wrapper function for the trx reading functionality in the gromacs library
  for the first frame of such a file
*****************************************************************************************/
int read_first_trxframe(tW_gmx_trxframe *frame, const char *trx_fnm)
{

    int in_flags = TRX_READ_X | TRX_READ_V | TRX_READ_F;

    frame->contents = emalloc(sizeof(t_trxframe));

    return read_first_frame(frame->oenv, &frame->status, trx_fnm, frame->contents, in_flags);
}


/*****************************************************************************************
read_next_trxframe(): Wrapper function for the trx reading functionality in the gromacs library
  for the first frame of such a file
*****************************************************************************************/
bool read_next_trxframe(tW_gmx_trxframe *frame)
{
    return read_next_frame(frame->oenv, frame->status, frame->contents);

}


/*****************************************************************************************
open_new_trxfile(): Opens a new trx file for writing
*****************************************************************************************/
void open_new_trxfile(tW_gmx_trxframe *frame, const char *trx_fnm)
{
    //frame->status = emalloc(sizeof(t_trxstatus));

    frame->status = open_trx(trx_fnm, "w");

}


/*****************************************************************************************
write_trxframe_to_file(): Writes a trx frame to disk
*****************************************************************************************/
void write_trxframe_to_file(tW_gmx_trxframe *frame)
{
    write_trxframe(frame->status, frame->contents, NULL);
}


/*****************************************************************************************
copy_trxframe_info(): copies the resolution-independent information from source to dest
*****************************************************************************************/
void copy_trxframe_info(tW_gmx_trxframe *source, tW_gmx_trxframe *dest)
{
    dest->contents->bTitle = FALSE;

    dest->contents->flags = source->contents->flags;

    dest->contents->bDouble = source->contents->bDouble;

    dest->contents->bTime = source->contents->bTime;
    if (dest->contents->bTime)
    {
        dest->contents->time = source->contents->time;
        dest->contents->t0 = source->contents->t0;
        dest->contents->tpf = source->contents->tpf;
        dest->contents->tppf = source->contents->tppf;
    }

    dest->contents->bStep = source->contents->bStep;
    if (dest->contents->bStep)
    {
	dest->contents->step = source->contents->step;
    }



    dest->contents->bLambda = source->contents->bLambda;
    if (source->contents->bLambda) {
	dest->contents->lambda = source->contents->lambda;
    }

    dest->contents->bPrec = source->contents->bPrec;
    if (source->contents->bPrec) {
	dest->contents->prec = source->contents->prec;
    }

    dest->contents->bBox = source->contents->bBox;
    if (source->contents->bBox) {
	copy_matrix(source->contents->box, dest->contents->box);
    }


    dest->contents->bX = source->contents->bX;
    dest->contents->bV = source->contents->bV;
    dest->contents->bF = source->contents->bF;

}


/*****************************************************************************************
map_trxframe(): maps source to dest according to the correspondence defined in map
*****************************************************************************************/
void map_trxframe(tW_gmx_trxframe *source, tW_gmx_trxframe *dest, tW_site_map *map)
{
    int I_site, i_atm, n_I, atm_no;

    dvec x, v, f, sc_x, sc_v;

    double c_Ii;

    for (I_site = 0; I_site < dest->contents->natoms; I_site++) {
	n_I = map[I_site].n_atms;

	if (dest->contents->bX) {
	    clear_dvec(dest->contents->x[I_site]);
	}
	if (dest->contents->bV) {
	    clear_dvec(dest->contents->v[I_site]);
	}
	if (dest->contents->bF) {
	    clear_dvec(dest->contents->f[I_site]);
	}



	for (i_atm = 0; i_atm < n_I; i_atm++) {
	    atm_no = map[I_site].i_atm[i_atm];

	    c_Ii = map[I_site].c_Ii[i_atm];

	    if (dest->contents->bX) {
		copy_vector(source->contents->x[atm_no], x);
		scal_times_vect(c_Ii, x, sc_x);
		vect_inc(dest->contents->x[I_site], sc_x);
	    }
	    if (dest->contents->bV) {
		copy_vector(source->contents->v[atm_no], v);
		scal_times_vect(c_Ii, v, sc_v);
		vect_inc(dest->contents->v[I_site], sc_v);
	    }
	    // NB at the moment we are assuming that no atom contributes 
	    // to more than one cg sites, s.t. c_Ii = d_Ii and the 
	    // \cal F_I = net atomistic force
	    if (dest->contents->bF) {
		copy_vector(source->contents->f[atm_no], f);
		vect_inc(dest->contents->f[I_site], f);
	    }
	}
    }
}

/*****************************************************************************************
get_filename(): returns the filename corresponding to the argument opt
*****************************************************************************************/
const char *get_filename(tW_gmx_cgmap_input *self, const char *opt)
{
    return opt2fn(opt, self->N_files, self->fnm);
}

/*****************************************************************************************
arg_is_set(): returns TRUE if the argument corresponding to opt has been set, and FALSE otherwise
*****************************************************************************************/
bool arg_is_set(tW_gmx_cgmap_input *self, const char *opt)
{
    return opt2bSet(opt, self->N_files, self->fnm);
}


/*****************************************************************************************
init_tW_gmx_cgmap_input(): Initializes a new instance of tW_gmx_cgmap_input
*****************************************************************************************/
tW_gmx_cgmap_input* init_tW_gmx_cgmap_input(int argc, char *argv[])
{
    tW_gmx_cgmap_input *input;
    const char *desc[] = { "select which atom you want to examine with the -n argument." };

    input = ecalloc(1,sizeof(tW_gmx_cgmap_input));
    input->oenv = emalloc(sizeof(output_env_t));

    input->get_filename = get_filename;
    input->arg_is_set = arg_is_set;

    t_filenm fnm[] = {
	{efTOP, "-p", "map", ffREAD},	/* input mapping */
	{efTRX, "-f", "aa", ffREAD},	/* input atomistic trajectory */
	{efTRX, "-o", "cg", ffWRITE},	/* output cg trajectory */
	{efTPR, "-s", "cg", ffOPTRD},   /* input CG topology (optional) */
	{efGRO, "-c", "cg", ffOPTWR}    /* output gro file (optional) */
    };

    #define NFILES asize(fnm)

    int i;

    for (i=0; i<NFILES; i++)
    {
	input->fnm[i].ftp = fnm[i].ftp;
	input->fnm[i].opt = strdup(fnm[i].opt);
	input->fnm[i].fn = strdup(fnm[i].fn);
	input->fnm[i].flag = fnm[i].flag;
     
    }

    input->N_files = NFILES;

    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW, input->N_files, input->fnm, 0, NULL, asize(desc), desc, 0, NULL, input->oenv);

    return input;
}


/*****************************************************************************************
init_tW_gmx_cgmap_input(): Initializes a new instance of tW_gmx_cgff_input
*****************************************************************************************/
//tW_gmx_cgff_input* init_tW_gmx_cgmap_input()
//{
//    tW_gmx_cgff_input *input;
//
//    input = emalloc(sizeof(tW_gmx_cgff_input));
//    input->oenv = emalloc(sizeof(output_env_t));
//
//    return input;
//}








/*****************************************************************************************
open_trr_file():
*****************************************************************************************/
//int open_trr_file( char *trr_fnm, int *status, tW_gmx_info *info, t_trxframe *fr ) // 4.0.7
//int open_trr_file(char *trr_fnm, output_env_t oenv, t_trxstatus ** status, tW_gmx_info * info, t_trxframe * fr)	//4.5.3
//{
//    int n_atms;
//    int flags = TRX_READ_X | TRX_READ_V | TRX_READ_F;
//
///*
// * header in statutil.h
// * source in src/gmxlib/trxio.c
// * 4.0.7: extern int read_first_frame(int *status,char *fn,t_trxframe *fr,int flags);
// * 4.5.3: int read_first_frame(const output_env_t oenv,t_trxstatus **status,
// *                   const char *fn,t_trxframe *fr,int flags)
// */
//    //n_atms = read_first_frame( status, trr_fnm, fr, flags ); //4.0.7
//    n_atms = read_first_frame(oenv, status, trr_fnm, fr, flags);	//4.5.3
//
//    update_info_trr(info, *fr);
//
//    return n_atms;
//}



/*****************************************************************************************
copy_trr_2_CGstruct():
*****************************************************************************************/
bool copy_trr_2_CGstruct(tW_gmx_trxframe *fr, tW_CG_site CG_struct[])
{
    int i;
    int n_atms = fr->contents->natoms;

    for (i = 0; i < n_atms; i++) {
	if (fr->contents->bX) {
	    copy_dvec(fr->contents->x[i], CG_struct[i].r);
	}
	if (fr->contents->bF) {
	    copy_dvec(fr->contents->f[i], CG_struct[i].f);
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

    for (i = 0; i < n_atms; i++) {
	if (fr->contents->bF) {
	    copy_dvec(fr->contents->f[i], CG_struct[i].ref_f);
	}
    }

    return fr->contents->bF;
}


/*****************************************************************************************
read_trr_2_CGstruct():
*****************************************************************************************/
//bool read_trr_2_CGstruct( int status, tW_gmx_info * info, t_trxframe * fr, tW_CG_site CG_struct[] ) //4.0.7
bool read_trr_2_CGstruct(tW_gmx_info * info, tW_gmx_trxframe *fr, tW_CG_site CG_struct[])	//4.5.3
{
    bool b_trr;

/*
 * header in statutil.h
 * source in src/gmxlib/trxio.c
 * 4.0.7: extern bool read_next_frame(int status,t_trxframe *fr);
 * 4.5.3: gmx_bool read_next_frame(const output_env_t oenv,t_trxstatus *status,t_trxframe *fr);
 *
 */

    //b_trr = read_next_frame( status, fr ); // 4.0.7
    //b_trr = read_next_frame(oenv, status, fr);	// 4.5.3
    b_trr = read_next_trxframe(fr);

    if (b_trr) {
	copy_trr_2_CGstruct(fr, CG_struct);
	update_info_trr(info, fr);
    }

    return b_trr;
}

/*****************************************************************************************
read_trr_2_CGstruct_ref(): JFR - 07.16.12: for reading in reference forces from a trr file
*****************************************************************************************/
bool read_trr_2_CGstruct_ref(tW_gmx_info * info, tW_gmx_trxframe *fr, tW_CG_site CG_struct[])	//4.5.3
{
    bool b_trr;

/*
 * header in statutil.h
 * source in src/gmxlib/trxio.c
 * 4.0.7: extern bool read_next_frame(int status,t_trxframe *fr);
 * 4.5.3: gmx_bool read_next_frame(const output_env_t oenv,t_trxstatus *status,t_trxframe *fr);
 *
 */

    //b_trr = read_next_frame( status, fr ); // 4.0.7
//  b_trr = read_next_frame(oenv, status, fr);	// 4.5.3
    b_trr = read_next_trxframe(fr);

    if (b_trr) {
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

    if (fr->contents->bF) {
	info->b_Forces_1 = TRUE;
    } else {
	info->b_Forces_N = FALSE;
    }

    if (fr->contents->bBox) {
	for (i = 0; i < DIM; i++) {
	    for (j = 0; j < DIM; j++) {
		info->box[i][j] = fr->contents->box[i][j];
	    }
	}
    }

    return 0;
}


















