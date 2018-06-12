/**
@file cgff_fn.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Helper functions for the cgff calculation
*/

//c library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cgff_types.h"
//#include "gmx-interface.h"
#include "gromacs_topology.h"
#include "cgff_fn.h"
#include "safe_mem.h"
#include "io_read.h"

/*****************************************************************************************
init_gmx_info(): Initializes flags to conditions so that b_Forces_1 can be set to TRUE
as soon as a frame is found with forces while b_Forces_N can be set to FALSE as soon as
a frame is found that does not have forces. This way we can determine wheter or not
b must be calculated from dg/dr by checking ( b_Forces1 && b_Forces_N ).
*****************************************************************************************/
void init_gmx_info(tW_gmx_info *info)
{
    info->b_Gromacs = FALSE;
    info->b_Forces = FALSE;
    info->b_PBC = FALSE;
    info->b_Forces_1 = FALSE;
    info->b_Forces_N = TRUE;
}

/*****************************************************************************************
reset_gmx_info(): Reset the flags for a new topology.
*****************************************************************************************/
void reset_gmx_info(tW_gmx_info *info)
{
    info->b_Gromacs = FALSE;
    info->b_Forces = FALSE;
    info->b_PBC = FALSE;
}




/*****************************************************************************************
setup_CG_struct(): Stores the information in top into sys and CG_struct. 
*****************************************************************************************/
int setup_CG_struct(tW_system *sys, tW_gmx_topology *top,
		    tW_CG_site *CG_struct,
		    tW_Bonded_Inter *Bonded_Inter_Types)
{
    /* Get information for all CG sites from top. */
    get_site_info(CG_struct, sys, top);

    /* Get bond information. */
    get_bonds(sys, top, top->get_natoms(top), CG_struct,
	      sys->N_Bond_Int_Types, Bonded_Inter_Types);

    /* Return no. of particles in top file. */
    return top->get_natoms(top);

}

/*****************************************************************************************
get_bonds(): Stores information from top into sys and CG_struct.
*****************************************************************************************/
void get_bonds(tW_system  *sys, 
               tW_gmx_topology *top, 
               int n_sites,
	       tW_CG_site *CG_struct, 
               int N_Bond_Int_Types,
	       tW_Bonded_Inter *Bond_Type_List)
{
    tW_bond_type_indices bond_type_indices;

    /* Allocate memory for list of indices for all bond interaction types. */
    allocate_memory_bond_type_indices(top, &bond_type_indices);

    /* Clear previous bond lists, if needed. */
    clear_bond_lists(N_Bond_Int_Types, Bond_Type_List);

    /* Clear information about bonds for each site from previous topology. */
    clear_site_bond_info(n_sites, CG_struct);

    /* Check that all bond interactions in top file check out with par.txt. */
    check_bond_interactions(top, CG_struct, N_Bond_Int_Types, Bond_Type_List, &bond_type_indices);

    /* Allocate memory for Bond_Type_ptr->Inter_List[]. */
    allocate_memory_for_Bond_Inter_List(N_Bond_Int_Types, Bond_Type_List);

    /* Allocate memory for bond info for each site. */
    allocate_memory_for_sites(n_sites, CG_struct);

    /* Stores bond information into CG_struct. */
    store_bond_interactions(top, CG_struct, N_Bond_Int_Types, Bond_Type_List, bond_type_indices, n_sites);
}


/*****************************************************************************************
allocate_memory_bond_type_indices(): Allocate memory to hold indices to 
tW_system.Bonded_Inter_Types[i] for each type of interaction found in the GROMACS 
topology.
*****************************************************************************************/
void allocate_memory_bond_type_indices(tW_gmx_topology *top,
				       tW_bond_type_indices *bond_type_indices)
{
    bond_type_indices->bond_type_i =  ecalloc(top->get_nbonds(top), sizeof(int));
    bond_type_indices->angle_type_i = (int *) ecalloc(top->get_nangles(top), sizeof(int));
    bond_type_indices->dihedral_type_i = (int *) ecalloc(top->get_ndihs(top), sizeof(int));
    bond_type_indices->pair_type_i = (int *) ecalloc(top->get_npairs(top), sizeof(int));
}

/*****************************************************************************************
clear_bond_lists(): Deletes previous bond interaction lists and initializes N_instances
to 0 for each interaction.
*****************************************************************************************/
void clear_bond_lists(int N_Bond_Int_Types,
		      tW_Bonded_Inter * Bond_Type_List)
{
    int i, j;
    tW_Bonded_Inter *Bond_Type_ptr;


    /* Loop over all bond interactions. */
    for (i = 0; i < N_Bond_Int_Types; i++) {
	/* Pointer to the ith bond interaction. */
	Bond_Type_ptr = &(Bond_Type_List[i]);

	/* N_instances is initialized to 0, so this loop is skipped the first time. */
	for (j = 0; j < Bond_Type_ptr->N_instances; j++) {
	    free(Bond_Type_ptr->Inter_List[j]);
	}

	/* Pointed to NULL initially. */
	if (Bond_Type_ptr->Inter_List != NULL) {
	    free(Bond_Type_List[i].Inter_List);
	}

	/* Initialize for new topology. */
	Bond_Type_List[i].N_instances = 0;
    }
}


/*****************************************************************************************
clear_site_bond_info(): Free lists from previous topology, and initialize nr_bonds and
nr_bond_coeffs.
*****************************************************************************************/
void clear_site_bond_info(int n_sites, tW_CG_site CG_struct[])
{
    int i;
    tW_CG_site *site_ptr;


    /* Loop over n_sites CG sites. */
    for (i = 0; i < n_sites; i++) {
	site_ptr = &(CG_struct[i]);

	/* Free memory from previous topology. */
	if (site_ptr->bond_site != NULL) {
	    free(site_ptr->bond_site);
	}

	/* Free memory from previous topology. */
	if (site_ptr->bond_type != NULL) {
	    free(site_ptr->bond_type);
	}

	/* Initialize for current topology. */
	CG_struct[i].nr_bonds = 0;
	CG_struct[i].nr_bond_coeffs = 0;
    }
}



/*****************************************************************************************
check_bond_interactions(): Store indices to Bond_Type_List[] for each bond interaction 
found in the GROMACS topology into the arrays in bond_type_indices. Keep track of the 
number of occurances for each bond interaction type: Bond_Type_List.N_instances. Keep 
track of the total number of bond interactions for each CG site: CG_struct[].nr_bonds
*****************************************************************************************/
void check_bond_interactions(tW_gmx_topology *top, 
                             tW_CG_site * CG_struct,
			     int N_Bond_Int_Types,
			     tW_Bonded_Inter * Bond_Type_List,
			     tW_bond_type_indices * bond_type_indices)
{
    check_bonds(top, CG_struct,
		bond_type_indices->bond_type_i, N_Bond_Int_Types,
		Bond_Type_List);
    check_angles(top, CG_struct,
		 bond_type_indices->angle_type_i, N_Bond_Int_Types,
		 Bond_Type_List);
    check_dihedrals(top, CG_struct,
		    bond_type_indices->dihedral_type_i, N_Bond_Int_Types,
		    Bond_Type_List);
    check_pairs(top, CG_struct, bond_type_indices->pair_type_i,
		N_Bond_Int_Types, Bond_Type_List);
}


/*****************************************************************************************
allocate_memory_for_Bond_Inter_List():
*****************************************************************************************/
void allocate_memory_for_Bond_Inter_List(int N_Bond_Int_Types,
					 tW_Bonded_Inter *Bond_Type_List)
{
    int i, j;
    int N_instances;
    tW_Bonded_Inter *Bond_Type_ptr;

    /* Loop over all bond interaction types. */
    for (i = 0; i < N_Bond_Int_Types; i++) {
	/* Local variable for the ith bond interaction type. */
	Bond_Type_ptr = &(Bond_Type_List[i]);

	/* Number of occurrences of ith bond interaction for current topology. */
	N_instances = Bond_Type_ptr->N_instances;

	/* Allocate memory for interaction list. */
	Bond_Type_List[i].Inter_List =
	    (int **) emalloc(N_instances * sizeof(int *));
	for (j = 0; j < N_instances; j++) {
	    Bond_Type_ptr->Inter_List[j] =
		(int *) emalloc(Bond_Type_ptr->N_Int_Sites * sizeof(int));
	}
    }
}


/*****************************************************************************************
allocate_memory_for_sites():
*****************************************************************************************/
void allocate_memory_for_sites(int n_sites, tW_CG_site *CG_struct)
{
    int i;
    tW_CG_site *site_ptr;

    /* Loop over all sites. */
    for (i = 0; i < n_sites; i++) {
	/* Local variable for ith site. */
	site_ptr = &(CG_struct[i]);

	/* Allocate the memory. */
	site_ptr->bond_type = (int *) emalloc(site_ptr->nr_bonds * sizeof(int));
	site_ptr->bond_site = (int **) emalloc(site_ptr->nr_bonds * sizeof(int *));
    }
}




/*****************************************************************************************
print_Bond_Types():
*****************************************************************************************/
void print_Bond_Types(FILE * fp, tW_system sys)
{
    int i, j, k;

    fprintf(fp, "\n\nIn print_Bond_Types().\n\n");

    fprintf(fp, "  There are %d bond interaction types.\n\n",
	    sys.N_Bond_Int_Types);

    for (i = 0; i < sys.N_Bond_Int_Types; i++) {
	fprintf(fp, "    %3d. %s    Sites:  ", i + 1,
		sys.Bonded_Inter_Types[i].name);

	for (j = 0; j < sys.Bonded_Inter_Types[i].N_Int_Sites; j++) {
	    fprintf(fp, "%s", sys.Bonded_Inter_Types[i].Site_Types[j]);
	    if (j != (sys.Bonded_Inter_Types[i].N_Int_Sites - 1)) {
		fprintf(fp, "-");
	    }
	}
	fprintf(fp, "\n");

	fprintf(fp, "      There are %d instances of this interaction.\n",
		sys.Bonded_Inter_Types[i].N_instances);
	for (j = 0; j < sys.Bonded_Inter_Types[i].N_instances; j++) {
	    fprintf(fp, "        %3d.  Sites:  ", j + 1);
	    for (k = 0; k < sys.Bonded_Inter_Types[i].N_Int_Sites; k++) {
		fprintf(fp, "%d",
			sys.Bonded_Inter_Types[i].Inter_List[j][k]);
		if (k != (sys.Bonded_Inter_Types[i].N_Int_Sites - 1)) {
		    fprintf(fp, "-");
		}
	    }
	    fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
    }

    fprintf(fp, "End print_Bond_Types().\n\n\n");
}

/*****************************************************************************************
store_bond_interactions():
*****************************************************************************************/
void store_bond_interactions(tW_gmx_topology *top, 
                             tW_CG_site * CG_struct,
			     int N_Bond_Int_Types,
			     tW_Bonded_Inter * Bond_Type_List,
			     tW_bond_type_indices bond_type_indices,
			     int n_sites)
{
    int i;
    int site_bond_ctr[n_sites];
    int bond_ctr[N_Bond_Int_Types];

    /* 1. bond_ctr[i] is incremented each time the interaction type Bond_Type_List[i] is found       */
    /*    in the GROMACS topology.                                                                   */
    /* 2. site_bond_ctr[i] is incremented each time CG site i is involved in a bond interaction      */
    /*    found in the GROMACS topology.                                                             */
    /* 3. Store CG site indices into Bond_Type_List[].Inter_List.                                    */
    /* 4. Store into CG_struct[].bond_type the indices to Bond_Type_List[] for each bond interaction */
    /*    that each CG site is involved.                                                             */
    /* 5. Count the number of coefficients to determine for bond interactions for each CG site.      */

    for (i = 0; i < N_Bond_Int_Types; i++) {
	bond_ctr[i] = 0;
    }
    for (i = 0; i < n_sites; i++) {
	site_bond_ctr[i] = 0;
    }

    store_bond_lists(top, CG_struct,
		     bond_type_indices.bond_type_i, N_Bond_Int_Types,
		     Bond_Type_List, site_bond_ctr, bond_ctr);

    store_angle_lists(top, CG_struct,
		      bond_type_indices.angle_type_i, N_Bond_Int_Types,
		      Bond_Type_List, site_bond_ctr, bond_ctr);

    store_dihedral_lists(top, CG_struct,
			 bond_type_indices.dihedral_type_i,
			 N_Bond_Int_Types, Bond_Type_List, site_bond_ctr,
			 bond_ctr);

    store_pair_lists(top, CG_struct,
		     bond_type_indices.pair_type_i, N_Bond_Int_Types,
		     Bond_Type_List, site_bond_ctr, bond_ctr);

}



/*****************************************************************************************
check_bonds(): Make sure bonds in topology match those in input file.
*****************************************************************************************/
void check_bonds(tW_gmx_topology *top, tW_CG_site * CG_struct, int *i_type,
		 int N_Bond_Int_Types, tW_Bonded_Inter * Bond_Type_List)
{
    int i = 0;
    int i_bond = 0;
    int i1, i2;
    tW_word site_names[2];
    int *bond_list;
    int n_bonds;

    bond_list = top->get_bond_list(top);
    n_bonds = top->get_nbonds(top);

    /* Loop over all bonds int GROMACS topology. */
    while (i < n_bonds) {

	/* Get site numbers from the GROMACS topology. */
	i1 = bond_list[i + 1];
	i2 = bond_list[i + 2];

	/* Copy site names into local variable. */
	strcpy(site_names[0], CG_struct[i1].name);
	strcpy(site_names[1], CG_struct[i2].name);

	/* Store indices to Bond_Type_List[] in order of interaction appreance in GROMACS topology. */
	i_type[i_bond] =
	    determine_bond_type(N_Bond_Int_Types, Bond_Type_List,
				site_names);

	/* Is there a match between site_names and Bond_Type_List? */
	if (i_type[i_bond] == -1) {
	    printf
		("\nERROR: BondStretch type in top file does not match with par.txt.\n");
	    printf("  Site_1: %s   Site_2: %s\n", site_names[0],
		   site_names[1]);
	    exit(EXIT_FAILURE);
	} else {
	    Bond_Type_List[i_type[i_bond]].N_instances++;
	}

	/* Count no. of bonds each site is involved in. */
	CG_struct[i1].nr_bonds++;
	CG_struct[i2].nr_bonds++;

	/* Increment counters. */
	i += 3;
	i_bond += 1;

    }				/* End list over bonds */

    free(bond_list);
}


/*****************************************************************************************
check_angles():
*****************************************************************************************/
void check_angles(tW_gmx_topology *top, tW_CG_site CG_struct[], int *i_type,
		  int N_Bond_Int_Types, tW_Bonded_Inter * Bond_Type_List)
{
    int i = 0;
    int i_bond = 0;
    int i1, i2, i3;
    tW_word site_names[3];
    int *angle_list;
    int n_angles;

    angle_list = top->get_angle_list(top);
    n_angles = top->get_nangles(top);

    /* Loop over all angles in GROMACS topology. */
    while (i < n_angles) {

	/* Get site numbers for bonds from the GROMACS bond_list.iatoms. */
	i1 = angle_list[i + 1];
	i2 = angle_list[i + 2];
	i3 = angle_list[i + 3];

	/* Copy site names into local variables. Names from get_site_info(). */
	strcpy(site_names[0], CG_struct[i1].name);
	strcpy(site_names[1], CG_struct[i2].name);
	strcpy(site_names[2], CG_struct[i3].name);

	/* This creates an array i_type of the bond types in order of top. */
	i_type[i_bond] =
	    determine_angle_type(N_Bond_Int_Types, Bond_Type_List,
				 site_names);

	/* Is there a match between site_names[] and Bond_Type_List[]? */
	if (i_type[i_bond] == -1) {
	    printf
		("\nERROR: Angle type in top file does not match with par.txt.\n");
	    printf("  Site_1: %s   Site_2: %s   Site_3: %s\n",
		   site_names[0], site_names[1], site_names[2]);
	    exit(EXIT_FAILURE);
	} else {
	    Bond_Type_List[i_type[i_bond]].N_instances++;
	}

	/* Count no. of angles each site is involved in. */
	CG_struct[i1].nr_bonds++;
	CG_struct[i2].nr_bonds++;
	CG_struct[i3].nr_bonds++;

	/* Increment counters. */
	i += 4;
	i_bond += 1;
    }

    free(angle_list);
}


/*****************************************************************************************
check_dihedrals(): 
*****************************************************************************************/
void check_dihedrals(tW_gmx_topology *top, tW_CG_site * CG_struct,
		     int *i_type, int N_Bond_Int_Types,
		     tW_Bonded_Inter * Bond_Type_List)
{
    int i = 0;
    int i_bond = 0;
    int i1, i2, i3, i4;
    tW_word site_names[4];
    int *dihedral_list;
    int n_dihs;

    dihedral_list = top->get_dih_list(top);
    n_dihs = top->get_ndihs(top);

    /* Loop over all dihedrals in GROMACS topology. */
    while (i < n_dihs) {

	/* Get site numbers for bonds from the GROMACS bond_list.iatoms. */
	i1 = dihedral_list[i + 1];
	i2 = dihedral_list[i + 2];
	i3 = dihedral_list[i + 3];
	i4 = dihedral_list[i + 4];

	/* Copy site names into local variables. Names from get_site_info(). */
	strcpy(site_names[0], CG_struct[i1].name);
	strcpy(site_names[1], CG_struct[i2].name);
	strcpy(site_names[2], CG_struct[i3].name);
	strcpy(site_names[3], CG_struct[i4].name);

	/* This creates an array i_type of the bond types in order of top. */
	i_type[i_bond] =
	    determine_dihedral_type(N_Bond_Int_Types, Bond_Type_List,
				    site_names);

	/* Is there a match between site_names[] and Bond_Type_List[]? */
	if (i_type[i_bond] == -1) {
	    printf
		("\nERROR: Dihedral type in top file does not match with par file.\n");
	    printf("  Site_1: %s   Site_2: %s   Site_3: %s   Site_4: %s\n",
		   site_names[0], site_names[1], site_names[2],
		   site_names[3]);
	    exit(0);
	} else {
	    Bond_Type_List[i_type[i_bond]].N_instances++;
	}

	/* Count no. of dihedral angles each site is involved in. */
	CG_struct[i1].nr_bonds++;
	CG_struct[i2].nr_bonds++;
	CG_struct[i3].nr_bonds++;
	CG_struct[i4].nr_bonds++;

	/* Increment counters. */
	i += 5;
	i_bond += 1;
    }

    free(dihedral_list);
}


/*****************************************************************************************
check_pairs():
*****************************************************************************************/
void check_pairs(tW_gmx_topology *top, tW_CG_site * CG_struct, int *i_type,
		 int N_Bond_Int_Types, tW_Bonded_Inter * Bond_Type_List)
{
    int i = 0;
    int i_bond = 0;
    int i1, i2;
    tW_word site_names[2];
    int *pair_list;
    int n_pairs;

    pair_list = top->get_pair_list(top);
    n_pairs = top->get_npairs(top);

    /* Loop over all pairs in GROMACS topology. */
    while (i < n_pairs) {
	/* Get site numbers for bonds from the GROMACS bond_list.iatoms. */
	i1 = pair_list[i + 1];
	i2 = pair_list[i + 2];

	/* Copy site names into local variables. Names from get_site_info(). */
	strcpy(site_names[0], CG_struct[i1].name);
	strcpy(site_names[1], CG_struct[i2].name);

	/* This creates an array i_type of the bond types in order of top. */
	i_type[i_bond] =
	    determine_pair_type(N_Bond_Int_Types, Bond_Type_List,
				site_names);

	/* Is there a match between site_names[] and Bond_Type_List[]? */
	if (i_type[i_bond] == -1) {
	    printf
		("\nERROR: Inter. Pair type in top file does not match with par file.\n");
	    printf("  Site_1: %s   Site_2: %s\n", site_names[0],
		   site_names[1]);
	    exit(0);
	} else {
	    Bond_Type_List[i_type[i_bond]].N_instances++;
	}

	/* Count no. of bonds each site is involved in */
	CG_struct[i1].nr_bonds++;
	CG_struct[i2].nr_bonds++;

	/* Increment counters. */
	i += 3;
	i_bond += 1;
    }

    free(pair_list);
}








/*****************************************************************************************
store_bond_lists():
*****************************************************************************************/
void store_bond_lists(tW_gmx_topology *top, 
		      tW_CG_site * CG_struct,
		      int *i_type, int N_Bond_Int_Types,
		      tW_Bonded_Inter * Bond_Type_List, 
                      int *site_bond_ctr,
		      int *bond_ctr)
{
    int i = 0;
    int i_bond = 0;
    int i1, i2;
    int i_bond_type;
    int bond_index;
    tW_Bonded_Inter *Bond_Type_ptr;
    int *bond_list;
    int n_bonds;

    bond_list = top->get_bond_list(top);
    n_bonds = top->get_nbonds(top);


    /* Loop over all BondStretch interactions in current topology. */
    while (i < n_bonds) {
	/* Index to Bond_Type_List for the ith interaction. */
	i_bond_type = i_type[i_bond];

	/* Pointer to the i_bond_type interaction. */
	Bond_Type_ptr = &(Bond_Type_List[i_bond_type]);

	/* Index for each CG site in the ith interaction. */
	i1 = bond_list[i + 1];
	i2 = bond_list[i + 2];

	/* Order the ith interaction occurs in Bond_Type_List[i_bond_type].Inter_List. */
	bond_index = bond_ctr[i_bond_type];

	/* Store CG site indices for the ith interaction. */
	Bond_Type_ptr->Inter_List[bond_index][0] = i1;
	Bond_Type_ptr->Inter_List[bond_index][1] = i2;

	/* Store i_bond_type in CG_struct[] for each site. */
	CG_struct[i1].bond_type[site_bond_ctr[i1]] = i_bond_type;
	CG_struct[i2].bond_type[site_bond_ctr[i2]] = i_bond_type;

	/*  Link lists for the ith interaction. */
	CG_struct[i1].bond_site[site_bond_ctr[i1]] =
	    Bond_Type_ptr->Inter_List[bond_index];
	CG_struct[i2].bond_site[site_bond_ctr[i2]] =
	    Bond_Type_ptr->Inter_List[bond_index];

	/* Keep track of the number of coeffs. for each site. */
	CG_struct[i1].nr_bond_coeffs += Bond_Type_ptr->N_coeff;
	CG_struct[i2].nr_bond_coeffs += Bond_Type_ptr->N_coeff;

	/* Update counters. */
	i += 3;
	i_bond += 1;
	site_bond_ctr[i1] += 1;
	site_bond_ctr[i2] += 1;
	bond_ctr[i_bond_type] += 1;
    }

    free(bond_list);
}


/*****************************************************************************************
store_angle_lists():
*****************************************************************************************/
void store_angle_lists(tW_gmx_topology *top, 
                       tW_CG_site * CG_struct,
		       int *i_type, int N_Bond_Int_Types,
		       tW_Bonded_Inter * Bond_Type_List,
		       int *site_bond_ctr, int *bond_ctr)
{
    int i = 0;
    int i_bond = 0;
    int i1, i2, i3;
    int i_bond_type;
    int bond_index;
    tW_Bonded_Inter *Bond_Type_ptr;
    int *angle_list;
    int n_angles;

    angle_list = top->get_angle_list(top);
    n_angles = top->get_nangles(top);

    /* Loop over all Angle interactions in current toplogy. */
    while (i < n_angles) {
	/* Index to Bond_Type_List for the ith interaction. */
	i_bond_type = i_type[i_bond];

	/* Pointer to i_bond_type interaction. */
	Bond_Type_ptr = &(Bond_Type_List[i_bond_type]);

	/* Index for each CG site in the ith interaction. */
	i1 = angle_list[i + 1];
	i2 = angle_list[i + 2];
	i3 = angle_list[i + 3];

	/* Order the ith interaction occurs in Bond_Type_List[i_bond_type].Inter_List. */
	bond_index = bond_ctr[i_bond_type];

	/* Store CG site indices for the ith interaction. */
	Bond_Type_ptr->Inter_List[bond_index][0] = i1;
	Bond_Type_ptr->Inter_List[bond_index][1] = i2;
	Bond_Type_ptr->Inter_List[bond_index][2] = i3;

	/* Store i_bond_type in CG_struct[] for each site. */
	CG_struct[i1].bond_type[site_bond_ctr[i1]] = i_bond_type;
	CG_struct[i2].bond_type[site_bond_ctr[i2]] = i_bond_type;
	CG_struct[i3].bond_type[site_bond_ctr[i3]] = i_bond_type;

	/* Link lists for the ith interaction. */
	CG_struct[i1].bond_site[site_bond_ctr[i1]] =
	    Bond_Type_ptr->Inter_List[bond_index];
	CG_struct[i2].bond_site[site_bond_ctr[i2]] =
	    Bond_Type_ptr->Inter_List[bond_index];
	CG_struct[i3].bond_site[site_bond_ctr[i3]] =
	    Bond_Type_ptr->Inter_List[bond_index];

	/* Keep track of the number of coeffs. for each site. */
	CG_struct[i1].nr_bond_coeffs += Bond_Type_ptr->N_coeff;
	CG_struct[i2].nr_bond_coeffs += Bond_Type_ptr->N_coeff;
	CG_struct[i3].nr_bond_coeffs += Bond_Type_ptr->N_coeff;

	/* Update counters. */
	i += 4;
	i_bond += 1;
	site_bond_ctr[i1] += 1;
	site_bond_ctr[i2] += 1;
	site_bond_ctr[i3] += 1;
	bond_ctr[i_bond_type] += 1;
    }

    free(angle_list);
}


/*****************************************************************************************
store_dihedral_lists():
*****************************************************************************************/
void store_dihedral_lists(tW_gmx_topology *top, 
		          tW_CG_site * CG_struct,
			  int *i_type, int N_Bond_Int_Types,
			  tW_Bonded_Inter * Bond_Type_List,
			  int *site_bond_ctr, int *bond_ctr)
{
    int i = 0;
    int i_bond = 0;
    int i1, i2, i3, i4;
    int i_bond_type;
    int bond_index;
    tW_Bonded_Inter *Bond_Type_ptr;
    int *dihedral_list;
    int n_dihs;

    dihedral_list = top->get_dih_list(top); 
    n_dihs = top->get_ndihs(top);

    /* Loop over all Dihedral interactions in current topology. */
    while (i < n_dihs) {
	/* Index to Bond_Type_List for the ith interaction. */
	i_bond_type = i_type[i_bond];

	/* Pointer to i_bond_type interaction. */
	Bond_Type_ptr = &(Bond_Type_List[i_bond_type]);

	/* Index for each CG site in the ith interaction. */
	i1 = dihedral_list[i + 1];
	i2 = dihedral_list[i + 2];
	i3 = dihedral_list[i + 3];
	i4 = dihedral_list[i + 4];

	/* Order the ith interaction occurs in Bond_Type_List[i_bond_type].Inter_List. */
	bond_index = bond_ctr[i_bond_type];

	/* Store CG site indices for the ith interaction. */
	Bond_Type_ptr->Inter_List[bond_index][0] = i1;
	Bond_Type_ptr->Inter_List[bond_index][1] = i2;
	Bond_Type_ptr->Inter_List[bond_index][2] = i3;
	Bond_Type_ptr->Inter_List[bond_index][3] = i4;

	/* Store i_bond_type in CG_struct[] for each site. */
	CG_struct[i1].bond_type[site_bond_ctr[i1]] = i_bond_type;
	CG_struct[i2].bond_type[site_bond_ctr[i2]] = i_bond_type;
	CG_struct[i3].bond_type[site_bond_ctr[i3]] = i_bond_type;
	CG_struct[i4].bond_type[site_bond_ctr[i4]] = i_bond_type;

	/* Link lists for the ith interaciton. */
	CG_struct[i1].bond_site[site_bond_ctr[i1]] =
	    Bond_Type_ptr->Inter_List[bond_index];
	CG_struct[i2].bond_site[site_bond_ctr[i2]] =
	    Bond_Type_ptr->Inter_List[bond_index];
	CG_struct[i3].bond_site[site_bond_ctr[i3]] =
	    Bond_Type_ptr->Inter_List[bond_index];
	CG_struct[i4].bond_site[site_bond_ctr[i4]] =
	    Bond_Type_ptr->Inter_List[bond_index];

	/* Keep track of the number of coeffs. for each site. */
	CG_struct[i1].nr_bond_coeffs += Bond_Type_ptr->N_coeff;
	CG_struct[i2].nr_bond_coeffs += Bond_Type_ptr->N_coeff;
	CG_struct[i3].nr_bond_coeffs += Bond_Type_ptr->N_coeff;
	CG_struct[i4].nr_bond_coeffs += Bond_Type_ptr->N_coeff;

	/* Update counters. */
	i += 5;
	i_bond += 1;
	site_bond_ctr[i1] += 1;
	site_bond_ctr[i2] += 1;
	site_bond_ctr[i3] += 1;
	site_bond_ctr[i4] += 1;
	bond_ctr[i_bond_type] += 1;
    }

    free(dihedral_list);
}


/*****************************************************************************************
store_pair_lists():
*****************************************************************************************/
void store_pair_lists(tW_gmx_topology *top, tW_CG_site * CG_struct,
		      int *i_type, int N_Bond_Int_Types,
		      tW_Bonded_Inter * Bond_Type_List, int *site_bond_ctr,
		      int *bond_ctr)
{
    int i = 0;
    int i_bond = 0;
    int i1, i2;
    int i_bond_type;
    int bond_index;
    tW_Bonded_Inter *Bond_Type_ptr;
    int *pair_list;
    int n_pairs;
 
    pair_list = top->get_pair_list(top);
    n_pairs = top->get_npairs(top);

    /* Loop over all InterMolec_NB_Pair interactions in current topology. */
    while (i < n_pairs) {
	/* Index to Bond_Type_List for the ith interaction. */
	i_bond_type = i_type[i_bond];

	/* Pointer to i_bond_type interaction. */
	Bond_Type_ptr = &(Bond_Type_List[i_bond_type]);

	/* Index for each CG site in the ith interaction. */
	i1 = pair_list[i + 1];
	i2 = pair_list[i + 2];

	/* Order the ith interaction occurs in Bond_Type_List[i_bond_type].Inter_List. */
	bond_index = bond_ctr[i_bond_type];

	/* Store CG site indices for the ith interaction. */
	Bond_Type_ptr->Inter_List[bond_index][0] = i1;
	Bond_Type_ptr->Inter_List[bond_index][1] = i2;

	/* Store i_bond_type in CG_struct[] for each site. */
	CG_struct[i1].bond_type[site_bond_ctr[i1]] = i_bond_type;
	CG_struct[i2].bond_type[site_bond_ctr[i2]] = i_bond_type;

	/* Link lists for the ith interaction. */
	CG_struct[i1].bond_site[site_bond_ctr[i1]] =
	    Bond_Type_ptr->Inter_List[bond_index];
	CG_struct[i2].bond_site[site_bond_ctr[i2]] =
	    Bond_Type_ptr->Inter_List[bond_index];

	/* Keep track of the number of coeffs. for each site. */
	CG_struct[i1].nr_bond_coeffs += Bond_Type_ptr->N_coeff;
	CG_struct[i2].nr_bond_coeffs += Bond_Type_ptr->N_coeff;

	/* Update counters. */
	i += 3;
	i_bond += 1;
	site_bond_ctr[i1] += 1;
	site_bond_ctr[i2] += 1;

	bond_ctr[i_bond_type] += 1;
    }

    free(pair_list);
}


/*****************************************************************************************
determine_bond_type():
*****************************************************************************************/
int determine_bond_type(int N_Bond_Int_Types,
			tW_Bonded_Inter * Bond_Type_List,
			tW_word * site_list)
{
    int i, j;
    tW_word Bond_Pair_List[2];

    /* Loop over all bond interaction types. */
    for (i = 0; i < N_Bond_Int_Types; i++) {
	/* Is the ith interaction a BondStretch? */
	if (strcmp(Bond_Type_List[i].name, B_BOND_STRETCH) == 0) {
	    /* Copy site pairs read from par.txt into local var. */
	    for (j = 0; j < 2; j++) {
		strcpy(Bond_Pair_List[j], Bond_Type_List[i].Site_Types[j]);
	    }

	    /* Do the two list agree? If so, return index of the ith interaction. */
	    if (check_word_list(2, site_list, Bond_Pair_List)) {
		return i;
	    }
	}
    }

    return -1;
}


/*****************************************************************************************
determine_angle_type():
*****************************************************************************************/
int determine_angle_type(int N_Bond_Int_Types,
			 tW_Bonded_Inter * Bond_Type_List,
			 tW_word * site_list)
{
    int i, j;
    tW_word Bond_Site_List[3];

    /* Loop over all bond interaction types. */
    for (i = 0; i < N_Bond_Int_Types; i++) {
	/* Is it an Angle interaction? */
	if (strcmp(Bond_Type_List[i].name, B_ANGLE) == 0) {
	    /* Copy sites from par.txt into local variable. */
	    for (j = 0; j < 3; j++) {
		strcpy(Bond_Site_List[j], Bond_Type_List[i].Site_Types[j]);
	    }

	    /* Note: site_list[] came from GROMACS topology, Bond_Site_list came from par.txt.       */
	    /* site_list[1] MUST be the SAME as Bond_Site_List[1], but other sites may be exchanged. */
	    if (strcmp(Bond_Type_List[i].Site_Types[1], site_list[1]) == 0) {
		/* Do the two list agree? If so, return index of the ith interaction. */
		if (check_word_list(3, site_list, Bond_Site_List)) {
		    return i;
		}
	    }
	}
    }

    return -1;
}


/*****************************************************************************************
determine_dihedral_type():
*****************************************************************************************/
int determine_dihedral_type(int N_Bond_Int_Types,
			    tW_Bonded_Inter * Bond_Type_List,
			    tW_word * site_list)
{
    int i, j;
    tW_word Bond_Site_List[4];

    /* Loop over all bond interaction types. */
    for (i = 0; i < N_Bond_Int_Types; i++) {
	/* Is is a dihedral interaction? */
	if (strcmp(Bond_Type_List[i].name, B_DIHEDRAL) == 0) {
	    /* Copy sites from par.txt into local variable. */
	    for (j = 0; j < 4; j++) {
		strcpy(Bond_Site_List[j], Bond_Type_List[i].Site_Types[j]);
	    }

	    /* The order of the two list must be the same or in reverse order for a dihedral angle. */
	    /* Do the two list agree? If so, return index of the ith interaction. */
	    if (match_word_list(4, site_list, Bond_Site_List)) {
		return i;
	    }
	}
    }

    return -1;
}


/*****************************************************************************************
determine_pair_type():
*****************************************************************************************/
int determine_pair_type(int N_Bond_Int_Types,
			tW_Bonded_Inter * Bond_Type_List,
			tW_word * site_list)
{
    int i, j;
    tW_word Bond_Pair_List[2];

    /* Loop over all bond interaction types. */
    for (i = 0; i < N_Bond_Int_Types; i++) {
	/* Is the ith interaction a bonded pair interaction. */
	if (strcmp(Bond_Type_List[i].name, B_NB_PAIR_BOND) == 0) {
	    /* Copy sites read from par.txt into local variable. */
	    for (j = 0; j < 2; j++) {
		strcpy(Bond_Pair_List[j], Bond_Type_List[i].Site_Types[j]);
	    }

	    /* Do the two list agree? If so, return index of the ith interaction. */
	    if (check_word_list(2, site_list, Bond_Pair_List)) {
		return i;
	    }
	}
    }

    return -1;
}

// MRD dumping everything
void dump_tW_system(char * fnm, tW_system * sys)
{
  FILE *fp;
  fp = fopen(fnm,"w");
  int i,j;
  fprintf(fp,"N_Site_Types: %d\n\n",sys->N_Site_Types);
if (sys->Site_Types)
{
  fprintf(fp,"Site_Types: ");
  for (i = 0; i < sys->N_Site_Types; ++i) { fprintf(fp," %s ",sys->Site_Types[i]); }
  fprintf(fp,"\n\n");
}
if (sys->Inter_Map)
{
  fprintf(fp,"Inter_Map: \n");
  for (i = 0; i < sys->N_Site_Types; ++i)
  {
    for (j = 0; j < sys->Inter_Map_Len[i]; ++j)
    {
      fprintf(fp," %s ",sys->Inter_Map[i][j]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
}

if (sys->Inter_iMap)
{ 
  fprintf(fp,"Inter_iMap: \n");
  for (i = 0; i < sys->N_Site_Types; ++i)
  {
    for (j = 0; j < sys->Inter_Map_Len[i]; ++j)
    {
      fprintf(fp," %d ",sys->Inter_iMap[i][j]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
}

if (sys->Inter_Map_Len)
{
  fprintf(fp,"Inter_Map_Len: ");
  for (i = 0; i < sys->N_Site_Types; ++i) { fprintf(fp," %d ",sys->Inter_Map_Len[i]); }
  fprintf(fp,"\n\n");
}

 
  fprintf(fp,"N_Inter_Types: %d\n\n",sys->N_Inter_Types);

if (sys->Inter_Types)
{
  int i;
  fprintf(fp,"idx        inter_name           inter_type       basis    i_basis     dr      R_min     R_max    N_pts   N_coeff   n_smooth     i0\n");
  fprintf(fp,"---   --------------------   ----------------   -------   -------   ------   -------   -------   -----   -------   --------   ------\n");
  for (i = 0; i < sys->N_Inter_Types; ++i)
  {
    fprintf(fp,"%3d   ",i);
    fprintf(fp,"%20s   ",sys->Inter_Types[i].inter_name);
    fprintf(fp,"%16s   ",sys->Inter_Types[i].inter_type);
    fprintf(fp,"%7s   ",sys->Inter_Types[i].basis);
    fprintf(fp,"%7d   ",sys->Inter_Types[i].i_basis);
    fprintf(fp,"%6.4f   ",sys->Inter_Types[i].dr);
    fprintf(fp,"%7.5f   ",sys->Inter_Types[i].R_min);
    fprintf(fp,"%7.5f   ",sys->Inter_Types[i].R_max);
    fprintf(fp,"%5d   ",sys->Inter_Types[i].N_pts);
    fprintf(fp,"%7d   ",sys->Inter_Types[i].N_coeff);
    fprintf(fp,"%8d   ",sys->Inter_Types[i].n_smooth);
    fprintf(fp,"%6d\n",sys->Inter_Types[i].i_0);
  }
  fprintf(fp,"\n");
}

  fprintf(fp,"N_Inter2_Types: %d\n\n",sys->N_Inter2_Types);

if (sys->Inter2_Type_List)
{
  fprintf(fp,"idx        inter_name         basis    name1   name2   N_inter   N_pts   N_coeff   i_basis    i_0      dr       R_0      R_max    n_smooth\n");
  fprintf(fp,"---   --------------------   -------   -----   -----   -------   -----   -------   -------   -----   ------   -------   -------   --------\n");
  for (i = 0; i < sys->N_Inter2_Types; ++i);
  {
    fprintf(fp,"%3d   ",i);
    fprintf(fp,"%20s   ",sys->Inter2_Type_List[i].inter_name);
    fprintf(fp,"%7s   ",sys->Inter2_Type_List[i].basis);
    fprintf(fp,"%5s   ",sys->Inter2_Type_List[i].name1);
    fprintf(fp,"%5s   ",sys->Inter2_Type_List[i].name2);
    fprintf(fp,"%7d   ",sys->Inter2_Type_List[i].N_inter);
    fprintf(fp,"%5d   ",sys->Inter2_Type_List[i].N_pts);
    fprintf(fp,"%7d   ",sys->Inter2_Type_List[i].N_coeff);
    fprintf(fp,"%5d   ",sys->Inter2_Type_List[i].i_basis);
    fprintf(fp,"%5d   ",sys->Inter2_Type_List[i].i_0);
    fprintf(fp,"%6.4f   ",sys->Inter2_Type_List[i].dr);
    fprintf(fp,"%7.5f   ",sys->Inter2_Type_List[i].R_0);
    fprintf(fp,"%7.5f   ",sys->Inter2_Type_List[i].R_max);
    fprintf(fp,"%8d    \n",sys->Inter2_Type_List[i].n_smooth);
  }
  fprintf(fp,"\n");
}


  fprintf(fp,"N_Bond_Int_Types: %d\n\n",sys->N_Bond_Int_Types);

if (sys->Bonded_Inter_Types)
{
  fprintf(fp,"idx          name               inter_name         basis    N_Int_Sites   site_i   site_j   site_k   site_l   i_basis   N_coeff    i_0      dr       R_0      R_max    N_pts   n_smooth   n_bonds_i...   N_instances\n");
  fprintf(fp,"---   ------------------   --------------------   -------   -----------   ------   ------   ------   ------   -------   -------   -----   ------   -------   -------   -----   --------   ------------   -----------\n");

  for (i = 0; i < sys->N_Bond_Int_Types; ++i)
  { 
    fprintf(fp,"%3d   ",i);
    fprintf(fp,"%18s   ",sys->Bonded_Inter_Types[i].name);
    fprintf(fp,"%20s   ",sys->Bonded_Inter_Types[i].inter_name);
    fprintf(fp,"%7s   ",sys->Bonded_Inter_Types[i].basis);
    fprintf(fp,"%11d   ",sys->Bonded_Inter_Types[i].N_Int_Sites);
    if (sys->Bonded_Inter_Types[i].N_Int_Sites == 2)
    {
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[0]);
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[1]);
      fprintf(fp,"%6s   ","N/A");
      fprintf(fp,"%6s   ","N/A");
    }
    else if (sys->Bonded_Inter_Types[i].N_Int_Sites == 3)
    {
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[0]);
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[1]);
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[2]);
      fprintf(fp,"%6s   ","N/A");
    }
    else if (sys->Bonded_Inter_Types[i].N_Int_Sites == 4)
    {
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[0]);
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[1]);
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[2]);
      fprintf(fp,"%6s   ",sys->Bonded_Inter_Types[i].Site_Types[3]);
    }
    else
    {
      fprintf(fp,"ERROR: Bonded_Inter_Types[%d].N_Int_Sites: %d\n",i,sys->Bonded_Inter_Types[i].N_Int_Sites);
    }
    fprintf(fp,"%7d   ",sys->Bonded_Inter_Types[i].i_basis);
    fprintf(fp,"%7d   ",sys->Bonded_Inter_Types[i].N_coeff);
    fprintf(fp,"%5d   ",sys->Bonded_Inter_Types[i].i_0);
    fprintf(fp,"%6.4f   ",sys->Bonded_Inter_Types[i].dr);
    fprintf(fp,"%7.5f   ",sys->Bonded_Inter_Types[i].R_0);
    fprintf(fp,"%7.5f   ",sys->Bonded_Inter_Types[i].R_max);
    fprintf(fp,"%5d   ",sys->Bonded_Inter_Types[i].N_pts);
    fprintf(fp,"%8d   ",sys->Bonded_Inter_Types[i].n_smooth);
    fprintf(fp,"%12d   ",sys->Bonded_Inter_Types[i].n_bonds_intramolec_pair_inter);
    fprintf(fp,"%11d   \n",sys->Bonded_Inter_Types[i].N_instances);
  }
  fprintf(fp,"\n");
}

  fprintf(fp,"N_coeff: %d\n",sys->N_coeff);
  fprintf(fp,"N_pack: %d\n",sys->N_pack);
  fprintf(fp,"nrexcl: %d\n",sys->nrexcl);
 
  fclose(fp);
}

void dump_tW_CG_struct(char *fnm, tW_CG_site * sites, int n_sites)
{
  FILE *fp;
  int i;
  fp = fopen(fnm,"w");
  fprintf(fp," idx     name    i_type   i_res   nr_excl   nr_bonds   nr_bond_coeffs\n");
  fprintf(fp,"-----   ------   ------   -----   -------   --------   --------------\n");
  for (i = 0; i < n_sites; ++i)
  {
    fprintf(fp,"%5d   ",i);
    fprintf(fp,"%6s   ",sites[i].name);
    fprintf(fp,"%6d   ",sites[i].i_type);
    fprintf(fp,"%5d   ",sites[i].i_res);
    fprintf(fp,"%7d   ",sites[i].nr_excl);
    fprintf(fp,"%8d   ",sites[i].nr_bonds);
    fprintf(fp,"%14d   \n",sites[i].nr_bond_coeffs);
    fprintf(fp,"\tr[0]: %g  %f   r[1]: %g  %f   r[2]: %g  %f   f[0]: %g  %f   f[1]: %g  %f   f[2]: %g  %f\n",sites[i].r[0],sites[i].r[0],sites[i].r[1],sites[i].r[1],sites[i].r[2],sites[i].r[2],sites[i].f[0],sites[i].f[0],sites[i].f[1],sites[i].f[1],sites[i].f[2],sites[i].f[2]);
  }
  fclose(fp);  
}
















