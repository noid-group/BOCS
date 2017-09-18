/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file cgff_fn.h 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
@brief Helper functions for the cgff calculation
*/

#ifndef CGFF_HELPER
#define CGFF_HELPER


void init_gmx_info(tW_gmx_info * info);

void reset_gmx_info(tW_gmx_info * info);



int setup_CG_struct(tW_system *sys, tW_gmx_topology *top, tW_CG_site *CG_struct, tW_Bonded_Inter *Bonded_Inter_Types);

void get_bonds(tW_system * sys, tW_gmx_topology *top, int n_sites,
	       tW_CG_site * CG_struct, int N_Bond_Int_Types,
	       tW_Bonded_Inter * Bond_Type_List);

void allocate_memory_bond_type_indices(tW_gmx_topology *top,
				       tW_bond_type_indices *bond_type_indices);

void clear_bond_lists(int N_Bond_Int_Types,
		      tW_Bonded_Inter * Bond_Type_List);

void clear_site_bond_info(int n_sites, tW_CG_site CG_struct[]);

void check_bond_interactions(tW_gmx_topology *top, tW_CG_site * CG_struct,
			     int N_Bond_Int_Types,
			     tW_Bonded_Inter * Bond_Type_List,
			     tW_bond_type_indices * bond_type_indices);

void allocate_memory_for_Bond_Inter_List(int N_Bond_Int_Types,
					 tW_Bonded_Inter * Bond_Type_List);

void allocate_memory_for_sites(int n_sites, tW_CG_site * CG_struct);

void check_bonds(tW_gmx_topology *top, tW_CG_site * CG_struct, int *i_type,
		 int N_Bond_Int_Types, tW_Bonded_Inter * Bond_Type_List);

void check_angles(tW_gmx_topology *top, tW_CG_site * CG_struct, int *i_type,
		  int N_Bond_Int_Types, tW_Bonded_Inter * Bond_Type_List);

void check_dihedrals(tW_gmx_topology *top, tW_CG_site * CG_struct,
		     int *i_type, int N_Bond_Int_Types,
		     tW_Bonded_Inter * Bond_Type_List);

void check_pairs(tW_gmx_topology *top, tW_CG_site * CG_struct, int *i_type,
		 int N_Bond_Int_Types, tW_Bonded_Inter * Bond_Type_List);

int determine_bond_type(int N_Bond_Int_Types,
			tW_Bonded_Inter * Bond_Type_List,
			tW_word * site_list);

int determine_angle_type(int N_Bond_Int_Types,
			 tW_Bonded_Inter * Bond_Type_List,
			 tW_word * site_list);

int determine_dihedral_type(int N_Bond_Int_Types,
			    tW_Bonded_Inter * Bond_Type_List,
			    tW_word * site_list);

int determine_pair_type(int N_Bond_Int_Types,
			tW_Bonded_Inter * Bond_Type_List,
			tW_word * site_list);

void store_bond_interactions(tW_gmx_topology *top, tW_CG_site * CG_struct,
			     int N_Bond_Int_Types,
			     tW_Bonded_Inter * Bond_Type_List,
			     tW_bond_type_indices bond_type_indices,
			     int n_sites);

void store_bond_lists(tW_gmx_topology *top, tW_CG_site * CG_struct,
		      int *i_type, int N_Bond_Int_Types,
		      tW_Bonded_Inter * Bond_Type_List, int *site_bond_ctr,
		      int *bond_ctr);

void store_angle_lists(tW_gmx_topology *top, tW_CG_site * CG_struct,
		       int *i_type, int N_Bond_Int_Types,
		       tW_Bonded_Inter * Bond_Type_List,
		       int *site_bond_ctr, int *bond_ctr);

void store_dihedral_lists(tW_gmx_topology *top, tW_CG_site * CG_struct,
			  int *i_type, int N_Bond_Int_Types,
			  tW_Bonded_Inter * Bond_Type_List,
			  int *site_bond_ctr, int *bond_ctr);

void store_pair_lists(tW_gmx_topology *top, tW_CG_site * CG_struct,
		      int *i_type, int N_Bond_Int_Types,
		      tW_Bonded_Inter * Bond_Type_List, int *site_bond_ctr,
		      int *bond_ctr);

void print_Bond_Types(FILE * fp, tW_system sys);

#endif
