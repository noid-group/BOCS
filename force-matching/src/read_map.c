/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file read_map.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael DeLyser
@brief Functions related to reading in a map.top file 
*/

//c library includes
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>

//local includes
#include "read_map.h"
#include "safe_mem.h"
#include "io_read.h"

static const char ATOMTYPE[] = "atomtypes";
static const char MOLTYPE[] = "moleculetype";
static const char SITETYPE[] = "sitetypes";
static const char MOLECULES[] = "molecules";

#define eps_norm 1.e-10

#define DEBUG_get_clean_line 0
#define DEBUG_get_map        0
#define DEBUG_get_map_1      0


void print_site_mapping(tW_site_map map)
{
    int i;

    printf("--- in print_site_mapping ---\n");

    if (map.cg_name) {
	printf("map.cg_name:  %s\n", map.cg_name);
    }
    if (map.mol_name) {
	printf("map.mol_name: %s\n", map.mol_name);
    }

    printf("map.res_no  : %d \n", map.res_no);
    printf("map.n_atms  : %d \n", map.n_atms);
    printf("map.map_type: %s \n", map.map_type);

    for (i = 0; i < map.n_atms; i++) {
	// NB reading and printing atm indices starting from 1, but stored
	// from 0
	printf("i_atm[%d]: %d \t ", i, map.i_atm[i] + 1);
	if (map.atm_name) {
	    printf("atm_name: %s \t", map.atm_name[i]);
	}
	printf("c_Ii: %lf \n", map.c_Ii[i]);
    }

}

void print_mapping_N(int N_sites, tW_site_map map[])
{
    int i;
    for (i = 0; i < N_sites; i++) {
	printf("--- site: %d ---\n", i);
	print_site_mapping(map[i]);
    }
}

bool is_blank(char *line)
{
	int i = 0;
	bool is_blank = TRUE;
	size_t len = strlen(line);

	for (i=0; i<len; i++)
	{
		if ( !isspace(line[i]) )
		{
			is_blank = FALSE;
		}
	}

	return is_blank;
}

int skip_space(char *line, size_t nbytes, int start)
{
    int end = start;

    while (isspace(line[end]) && end <= nbytes) {
	end++;
    }

    return end;
}


char *get_directive(char *line, char *directive, bool * have_directive,
		    int start)
{
    int j, k;
    int len_directive;


    j = start + 1;
    len_directive = 0;
    while (line[j] != ']') {
	if (!isspace(line[j])) {
	    len_directive++;
	}
	j++;
    }

    directive = ecalloc(len_directive + 1, sizeof(char));

    j = start + 1;
    k = 0;
    while (line[j] != ']') {
	if (!isspace(line[j])) {
	    directive[k] = line[j];
	    k++;
	}
	j++;
    }
    if (strcmp(directive, MOLTYPE) == 0 || strcmp(directive, ATOMTYPE) == 0
	|| strcmp(directive, SITETYPE) == 0
	|| strcmp(directive, MOLECULES) == 0) {
	*have_directive = true;
    } else {
	fprintf(stderr,
		"ERROR: Directive '%s' not recognized. Terminating.\n",
		directive);
	exit(1);
    }

    return directive;
}

/*
MRD 1.12.2018
We have recently had a problem with the following function (get_CG_map). That problem has been patched by 
moving an error check to a different block of code (look for the comments, I pointed it out). I also added
a new error check to make sure all CG sites have the expected number of atoms listed under the [ atomtypes ]
directive (although the regular checks should catch this issue for sites having too many atoms. Originally,
it would catch if any site except the last site was atom deficient, but I removed this check in order to allow
"interleaved" mappings (look for the comment explanation later)). 
*/
tW_site_map *get_CG_map(FILE * map_top, int *N_sites)
{
    int h, i, j, k, l;		// generic indices
    int ATOM, SITE, N_MOL, TOTAL_ATOMS, *ATOMS, moltype;	// named indices for navigating input file
    int count = 0;
    int atm_count = 0;
    int n_atm_chk = 0;
    int n_site_chk = 0;
    int res_no = 0;
    int n_sites = 0;
    int total_atoms = 0;
    int total_sites = 0;

    bool have_directive = false;
    bool not_found = true;

    size_t nbytes = 100;

    mapping_op map;
    tW_site_map *CG_map;

    char *site_type = NULL;
    char *molname = NULL;
    char *directive = NULL;
    char *line;

    map.n_types = 0;
    TOTAL_ATOMS = 0;
    ATOM = 0;
    SITE = 0;
    N_MOL = 0;

    line = ecalloc(nbytes, sizeof(char));

   /* counts molecule types by finding [ moleculetype ] directives so memory can be allocated */
    while (getline(&line, &nbytes, map_top) != -1) {
	i = 0;
	i = skip_space(line, nbytes, i);

	if (line[i] != ';' && nbytes > i) {
	    if (line[i] == '[') {
		directive = get_directive(line, directive, &have_directive, i);

		if (strcmp(directive, MOLTYPE) == 0) {
		    map.n_types++;
		}
	    } else if (is_blank(line) && have_directive) {
		free(directive);
		have_directive = false;
	    }
	}

    }
    rewind(map_top);

    if (map.n_types == 0)
    {
      fprintf(stderr,"ERROR: map.n_types is 0!\n");
      fprintf(stderr,"\tDid you include the [ moleculetype ] directive for each type of CG molecule?\n");
      exit(1);
    }

    map.n_molecules = ecalloc(map.n_types, sizeof(int));
    map.molecule_types = ecalloc(map.n_types, sizeof(mapping_type));

    for (i = 0; i < map.n_types; i++) {
	map.molecule_types[i].n_sites = 0;
    }


    // reads in data for each molecule type
    moltype = -1;
    while (getline(&line, &nbytes, map_top) != -1) {

//	line = strip_leading_space(line, nbytes);

	i = 0;
	i = skip_space(line, nbytes, i);

	// check for array overrun before it occurs
	if (moltype > map.n_types)
	{
		fprintf(stderr, "ERROR: Incorrect number of molecule types in your map.top file.\n");
		exit(1);
	} 

/* MRD Moved this down a bit 1.12.2018
        else if (SITE > map.molecule_types[moltype].n_sites)
	{
		fprintf(stderr, "ERROR: Found more site types than expected under [ sitetypes ] directory.\n");
		exit(1);
	}*/


	if (line[i] != ';' && nbytes > i) 
        {
	    if (line[i] == '[') 
            {
		directive = get_directive(line, directive, &have_directive, i);

		if (strcmp(directive, SITETYPE) == 0) {
		    // read in # of sites and allocate memory
		    sscanf(line, "%*[[ ]%*s%*[] ]%d", &map.molecule_types[moltype].n_sites);

		    ATOMS = ecalloc(map.molecule_types[moltype].n_sites, sizeof(int));

		    // Allocate memory for variables of size (n_sites) 
		    map.molecule_types[moltype].site_names =
			ecalloc(map.molecule_types[moltype].n_sites,
			       sizeof(char *));
		    map.molecule_types[moltype].site_masses =
			ecalloc(map.molecule_types[moltype].n_sites,
			       sizeof(char *));
		    map.molecule_types[moltype].atom_names =
			ecalloc(map.molecule_types[moltype].n_sites,
			       sizeof(char **));
		    map.molecule_types[moltype].n_atoms =
			ecalloc(map.molecule_types[moltype].n_sites,
			       sizeof(int));
		    map.molecule_types[moltype].indices =
			ecalloc(map.molecule_types[moltype].n_sites,
			       sizeof(int *));
		    map.molecule_types[moltype].weights =
			ecalloc(map.molecule_types[moltype].n_sites,
			       sizeof(double *));
		    map.molecule_types[moltype].map_type =
			ecalloc(map.molecule_types[moltype].n_sites,
			       sizeof(char *));

		}
	    } 
            else if (is_blank(line) && have_directive) 
            {
		if (strcmp(directive, SITETYPE) == 0 && SITE != map.molecule_types[moltype].n_sites )
		{
		    fprintf(stderr, "ERROR: Wrong # of site types for %s. Expected %d, but found %d.\n", 
			map.molecule_types[moltype].molname, map.molecule_types[moltype].n_sites, SITE);
		    exit(1);
                } else if (strcmp(directive, ATOMTYPE) == 0) {
		      if ( TOTAL_ATOMS != map.molecule_types[moltype].total_atoms ) {
		 	   fprintf(stderr, "ERROR: Wrong # of atom types for %s. Expected %d, but found %d.\n", 
				map.molecule_types[moltype].molname, map.molecule_types[moltype].total_atoms, TOTAL_ATOMS);
		   	    exit(1);
		      }
		}

		free(directive);
		have_directive = false;

	    } 
            else if (have_directive) 
            {
		if (strcmp(directive, MOLTYPE) == 0) {

		    if (moltype >=0 && SITE > map.molecule_types[moltype].n_sites)
		    {
			fprintf(stderr, "ERROR: Found more site types than expected under [ sitetypes ] directory.\n");
		   	exit(1);
		    }

		    moltype++;
		    map.molecule_types[moltype].molname = ecalloc(strlen(line) + 1, sizeof(char));
		    sscanf(line, " %s %d ", map.molecule_types[moltype].molname, &map.molecule_types[moltype].total_atoms);
// MRD 1.12.2018
                    map.molecule_types[moltype].atom_sites = ecalloc(map.molecule_types[moltype].total_atoms,sizeof(char *));
                    for (j = 0; j < map.molecule_types[moltype].total_atoms; ++j)
                    {
                        map.molecule_types[moltype].atom_sites[j] = ecalloc(20,sizeof(char));
                    }

		    SITE = 0;
		    ATOM = 0;
		    TOTAL_ATOMS = 0;

		} else if (strcmp(directive, SITETYPE) == 0) {

                    if (SITE >= map.molecule_types[moltype].n_sites) // MRD moved this here 1.12.2018
                    {
                        fprintf(stderr, "ERROR: Found more site types than expected under [ sitetypes ] directory.\n");
                        exit(1);
                    }

		    n_site_chk++;
		    // parse line to find how many atoms per site there are
		    map.molecule_types[moltype].site_names[SITE] = ecalloc(strlen(line) + 1, sizeof(char));
		    map.molecule_types[moltype].map_type[SITE] = ecalloc(strlen(line) + 1, sizeof(char));

		    sscanf(line, " %s %s %d ",
			   map.molecule_types[moltype].site_names[SITE],
			   map.molecule_types[moltype].map_type[SITE],
			   &map.molecule_types[moltype].n_atoms[SITE]);

// MRD 1.12.2018    We generally don't consider models that have CG sites with no atoms mapped to them.
//                  HOWEVER, IIRC, in the literature, people have coupled the "actual" cg sites with
//                  virtual sites in order to try to fix the dynamics of the cg sites. SO, I will allow
//                  CG sites with 0 atoms mapped to them, but I'm going to print a noisy warning message.
                    if (map.molecule_types[moltype].n_atoms[SITE] == 0)
                    {
                        fprintf(stderr,"\n\nWARNING: You have included cg site type %s that has 0 atoms that map to it.\n",
                                                                                map.molecule_types[moltype].site_names[SITE]);
                        fprintf(stderr,"\t This is not part of the \"traditional\" MSCG method, but other CG methods have\n");
                        fprintf(stderr,"\t since been developed that would use these sites. Note, all CG sites with 0 atoms\n");
                        fprintf(stderr,"\t that map to them will have 0s for all components of the x, v, and f vectors\n\n\n");
                    }

		    map.molecule_types[moltype].indices[SITE] =
			ecalloc(map.molecule_types[moltype].n_atoms[SITE],
			       sizeof(int));
		    map.molecule_types[moltype].weights[SITE] =
			ecalloc(map.molecule_types[moltype].n_atoms[SITE],
			       sizeof(double));
		    map.molecule_types[moltype].atom_names[SITE] =
			ecalloc(map.molecule_types[moltype].n_atoms[SITE],
			       sizeof(char *));

		    for (j = 0; j < map.molecule_types[moltype].n_atoms[SITE]; j++) {
			map.molecule_types[moltype].atom_names[SITE][j] = ecalloc(strlen(line) + 1, sizeof(char));
		    }

		    if (strcmp("com", map.molecule_types[moltype].map_type[SITE]) == 0)
		    {
			map.molecule_types[moltype].site_masses[SITE] = 0;

		    } else if (strcmp("cog",map.molecule_types[moltype].map_type[SITE]) == 0) {
			map.molecule_types[moltype].site_masses[SITE] =
			    map.molecule_types[moltype].n_atoms[SITE];

		    } else if (strcmp("usr",map.molecule_types[moltype].map_type[SITE]) == 0) {
			// Currently, usr is a com mode where the user is
			// forced to pre-normalize the weights or the program
			// will exit
			map.molecule_types[moltype].site_masses[SITE] = 0;

		    } else {
			fprintf(stderr,
				"ERROR: Mapping type '%s' not recognized. Terminating.\n",
				map.molecule_types[moltype].
				map_type[SITE]);
			exit(1);
		    }


		    SITE++;

		} else if (strcmp(directive, ATOMTYPE) == 0) {
		    TOTAL_ATOMS++;
		    /* parse line to find out what site the atom belongs to,  if any */
		    site_type = ecalloc(strlen(line) + 1, sizeof(char));
		    sscanf(line, " %*d %*s %s %*d", site_type);

// MRD 1.12.2018
                    strcpy(map.molecule_types[moltype].atom_sites[ATOM],site_type);

		    n_atm_chk += 1;

		    not_found = true;

		    /*loop over sites to find which this atom belongs to */ 
		    for (j = 0; j < map.molecule_types[moltype].n_sites; j++)
		    {

			if (strcmp (site_type, map.molecule_types[moltype].site_names[j]) == 0) {

			    ATOMS[j]++;

			    // Did we find the right number of atoms for the previous site? 
/* MRD 1.12.2018 I commented this out because this case will be handled by my new error check.
 * However, this check does not allow for "interleaved" CG mappings, e.g.
 * 
 * [ atomtypes ]
 * 1 C Site1 12.011
 * 2 C Site2 12.011
 * 3 C Site1 12.011
 * ...
 *
 * wouldn't work because this error check assumes that if we've found Site2, then we all atoms that map to Site1 must
 * have already been given.
			    if (j>0 && map.molecule_types[moltype].n_atoms[j-1] != ATOMS[j-1])
			    {
		                    fprintf(stderr, "ERROR: Too few atoms for site %s in %s. Expected %d, but found %d.\n", 
		               	        map.molecule_types[moltype].site_names[j-1], map.molecule_types[moltype].molname, 
		               	        map.molecule_types[moltype].n_atoms[j-1], ATOMS[j-1]);
		                   exit(1);
			    }*/

			    // Are there too many atoms in this site?
			    if (ATOMS[j] > map.molecule_types[moltype].n_atoms[j])
			    {
			        fprintf(stderr, "ERROR: Too many atoms for site %s in %s. Expected %d, but found at least %d.\n", 
				   map.molecule_types[moltype].site_names[j], map.molecule_types[moltype].molname, 
				   map.molecule_types[moltype].n_atoms[j], ATOMS[j]);
			        exit(1);
			    }

			    // find the index l of the atom within its site
//			    l = ATOM;
//			    for (k = 0; k < j; k++) {
//				l = l - map.molecule_types[moltype].n_atoms[k];
//			    }
//			    l = l % map.molecule_types[moltype].n_atoms[j];
			    l = ATOMS[j]-1;

			    if (strcmp("com", map.molecule_types[moltype].map_type[j]) == 0) {
				sscanf(line, " %d %s %*s %lf",
				       &map.molecule_types[moltype].indices[j][l],
				       map.molecule_types[moltype].atom_names[j][l],
				       &map.molecule_types[moltype].weights[j][l]);

				map.molecule_types[moltype].site_masses[j]
				    += map.molecule_types[moltype].weights[j][l];

			    } else if (strcmp("cog", map.molecule_types[moltype].map_type[j]) == 0) {
				sscanf(line, " %d %s %*s ",
				       &map.molecule_types[moltype].indices[j][l],
				       map.molecule_types[moltype].atom_names[j][l]);

				map.molecule_types[moltype].weights[j][l] =  1;

			    } else if (strcmp("usr", map.molecule_types[moltype].map_type[j]) == 0) {
				sscanf(line, " %d %s %*s %lf", &map.molecule_types[moltype].indices[j][l],
				       map.molecule_types[moltype].atom_names[j][l],
				       &map.molecule_types[moltype].weights[j][l]);

				map.molecule_types[moltype].site_masses[j]
				   += map.molecule_types[moltype].weights[j][l];
			    }

			    not_found = false;
			    ATOM++;
			} else if (strcmp("none", site_type) == 0) {
			    not_found = false;
			}

		    }

		    if (not_found) {
			fprintf(stderr,	"ERROR: Site type %s not found in site list. Terminating.\n", site_type);
			exit(1);
		    }

		    free(site_type);
		} else if (strcmp(directive, MOLECULES) == 0) {
		    // get final counts for molecule_types in the system
		    molname = ecalloc(strlen(line) + 1, sizeof(char));
		    sscanf(line, " %s %d ", molname, &N_MOL);
                    bool mol_found = FALSE;
		    for (j = 0; j < map.n_types; j++)	
		    {
			if (strcmp(molname, map.molecule_types[j].molname) == 0) 
                        {
			    map.n_molecules[j] = N_MOL;
                            mol_found = TRUE;
			}
		    }
                    if (! mol_found) // MRD 03.28.2019
                    {
                       fprintf(stderr,"ERROR: found site type %s under directive [molecules]\n",molname);
                       fprintf(stderr,"\tHowever, we only found the following %d types of molecules: \n",map.n_types);
                       for (j = 0; j < map.n_types; ++j)
                       {
                           fprintf(stderr,"\t%s \n",map.molecule_types[j].molname);
                       }
                       exit(EXIT_FAILURE); 
                    }
		    free(molname);
		}

	    }
            else
            {
              fprintf(stderr,"ERROR: found line: %s",line);
              fprintf(stderr,"\toutside of a directive. Directive lines should look like, e.g., \n[ moleculetypes ]\n");
              exit(1);
            }
	}
    }


/* MRD 1.12.2018 New Error Checks */
    for (i = 0; i < map.n_types; ++i)
    {
        for (j = 0; j < map.molecule_types[i].n_sites; ++j)
        {
            int n_atoms_for_site = 0;
            for (k = 0; k < map.molecule_types[i].total_atoms; ++k)
            {
                if (strcmp(map.molecule_types[i].site_names[j],map.molecule_types[i].atom_sites[k]) == 0)
                {
                    ++n_atoms_for_site;
                }
            }
            if (n_atoms_for_site != map.molecule_types[i].n_atoms[j])
            {
                fprintf(stderr,"ERROR:  Under the [ sitetypes ] directive, site type %s is said to have %d atoms map to it\n",
                                        map.molecule_types[i].site_names[j],map.molecule_types[i].n_atoms[j]);
                fprintf(stderr,"\tHowever, we found %d atoms under the [ atomtypes ] directive with site %s given\n",
                                        n_atoms_for_site,map.molecule_types[i].site_names[j]);
                exit(1);
            }
        }
    }

//    // check for common input errors - wrong # of sites or atoms in map.top
//    for (j = 0; j < map.n_types; j++) {
//	total_atoms += map.molecule_types[j].total_atoms;
//	total_sites += map.molecule_types[j].n_sites;
//    }
//
//    if (n_site_chk != total_sites) {
//	fprintf(stderr, "ERROR: Incorrect site count (%d). Expected %d sites. Terminating.\n", n_site_chk, total_sites);
//	exit(1);
//    }
//
//    if (n_atm_chk != total_atoms) {
//	fprintf(stderr,
//		"ERROR: Incorrect atom count (%d). Expected %d atoms. Terminating.\n",
//		n_atm_chk, total_atoms);
//	exit(1);
//    }

    for (i = 0; i < map.n_types; i++) {
	n_sites += map.n_molecules[i] * map.molecule_types[i].n_sites;
    }

    CG_map = ecalloc(n_sites, sizeof(tW_site_map));
    (*N_sites) = n_sites;

    for (i = 0; i < map.n_types; i++) {
	for (h = 0; h < map.n_molecules[i]; h++) {
	    for (j = 0; j < map.molecule_types[i].n_sites; j++) {
		CG_map[count].n_atms = map.molecule_types[i].n_atoms[j];
		CG_map[count].res_no = res_no;

		CG_map[count].cg_name =  ecalloc(strlen(map.molecule_types[i].site_names[j]) + 1, sizeof(char));
		CG_map[count].mol_name = ecalloc(strlen(map.molecule_types[i].molname) + 1,  sizeof(char));
		CG_map[count].i_atm =    ecalloc(map.molecule_types[i].n_atoms[j], sizeof(int));
		CG_map[count].atm_name = ecalloc(map.molecule_types[i].n_atoms[j],  sizeof(char *));
		CG_map[count].c_Ii =  ecalloc(map.molecule_types[i].n_atoms[j],  sizeof(double));


		if (strcmp("usr", map.molecule_types[i].map_type[j]) == 0) {
		    if (map.molecule_types[i].site_masses[j] - 1 > eps_norm) {
			fprintf(stderr,
				"ERROR: Mapping coefficients not normalized. (Norm = %lg)\n",
				map.molecule_types[i].site_masses[j]);
			exit(1);
		    }
		}

		for (k = 0; k < map.molecule_types[i].n_atoms[j]; k++) {
		    CG_map[count].atm_name[k] = ecalloc(strlen(map.molecule_types[i].atom_names[j][k]) + 1, sizeof(char));
		    strcpy(CG_map[count].atm_name[k], map.molecule_types[i].atom_names[j][k]);
		    CG_map[count].i_atm[k] = map.molecule_types[i].indices[j][k] + atm_count - 1;
		    CG_map[count].c_Ii[k] = map.molecule_types[i].weights[j][k] / map.molecule_types[i].site_masses[j];
		}

		strcpy(CG_map[count].cg_name, map.molecule_types[i].site_names[j]);
		strcpy(CG_map[count].mol_name, map.molecule_types[i].molname);
		strcpy(CG_map[count].map_type, map.molecule_types[i].map_type[j]);
		count++;

		if (j == map.molecule_types[i].n_sites - 1) {
		    res_no++;
		    atm_count += map.molecule_types[i].total_atoms;
		}
	    }
	}
    }

    free(line);
    free(map.n_molecules);

    return CG_map;
}

tW_site_map *NEW_get_CG_map(FILE * map_top, int *N_sites) // THIS IS THE NEW VERSION. MRD 1/9/2018
{
  tW_line inp_line;
  int n_molecule_types = 0;

  while (get_next_line(map_top,inp_line) != -1)
  {
    if (strstr(inp_line,"moleculetype") != NULL) { ++n_molecule_types; }
  }
  rewind(map_top);

  char **molecule_names = (char **) ecalloc(n_molecule_types,sizeof(char *));
  int *n_atoms_per_molecule = (int *) ecalloc(n_molecule_types,sizeof(int));
  int *n_sites_per_molecule = (int *) ecalloc(n_molecule_types,sizeof(int));
  int *n_molecules = (int *) ecalloc(n_molecule_types,sizeof(int));
  char ***site_names = (char ***) ecalloc(n_molecule_types,sizeof(char **));
  char ***mapping_types = (char ***) ecalloc(n_molecule_types,sizeof(char **));
  int **n_atoms_per_site = (int **) ecalloc(n_molecule_types,sizeof(int *));
  int **n_atom_in_molecule = (int **) ecalloc(n_molecule_types,sizeof(int *));
  char ***atom_type = (char ***) ecalloc(n_molecule_types,sizeof(char **));
  char ***atom_n_site_type = (char ***) ecalloc(n_molecule_types,sizeof(char **));
  float **atom_weight = (float **) ecalloc(n_molecule_types,sizeof(float *));
  float **site_weight = (float **) ecalloc(n_molecule_types,sizeof(float *));

  int mol_type_idx = 0;
  int i, j, k, l;
  int test_sscanf;
  char *p;
  char *inp_word = (char *) ecalloc(50,sizeof(char));

  for (i = 0; i < n_molecule_types; ++i)
  {
    molecule_names[i] = (char *) ecalloc(50,sizeof(char));
  }

  while (mol_type_idx < n_molecule_types)
  {
    l = mol_type_idx;
    get_next_line(map_top,inp_line);
    if (strstr(inp_line,"moleculetype") != NULL)
    {
      get_next_line(map_top,inp_line);
      test_sscanf = sscanf(inp_line," %s %d ",(molecule_names[mol_type_idx]),&(n_atoms_per_molecule[mol_type_idx]));
      if (test_sscanf != 2)
      {
        fprintf(stderr,"ERROR: expected to find name of moleculetype %d and number of atoms on line after directive [ moleculetype ]\n",mol_type_idx);
        fprintf(stderr,"\tinp_line: %s",inp_line);
        exit(1);
      }

      n_atom_in_molecule[mol_type_idx] = (int *) ecalloc(n_atoms_per_molecule[mol_type_idx],sizeof(int));
      atom_type[mol_type_idx] = (char **) ecalloc(n_atoms_per_molecule[mol_type_idx],sizeof(char *));
      atom_n_site_type[mol_type_idx] = (char **) ecalloc(n_atoms_per_molecule[mol_type_idx],sizeof(char *));
      atom_weight[mol_type_idx] = (float *) ecalloc(n_atoms_per_molecule[mol_type_idx],sizeof(float));
      site_weight[mol_type_idx] = (float *) ecalloc(n_sites_per_molecule[mol_type_idx],sizeof(float));
      for (i = 0; i < n_atoms_per_molecule[mol_type_idx]; ++i)
      {
        atom_type[l][i] = (char *) ecalloc(50,sizeof(char));
        atom_n_site_type[l][i] = (char *) ecalloc(50,sizeof(char));
        site_weight[l][i] = 0.0;
      }
    }
    else
    {
      fprintf(stderr,"ERROR: expected directive [ moleculetype ] for molecule type %d\n",mol_type_idx+1);
      fprintf(stderr,"\tinp_line: %s",inp_line);
      exit(1);
    }

    get_next_line(map_top,inp_line);
    if (strstr(inp_line,"sitetypes") != NULL)
    {
      p = strstr(inp_line,"]");
      ++p;
      test_sscanf = sscanf(p," %d ", &(n_sites_per_molecule[mol_type_idx]));
      if (test_sscanf != 1)
      {
        fprintf(stderr,"ERROR: expected to find number of sites in moleculetype %s on line with directive [ sitetypes ]\n",molecule_names[mol_type_idx]);
        fprintf(stderr,"\tinp_line: %s",inp_line);
        exit(1);
      }
      site_names[mol_type_idx] = (char **) ecalloc(n_sites_per_molecule[mol_type_idx],sizeof(char *));
      mapping_types[mol_type_idx] = (char **) ecalloc(n_sites_per_molecule[mol_type_idx],sizeof(char *));
      n_atoms_per_site[mol_type_idx] = (int *) ecalloc(n_sites_per_molecule[mol_type_idx],sizeof(int));
      for (i = 0; i < n_sites_per_molecule[mol_type_idx]; ++i)
      {
        site_names[mol_type_idx][i] = (char *) ecalloc(50,sizeof(char));
        mapping_types[mol_type_idx][i] = (char *) ecalloc(50,sizeof(char));
        get_next_line(map_top,inp_line);
        test_sscanf = sscanf(inp_line," %s %s %d ",(site_names[mol_type_idx][i]),(mapping_types[mol_type_idx][i]),&(n_atoms_per_site[mol_type_idx][i]));
        if (test_sscanf != 3)
        {
          fprintf(stderr,"ERROR: expected site_name, mapping_type, and n_atoms_per_site on line %d after [ sitetypes ] directive!\n",i);
          fprintf(stderr,"\tinp_line: %s",inp_line);
          exit(1);
        }
      }
    }
    else
    {
      fprintf(stderr,"ERROR: expected to find directive [ sitetypes ] # after [ moleculetype ] directive!\n");
      fprintf(stderr,"\tinp_line: %s",inp_line);
      exit(1);
    }

    get_next_line(map_top,inp_line);
    if (strstr(inp_line,"atomtypes") != NULL)
    {
      for (i = 0; i < n_atoms_per_molecule[mol_type_idx]; ++i)
      {
        get_next_line(map_top,inp_line);
        test_sscanf = sscanf(inp_line," %d %s %s %f ",&(n_atom_in_molecule[mol_type_idx][i]),
                                                      (atom_type[mol_type_idx][i]),
                                                      (atom_n_site_type[mol_type_idx][i]),
                                                      &(atom_weight[mol_type_idx][i]));
        if (test_sscanf != 4)
        {
          fprintf(stderr,"ERROR: expected to find atom_number, atom_type, site_type, atom_weighting on line %d after directive [ atomtypes ]\n",i);
          fprintf(stderr,"\tinp_line: %s",inp_line);
          fprintf(stderr,"Make sure the number of lines following the [ atomtypes ] is equal to the number\n");
          fprintf(stderr,"of atoms specified under the most recent [ moleculetype ] directive\n");
          exit(1);
        }
        // increment site weights
        for (j = 0; j < n_sites_per_molecule[mol_type_idx]; ++j)
        {
          if (strcmp(atom_n_site_type[mol_type_idx][i],site_names[mol_type_idx][j]) == 0)
          {
            site_weight[mol_type_idx][j] += atom_weight[mol_type_idx][i];
          }
        }
      }
    }
    else
    {
      fprintf(stderr,"ERROR: expected to find directive [ atomtypes ] after [ sitetypes ] directive!\n");
      fprintf(stderr,"\tinp_line: %s",inp_line);
      fprintf(stderr,"Be sure the number of (non comment and non-blank) lines after [ sitetypes ]\n");
      fprintf(stderr,"is equal to the number of site types specified on the line with directive [ sitetypes ]\n");
      fprintf(stderr,"expected n_sitetypes: %d\n",n_sites_per_molecule[mol_type_idx]);
      exit(1);
    }
    ++mol_type_idx;
  }


  get_next_line(map_top,inp_line);
  if (strstr(inp_line,"molecules") == NULL)
  {
    fprintf(stderr,"ERROR: expected to find directive [ molecules ]\n");
    fprintf(stderr,"\tinp_line: %s",inp_line);
    exit(1);
  }

  for (mol_type_idx = 0; mol_type_idx < n_molecule_types; ++mol_type_idx)
  {
    get_next_line(map_top,inp_line);
    test_sscanf = sscanf(inp_line," %s %d ",inp_word,&(n_molecules[mol_type_idx]));
    if (test_sscanf != 2)
    {
      fprintf(stderr,"ERROR: expected to find molecule_type and n_molecules for molecule idx %d\n",mol_type_idx);
      fprintf(stderr,"\tinp_line: %s",inp_line);
      exit(1);
    }
  }
  if (get_next_line(map_top,inp_line) != -1)
  {
    fprintf(stderr,"WARNING: expected to be out of lines for mapping topology file\n");
    fprintf(stderr,"\tinstead, found inp_line: %s",inp_line);
  }
  fclose(map_top);
  int n_CG_sites = 0;
  for (i = 0; i < n_molecule_types; ++i) { n_CG_sites += n_molecules[i] * n_sites_per_molecule[i]; }

  *N_sites = n_CG_sites;

  tW_site_map *CG_map = (tW_site_map *) ecalloc(n_CG_sites,sizeof(tW_site_map));

  int site_idx = 0;
  int res_no = 0;
  int abs_at_num = 0;
  for (i = 0; i < n_molecule_types; ++i)
  {
    for (j = 0; j < n_molecules[i]; ++j)
    {
      int at_idx_in_mol = 0;
      for (k = 0; k < n_sites_per_molecule[i]; ++k)
      {
        CG_map[site_idx].res_no = res_no;
        CG_map[site_idx].n_atms = n_atoms_per_site[i][k];
        CG_map[site_idx].i_atm = (int *) ecalloc(CG_map[site_idx].n_atms,sizeof(int));
        CG_map[site_idx].c_Ii = (double *) ecalloc(CG_map[site_idx].n_atms,sizeof(double));
        CG_map[site_idx].cg_name = (char *) ecalloc(50,sizeof(char));
        CG_map[site_idx].mol_name = (char *) ecalloc(50,sizeof(char));

        CG_map[site_idx].atm_name = (char **) ecalloc(CG_map[site_idx].n_atms,sizeof(char *));
        for (l = 0; l < CG_map[site_idx].n_atms; ++l)
        {
          CG_map[site_idx].atm_name[l] = (char *) ecalloc(50,sizeof(char));
        }

        for (l = 0; l < n_atoms_per_site[i][k]; ++l)
        {
          CG_map[site_idx].c_Ii[l] = (double) ((double)atom_weight[i][at_idx_in_mol] / (double)site_weight[i][k]);
          CG_map[site_idx].i_atm[l] = abs_at_num;
          strcpy(CG_map[site_idx].atm_name[l],atom_type[i][at_idx_in_mol]);

          ++at_idx_in_mol;
          ++abs_at_num;
        }

        strcpy(CG_map[site_idx].cg_name,site_names[i][k]);
        strcpy(CG_map[site_idx].mol_name,molecule_names[i]);
        strcpy(CG_map[site_idx].map_type,mapping_types[i][k]);

        ++site_idx;
      }
      ++res_no;
    }
  }
  return CG_map;
}

