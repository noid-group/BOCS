/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

/**
@file read_map.c 
@authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn
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
	} else if (SITE > map.molecule_types[moltype].n_sites)
	{
		fprintf(stderr, "ERROR: Found more site types than expected under [ sitetypes ] directory.\n");
		exit(1);
	}


	if (line[i] != ';' && nbytes > i) {
	    if (line[i] == '[') {
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
	    } else if (is_blank(line) && have_directive) {
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

	    } else if (have_directive) {
		if (strcmp(directive, MOLTYPE) == 0) {

		    if (moltype >=0 && SITE > map.molecule_types[moltype].n_sites)
		    {
			fprintf(stderr, "ERROR: Found more site types than expected under [ sitetypes ] directory.\n");
		   	exit(1);
		    }

		    moltype++;
		    map.molecule_types[moltype].molname = ecalloc(strlen(line) + 1, sizeof(char));
		    sscanf(line, " %s %d ", map.molecule_types[moltype].molname, &map.molecule_types[moltype].total_atoms);
		    SITE = 0;
		    ATOM = 0;
		    TOTAL_ATOMS = 0;

		} else if (strcmp(directive, SITETYPE) == 0) {
		    n_site_chk++;
		    // parse line to find how many atoms per site there are
		    map.molecule_types[moltype].site_names[SITE] = ecalloc(strlen(line) + 1, sizeof(char));
		    map.molecule_types[moltype].map_type[SITE] = ecalloc(strlen(line) + 1, sizeof(char));

		    sscanf(line, " %s %s %d ",
			   map.molecule_types[moltype].site_names[SITE],
			   map.molecule_types[moltype].map_type[SITE],
			   &map.molecule_types[moltype].n_atoms[SITE]);

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

		    n_atm_chk += 1;

		    not_found = true;

		    /*loop over sites to find which this atom belongs to */ 
		    for (j = 0; j < map.molecule_types[moltype].n_sites; j++)
		    {

			if (strcmp (site_type, map.molecule_types[moltype].site_names[j]) == 0) {

			    ATOMS[j]++;

			    // Did we find the right number of atoms for the previous site?
			    if (j>0 && map.molecule_types[moltype].n_atoms[j-1] != ATOMS[j-1])
			    {
		                    fprintf(stderr, "ERROR: Too few atoms for site %s in %s. Expected %d, but found %d.\n", 
		               	        map.molecule_types[moltype].site_names[j-1], map.molecule_types[moltype].molname, 
		               	        map.molecule_types[moltype].n_atoms[j-1], ATOMS[j-1]);
		                   exit(1);


			    }

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

		    for (j = 0; j < map.n_types; j++)	
		    {
			if (strcmp(molname, map.molecule_types[j].molname) == 0) {
			    map.n_molecules[j] = N_MOL;
			}

		    }
		    free(molname);
		}

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
