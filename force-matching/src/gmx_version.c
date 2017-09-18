/* ----------------------------------------------------------------------
 *    BOCS - Bottom-up Open-source Coarse-graining Software
 *    http://github.org/noid-group/bocs
 *    Dr. William Noid, wgn1@psu.edu
 *
 *    This software is distributed under the GNU General Public License v3.
 *    ------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>

#if GMX
#include "copyrite.h"
#endif


int main()
{
	int i;
	int N_digit=0;

	const char *version;
	char version_id[2];

	#ifdef GMX
	version = GromacsVersion();

	for (i=0; i<strlen(version); i++)
	{
		if ( isdigit(version[i]) && (N_digit < 2) )
		{
			version_id[N_digit] = version[i];
			fprintf(stderr, "%c", version[i]);
			N_digit++;
		}
	}

	#else
	fprintf(stderr, "Gromacs not found by version checking tool. Usually this is caught by cmake - have you edited the build files?\n");
	#endif

	return 0;
}
