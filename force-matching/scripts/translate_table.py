#!/usr/bin/python

import sys

if (len(sys.argv) < 4  ):
	print "ERROR: Too few arguments.  Accepted modes are:"
	print "'translate_table.py [gromacs_table] [table_type] [lammps_table] [interaction_name]'"
	print "[table_type] = {bond, nb, angle, dih}"
	exit(1)

gromacs_file = open(sys.argv[1],"r")
lammps_file = open(sys.argv[3],"w")
table_type = sys.argv[2]


lammps_file.write("#Converted from {0}\n".format(sys.argv[1]))

lammps_file.write("%s\n" %(sys.argv[4]))

N_lines = 0
for line in gromacs_file:
	line = line.strip()
	N_lines += 1
	data = line.split()



if ( table_type== "nb" ):
	lammps_file.write("N %d\n\n" %(N_lines-1))
else:
	lammps_file.write("N %d\n\n" %(N_lines))


lindex= 1
gromacs_file.seek(0)
for line in gromacs_file:
	line = line.strip()
	data = line.split()

	if (table_type == "nb"):
		r = float( data[0].strip() )
		u = float( data[5].strip() )
		f = float( data[6].strip() )

	else:
		r = float( data[0].strip() )
		u = float( data[1].strip() )
		f = float( data[2].strip() )


	if (table_type == "nb" or table_type == "bond"):
		r = r * 10.

	u = u / 4.184

	f = f / ( 4.184 * 10. ) 

	print r

	if ( !(table_type=="nb" and r==0)  ):
		lammps_file.write("%d     %f  %f  %f\n" %(lindex, r, u, f) )

		lindex += 1
