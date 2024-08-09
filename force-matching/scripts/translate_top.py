#!/usr/bin/python
# Filename: translate_top,py

import translate_top_fn
import sys

if (len(sys.argv) < 5  ):
	print("ERROR: Too few arguments.  Accepted modes are:")
	print("'translate_top.py [gromacs_gro] [gromacs_top] [gromacs_ndx] [lammps_data]'")
	exit(1)

input_gro = open(sys.argv[1],"r")
input_top = open(sys.argv[2],"r")
input_ndx = open(sys.argv[3],"r")
output_lammps = open(sys.argv[4],"w")


lammps_top = translate_top_fn.lammps_datafile()

# Read N_atoms from gro file



# Read info from .ndx file
lammps_top.read_ndx_file(input_ndx)

lammps_top.parse_ndx_structure()

# Read info from .top file
lammps_top.read_top_file(input_top)

lammps_top.parse_top_structure()

for atom in lammps_top.atom_types:
	print("atom", atom)

for bond in lammps_top.bond_types:
	print("bond", bond[0], bond[1])

for angle in lammps_top.angle_types:
	print("angle", angle[0], angle[1], angle[2])

# Read box info from .gro file
lammps_top.read_gro_file(input_gro)

lammps_top.parse_gro_box()

lammps_top.parse_gro_atoms()

print(lammps_top.N_atoms)


# Parse bonds and angles
lammps_top.parse_bonds()
lammps_top.parse_angles()




# Print output file
lammps_top.print_output(output_lammps)



