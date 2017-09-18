#!/usr/bin/python
# Filename: translate_top_fn.py

class lammps_datafile:
	def __init__(self):
		self.N_atoms = 0
		self.N_bonds = 0
		self.N_angles = 0
		self.N_dihedrals = 0
		self.N_impropers = 0
		self.atom_types = []
		self.bond_types = []
		self.angle_types = []
		self.dihedral_types = []
		self.improper_types = []
		self.box_dim = []

		self.gro_data = []
		self.ndx_data = []
		self.top_data = []

		self.ndx_bonds = []
		self.ndx_angles = []


		self.structure = []
		self.atoms = []
		self.bonds = []
		self.angles = []


	def read_gro_file(self, input_gro):
		for line in input_gro:
			line = line.strip()
			data = line.split()
			self.gro_data.append(data)
		self.N_atoms = int(self.gro_data[1][0])

	def parse_gro_box(self):
		box_vec = self.gro_data[-1]
		for coord in box_vec:
			self.box_dim.append( [0, (float(coord) * 10.)] )

	def parse_gro_atoms(self):
		i = 2
		while len(self.atoms) < self.N_atoms:
			if self.gro_data[i][1] == "CTA":
				atom_id = "CT"
			elif self.gro_data[i][1] == "CTB":
				atom_id = "CT"
			elif self.gro_data[i][1] == "CBA":
				atom_id = "CB"
			elif self.gro_data[i][1] == "CBB":
				atom_id = "CB"
			else:
				atom_id = self.gro_data[i][1]



			atom_id_index = -2
			for atom in self.atom_types:
				if atom_id == atom[0]:
					atom_id_index = self.atom_types.index(atom)
					

			self.atoms.append([ i-1, self.gro_data[i][0], ( atom_id_index + 1 ), self.gro_data[i][3], self.gro_data[i][4], self.gro_data[i][5], atom_id ] )
			i += 1


	def read_ndx_file(self, input_ndx):
		for line in input_ndx:
			line = line.strip()
			data = line.split()
			self.ndx_data.append(data)


	def read_top_file(self, input_top):
		for line in input_top:
			line = line.strip()
			self.top_data.append(line)

	def __parse_top(self, target_array, start_key):
		start_index = -1
		keep_going = 0
		for line in self.top_data:
			if len(line) >= 2:
				if start_key in line:
					start_index = self.top_data.index(line)
					keep_going = 1
					#print line, start_index
				if start_index != -1 and start_key not in line and line[0][0] != ";" and keep_going:
					if line[0][0] != "[":
						#print line
						line = line.split()
						target_array.append([line[0], line[1], line[2]])
					else:
						keep_going = 0


	def parse_top_atoms(self):
		start_key = "atomtypes"
		target_array = self.atom_types
		self.__parse_top(target_array, start_key)
				

	def parse_top_bonds(self):
		start_key = "bondtypes"
		target_array = self.bond_types
		self.__parse_top(target_array, start_key)


	def parse_top_angles(self):
		start_key = "angletypes"
		target_array = self.angle_types
		self.__parse_top(target_array, start_key)

	def parse_top_structure(self):
		self.parse_top_angles()
		self.parse_top_bonds()
		self.parse_top_atoms()

	def parse_ndx_structure(self):
		ndx_index = 0
		header_list = []
		for entry in self.ndx_data:
			if (len(entry) > 0):
				if entry[0] == "[":
					print ndx_index, self.ndx_data.index(entry), entry[0], entry[1], entry[2]
					header_list.append(self.ndx_data.index(entry))
					ndx_index += 1

		i = 0
		while i < len(header_list):
			header = header_list[i]
			header_name = self.ndx_data[header][0] + self.ndx_data[header][1] + self.ndx_data[header][2]
			self.structure.append([ [header_name], [] ])

			if i+1 < len(header_list):
				next_header = header_list[i+1]
			else:
				next_header = len(self.ndx_data)

			j = header + 1
			while j < next_header:
				for atom in self.ndx_data[j]:
					atom = int( atom.strip() )
					self.structure[i][1].append( atom )
				j += 1
			i += 1

		for entry in self.structure:
			if "bond" in entry[0][0] or "Bond" in entry[0][0]:
				self.N_bonds += len(entry[1])/2
				self.ndx_bonds.append(entry[1])
				print entry[0][0], "is a bond entry with", len(entry[1])/2, "members."
			if "angle" in entry[0][0] or "Angle" in entry[0][0]:
				self.N_angles += len(entry[1])/3
				self.ndx_angles.append(entry[1])
				print entry[0][0], "is an angle entry with", len(entry[1])/3, "members."


	def parse_bonds(self):
		bond_id = 1
		for bond_type in self.ndx_bonds:
			i = 0
			#print self.atoms[bond_type[0]-1][2], self.atoms[bond_type[1]-1][2]
			while i < len(bond_type) - 1:
				# FIXME - write a function that can get the bond_id# from the atoms involved - will need for TOL and for mixed systems
				self.bonds.append( [ bond_id, self.atoms[bond_type[i]-1][0], self.atoms[bond_type[i+1]-1][0]] )
				i += 2
			bond_id += 1

	def parse_angles(self):
		i = 0;
		for angle_type in self.ndx_angles:
			while i < len(angle_type) - 1:
				angle_id = 1 # FIXME - write a function that can get the angle_id# from the atoms involved - will need for TOL and for mixed systems
				self.angles.append( [ angle_id, self.atoms[angle_type[i]-1][0], self.atoms[angle_type[i+1]-1][0], self.atoms[angle_type[i+2]-1][0] ] )
				i += 3



	def print_output(self, outfile):
		outfile.write("Topology of CG heptane translated from gromacs input\n\n")
		outfile.write(" %d atoms\n" %self.N_atoms)
		outfile.write(" %d bonds\n" %self.N_bonds)
		outfile.write(" %d angles\n" %self.N_angles)
		outfile.write(" %d dihedrals\n" %self.N_dihedrals)
		outfile.write(" %d impropers\n\n" %self.N_impropers)

		outfile.write(" %d atom types\n" %(len(self.atom_types)))
		outfile.write(" %d bond types\n" %(len(self.bond_types)))
		outfile.write(" %d angle types\n" %(len(self.angle_types)))
		outfile.write(" %d dihedral types\n" %(len(self.dihedral_types)))
		outfile.write(" %d improper types\n\n" %(len(self.improper_types)))

		outfile.write(" %f %f xlo xhi\n" %(self.box_dim[0][0], self.box_dim[0][1]))
		outfile.write(" %f %f ylo yhi\n" %(self.box_dim[1][0], self.box_dim[1][1]))
		outfile.write(" %f %f zlo zhi\n\n\n\n" %(self.box_dim[2][0], self.box_dim[2][1]))

		#outfile.write("Pair Coeffs\n\n")

		N_types = len(self.atom_types)
#		i = 1
#		while i <= N_types:
#			j = i
#			while j <= N_types:
#				outfile.write(" %d %d   # %s %s\n" %(i, j, self.atom_types[i-1][0], self.atom_types[j-1][0]) )
#				j += 1
#			i += 1

		outfile.write("\n\n\nBond Coeffs\n\n")
		N_bonds = len(self.bond_types)
		i = 1
		while i <= N_bonds:
			outfile.write( " %d # %s-%s\n " %(i, self.bond_types[i-1][0], self.bond_types[i-1][1]) )
			i += 1


		outfile.write("\n\n\nAngle Coeffs\n\n")
		N_angles = len(self.angle_types)
		i = 1
		while i <= N_angles:
			outfile.write(" %d # %s-%s-%s\n" %(i, self.angle_types[i-1][0], self.angle_types[i-1][1], self.angle_types[i-1][2] ))
			i += 1


		outfile.write("\n\n\nMasses\n\n")
		i = 1
		while i <= N_types:
			outfile.write(" %d %s # %s\n" %(i, self.atom_types[i-1][1], self.atom_types[i-1][0]))
			i += 1

		outfile.write("\n\n\nAtoms\n\n")
		i = 1		
		for atom in self.atoms:
			outfile.write(" %d 1 %d 0.0 %f %f %f # %s\n" %(i, atom[2], float(atom[3])*10, float(atom[4])*10, float(atom[5])*10, atom[6]))
			i += 1

		outfile.write("\n\n\nBonds\n\n")
		i = 1
		for bond in self.bonds:
			outfile.write(" %d %d %d %d\n" %(i, bond[0], bond[1], bond[2]))
			i += 1

		outfile.write("\n\n\nAngles\n\n")
		i = 1
		for angle in self.angles:
			outfile.write(" %d %d %d %d %d\n" %(i, angle[0], angle[1], angle[2], angle[3]))
			i += 1







		


version = '0.1'
