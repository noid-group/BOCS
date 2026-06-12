This folder contains examples for mapping a 50:50 mixture of Heptane(23 atoms) and Toluene(15 atoms) to multiple CG sites. 
It is a subset of the extended-ensemble/pressure matching tutorial included in the main folders. 
(See N.J.H. Dunn's work on Heptane-Toluene Mixtures)


For the files in this folder: 
Heptane is mapped to 3 sites at the center of mass of composite atoms, sites include: a terminal, a middle, and another terminal site. 
(Terminal sites will be treated as equivalent in the cg force-field indicated by the top. See the 1-component example for marked up syntax)
Toluene is also mapped to 3 sites, but not all atoms in the AA topology influence the CG force-field.
site = "none" -> Atom is decimated

Gromacs files are formatted for gromacs version 4.5, but contain the same readable topology as other folders.
