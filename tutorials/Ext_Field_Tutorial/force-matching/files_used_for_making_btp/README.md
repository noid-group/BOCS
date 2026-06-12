Though not listed as a necessary step for this tutorial (see prior tutorials 1-4 for practice), a topology file is always needed to force-match in BOCS. 
This external fields tutorial already provides a BOCS .btp file for convenience, since the point of the tutorial is more to learn the cgff end of BOCS with external fields. 
The files herein are included for reference/transparency, in case one wants to see the gromacs files used to generate the tpr that BOCS's translator generated .btp from


i.e. 
# Gromacs
$ gmx_mpi grompp -p simpleABC.top -c out.gro -f general.mdp -maxwarn 2 -o 1200.tpr
$ gmx_mpi dump -s 1200.tpr > 1200.dump
# BOCS
$ translator -f 1200.dump -o 1200.btp
