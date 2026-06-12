# Key Files for cgmap 
map.\*.top


# Typical CL Usage
$ cgmap -f wrapped_trajectory.trr -p map.\*.top -o mapped_trajectory.trr -c mapped_first_frame.gro
$ cgmap -f wrapped_trajectory.lmp -p map.\*.top -c mapped_first_frame.gro

# In these folders
- Mapping files for different complexities of system components / sites (map.<systemname>.top)
- Reference AA / CG topology files in gromacs format (AA.top, CG.top) (For comparing before/after cgmap, these are usually made by the user during AA/CG simulations) 


Please note, for systems including intermolecular interactions AA and CG trajectories should be wrapped before / after mapping to avoid splitting atomistic / CG bonds across pbc in the usual force-matching procedure.

See the cgmap section of the user manual for a detailed description of file syntax and complete options
