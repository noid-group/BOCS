This folder contains examples files for using the ldcalc tool
on a system with many CG types / branched topology. 

FILE: 
info.txt (Structure is described in the BOCS user manual - See LD tutorial) 

PROCEDURE: 
- Make sure that the CG trajectory is wrapped such that the CG molecules are whole. (No bond Splitting. Ideally hit the CG trajectory with -pbc whole in gromacs) 
- ldcalc -f wrapped_trajectory.trr -p branched.btp -i info.txt -self -o LJ_w_self.btj # self if using self term, do not flag if not. Trajectory name must end in .btj
- ^ The above is a fairly intensive calculation, best run inside of a compute node. The BOCS user manual describes using ldcalc on parallel threads (see LD tutorial) 
