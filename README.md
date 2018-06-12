# BOCS

**B**ottom-up **O**pen-source **C**oarse-graining **S**oftware (BOCS, pronounced 'box') is 
software for  parameterizing  thermodynamically  accurate  and transferable bottom-up 
coarse-grained (CG) force fields.  This software supports the calculation of force fields 
using both force and structural information from a reference all-atom (AA) simulation.  
Additionally, tools are provided for determining transferable force fields from an extended 
ensemble of AA systems, and for determining and simulating with a pressure correction that 
allows CG models to reproduce the density and compressibility of the reference AA model.

## Dependencies:
* gcc/g++ >= 4.9.2 OR icc/icpc >= 13.1
* CMake >= 3.0
* OpenMPI >= 1.9.1
* A linear algebra library, such as LAPACK/BLAS or Intel MKL

## Building and Installation
BOCS is a loosely coupled set of tools for working with coarse-grained models of molecular
systems. The components of BOCS, force-matching and pressure-matching, are each compiled separately. 
Furthermore, to run pressure-matched models, you will need to utilize the USER-BOCS package that is
distributed with the official LAMMPS distribution. Installation instructions can be found in the User 
Manual in the /docs folder. Briefly, the force-matching and pressure-matching components both utilize 
CMake. The force-matching component NO LONGER relies on an existing GROMACS installation. 
