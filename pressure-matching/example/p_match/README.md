# p_match exec
This is the executable for determining U(V) (the volume potential based pressure correction) based on reference atomistic and cg data. 
The BOCS user manual provides a tutorial example of a complete pressure matching work through / options, but below we include a quick overview.
 
All usage boils down to

$p_match config.ini

The contents of config.ini are different dependent on the iteration of pressure matching being done. 

# Usage Overview
Using the Dunn-Noid iterative pressure correction usually means going through multiple rounds of pressure matching. 
(See BOCS/tutorials/heptol_xn for a complete tutorial of each step using pressure matching. Select one composition from the AA/CG folder for a simpler tutorial.)

See the BOCS User Manual for all Options / File Syntax

## First Round 
Usually the user will have a CG force-field obtained from mapped atomistic data. When this is the case, we recommend users rerun the mapped atomistic trajectory to obtain CG pressures as a function of volume.
In this case BOCS will require that the PV data passed for the AA(original) / CG(mapped rerun) simulations are the same length.
Setting the "Iterative" directive to 0 enforces this constraint. 

## Future Rounds 
The user can come up with a pressure correction based on the differences between any two sets of PV data. Usually the atomistic PV data and the PV data generated from CG simulation..
Setting the "Iterative" directive to 1 relaxes the constraint that PV data is the same length. 


## Output:
When using the analytic Das-Anderson basis set by the "Basis_type: das_andersen" option in config.ini, psi.dat will typically contain the coefficients of the pressure correction basis in bar nm^3. 
PKG-BOCS in LAMMPS (see https://docs.lammps.org/fix_bocs.html) can simulate these coefficients and also tabulated pressure corrections. 
Typically, "real" units are employed for CG simulations in LAMMPS. Thus the conversion factor to atm A^3 is: LAMMPS_INPUT = BOCS_PSI_OUTPUT * 986.9232667

# Included in this folder:
Pressure matching files for a 5000 site (1-site water) system pressure-matched with the analytic basis.
- config.ini0 - Round0 - Pressure matching from mapped reference data
- config.ini - RoundX - Pressure matching from CG sim data 
- cg.water.lmp - LAMMPS simulation file example with "UNO" and "DOS" replacing psi_1 and psi_2 for the analytic basis
- pgo.sh - Shell script demonstrating an automated pressure-matching workflow given template files and AA data. 

