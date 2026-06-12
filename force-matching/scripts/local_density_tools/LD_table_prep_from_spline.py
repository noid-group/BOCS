# @author Maria Lesniewski
# @desc this is a file to quickly convert from BOCS LD units to LAMMPS
#        and to print the table you need for lammps. This version of the script is meant for Bspline output / f_force_coeff. 
#        It will extrapolate naively, using the last coefficients of the last solved spline to pad the table.
# @desc It uses the spline degree passed by the user to create table files for lammps with the same coarse-ness as the original spline


import sys 
import numpy as np
from scipy.interpolate import splrep, BSpline

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 4:
    print("Usage: <f_forces_coeff.Local_Density.*.dat> <outfilename> <spline degree>")
    sys.exit(1)

print("LD_table_prep_from_spline.py finds the value of the recovered force function F(rho)")
print("with the same table granularity as the coefficients file (f_forces_coeff.*)")
print("This result will be stored in the output filename passed in LAMMPS table format")
print("We pad this file naively using scipy's Bspline extrapolation and save a padded LAMMPS table file too")
print("It is designed for cases where the \"Bspline\" option is selected in par.txt")
print("Check results against the f_forces.*.dat file, which reports the function on a 10x granularity when splines are used.")

spline_degree = int(sys.argv[3]) # This is the degree passed by the user

print("Received args")
print("spline_degree = " + str(spline_degree))
print("f_forces_coeff file = " + str(sys.argv[1]))
print("output_file = " + str(sys.argv[2]))
print("padded version will be " + str(sys.argv[2] + "PADDED_TABLE_VERSION"))


### Interp function from coeffs on Domain Passed #####
k=spline_degree + 1 # k is the order, passed into bocs
LDX, c = np.loadtxt(sys.argv[1], unpack=True) # f_forces_coeff.Local_Density.*.dat

dx = LDX[1] - LDX[0]
t = LDX[0] - np.floor(k/2)*dx + dx*np.arange(0, len(c) + 1 + spline_degree)


interp = BSpline(t,c, k=spline_degree)
uinterp = interp.antiderivative(nu=1)


LDF = interp(t)
LDU = -uinterp(t)

LDXLMP = t /1000

LDULMP = LDU/(4.184)

LDFLMP = LDF/(4.184)*1e3

np.savetxt(sys.argv[2], np.column_stack((LDXLMP, LDULMP, LDFLMP)), fmt = "%.16f")


### Extrapolative function beyond domain passed ### 
# Make sure we have table values down to the lowest possible LD (0) 
# Add (arbitrarily) 200 bins of extrap on the right edge)
while (t[0] > 0):
    t = np.append(t[0]-dx, t)

for i in range(1,200):
    t = np.append(t, t[-1] + dx*i)

LDU = -uinterp(t)
LDF = interp(t)

LDXLMP = t /1000

LDULMP = LDU/(4.184)

LDFLMP = LDF/(4.184)*1e3

np.savetxt(sys.argv[2] + "PADDED_TABLE_VERSION", np.column_stack((LDXLMP, LDULMP, LDFLMP)), fmt = "%.16f")

