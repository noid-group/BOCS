# @author Maria Lesniewski
# @desc this is a file to quickly convert from BOCS SG units to LAMMPS
#        and to print the table you need for lammps. This version of the script is meant for Bspline output / f_force_coeff. 
#        It will extrapolate naively, using the last coefficients of the last solved spline to pad the table.
# @desc It uses the spline degree passed by the user to create table files for lammps with the same coarse-ness as the original spline


import sys 
import numpy as np
from scipy.interpolate import splrep, BSpline

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 4:
    print("Usage: <f_forces_coeff.LDGrad> <outfilename> <spline degree>")
    sys.exit(1)


spline_degree = int(sys.argv[3])
LDX, c = np.loadtxt(sys.argv[1], unpack=True)


print("SG_table_prep_from_spline.py finds the value of the recovered force function F(rho)")
print("with the same table granularity as the coefficients file (f_forces_coeff.*)")
print("This result will be stored in the output filename passed in LAMMPS table format")
print("We pad this file naively using scipy's Bspline extrapolation and save a padded LAMMPS table file too")
print("It is designed for cases where the \"Bspline\" option is selected in par.txt")
print("Check results against the f_forces.*.dat file, which reports the function on a 10x granularity when splines are used.")



print("Received args")
print("spline_degree = " + str(spline_degree))
print("f_forces_coeff file = " + str(sys.argv[1]))
print("output_file = " + str(sys.argv[2]))
print("padded version will be " + str(sys.argv[2] + "PADDED_TABLE_VERSION"))


### Interp function from coeffs on Domain Passed #####
k=spline_degree+1
dx = LDX[1] - LDX[0]
x = LDX[0] - np.floor(k/2)*dx + dx*np.arange(0, len(c) + 1 + spline_degree)


interp = BSpline(x,c, k=spline_degree)
dinterp = interp.derivative(nu=1)


SGU = interp(x)
SGF = -dinterp(x)

LDXLMP = x/1000

SGULMP = SGU/(4.184)*1e8

FLMP = SGF/(4.184)*1e11

np.savetxt(sys.argv[2], np.column_stack((LDXLMP, SGULMP, FLMP)), fmt = "%.16f")

### Extrapolative Region ### 
# Make sure we have table values down to the lowest possible LD (0) 
# Add (arbitrarily) 200 bins of extrap on the right edge)
while (x[0] > 0):
    x = np.append(x[0]-dx, x)

for i in range(1,200):
    x = np.append(x, x[-1] + dx*i)

SGU = interp(x)
SGF = -dinterp(x)

LDXLMP = x/1000

SGULMP = SGU/(4.184)*1e8

FLMP = SGF/(4.184)*1e11

np.savetxt(sys.argv[2] + "PADDED_TABLE_VERSION", np.column_stack((LDXLMP, SGULMP, FLMP)), fmt = "%.16f")

