# This is a script written for the 1-site LD tutorial 
# Its purpose is to integrate the f_forces file from BOCS quickly
# And translate all potentials / forces into LAMMPS ldd format
# BOCS outputs energies in kJ/mol and local densities in nm^-3
# Lammps requires energies in kcal/mol and local densities in A^-3

# Note the reason for this disclaimer is because this script naively uses np.integrate 
# and flat extrapolations in force at the begin/end tails of the passed file 
# You have to think harder about interpolation/extrapolation in the event that you want to 
# simulate a slab or that your training set samples the domain unevenly. (e.g. drho = 0.1 and you have data at 0.5 0.6 0.8)


import sys
import numpy as np

# Tell the user how to use
if (len(sys.argv) != 2):
    print("ERROR: pass in f_force.Local_Dentity.total._____.dat file as command line arg!")
    exit()

force_file = sys.argv[1]

# Read x = rho y = LD force
data = np.loadtxt(force_file)
rho = data[:,0]
f = data[:,1]

# Convert kJ/(mol nm^3) to kcal / (mol A^3)
f_lammps = f * 1000/4.184
rho_lammps = rho * 0.001

# Figure out how far to pad domain
#I'll just go to zero and add the same number of bins to the end

# Get the padding domain
drho_lammps = rho_lammps[-1] - rho_lammps[-2]
extrap_len = rho_lammps[0]/drho_lammps
low_extrap_domain = np.linspace(0,rho_lammps[0]-drho_lammps, int(extrap_len))
high_extrap_domain = np.linspace(rho_lammps[-1] + drho_lammps, rho_lammps[-1]+drho_lammps*extrap_len, int(extrap_len))

# Get the flat padded forces (repeat the last ones we found)
low_forces = f_lammps[0]*np.ones_like(low_extrap_domain)
high_forces = f_lammps[-1]*np.ones_like(high_extrap_domain)

# Stitch it all together
final_rho_lammps = np.hstack((low_extrap_domain, rho_lammps, high_extrap_domain))
final_f_lammps = np.hstack((low_forces, f_lammps, high_forces))

# Find the potential
def integrate(x,y): # Define rather than require they also have scipy installed
  '''
  This function integrates a function the same way xmgrace does
  '''
  if (len(x) != len(y)):
    print("ERROR: unable to integrate bc x and y are different lengths!")
    print("len(x): %d  len(y): %d" % (len(x),len(y)))
    exit()
  integral = [0.0]
  for i in range(0,len(x)-2):
    integral += [integral[-1] + (y[i]+y[i+1])/(2.0) * (x[i+1] - x[i])]
  integral += [integral[-1] + (y[-1] / 2.0 * (x[-1] - x[-2]))]
  return integral

u_lammps = integrate(final_rho_lammps, -1*final_f_lammps)

out = np.column_stack((final_rho_lammps, u_lammps, final_f_lammps))
outfilename = "LD_table.MEO.MEO.dat" #str(force_file + "_lammps.table")
np.savetxt(outfilename, out[1:])
print("Saved lammps compatible LD file to " + outfilename)

