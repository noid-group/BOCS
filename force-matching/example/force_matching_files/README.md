This folder contains examples of files used for the cgff application included in BOCS: 

# KEY cgff FILES
- inp.txt - Tells BOCS where the trajectory data is. In practice, the name of this file is specified in par.txt
- par.txt - Tells BOCS which kinds of spline parameters / order parameters to include in the force-field recovered. (Must be named par.txt)
- cg.btp - Tells BOCS the kind of molecular topology in the system / how many / the expected read order from the trajectory 


# cgff and Many Body (LD/SG) interactions
For Many Body Interactions, the typical workflow is:
- Generate an atomistic trajectory 
- Map that trajectory (cgmap)
- Calculate Local Densities/Gradients of the Local Density (ldcalc)
- Force-Match the resulting trajectory (cgff) 
Thus for examples including many body interactions, we have included sample files for the requisite ldcalc command (info.txt) 

# This Folder
Contains variations on par.txt and cg.btp for more/less complex topologies / force-fields (Multi-site, Mixtures, Many Body Interactions)
Contains variations on inp.txt / par.txt for extended ensemble sampling (multiple topologies/trajectories) 
