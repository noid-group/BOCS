# 1) Automate finding BLAS/LAPACK or MKL with provided modules
# 2) Can I automate a find gromacs, or have the path location to gmxdir
# available as a command line argument?

add_executable( cgmap   cgmap.c 
			read_map.c
			transf_map.c
			safe_mem.c
			io_read.c
			gromacs_topology.c
			wnoid_math.c )

add_executable( cgff    cgff.c
			cgff_fn.c
 			calc_grids.c
			gromacs_topology.c
			io_read.c
			read_parameters.c
			wnoid_math.c
			solv_lin_eqns.c
			io_output.c
			safe_mem.c
			calc_ref_potential.c )



add_executable( tables  tables.c
			safe_mem.c
			table_functions.c
			read_parameters.c
			wnoid_math.c
			io_read.c )

add_executable( translator	translator.c
				gromacs_topology.c 
				io_read.c
				wnoid_math.c
				safe_mem.c )

target_link_libraries( tables ${TLIBS} )

target_link_libraries( cgmap ${TLIBS} )

target_link_libraries( cgff  ${MLIBS} )

target_link_libraries( translator ${TLIBS} )
