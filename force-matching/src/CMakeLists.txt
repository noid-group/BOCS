# 1) Automate finding BLAS/LAPACK or MKL with provided modules
# 2) Can I automate a find gromacs, or have the path location to gmxdir
# available as a command line argument?

if ( "${GMX_VERSION}" STREQUAL "45" ) 

    add_executable( cgmap   cgmap.c 
			    read_map.c
			    transf_map.c
			    safe_mem.c
			    io_read.c
			    gmx-interface_45.c
			    wnoid_math.c )


  add_executable( cgff  cgff.c
			cgff_fn.c
 			calc_grids.c
			gmx-interface_45.c
			io_read.c
			read_parameters.c
			wnoid_math.c
			solv_lin_eqns.c
			io_output.c
			safe_mem.c
			calc_ref_potential.c )

elseif ( "${GMX_VERSION}" STREQUAL "46" ) 

    add_executable( cgmap   cgmap.c 
			    read_map.c
			    transf_map.c
			    safe_mem.c
			    io_read.c
			    gmx-interface_46.c
			    wnoid_math.c )


  add_executable( cgff  cgff.c
			cgff_fn.c
 			calc_grids.c
			gmx-interface_46.c
			io_read.c
			read_parameters.c
			wnoid_math.c
			solv_lin_eqns.c
			io_output.c
			safe_mem.c
			calc_ref_potential.c )

elseif ( "${GMX_VERSION}" STREQUAL "50" )

    add_executable( cgmap   cgmap.c 
			    read_map.c
			    transf_map.c
			    safe_mem.c
			    io_read.c
			    gmx-interface_50.c
			    wnoid_math.c )

  add_executable( cgff  cgff.c
			cgff_fn.c
 			calc_grids.c
			gmx-interface_50.c
			io_read.c
			read_parameters.c
			wnoid_math.c
			solv_lin_eqns.c
			io_output.c
			safe_mem.c
			calc_ref_potential.c )

elseif( "${GMX_VERSION}" STREQUAL "51" ) 

    add_executable( cgmap   cgmap.c 
			    read_map.c
			    transf_map.c
			    safe_mem.c
			    io_read.c
			    gmx-interface_51.c
			    wnoid_math.c )

  add_executable( cgff  cgff.c
			cgff_fn.c
 			calc_grids.c
			gmx-interface_51.c
			io_read.c
			read_parameters.c
			wnoid_math.c
			solv_lin_eqns.c
			io_output.c
			safe_mem.c
			calc_ref_potential.c )

endif ( "${GMX_VERSION}" STREQUAL "45" )


add_executable( tables  tables.c
			safe_mem.c
			table_functions.c
			read_parameters.c
			wnoid_math.c
			io_read.c)


target_link_libraries( tables ${TLIBS} )

target_link_libraries( cgmap ${GLIBS} )

target_link_libraries( cgff  ${MLIBS} )
