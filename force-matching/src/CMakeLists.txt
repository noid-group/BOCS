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
                        LDD.c
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
                        LDD.c
                        gromacs_topology.c
			table_functions.c
			read_parameters.c
			wnoid_math.c
			io_read.c )

add_executable( translator	translator.c
                                LDD.c
				gromacs_topology.c 
				io_read.c
				wnoid_math.c
				safe_mem.c )

add_executable( ldcalc          LDD_calc.c
                                safe_mem.c
                                wnoid_math.c
                                io_read.c
                                gromacs_topology.c
                                LDD.c )

add_executable (LDPDF make_LDPDF.c
                       safe_mem.c
                       wnoid_math.c
                       io_read.c
                       gromacs_topology.c
                       LDD.c  )


IF(ENABLE_DEBUG_BUILD)
	# Then we'll enable these flags
	target_compile_definitions(cgmap PRIVATE DEBUG_BUILD=1)
	target_compile_definitions(cgff PRIVATE DEBUG_BUILD=1)
	target_compile_definitions(tables PRIVATE DEBUG_BUILD=1)
	target_compile_definitions(translator PRIVATE DEBUG_BUILD=1)
        target_compile_definitions(ldcalc PRIVATE DEBUG_BUILD=1)
	target_compile_definitions(LDPDF PRIVATE DEBUG_BUILD=1)
ELSE()
	# Then we'll turn them off
	target_compile_definitions(cgmap PRIVATE DEBUG_BUILD=0)
	target_compile_definitions(cgff PRIVATE DEBUG_BUILD=0)
        target_compile_definitions(tables PRIVATE DEBUG_BUILD=0)
        target_compile_definitions(translator PRIVATE DEBUG_BUILD=0)
        target_compile_definitions(ldcalc PRIVATE DEBUG_BUILD=0)
        target_compile_definitions(LDPDF PRIVATE DEBUG_BUILD=0)
ENDIF()


target_link_libraries( tables ${TLIBS} )

target_link_libraries( cgmap ${TLIBS} )

target_link_libraries( cgff  ${MLIBS} )

target_link_libraries( translator ${TLIBS} )

target_link_libraries( ldcalc ${MLIBS} )

target_link_libraries( LDPDF ${TLIBS} )
