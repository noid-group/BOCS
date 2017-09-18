

add_executable( p_match pressure_matching.c
			distro.c
			io.c
			basis.c
			pressure_matching_fn.c
			solv_lin_eqns.c
			safe_mem.c )


add_executable( plot_basis 	plot_basis.c
				distro.c
				io.c
				basis.c
				pressure_matching_fn.c
				solv_lin_eqns.c
				safe_mem.c )


target_link_libraries( p_match ${PLIBS} )
target_link_libraries( plot_basis ${PLIBS} )

