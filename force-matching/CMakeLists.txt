cmake_minimum_required( VERSION 3.0 )
project( "CGFF_toolkit" )

set ( CGFF_toolkit_VERSION_MAJOR 1 )
set ( CGFF_toolkit_VERSION_MINOR 0 )

if( COMMAND cmake_policy )
  cmake_policy(SET CMP0003 NEW)
endif( COMMAND cmake_policy )

set ( CMAKE_BUILD_TYPE Debug )

# Change default output directories so not everything lives in src
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

# Set compiler flags
set ( CMAKE_C_FLAGS_DEBUG "-g -O2" )

#############################################################################################
# Set path for MKL / LAPACK/BLAS libraries
##############################################################################################
set ( MKLDIR ${MKLDIR} /soft/intel/x86_64/2013/composer_xe_2013/composer_xe_2013.5.192/mkl/lib/intel64 )

# Find the libraries we need
find_library( MKLRT mkl_rt ${MKLDIR} )
find_library( PTHREAD pthread ${MKLDIR} )

##############################################################################################
# Determine gromacs version we are compiling against
##############################################################################################

# Gromacs compilation definitions
add_definitions( -DGMX_LIB_MPI -DGMX_DOUBLE )

# Set path for gromacs 4.5.3 installation and relevant subdirectories
set ( GMXDIR ${GMXDIR} //soft/gromacs/4.6.3-ompi-double/ )
set ( GMXLIB ${GMXLIB} ${GMXDIR}/lib )

##############################################################################################


# Try first to find the gromacs 4.5 (double precision MPI) libraries
find_library( GMX_LIBRARY gmx_mpi_d ${GMXLIB} REQUIRED NO_DEFAULT_PATH)
find_library( GMX_MD_LIBRARY md_mpi_d ${GMXLIB} REQUIRED NO_DEFAULT_PATH)

if (GMX_LIBRARY AND GMX_MD_LIBRARY)
	set (GMXINC ${GMXINC} ${GMXDIR}/include/gromacs )

	set( GLIBS  ${GMX_MD_LIBRARY} ${GMX_LIBRARY} m dl )

	add_definitions( -DGMX_LIB_MPI -DGMX_DOUBLE )
	message("Detected libgmx_mpi_d and libmd_mpi_d. Attempting compilation with Gromacs 4.5")

	include_directories(${GMXINC})

	# Run a small executable that uses a gromacs library function to print
	# the version of the gromacs library
	try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
		${CMAKE_BINARY_DIR}/version_check ${CMAKE_CURRENT_SOURCE_DIR}/src/gmx_version.c
		CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${GMXINC}
		COMPILE_DEFINITIONS -DGMX=4
		LINK_LIBRARIES ${GLIBS}
		RUN_OUTPUT_VARIABLE GMX_VERSION
	)


        message("GLIBS: ${GLIBS}")
	message("GMXLIB: ${GMXLIB}")
	message("GMXINC: ${GMXINC}")
	message("GMX_LIBRARY: ${GMX_LIBRARY}")
	message("GMX_MD_LIBRARY: ${GMX_MD_LIBRARY}")



else(GMX_LIBRARY AND GMX_MD_LIBRARY) # If we didn't find the previous libraries, look for gromacs 5.x libraries.
	set (GMXINC ${GMXINC} ${GMXDIR}/include )
	set (GMXLIB ${GMXLIB} ${GMXDIR}/lib64 )

	find_library( GMX_LIBRARY gromacs_mpi_d ${GMXLIB} REQUIRED NO_DEFAULT_PATH)

	if(NOT GMX_LIBRARY)
		message(FATAL_ERROR "Gromacs library not found. Check that you are attempting to compile against v4.5 or v5.x, which themselves were compiled in double precision with openMPI. We look for the default \"_mpi_d\" suffixes on the library names.")
	endif()

	set( GLIBS ${GMX_LIBRARY} m dl )

	add_definitions( -DGMX_LIB_MPI -DGMX_DOUBLE )
	message("Detected libgromacs_mpi_d. Attempting compilation with Gromacs 5.x")

	include_directories(${GMXINC})
	include_directories(${GMXINC}/gromacs)
	include_directories(${GMXINC}/gromacs/legacyheaders)

	# Run a small executable that uses a gromacs library function to print
	# the version of the gromacs library
	try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
		${CMAKE_BINARY_DIR}/version_check ${CMAKE_CURRENT_SOURCE_DIR}/src/gmx_version.c
		CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${GMXINC}/gromacs/legacyheaders
		COMPILE_DEFINITIONS -DGMX=5
		LINK_LIBRARIES ${GLIBS}
 		COMPILE_OUTPUT_VARIABLE COMPILE 
		RUN_OUTPUT_VARIABLE GMX_VERSION
	)

	message("${COMPILE}")

endif(GMX_LIBRARY AND GMX_MD_LIBRARY)

message( "Detected GMX Version: ${GMX_VERSION}" )
if (    NOT "${GMX_VERSION}" STREQUAL "45" AND 
	NOT "${GMX_VERSION}" STREQUAL "46" AND 
	NOT "${GMX_VERSION}" STREQUAL "50" AND 
	NOT "${GMX_VERSION}" STREQUAL "51"     )
	message(FATAL_ERROR, "Incompatible gromacs version detected. v4.5 and
v5.x are the only currently supported versions.")
endif()

add_definitions(-DGMX=${GMX_VERSION})

######################################################################################################################

# Set the library combinations we need for compiling different executables
set ( TLIBS ${TLIBS} m dl )
set ( MLIBS ${MLIBS} ${GLIBS} ${MKLRT} ${PTHREAD} )

add_subdirectory( src )

#file( COPY doc DESTINATION ${CMAKE_BINARY_DIR} ) 
file( COPY src DESTINATION ${CMAKE_BINARY_DIR} )

