cmake_minimum_required( VERSION 2.8 )
project( "Pressure matching" )

if( COMMAND cmake_policy )
 cmake_policy(SET CMP0003 NEW)
endif( COMMAND cmake_policy)

set( CMAKE_BUILD_TYPE Debug )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_C_FLAGS_DEBUG "-g -O2" )

set ( MKLDIR ${MKLDIR} /soft/intel/x86_64/2013/composer_xe_2013/composer_xe_2013.5.192/mkl/lib/intel64 )

find_library(MKLRT mkl_rt ${MKLDIR})

set ( PLIBS ${PLIBS} m dl ${MKLRT} )

add_subdirectory(src) 
