cmake_minimum_required (VERSION 2.6)
project (migrated)

set (Dataio_VERSION_MAJOR 0)
set (Dataio_VERSION_MINOR 0)

include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/")
SET(ENV{LD_LIBRARY_PATH} "${PROJECT_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}")
SET(Boost_USE_STATIC_LIBS ON) 
find_package(Boost 1.66.0 REQUIRED COMPONENTS system filesystem timer)
add_definitions(-D_LAPACK_COMPLEX_STRUCTURE -D_HAVE_LAPACK_CONFIG_H -D_BOOST_NO_CXX11_SCOPED_ENUMS -D_Boost_USE_STATIC_LIBS=ON)
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O4 -static -fopenmp -std=c++11 -Wl,--no-as-needed -L${PROJECT_BINARY_DIR}/lib/libloop.a -L/usr/local/lib64/libgomp.a -L/usr/lib64/liblapacke.a -L/usr/lib64/liblapack.a -L/usr/lib64/libblas.a -L/usr/local/lib/libboost_system.a -L/usr/local/lib/libboost_filesystem.a -L/usr/local/lib/libboost_timer.a -L/usr/lib64/libpthread_nonshared.a -L/usr/lib64/libc.a -L/usr/lib64/libm.a")
SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/test" )
SET( LIBRARY_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/lib" )
add_library(loop STATIC IMPORTED GLOBAL)
set(C_LIB ${CMAKE_BINARY_DIR}/lib/libloop.a)
set_target_properties(loop
        PROPERTIES
        IMPORTED_LOCATION ${C_LIB}
        )
add_dependencies(loop combinelib)
add_custom_target(combinelib 
        rm -f libloop.a
	COMMAND ar rcT libloop.a libnsppdbstatistics.a libbackbone.a libnsploopclosure.a libnspdataio.a libnspgeometry.a libnspdstl.a libsd.a libnspproteinrep.a libfullsite.a libnspdesignseq.a
	COMMAND ranlib libloop.a
	DEPENDS nsppdbstatistics backbone nsploopclosure nspdataio nspgeometry nspdstl nspproteinrep sd fullsite nspdesignseq
        WORKING_DIRECTORY ${LIBRARY_OUTPUT_PATH})
add_subdirectory (dstl)
add_subdirectory (dataio)
add_subdirectory (geometry)
add_subdirectory (proteinrep)
add_subdirectory (pdbstatistics)
add_subdirectory (backbone)
add_subdirectory (loopclosure)
add_subdirectory (sd)
add_subdirectory (fullsite)
add_subdirectory (designseq)
