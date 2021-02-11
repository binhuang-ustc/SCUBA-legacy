SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBS OFF)
SET(CMAKE_EXE_LINKER_FLAGS "-static")

link_directories(/usr/local/lib)
link_directories(/usr/local/lib64)
link_directories(/usr/lib)
link_directories(/usr/lib64)

link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(ABACUS-Score TestScoring.cpp)
add_executable(ABACUS-DesignSeq TestDesignSeq.cpp)

target_link_libraries (ABACUS-Score loop lapacke lapack blas gfortran quadmath boost_system boost_filesystem boost_timer pthread_nonshared m)
target_link_libraries (ABACUS-DesignSeq loop lapacke lapack blas gfortran quadmath boost_system boost_filesystem boost_timer pthread_nonshared m)