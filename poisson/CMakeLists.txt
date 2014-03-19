project(poisson)
cmake_minimum_required(VERSION 2.6)

enable_language(Fortran)

option(ENABLE_OPENMP "Enable OpenMP support?" ON)
option(ENABLE_MPI    "Enable MPI support?" ON)


include_directories(${COMMON_INCLUDE_PATHS})

add_library(poisson_common fst.f pcommon.f90)

list(APPEND COMMON_LIBRARIES poisson_common)

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    list(APPEND COMMON_DEFINITIONS -DHAVE_OPENMP)
    set(COMMON_Fortran_FLAGS "${COMMON_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  endif()
endif()

if(ENABLE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    list(APPEND COMMON_DEFINITIONS -DHAVE_MPI)
    list(APPEND COMMON_Fortran_FLAGS ${MPI_Fortran_FLAGS})
    set(COMMON_Fortran_FLAGS "${COMMON_Fortran_FLAGS} ${MPI_Fortran_FLAGS}")
    list(APPEND COMMON_LIBRARIES ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    list(APPEND COMMON_INCLUDE_PATHS ${MPI_Fortran_INCLUDE_PATH})
  endif()
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${COMMON_Fortran_FLAGS} -cpp")
add_definitions(${COMMON_DEFINITIONS})
include_directories(${COMMON_INCLUDE_PATHS})
add_executable(poisson poisson.f90)
target_link_libraries(poisson poisson_common ${COMMON_LIBRARIES})
