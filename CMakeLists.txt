cmake_minimum_required(VERSION 3.2)
set(CMAKE_CXX_STANDARD 11)

project(SAMRAITEST)
#add_definitions(-DHAVE_HDF5)
#set(CMAKE_EXE_LINKER_FLAGS "-no-pie")
#set(POSITION_INDEPENDENT_CODE = True)

set(SAMRAI_DIR ~/thirdParty/SAMRAI/share/SAMRAI/cmake)
set(HDF5_DIR ~/thirdParty/HDF5/cmake)

# set(CMAKE_C_FLAGS"${CMAKE_C_FLAGS} -fPIE")
# set(CMAKE_CXX_FLAGS"${CMAKE_CXX_FLAGS} -fPIE")
set(INCLUDE_DIR 
${CMAKE_CURRENT_SOURCE_DIR}/Depend/include)
set(LIB_DIR
${CMAKE_CURRENT_SOURCE_DIR}/Depend/lib)


set(METIS_INCLUDE_DIR
${INCLUDE_DIR}/Metis)

set(VENDOR_INCLUDE_DIR
${CMAKE_CURRENT_SOURCE_DIR}/vendor/include
)

find_package(SAMRAI REQUIRED)
message("SAMRAI found at ${SAMRAI_DIR}")
message("SAMRAI_INCLUDE_DIRS = ${SAMRAI_INCLUDE_DIRS}")
message("SAMRAI_LIBRARIES = ${SAMRAI_LIBRARIES}")

find_package(HDF5 REQUIRED)
message("HDF5 found at ${HDF5_DIR}")
message("HDF5_INCLUDE_DIR = ${HDF5_INCLUDE_DIR}")
message("HDF5_LIBRARIES = ${HDF5_LIBRARIES}")

find_package(MPI REQUIRED)
message("MPI found at ${MPI_FOUND}")
message("MPI_INCLUDE_DIRS = ${MPI_INCLUDE_DIRS}")
message("MPI_CXX_INCLUDE_DIRS = ${MPI_CXX_INCLUDE_DIRS}")
message("MPI_LIBRARIES = ${MPI_LIBRARIES}")


set(METIS_LIBS
${LIB_DIR}/Metis/libmetis.so
)



add_subdirectory(src)

