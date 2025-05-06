cmake_minimum_required(VERSION 3.10)
project(sssp_mpi_metis)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH} include)

add_executable(sssp-updater
    src/main.cpp
    src/graph.cpp
    src/metis_partition.cpp
    src/sssp_local.cpp
    src/mpi_utils.cpp
    src/dynamic_update.cpp 
)

target_link_libraries(sssp-updater ${MPI_CXX_LIBRARIES} metis)
set_target_properties(sssp-updater PROPERTIES CXX_STANDARD 17)
