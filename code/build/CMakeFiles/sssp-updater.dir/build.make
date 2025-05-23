# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build

# Include any dependencies generated for this target.
include CMakeFiles/sssp-updater.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/sssp-updater.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/sssp-updater.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/sssp-updater.dir/flags.make

CMakeFiles/sssp-updater.dir/src/main.cpp.o: CMakeFiles/sssp-updater.dir/flags.make
CMakeFiles/sssp-updater.dir/src/main.cpp.o: ../src/main.cpp
CMakeFiles/sssp-updater.dir/src/main.cpp.o: CMakeFiles/sssp-updater.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/sssp-updater.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/sssp-updater.dir/src/main.cpp.o -MF CMakeFiles/sssp-updater.dir/src/main.cpp.o.d -o CMakeFiles/sssp-updater.dir/src/main.cpp.o -c /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/main.cpp

CMakeFiles/sssp-updater.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sssp-updater.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/main.cpp > CMakeFiles/sssp-updater.dir/src/main.cpp.i

CMakeFiles/sssp-updater.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sssp-updater.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/main.cpp -o CMakeFiles/sssp-updater.dir/src/main.cpp.s

CMakeFiles/sssp-updater.dir/src/graph.cpp.o: CMakeFiles/sssp-updater.dir/flags.make
CMakeFiles/sssp-updater.dir/src/graph.cpp.o: ../src/graph.cpp
CMakeFiles/sssp-updater.dir/src/graph.cpp.o: CMakeFiles/sssp-updater.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/sssp-updater.dir/src/graph.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/sssp-updater.dir/src/graph.cpp.o -MF CMakeFiles/sssp-updater.dir/src/graph.cpp.o.d -o CMakeFiles/sssp-updater.dir/src/graph.cpp.o -c /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/graph.cpp

CMakeFiles/sssp-updater.dir/src/graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sssp-updater.dir/src/graph.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/graph.cpp > CMakeFiles/sssp-updater.dir/src/graph.cpp.i

CMakeFiles/sssp-updater.dir/src/graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sssp-updater.dir/src/graph.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/graph.cpp -o CMakeFiles/sssp-updater.dir/src/graph.cpp.s

CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o: CMakeFiles/sssp-updater.dir/flags.make
CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o: ../src/metis_partition.cpp
CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o: CMakeFiles/sssp-updater.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o -MF CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o.d -o CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o -c /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/metis_partition.cpp

CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/metis_partition.cpp > CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.i

CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/metis_partition.cpp -o CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.s

CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o: CMakeFiles/sssp-updater.dir/flags.make
CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o: ../src/sssp_local.cpp
CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o: CMakeFiles/sssp-updater.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o -MF CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o.d -o CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o -c /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/sssp_local.cpp

CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/sssp_local.cpp > CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.i

CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/sssp_local.cpp -o CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.s

CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o: CMakeFiles/sssp-updater.dir/flags.make
CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o: ../src/mpi_utils.cpp
CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o: CMakeFiles/sssp-updater.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o -MF CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o.d -o CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o -c /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/mpi_utils.cpp

CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/mpi_utils.cpp > CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.i

CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/mpi_utils.cpp -o CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.s

CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o: CMakeFiles/sssp-updater.dir/flags.make
CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o: ../src/dynamic_update.cpp
CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o: CMakeFiles/sssp-updater.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o -MF CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o.d -o CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o -c /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/dynamic_update.cpp

CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/dynamic_update.cpp > CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.i

CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/src/dynamic_update.cpp -o CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.s

# Object files for target sssp-updater
sssp__updater_OBJECTS = \
"CMakeFiles/sssp-updater.dir/src/main.cpp.o" \
"CMakeFiles/sssp-updater.dir/src/graph.cpp.o" \
"CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o" \
"CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o" \
"CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o" \
"CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o"

# External object files for target sssp-updater
sssp__updater_EXTERNAL_OBJECTS =

sssp-updater: CMakeFiles/sssp-updater.dir/src/main.cpp.o
sssp-updater: CMakeFiles/sssp-updater.dir/src/graph.cpp.o
sssp-updater: CMakeFiles/sssp-updater.dir/src/metis_partition.cpp.o
sssp-updater: CMakeFiles/sssp-updater.dir/src/sssp_local.cpp.o
sssp-updater: CMakeFiles/sssp-updater.dir/src/mpi_utils.cpp.o
sssp-updater: CMakeFiles/sssp-updater.dir/src/dynamic_update.cpp.o
sssp-updater: CMakeFiles/sssp-updater.dir/build.make
sssp-updater: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
sssp-updater: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
sssp-updater: CMakeFiles/sssp-updater.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable sssp-updater"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sssp-updater.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/sssp-updater.dir/build: sssp-updater
.PHONY : CMakeFiles/sssp-updater.dir/build

CMakeFiles/sssp-updater.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sssp-updater.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sssp-updater.dir/clean

CMakeFiles/sssp-updater.dir/depend:
	cd /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build /home/munjenko/PDC/PDC-SSSP-MPI-METIS/code/build/CMakeFiles/sssp-updater.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/sssp-updater.dir/depend

