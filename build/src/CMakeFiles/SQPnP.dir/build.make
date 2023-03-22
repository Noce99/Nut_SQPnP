# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_SOURCE_DIR = /home/enrico/Progetti/NutPnP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/enrico/Progetti/NutPnP/build

# Include any dependencies generated for this target.
include src/CMakeFiles/SQPnP.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/SQPnP.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/SQPnP.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/SQPnP.dir/flags.make

src/CMakeFiles/SQPnP.dir/pnp_solver.cpp.o: src/CMakeFiles/SQPnP.dir/flags.make
src/CMakeFiles/SQPnP.dir/pnp_solver.cpp.o: /home/enrico/Progetti/NutPnP/src/pnp_solver.cpp
src/CMakeFiles/SQPnP.dir/pnp_solver.cpp.o: src/CMakeFiles/SQPnP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrico/Progetti/NutPnP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/SQPnP.dir/pnp_solver.cpp.o"
	cd /home/enrico/Progetti/NutPnP/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/SQPnP.dir/pnp_solver.cpp.o -MF CMakeFiles/SQPnP.dir/pnp_solver.cpp.o.d -o CMakeFiles/SQPnP.dir/pnp_solver.cpp.o -c /home/enrico/Progetti/NutPnP/src/pnp_solver.cpp

src/CMakeFiles/SQPnP.dir/pnp_solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SQPnP.dir/pnp_solver.cpp.i"
	cd /home/enrico/Progetti/NutPnP/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrico/Progetti/NutPnP/src/pnp_solver.cpp > CMakeFiles/SQPnP.dir/pnp_solver.cpp.i

src/CMakeFiles/SQPnP.dir/pnp_solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SQPnP.dir/pnp_solver.cpp.s"
	cd /home/enrico/Progetti/NutPnP/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrico/Progetti/NutPnP/src/pnp_solver.cpp -o CMakeFiles/SQPnP.dir/pnp_solver.cpp.s

# Object files for target SQPnP
SQPnP_OBJECTS = \
"CMakeFiles/SQPnP.dir/pnp_solver.cpp.o"

# External object files for target SQPnP
SQPnP_EXTERNAL_OBJECTS =

src/libSQPnP.a: src/CMakeFiles/SQPnP.dir/pnp_solver.cpp.o
src/libSQPnP.a: src/CMakeFiles/SQPnP.dir/build.make
src/libSQPnP.a: src/CMakeFiles/SQPnP.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enrico/Progetti/NutPnP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libSQPnP.a"
	cd /home/enrico/Progetti/NutPnP/build/src && $(CMAKE_COMMAND) -P CMakeFiles/SQPnP.dir/cmake_clean_target.cmake
	cd /home/enrico/Progetti/NutPnP/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SQPnP.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/SQPnP.dir/build: src/libSQPnP.a
.PHONY : src/CMakeFiles/SQPnP.dir/build

src/CMakeFiles/SQPnP.dir/clean:
	cd /home/enrico/Progetti/NutPnP/build/src && $(CMAKE_COMMAND) -P CMakeFiles/SQPnP.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/SQPnP.dir/clean

src/CMakeFiles/SQPnP.dir/depend:
	cd /home/enrico/Progetti/NutPnP/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enrico/Progetti/NutPnP /home/enrico/Progetti/NutPnP/src /home/enrico/Progetti/NutPnP/build /home/enrico/Progetti/NutPnP/build/src /home/enrico/Progetti/NutPnP/build/src/CMakeFiles/SQPnP.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/SQPnP.dir/depend

