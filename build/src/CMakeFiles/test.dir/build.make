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
CMAKE_SOURCE_DIR = /home/enrico/Progetti/Nut_SQPnP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/enrico/Progetti/Nut_SQPnP/build

# Include any dependencies generated for this target.
include src/CMakeFiles/test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/test.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/test.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/test.dir/flags.make

src/CMakeFiles/test.dir/test.cpp.o: src/CMakeFiles/test.dir/flags.make
src/CMakeFiles/test.dir/test.cpp.o: /home/enrico/Progetti/Nut_SQPnP/src/test.cpp
src/CMakeFiles/test.dir/test.cpp.o: src/CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrico/Progetti/Nut_SQPnP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/test.dir/test.cpp.o"
	cd /home/enrico/Progetti/Nut_SQPnP/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/test.dir/test.cpp.o -MF CMakeFiles/test.dir/test.cpp.o.d -o CMakeFiles/test.dir/test.cpp.o -c /home/enrico/Progetti/Nut_SQPnP/src/test.cpp

src/CMakeFiles/test.dir/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/test.cpp.i"
	cd /home/enrico/Progetti/Nut_SQPnP/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enrico/Progetti/Nut_SQPnP/src/test.cpp > CMakeFiles/test.dir/test.cpp.i

src/CMakeFiles/test.dir/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/test.cpp.s"
	cd /home/enrico/Progetti/Nut_SQPnP/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enrico/Progetti/Nut_SQPnP/src/test.cpp -o CMakeFiles/test.dir/test.cpp.s

# Object files for target test
test_OBJECTS = \
"CMakeFiles/test.dir/test.cpp.o"

# External object files for target test
test_EXTERNAL_OBJECTS =

src/test: src/CMakeFiles/test.dir/test.cpp.o
src/test: src/CMakeFiles/test.dir/build.make
src/test: src/CMakeFiles/test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enrico/Progetti/Nut_SQPnP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test"
	cd /home/enrico/Progetti/Nut_SQPnP/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/test.dir/build: src/test
.PHONY : src/CMakeFiles/test.dir/build

src/CMakeFiles/test.dir/clean:
	cd /home/enrico/Progetti/Nut_SQPnP/build/src && $(CMAKE_COMMAND) -P CMakeFiles/test.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/test.dir/clean

src/CMakeFiles/test.dir/depend:
	cd /home/enrico/Progetti/Nut_SQPnP/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enrico/Progetti/Nut_SQPnP /home/enrico/Progetti/Nut_SQPnP/src /home/enrico/Progetti/Nut_SQPnP/build /home/enrico/Progetti/Nut_SQPnP/build/src /home/enrico/Progetti/Nut_SQPnP/build/src/CMakeFiles/test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/test.dir/depend

