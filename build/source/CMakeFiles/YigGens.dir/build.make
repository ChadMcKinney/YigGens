# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/octopian/Documents/source/QtDev/YigGens

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/octopian/Documents/source/QtDev/YigGens/build

# Include any dependencies generated for this target.
include source/CMakeFiles/YigGens.dir/depend.make

# Include the progress variables for this target.
include source/CMakeFiles/YigGens.dir/progress.make

# Include the compile flags for this target's objects.
include source/CMakeFiles/YigGens.dir/flags.make

source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o: source/CMakeFiles/YigGens.dir/flags.make
source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o: ../source/YigGens/YigGens.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/octopian/Documents/source/QtDev/YigGens/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o"
	cd /home/octopian/Documents/source/QtDev/YigGens/build/source && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o -c /home/octopian/Documents/source/QtDev/YigGens/source/YigGens/YigGens.cpp

source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.i"
	cd /home/octopian/Documents/source/QtDev/YigGens/build/source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/octopian/Documents/source/QtDev/YigGens/source/YigGens/YigGens.cpp > CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.i

source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.s"
	cd /home/octopian/Documents/source/QtDev/YigGens/build/source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/octopian/Documents/source/QtDev/YigGens/source/YigGens/YigGens.cpp -o CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.s

source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o.requires:
.PHONY : source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o.requires

source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o.provides: source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o.requires
	$(MAKE) -f source/CMakeFiles/YigGens.dir/build.make source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o.provides.build
.PHONY : source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o.provides

source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o.provides.build: source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o

# Object files for target YigGens
YigGens_OBJECTS = \
"CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o"

# External object files for target YigGens
YigGens_EXTERNAL_OBJECTS =

source/YigGens.so: source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o
source/YigGens.so: source/CMakeFiles/YigGens.dir/build.make
source/YigGens.so: source/CMakeFiles/YigGens.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared module YigGens.so"
	cd /home/octopian/Documents/source/QtDev/YigGens/build/source && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/YigGens.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
source/CMakeFiles/YigGens.dir/build: source/YigGens.so
.PHONY : source/CMakeFiles/YigGens.dir/build

source/CMakeFiles/YigGens.dir/requires: source/CMakeFiles/YigGens.dir/YigGens/YigGens.cpp.o.requires
.PHONY : source/CMakeFiles/YigGens.dir/requires

source/CMakeFiles/YigGens.dir/clean:
	cd /home/octopian/Documents/source/QtDev/YigGens/build/source && $(CMAKE_COMMAND) -P CMakeFiles/YigGens.dir/cmake_clean.cmake
.PHONY : source/CMakeFiles/YigGens.dir/clean

source/CMakeFiles/YigGens.dir/depend:
	cd /home/octopian/Documents/source/QtDev/YigGens/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/octopian/Documents/source/QtDev/YigGens /home/octopian/Documents/source/QtDev/YigGens/source /home/octopian/Documents/source/QtDev/YigGens/build /home/octopian/Documents/source/QtDev/YigGens/build/source /home/octopian/Documents/source/QtDev/YigGens/build/source/CMakeFiles/YigGens.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : source/CMakeFiles/YigGens.dir/depend

