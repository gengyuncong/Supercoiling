# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.13.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.13.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/apple/test_lmes_v9/lmes

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/apple/test_lmes_v9/lmes/build

# Include any dependencies generated for this target.
include gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/depend.make

# Include the progress variables for this target.
include gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/progress.make

# Include the compile flags for this target's objects.
include gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/flags.make

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o: gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/flags.make
gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o: ../gtest/c/gmock-1.7.0/gtest/src/gtest-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o -c /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/gtest/src/gtest-all.cc

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.i"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/gtest/src/gtest-all.cc > CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.i

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.s"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/gtest/src/gtest-all.cc -o CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.s

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o: gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/flags.make
gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o: ../gtest/c/gmock-1.7.0/src/gmock-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock_main.dir/src/gmock-all.cc.o -c /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/src/gmock-all.cc

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock_main.dir/src/gmock-all.cc.i"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/src/gmock-all.cc > CMakeFiles/gmock_main.dir/src/gmock-all.cc.i

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock_main.dir/src/gmock-all.cc.s"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/src/gmock-all.cc -o CMakeFiles/gmock_main.dir/src/gmock-all.cc.s

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/flags.make
gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: ../gtest/c/gmock-1.7.0/src/gmock_main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock_main.dir/src/gmock_main.cc.o -c /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/src/gmock_main.cc

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock_main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock_main.dir/src/gmock_main.cc.i"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/src/gmock_main.cc > CMakeFiles/gmock_main.dir/src/gmock_main.cc.i

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock_main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock_main.dir/src/gmock_main.cc.s"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0/src/gmock_main.cc -o CMakeFiles/gmock_main.dir/src/gmock_main.cc.s

# Object files for target gmock_main
gmock_main_OBJECTS = \
"CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o" \
"CMakeFiles/gmock_main.dir/src/gmock-all.cc.o" \
"CMakeFiles/gmock_main.dir/src/gmock_main.cc.o"

# External object files for target gmock_main
gmock_main_EXTERNAL_OBJECTS =

gtest/c/gmock-1.7.0/libgmock_main.a: gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/gtest/src/gtest-all.cc.o
gtest/c/gmock-1.7.0/libgmock_main.a: gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock-all.cc.o
gtest/c/gmock-1.7.0/libgmock_main.a: gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o
gtest/c/gmock-1.7.0/libgmock_main.a: gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/build.make
gtest/c/gmock-1.7.0/libgmock_main.a: gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libgmock_main.a"
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && $(CMAKE_COMMAND) -P CMakeFiles/gmock_main.dir/cmake_clean_target.cmake
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gmock_main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/build: gtest/c/gmock-1.7.0/libgmock_main.a

.PHONY : gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/build

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/clean:
	cd /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 && $(CMAKE_COMMAND) -P CMakeFiles/gmock_main.dir/cmake_clean.cmake
.PHONY : gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/clean

gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/depend:
	cd /Users/apple/test_lmes_v9/lmes/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/apple/test_lmes_v9/lmes /Users/apple/test_lmes_v9/lmes/gtest/c/gmock-1.7.0 /Users/apple/test_lmes_v9/lmes/build /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0 /Users/apple/test_lmes_v9/lmes/build/gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : gtest/c/gmock-1.7.0/CMakeFiles/gmock_main.dir/depend

