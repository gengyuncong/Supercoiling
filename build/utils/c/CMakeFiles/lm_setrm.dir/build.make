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
include utils/c/CMakeFiles/lm_setrm.dir/depend.make

# Include the progress variables for this target.
include utils/c/CMakeFiles/lm_setrm.dir/progress.make

# Include the compile flags for this target's objects.
include utils/c/CMakeFiles/lm_setrm.dir/flags.make

utils/c/CMakeFiles/lm_setrm.dir/lm_setrm.cpp.o: utils/c/CMakeFiles/lm_setrm.dir/flags.make
utils/c/CMakeFiles/lm_setrm.dir/lm_setrm.cpp.o: ../utils/c/lm_setrm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object utils/c/CMakeFiles/lm_setrm.dir/lm_setrm.cpp.o"
	cd /Users/apple/test_lmes_v9/lmes/build/utils/c && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lm_setrm.dir/lm_setrm.cpp.o -c /Users/apple/test_lmes_v9/lmes/utils/c/lm_setrm.cpp

utils/c/CMakeFiles/lm_setrm.dir/lm_setrm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lm_setrm.dir/lm_setrm.cpp.i"
	cd /Users/apple/test_lmes_v9/lmes/build/utils/c && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/test_lmes_v9/lmes/utils/c/lm_setrm.cpp > CMakeFiles/lm_setrm.dir/lm_setrm.cpp.i

utils/c/CMakeFiles/lm_setrm.dir/lm_setrm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lm_setrm.dir/lm_setrm.cpp.s"
	cd /Users/apple/test_lmes_v9/lmes/build/utils/c && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/test_lmes_v9/lmes/utils/c/lm_setrm.cpp -o CMakeFiles/lm_setrm.dir/lm_setrm.cpp.s

utils/c/CMakeFiles/lm_setrm.dir/util.cpp.o: utils/c/CMakeFiles/lm_setrm.dir/flags.make
utils/c/CMakeFiles/lm_setrm.dir/util.cpp.o: ../utils/c/util.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object utils/c/CMakeFiles/lm_setrm.dir/util.cpp.o"
	cd /Users/apple/test_lmes_v9/lmes/build/utils/c && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lm_setrm.dir/util.cpp.o -c /Users/apple/test_lmes_v9/lmes/utils/c/util.cpp

utils/c/CMakeFiles/lm_setrm.dir/util.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lm_setrm.dir/util.cpp.i"
	cd /Users/apple/test_lmes_v9/lmes/build/utils/c && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/test_lmes_v9/lmes/utils/c/util.cpp > CMakeFiles/lm_setrm.dir/util.cpp.i

utils/c/CMakeFiles/lm_setrm.dir/util.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lm_setrm.dir/util.cpp.s"
	cd /Users/apple/test_lmes_v9/lmes/build/utils/c && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/test_lmes_v9/lmes/utils/c/util.cpp -o CMakeFiles/lm_setrm.dir/util.cpp.s

# Object files for target lm_setrm
lm_setrm_OBJECTS = \
"CMakeFiles/lm_setrm.dir/lm_setrm.cpp.o" \
"CMakeFiles/lm_setrm.dir/util.cpp.o"

# External object files for target lm_setrm
lm_setrm_EXTERNAL_OBJECTS =

utils/c/lm_setrm: utils/c/CMakeFiles/lm_setrm.dir/lm_setrm.cpp.o
utils/c/lm_setrm: utils/c/CMakeFiles/lm_setrm.dir/util.cpp.o
utils/c/lm_setrm: utils/c/CMakeFiles/lm_setrm.dir/build.make
utils/c/lm_setrm: src/c/liblm_c_lib.a
utils/c/lm_setrm: /usr/local/lib/libprotobuf.dylib
utils/c/lm_setrm: /usr/local/lib/libhdf5.dylib
utils/c/lm_setrm: /usr/local/lib/libsz.dylib
utils/c/lm_setrm: /usr/lib/libz.dylib
utils/c/lm_setrm: /usr/lib/libdl.dylib
utils/c/lm_setrm: /usr/lib/libm.dylib
utils/c/lm_setrm: /usr/local/lib/libhdf5_hl.dylib
utils/c/lm_setrm: /usr/lib/libz.dylib
utils/c/lm_setrm: /usr/lib/libdl.dylib
utils/c/lm_setrm: /usr/lib/libm.dylib
utils/c/lm_setrm: /usr/local/lib/libhdf5_hl.dylib
utils/c/lm_setrm: utils/c/CMakeFiles/lm_setrm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable lm_setrm"
	cd /Users/apple/test_lmes_v9/lmes/build/utils/c && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lm_setrm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
utils/c/CMakeFiles/lm_setrm.dir/build: utils/c/lm_setrm

.PHONY : utils/c/CMakeFiles/lm_setrm.dir/build

utils/c/CMakeFiles/lm_setrm.dir/clean:
	cd /Users/apple/test_lmes_v9/lmes/build/utils/c && $(CMAKE_COMMAND) -P CMakeFiles/lm_setrm.dir/cmake_clean.cmake
.PHONY : utils/c/CMakeFiles/lm_setrm.dir/clean

utils/c/CMakeFiles/lm_setrm.dir/depend:
	cd /Users/apple/test_lmes_v9/lmes/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/apple/test_lmes_v9/lmes /Users/apple/test_lmes_v9/lmes/utils/c /Users/apple/test_lmes_v9/lmes/build /Users/apple/test_lmes_v9/lmes/build/utils/c /Users/apple/test_lmes_v9/lmes/build/utils/c/CMakeFiles/lm_setrm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : utils/c/CMakeFiles/lm_setrm.dir/depend

