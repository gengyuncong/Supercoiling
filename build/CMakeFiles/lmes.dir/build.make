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
include CMakeFiles/lmes.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lmes.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lmes.dir/flags.make

CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.o: CMakeFiles/lmes.dir/flags.make
CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.o: ../src/c/lm/main/Main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.o -c /Users/apple/test_lmes_v9/lmes/src/c/lm/main/Main.cpp

CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/test_lmes_v9/lmes/src/c/lm/main/Main.cpp > CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.i

CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/test_lmes_v9/lmes/src/c/lm/main/Main.cpp -o CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.s

CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.o: CMakeFiles/lmes.dir/flags.make
CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.o: ../src/c/lm/main/MainArgs.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.o -c /Users/apple/test_lmes_v9/lmes/src/c/lm/main/MainArgs.cpp

CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/test_lmes_v9/lmes/src/c/lm/main/MainArgs.cpp > CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.i

CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/test_lmes_v9/lmes/src/c/lm/main/MainArgs.cpp -o CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.s

# Object files for target lmes
lmes_OBJECTS = \
"CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.o" \
"CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.o"

# External object files for target lmes
lmes_EXTERNAL_OBJECTS =

lmes: CMakeFiles/lmes.dir/src/c/lm/main/Main.cpp.o
lmes: CMakeFiles/lmes.dir/src/c/lm/main/MainArgs.cpp.o
lmes: CMakeFiles/lmes.dir/build.make
lmes: src/c/liblm_c_lib.a
lmes: /usr/local/lib/libprotobuf.dylib
lmes: /usr/local/lib/libhdf5.dylib
lmes: /usr/local/lib/libsz.dylib
lmes: /usr/lib/libz.dylib
lmes: /usr/lib/libdl.dylib
lmes: /usr/lib/libm.dylib
lmes: /usr/local/lib/libhdf5_hl.dylib
lmes: /usr/lib/libz.dylib
lmes: /usr/lib/libdl.dylib
lmes: /usr/lib/libm.dylib
lmes: /usr/local/lib/libhdf5_hl.dylib
lmes: CMakeFiles/lmes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/apple/test_lmes_v9/lmes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable lmes"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lmes.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lmes.dir/build: lmes

.PHONY : CMakeFiles/lmes.dir/build

CMakeFiles/lmes.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lmes.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lmes.dir/clean

CMakeFiles/lmes.dir/depend:
	cd /Users/apple/test_lmes_v9/lmes/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/apple/test_lmes_v9/lmes /Users/apple/test_lmes_v9/lmes /Users/apple/test_lmes_v9/lmes/build /Users/apple/test_lmes_v9/lmes/build /Users/apple/test_lmes_v9/lmes/build/CMakeFiles/lmes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lmes.dir/depend

