############################
# general options
############################

# Uncomment to set the build type. Can be one of Debug, Release, RelWithDebInfo or MinSizeRel
SET(CMAKE_BUILD_TYPE "Debug")

# Set the prefix dir for installation of lmes (i.e. when running "make install").
SET(CMAKE_INSTALL_PREFIX "/usr/local")

# USE_MPI can be set to yes, no, or optional, in which case MPI will be used if cmake can find it on your system.
SET(USE_MPI no)

# USE_CUDA can be set to yes, no, or optional, in which case CUDA will be used if cmake can find it on your system.
SET(USE_CUDA no)

# Set the verbosity level, as described in the lm::Print function. 10 is complete, 4 is normal-ish.
SET(VERBOSITY_LEVEL 4)

# Set to yes to turn on profiling.
SET(USE_PROFILING no)

# Uncomment and change the last part of the next line if some of your libraries/header files are on a nonstandard path.
#SET(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" "/usr/local/custom/path")

# Add any extra flags for the compiler.
SET(EXTRA_C_FLAGS "")
SET(EXTRA_CXX_FLAGS "")


############################
# cpu feature flags
############################

# USE_AVX can be set to yes, no, or optional, in which case AVX will be used if cmake can build a test AVX program with your compiler.
SET(USE_AVX optional)

# USE_FMA can be set to yes, no, or optional, in which case FMA will be used if cmake can build a test FMA program with your compiler.
SET(USE_FMA no)

# USE_SVML can be set to yes, no, or optional, in which case SVML will be used if cmake can build a test SVML program with your compiler.
SET(USE_SVML no)

############################
# compiler flags
############################

# uncomment to set a custom machine architecture flag. Defaults to "-march=native" for clang and gcc, and to "-xHost" for icc
#SET(CUSTOM_ARCH_FLAG "-march=native")

# uncomment to set any custom compiler flags for the default build
#SET(CUSTOM_C_FLAGS "-m64 -g -fPIC -Wall -fmessage-length=0")
#SET(CUSTOM_CXX_FLAGS "-m64 -g -fPIC -Wall -fmessage-length=0")

# uncomment to set any named-build specific (ie Debug, Release, etc) compiler flags.
#SET(CUSTOM_C_FLAGS_DEBUG "-m64 -g -fPIC -Wall -fmessage-length=0")
#SET(CUSTOM_CXX_FLAGS_DEBUG "-m64 -g -fPIC -Wall -fmessage-length=0")

#SET(CUSTOM_C_FLAGS_RELEASE "-m64 -O3 -DNDEBUG -fPIC -Wall -fmessage-length=0")
#SET(CUSTOM_CXX_FLAGS_RELEASE "-m64 -O3 -DNDEBUG -fPIC -Wall -fmessage-length=0")

#SET(CUSTOM_C_FLAGS_RELWITHDEBINFO "-m64 -g -fPIC -Wall -fmessage-length=0")
#SET(CUSTOM_CXX_FLAGS_RELWITHDEBINFO "-m64 -g -fPIC -Wall -fmessage-length=0")

#SET(CUSTOM_C_FLAGS_MINSIZEREL "-m64 -O3 -DNDEBUG -fPIC -Wall -fmessage-length=0")
#SET(CUSTOM_CXX_FLAGS_MINSIZEREL "-m64 -O3 -DNDEBUG -fPIC -Wall -fmessage-length=0")

############################
# paths to libraries
############################
# uncomment and set if CMake is having trouble finding one of your libraries during building.

#SET(CUDA_TOOLKIT_ROOT_DIR $CUDA_TOOLKIT_ROOT_DIR)
#SET(CUDA_SDK_ROOT_DIR $CUDA_SDK_ROOT_DIR)

#SET(ENV{HDF5_ROOT} $HDF5_ROOT)

#SET(MPI_C_COMPILER $MPI_C_COMPILER)
#SET(MPI_C_EXEC $MPI_C_EXEC)
#SET(MPI_CXX_COMPILER $MPI_CXX_COMPILER)
#SET(MPI_CXX_EXEC $MPI_CXX_EXEC)

#SET(PYTHON_LIBRARY $PYTHON_LIBRARY)
#SET(PYTHON_INCLUDE_DIR $PYTHON_INCLUDE_DIR)

#SET(PROTOBUF_ROOT $PROTOBUF_ROOT)

#SET(SBML_ROOT $SBML_ROOT)

############################
# trouble-shooting options
############################
# if you're having trouble compiling or running Lattice Microbes, you can try tweaking one of the options in this section

############################
# trouble-shooting options - GCC
############################
# if you're using GCC to compile Lattice Microbes and you're having trouble, you can try tweaking one of the options in this section

# ...and you get an error at link time that includes the phrase "cxx11:abi" somewhere, try uncommenting the following line
#add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

# ...and you're having trouble getting AVX support to work or get an error message durring compilation like "no such instruction: `vmovd %xmm0, %eax'", try uncommenting the following line
# on a mac, this will tell GCC to use the built-in clang assembler rather than the gnu as
#set(TROUBLESHOOTING_FLAGS "-Wa,-q ${TROUBLESHOOTING_FLAGS}")

############################
# deprecated library specific path settings
############################

#SET(BOOST_ROOT $BOOST_ROOT)
#SET(SWIG_EXECUTABLE $SWIG_EXECUTABLE)
