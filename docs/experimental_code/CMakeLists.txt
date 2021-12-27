# experimental_code c directory

include(SourceHelper)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

################
#TARGET: experimental_code
################
add_executable_standalone(TARGET MapMacro SRCS MapMacro.cpp)
add_executable_standalone(TARGET ResolveMacro SRCS ResolveMacro.cpp)
add_executable_standalone(TARGET Sparse SRCS Sparse.cpp)
add_executable_standalone(TARGET Sparse11 SRCS Sparse.cpp BUILD_WITH_CPP11)
add_executable_standalone(TARGET standaloneLimitChecking SRCS standaloneLimitChecking.cpp)
add_executable_standalone(TARGET standaloneOptionParsing SRCS standaloneOptionParsing.cpp)
add_executable_standalone(TARGET TypeTraits SRCS TypeTraits.cpp)
add_executable_standalone(TARGET TypeTraits11 SRCS TypeTraits.cpp BUILD_WITH_CPP11)
