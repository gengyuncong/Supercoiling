# a set of functions to test if the fully configured compiler supports various features
INCLUDE(CheckCXXSourceRuns)

########################
# test functions
########################
# All the functions named test<some-feature> take two arguments, varName and required
# The result of the test will be stored at the parent scope in a variable named varName
# If required is yes, set varName to true.
# If required is no, set varName to false.
# If required is optional, directly check if the compiler can build avx code and set varName to the test result

# test if the compiler supports avx
function(testAVX varName required)
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/TestExecutionSnippets/avxTestSnippet.cpp avxTestSnippet)
    featureTest(${varName} ${required} AVX "${avxTestSnippet}" ${ARGV2})
    set(${varName} ${${varName}} PARENT_SCOPE)
endfunction()

# test if the compiler supports fma
function(testFMA varName required)
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/TestExecutionSnippets/fmaTestSnippet.cpp fmaTestSnippet)
    featureTest(${varName} ${required} FMA "${fmaTestSnippet}" ${ARGV2})
    set(${varName} ${${varName}} PARENT_SCOPE)
endfunction()

# test if the compiler supports svml. currently this is just a stub that tests if the compiler is ICC
function(testSVML varName required)
    set(featureName SVML)

    if(required STREQUAL yes)
        message(STATUS "${featureName} support required by user option, setting ${varName} to true")
        set(${varName} true)

    elseif(required STREQUAL no)
        message(STATUS "${featureName} support disabled by user option, setting ${varName} to false")
        set(${varName} false)

    else(required STREQUAL yes)
        if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
            set(${varName} true)
        else(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
            set(${varName} false)
        endif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")

        if(${varName})
            message(STATUS "${featureName} is supported by the compiler, setting ${varName} to true")
        else(${varName})
            message(STATUS "${featureName} is not supported by the compiler, setting ${varName} to false")
        endif(${varName})

    endif(required STREQUAL yes)
    set(${varName} ${${varName}} PARENT_SCOPE)
endfunction()

#################################
# the remaining functions are meant to be internal to the module
#################################
function(featureTest varName required featureName testSnippet)
    if(required STREQUAL yes OR required STREQUAL Yes OR required STREQUAL true OR required STREQUAL True)
        message(STATUS "${featureName} support required by user option, setting ${varName} to true")
        set(${varName} true)

    elseif(required STREQUAL no OR required STREQUAL No OR required STREQUAL false OR required STREQUAL False)
        message(STATUS "${featureName} support disabled by user option, setting ${varName} to false")
        set(${varName} false)

    else(required STREQUAL yes OR required STREQUAL Yes OR required STREQUAL true OR required STREQUAL True)
        # set flags for the compile/run attempt, if provided
        if(ARGV4)
            set(CMAKE_REQUIRED_FLAGS ${ARGV4})
        endif(ARGV4)

        # intel compilers may need a hand finding the c++ std lib
        if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
            find_library(libstdcxx_full_path stdc++)
            set(CMAKE_REQUIRED_LIBRARIES ${libstdcxx_full_path})
        endif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")

        CHECK_CXX_SOURCE_RUNS("${testSnippet}" ${varName})

        if(${varName})
            message(STATUS "${featureName} is supported by the compiler, setting ${varName} to true")
        else(${varName})
            message(STATUS "${featureName} is not supported by the compiler, setting ${varName} to false")
        endif(${varName})

    endif(required STREQUAL yes OR required STREQUAL Yes OR required STREQUAL true OR required STREQUAL True)
    set(${varName} ${${varName}} PARENT_SCOPE)
endfunction()
