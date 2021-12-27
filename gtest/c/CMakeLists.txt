# gtest c src directory
add_subdirectory(gmock-1.7.0)
include_directories(${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR}
                    ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR})

##############################################
#### setup the test source and data files ####
##############################################

## set the testdata directory var. This will be directly substituted into the gtest source files
#set(TESTDATA_ROOT ${CMAKE_CURRENT_BINARY_DIR}/testdata)
#
## copy over the test data files
#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/testdata DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# set the testdata directory var. This will be directly substituted into the gtest source files
set(TESTDATA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/testdata)

# add a compiler definition with the path to the root testdata dir
add_definitions(-DTESTDATA_ROOT="${TESTDATA_ROOT}")

############################
#### TARGET: lmes_gtest ####
############################

# get the source files
file(GLOB_RECURSE lmes_gtest_c_source lm/*.cpp lm/*.cc lm/*.h)

# remove all files with certain prefixes
file(GLOB_RECURSE */to_remove */off_* */old_*)
if(NOT "${to_remove}" STREQUAL "")
    # REMOVE_ITEM causes problems if to_remove is empty, so guard against that
    list(REMOVE_ITEM lmes_gtest_c_source ${to_remove})
endif(NOT "${to_remove}" STREQUAL "")

add_executable(lmes_gtest ${lmes_gtest_c_source})
set_target_properties(gmock gmock_main gtest gtest_main lmes_gtest PROPERTIES EXCLUDE_FROM_ALL TRUE)

SET(LM_C_LIB lm_c_lib)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(lmes_gtest -Wl,--whole-archive)
    target_link_libraries(lmes_gtest ${LM_C_LIB})
    target_link_libraries(lmes_gtest -Wl,--no-whole-archive)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(lmes_gtest -Wl,-force_load)
    target_link_libraries(lmes_gtest ${LM_C_LIB})
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
if(BUILD_WITH_CUDA)
    target_link_libraries(lmes_gtest lm_cuda_lib)
endif(BUILD_WITH_CUDA)

target_link_libraries(lmes_gtest ${lm_dependency_libs}
                                 gmock
                                 gmock_main
                                 gtest
                                 gtest_main)

## add the gtest unittest to ctest, so it will run with `make test`
#add_test(
#    NAME lmes_gtest
#    COMMAND lmes_gtest
#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#)