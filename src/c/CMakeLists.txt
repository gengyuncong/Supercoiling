# c src directory

#################
#TARGET: lm_c_lib
#################
# list of files/dirs to exclude from the build of lm_c_lib
set(to_remove_globs lm/main/Main.cpp lm/main/MainArgs.cpp lm/neus/*)
set(to_remove_recurse_globs old_*)

# get the cpp source files
file(GLOB_RECURSE lm_c_source *.cpp *.cc *.h)

# find any files we need to remove from the cpp sources
file(GLOB to_remove_srcs ${to_remove_globs})
file(GLOB_RECURSE to_remove_recurse_srcs ${to_remove_recurse_globs})

# remove files from the cpp sources as appropriate
list(REMOVE_ITEM lm_c_source ${to_remove_srcs} ${to_remove_recurse_srcs})

# get the proto source files
file(GLOB_RECURSE lm_proto_src ../protobuf/*.proto)

#include_directories(${CMAKE_CURRENT_BINARY_DIR})
PROTOBUF_GENERATE_CPP(proto_srcs proto_hdrs ${lm_proto_src})
PROTOBUF_GENERATE_PYTHON(proto_pb2s ${lm_proto_src})

set_source_files_properties(${proto_srcs} PROPERTIES GENERATED TRUE)
set_source_files_properties(${proto_hdrs} PROPERTIES GENERATED TRUE)

################
# Add a custom command that generates a Version.cpp file.
################

add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/src/c/lm/Version.cpp
  COMMAND echo "#include \"lm/Version.h\"" > ${CMAKE_BINARY_DIR}/src/c/lm/Version.cpp
  COMMAND date "+const char BUILD_TIMESTAMP[] = \"%Y-%m-%d-%H-%M-%S\";" >> ${CMAKE_BINARY_DIR}/src/c/lm/Version.cpp
  VERBATIM
)

add_library(lm_c_lib STATIC ${lm_c_source} ${proto_hdrs} ${proto_srcs} ${CMAKE_BINARY_DIR}/src/c/lm/Version.cpp)

################
# Add a custom command to regenerate Version.cpp whenever lmes is rebuilt.
################

add_custom_command(
  TARGET lm_c_lib
  PRE_BUILD
  COMMAND rm ${CMAKE_BINARY_DIR}/src/c/lm/Version.cpp
)


#set(proto_hdrs ${proto_hdrs} PARENT_SCOPE)
