# a module with functions that extend/fix issues with CMake integration if you're using the CLion IDE

function(include_directories_clion_generated)
    message("lm_BINARY_DIR: ${lm_BINARY_DIR}")
    string(REGEX MATCH "/(cmake-build-[a-zA-Z]+)" BUILT_WITH_CLION "${lm_BINARY_DIR}")
    if(BUILT_WITH_CLION)
        include_directories(${CMAKE_BINARY_DIR}/src/c/)
    endif(BUILT_WITH_CLION)
endfunction()