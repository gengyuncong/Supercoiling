# discover the available regression tests
execute_process(COMMAND ./test.py -ls
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_VARIABLE testIDs
                ERROR_VARIABLE testPyErr
                RESULT_VARIABLE testPyReturn)

# split the discovery results into a list
string(REGEX REPLACE "\n" ";" testIDs "${testIDs}")

foreach(testID IN LISTS testIDs)
    string(REGEX REPLACE ".*:([a-zA-Z0-9_]+)" "\\1" testName "${testID}")
    string(REGEX REPLACE ".*/([a-zA-Z0-9_]+)/[a-zA-Z0-9_]*\\.ipynb:.*" "\\1" testDir "/${testID}")

    # non AVX version of the test
    string(CONCAT testNameVisible ${testDir} ":" ${testName})
    message(STATUS "adding test => ${testNameVisible}")
    add_test(NAME ${testNameVisible}
             COMMAND ./test.py -tid ${testID} -out ${CMAKE_CURRENT_BINARY_DIR} -lm_sbml_import $<TARGET_FILE:lm_sbml_import> -lm_setp $<TARGET_FILE:lm_setp> -lmes $<TARGET_FILE:lmes>
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # AVX version of the test
    string(CONCAT testNameVisibleAVX ${testNameVisible} "_AVX")
    message(STATUS "adding test => ${testNameVisibleAVX}")
    add_test(NAME ${testNameVisibleAVX}
             COMMAND ./test.py -tid ${testID} -useavx -out ${CMAKE_CURRENT_BINARY_DIR} -lm_sbml_import $<TARGET_FILE:lm_sbml_import> -lm_setp $<TARGET_FILE:lm_setp> -lmes $<TARGET_FILE:lmes>
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach(testID)
