function (add_coverage TNAME)
	# only add coverage for target, if coverage is support and enabled.
    foreach (TNAME ${ARGV})
        add_coverage_target(${TNAME})
    endforeach ()
endfunction (add_coverage)


# Add coverage support for target ${TNAME} and register target for coverage
# evaluation.
function(add_coverage_target TNAME)
    target_compile_options(${TNAME} 
        PUBLIC
            ${COVERAGE_FLAGS}
    )

    target_link_libraries(${TNAME}
        PUBLIC
            gcov
    )
endfunction(add_coverage_target)

if(ENABLE_COVERAGE)
    # Coverage flags
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang")
        message(FATAL_ERROR "Not Implemented! Yet. Use GCC for coverage")
    elseif(CMAKE_COMPILER_IS_GNUCXX)

        find_program(GCOV_EXECUTABLE gcov)
        find_program(LCOV_EXECUTABLE lcov)

        if(NOT GCOV_EXECUTABLE)
            message(FATAL_ERROR "Gcov not found.")
        elseif(NOT LCOV_EXECUTABLE)
            message(FATAL_ERROR "Lcov not found.")
        endif()

        list(APPEND COVERAGE_FLAGS 
            "--coverage"
            "-fprofile-arcs"
            "-ftest-coverage"
            )

        add_custom_target(coverage
            COMMAND ${LCOV_EXECUTABLE} --gcov-tool ${GCOV_EXECUTABLE}
                -c -d . -o "${PROJECT_NAME}.info"

            COMMAND ${LCOV_EXECUTABLE} --gcov-tool ${GCOV_EXECUTABLE}
                --remove "${PROJECT_NAME}.info" ${LCOV_REMOVE_PATTERNS} 
                -o "${PROJECT_NAME}.info"

            COMMAND ${LCOV_EXECUTABLE} --gcov-tool ${GCOV_EXECUTABLE}
                -l ${PROJECT_NAME}.info

            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )

    else()
        message(FATAL_ERROR "Coverage needs GCC or Clang")
    endif()
endif(ENABLE_COVERAGE)

