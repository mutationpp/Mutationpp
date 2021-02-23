cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0022 NEW)

# set(MPP_INSTALL_DIR  ../../install)
# include_directories(${MPP_INSTALL_DIR}/include/mutation++)
# include_directories(${MPP_INSTALL_DIR}/include/)
# link_directories(${MPP_INSTALL_DIR}/lib)

# include_directories(../../thirdparty/eigen/)




pybind11_add_module(mutationpp
        src/mutationpp_python.cpp
        src/python_MixtureOptions.cpp
        src/python_Mixture.cpp
)


target_link_libraries(mutationpp PRIVATE mutation++)

find_program(PYTHON "python")

if (PYTHON)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    # set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    # add_custom_command(
    # OUTPUT ${OUTPUT}/timestamp
    # COMMAND ${PYTHON} setup.py build
    # COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}/timestamp
    # )

    # add_custom_target(target ALL ${OUTPUT}/timestamp)

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
endif()
