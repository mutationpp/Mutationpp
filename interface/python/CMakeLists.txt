cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0022 NEW)

pybind11_add_module(mutationpp
        src/mutationpp_python.cpp
        src/pyMixtureOptions.cpp
        src/pyMixture.cpp
)


target_link_libraries(mutationpp PRIVATE mutation++)

find_program(PYTHON "python")

if (PYTHON)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
endif()
