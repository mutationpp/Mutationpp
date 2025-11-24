if (CMAKE_VERSION VERSION_LESS 3.18)
    set(PY_DEV_MODULE Development)
else()
    set(PY_DEV_MODULE Development.Module)
endif()

find_package(Python COMPONENTS Interpreter ${PY_DEV_MODULE} REQUIRED)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE pybind11_ROOT)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_mutationpp
    src/mutationpp_python.cpp
    src/pyMixtureOptions.cpp
    src/pyMixture.cpp
)

target_link_libraries(_mutationpp PRIVATE mutation++)

install(TARGETS _mutationpp DESTINATION .)
