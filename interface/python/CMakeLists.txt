find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_mutationpp
    src/mutationpp_python.cpp
    src/pyGlobalOptions.cpp
    src/pyMixtureOptions.cpp
    src/pyMixture.cpp
)

target_link_libraries(_mutationpp PRIVATE mutation++)

install(DIRECTORY mutationpp DESTINATION .)

install(TARGETS _mutationpp DESTINATION mutationpp)
