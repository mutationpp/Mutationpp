pybind11_add_module(_mutationpp
        src/mutationpp_python.cpp
        src/pyMixtureOptions.cpp
        src/pyMixture.cpp
)

target_link_libraries(_mutationpp PRIVATE mutation++)

install(TARGETS _mutationpp DESTINATION .)
