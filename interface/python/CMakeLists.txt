find_package(Python COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT} REQUIRED)

find_package(nanobind CONFIG REQUIRED)

if(NOT "${SKBUILD_SABI_VERSION}" STREQUAL "")
  set(STABLE_ABI STABLE_ABI)
endif()

nanobind_add_module(_mutationpp
    ${STABLE_ABI}
    src/mutationpp_python.cpp
    src/pyGlobalOptions.cpp
    src/pyMixtureOptions.cpp
    src/pyMixture.cpp
)

target_link_libraries(_mutationpp PRIVATE mutation++)

install(DIRECTORY mutationpp DESTINATION .)

install(TARGETS _mutationpp DESTINATION mutationpp)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/data" DESTINATION mutationpp)
