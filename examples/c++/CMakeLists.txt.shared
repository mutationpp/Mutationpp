cmake_policy(SET CMP0048 NEW)
cmake_minimum_required(VERSION 3.0)

get_filename_component(example_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(${example_name})

# Let's require mutation++. It should export a suitable config that also
# exports include directories for Eigen. The `if` is needed to allow the
# examples to build during the process of building the entire mutation++
# tree.
if(NOT TARGET mutation++)
    find_package(mutation++ REQUIRED)
endif()

add_executable(${example_name} ${example_name}.cpp)
target_link_libraries(${example_name}
    PRIVATE
        mutation++
)
