cmake_minimum_required(VERSION 3.16)
project(PROJECT)
set(CMAKE_CXX_STANDARD 17)


ADD_EXECUTABLE(PROJECT
        source/whitted.cpp
        source/utils.cpp
        source/utils.h
        source/Hittable.h
        source/Hittable.cpp
        source/Sphere.h
        source/Sphere.cpp
        source/MeshTriangle.h
        source/MeshTriangle.cpp
        source/Light.h
        source/Light.cpp
)

target_link_libraries(
        PROJECT LINK_PUBLIC -lyaml-cpp
#        -lGLEW -lGL -lGLU -lglut -lpthread -lm
        )
# FIND ALL YAML FILES
file(GLOB_RECURSE YAML_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.yaml")

# copy config.yaml and input.yaml to build directory
foreach(YAML_FILE ${YAML_FILES})
    configure_file(${YAML_FILE} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()







