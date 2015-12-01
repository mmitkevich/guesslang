cmake_minimum_required(VERSION 3.3)
project(guesslang)
SET(CMAKE_BUILD_TYPE debug)
#SET(CMAKE_BUILD_TYPE release)
add_definitions("-std=c++0x -ggdb -O0 -fkeep-inline-functions -fpermissive")
#add_definitions("-std=c++0x -O3")

aux_source_directory(. SRC_LIST)

include(../deps/dlib/dlib/cmake)

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} dlib)
#target_link_libraries(${PROJECT_NAME} utf8proc)
