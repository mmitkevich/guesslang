cmake_minimum_required(VERSION 3.3)
project(guesslang)
aux_source_directory(. SRC_LIST)
ADD_DEFINITIONS("-std=c++0x -g -O0 -fpermissive")

include(deps/dlib/dlib/cmake)

add_subdirectory(deps/utfcpp)
add_subdirectory(deps/utf8proc)
include_directories(deps/utfcpp/source)
include_directories(deps/utf8proc)

add_executable(${PROJECT_NAME} ${SRC_LIST} utf8pp.h)
target_link_libraries(${PROJECT_NAME} dlib)
target_link_libraries(${PROJECT_NAME} utf8proc)
