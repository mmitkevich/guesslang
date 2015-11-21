cmake_minimum_required(VERSION 3.3)
project(guesslang)
aux_source_directory(guesslang SRC_LIST)
include_directories(guesslang)

ADD_DEFINITIONS("-std=c++0x -g -O0 -fpermissive")

include(deps/dlib/dlib/cmake)

add_subdirectory(deps/utf8proc)
include_directories(deps/utfcpp/source)
include_directories(deps/utf8proc)

add_executable(${PROJECT_NAME} ${SRC_LIST} guesslang/config.h)
target_link_libraries(${PROJECT_NAME} dlib)
target_link_libraries(${PROJECT_NAME} utf8proc)
