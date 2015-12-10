cmake_minimum_required(VERSION 3.3)
project(guesslang)
SET(CMAKE_BUILD_TYPE debug)
#SET(CMAKE_BUILD_TYPE release)
add_definitions("-std=c++0x -ggdb -O0 -fkeep-inline-functions  -DQT_FORCE_ASSERTS")
#-fpermissive
#add_definitions("-std=c++0x -O3")

#aux_source_directory(. SRC_LIST)

file(GLOB SRC "*.h" "*.cpp")
file(GLOB QQ_SRC "QQ/*.h" "QQ/*.cpp")

#set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
#set(CMAKE_AUTOMOC ON)
#include_directories(.)

find_package(Qt5Core)




#ADD_DEFINITIONS("-std=c++0x -ggdb -g3 -fpermissive")

#-O0

#include(../deps/dlib/dlib/cmake)


#add_subdirectory(../deps/utf8proc "${CMAKE_CURRENT_BINARY_DIR}/deps")
#include_directories(../deps/utfcpp/source)
#include_directories(../deps/utf8proc)
#include_directories(../deps/cxxopt/src)
#include_directories(../deps/format)

add_executable(${PROJECT_NAME} ${SRC} ${QQ_SRC})
#target_link_libraries(${PROJECT_NAME} dlib)
#target_link_libraries(${PROJECT_NAME} utf8proc)
target_link_libraries(${PROJECT_NAME} Qt5::Core)

add_executable(test_qq ${QQ_SRC} test/test_qq.cpp)
target_link_libraries(test_qq Qt5::Core)
