cmake_minimum_required(VERSION 2.8)
project(guesslang)
#SET(CMAKE_BUILD_TYPE debug)
#SET(CMAKE_BUILD_TYPE release)
SET(CMAKE_BUILD_TYPE relwithdebinfo)
add_definitions("-std=c++0x -ggdb") #-O3 -fkeep-inline-functions
#ADD_DEFINITIONS("-std=c++0x -ggdb -g3 -fpermissive")

aux_source_directory(. SRC_LIST)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_AUTOMOC ON)
#include_directories(.)

find_package(Qt5Core)

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} Qt5::Core)
