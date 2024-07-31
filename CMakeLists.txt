cmake_minimum_required(VERSION 3.19.6 FATAL_ERROR)
project(LArProfMaker)

# require C++11
set (CMAKE_CXX_STANDARD 17)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT 6.28.12 CONFIG REQUIRED)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT package not found.")
  return()
endif()

###### Compiler set up
# Add ROOT include directories
include_directories(${ROOT_INCLUDE_DIRS})

# add sub directories
add_subdirectory(ProfileMakers)
