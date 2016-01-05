cmake_minimum_required(VERSION 2.8)
project(rlShaders)

set(PROJECT_VERSION_MAJOR  0)
set(PROJECT_VERSION_MINOR  1)
set(PROJECT_VERSION_PATCH  0)
set(PROJECT_VERSION
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
set(VERSION ${PROJECT_VERSION})

if(DEFINED ENV{ARNOLD_ROOT})
    set(ARNOLD_ROOT $ENV{ARNOLD_ROOT})
endif()

if(NOT ARNOLD_ROOT)
    message(FATAL_ERROR "Can't find ARNOLD_ROOT!")
endif()

find_path(ARNOLD_INCLUDE_DIRS ai.h
          PATHS
          "${ARNOLD_ROOT}/include"
          NO_DEFAULT_PATH
          NO_CMAKE_ENVIRONMENT_PATH
          NO_CMAKE_PATH
          NO_SYSTEM_ENVIRONMENT_PATH
          NO_CMAKE_SYSTEM_PATH)

find_library(ARNOLD_LIBRARIES ai
             PATHS
             "${ARNOLD_ROOT}/lib/"
             "${ARNOLD_ROOT}/bin/"
             NO_DEFAULT_PATH
             NO_CMAKE_ENVIRONMENT_PATH
             NO_CMAKE_PATH
             NO_SYSTEM_ENVIRONMENT_PATH
             NO_CMAKE_SYSTEM_PATH)

if(NOT ARNOLD_INCLUDE_DIRS)
    message(FATAL_ERROR "Can't find ARNOLD SDK!")
else()
    message(STATUS "Found include directory of ARNOLD SDK: ${ARNOLD_INCLUDE_DIRS}")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(ORIGINAL_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CMAKE_INSTALL_PREFIX ${ORIGINAL_INSTALL_PREFIX}/rlShaders-${VERSION})
message(STATUS "rlShaders will install to: ${CMAKE_INSTALL_PREFIX}")

if(MSVC)
    add_definitions (-D_CRT_SECURE_NO_DEPRECATE)
    add_definitions (-D_CRT_SECURE_NO_WARNINGS)
    add_definitions (-D_CRT_NONSTDC_NO_WARNINGS)
    add_definitions (-D_SCL_SECURE_NO_WARNINGS)
    set(ARNOLD_COMPILE_FLAGS "/EHsc")
    set(ARNOLD_LINK_FLAGS "/nologo")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(ARNOLD_COMPILE_FLAGS "-c -fPIC -D_LINUX -std=c++14")
    set(ARNOLD_LINK_FLAGS "-shared")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(ARNOLD_COMPILE_FLAGS "-c -std=c++14")
    set(ARNOLD_LINK_FLAGS "-bundle -undefined dynamic_lookup")
endif()

add_subdirectory(scripts)
add_subdirectory(src)
