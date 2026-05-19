cmake_minimum_required(VERSION 3.20)
project(EVCSP_Application LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(APPLE)
    set(DEFAULT_GUROBI_HOME "/Library/gurobi1302/macos_universal2")
    set(GUROBI_LIB_NAME "gurobi130")
    set(GUROBI_CPP_LIB_NAME "gurobi_c++")
else()
    set(DEFAULT_GUROBI_HOME "/opt/gurobi1101/linux64")
    set(GUROBI_LIB_NAME "gurobi110")
    set(GUROBI_CPP_LIB_NAME "gurobi_c++")
endif()

set(GUROBI_HOME "${DEFAULT_GUROBI_HOME}" CACHE PATH "Gurobi Installation Path")

find_path(GUROBI_INCLUDE_DIR NAMES gurobi_c++.h PATHS "${GUROBI_HOME}/include")
find_library(GUROBI_LIBRARY NAMES ${GUROBI_LIB_NAME} PATHS "${GUROBI_HOME}/lib")
find_library(GUROBI_CPP_LIBRARY NAMES ${GUROBI_CPP_LIB_NAME} PATHS "${GUROBI_HOME}/lib")

if(NOT GUROBI_INCLUDE_DIR OR NOT GUROBI_LIBRARY OR NOT GUROBI_CPP_LIBRARY)
    message(FATAL_ERROR "Gurobi not found. Check GUROBI_HOME path: ${GUROBI_HOME}")
endif()

find_package(Threads REQUIRED)

add_executable(${PROJECT_NAME} src/main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE src/include ${GUROBI_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY} Threads::Threads)

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -O3)
endif()

enable_testing()

file(GLOB TEST_SOURCES "tests/*.cpp")

foreach(test_file ${TEST_SOURCES})
    get_filename_component(test_name ${test_file} NAME_WE)
    
    add_executable(${test_name} ${test_file})
    target_include_directories(${test_name} PRIVATE src/include ${GUROBI_INCLUDE_DIR})
    target_link_libraries(${test_name} PRIVATE ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY} Threads::Threads)
    
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()