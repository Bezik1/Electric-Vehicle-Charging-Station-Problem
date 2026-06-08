cmake_minimum_required(VERSION 3.20)

set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
cmake_policy(SET CMP0135 NEW)

project(EVCSP_Application LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

FetchContent_Declare(
  json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

FetchContent_Declare(
  dlib
  GIT_REPOSITORY https://github.com/davisking/dlib.git
  GIT_TAG        master
)

FetchContent_MakeAvailable(json googletest dlib)

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

set(COMMON_SOURCES
    src/LinearProgrammingModel.cpp
    src/MetaheuristicModel.cpp
    src/Evaluator.cpp
    src/Bacterium.cpp
    src/Maps.cpp
    src/Validator.cpp
    src/SolutionPrinter.cpp
)

add_executable(${PROJECT_NAME} 
    src/main.cpp
    ${COMMON_SOURCES}
)

target_include_directories(${PROJECT_NAME} PRIVATE src/include ${GUROBI_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${GUROBI_CPP_LIBRARY} 
    ${GUROBI_LIBRARY} 
    Threads::Threads
    nlohmann_json::nlohmann_json
    dlib::dlib
)

add_executable(optimize 
    src/optimize.cpp
    ${COMMON_SOURCES}
)

target_include_directories(optimize PRIVATE src/include ${GUROBI_INCLUDE_DIR})
target_link_libraries(optimize PRIVATE 
    ${GUROBI_CPP_LIBRARY} 
    ${GUROBI_LIBRARY} 
    Threads::Threads
    nlohmann_json::nlohmann_json
    dlib::dlib
)

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
    target_compile_options(optimize PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -O3)
    target_compile_options(optimize PRIVATE -Wall -Wextra -O3)
endif()

enable_testing()

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

file(GLOB TEST_SOURCES "tests/*.cpp")

foreach(test_file ${TEST_SOURCES})
    get_filename_component(test_name ${test_file} NAME_WE)
    
    add_executable(${test_name} 
        ${test_file}
        ${COMMON_SOURCES}
    )
    
    target_include_directories(${test_name} PRIVATE src/include ${GUROBI_INCLUDE_DIR})
    target_link_libraries(${test_name} 
        PRIVATE 
            ${GUROBI_CPP_LIBRARY} 
            ${GUROBI_LIBRARY} 
            Threads::Threads
            GTest::gtest
            GTest::gtest_main
            nlohmann_json::nlohmann_json
            dlib::dlib
    )
    
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()