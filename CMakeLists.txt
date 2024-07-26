cmake_minimum_required(VERSION 3.10)
project(radiation_cpp)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the directories containing the include files
include_directories(/usr/local/include/OpenXLSX)

# Add the directory containing the OpenXLSX library
link_directories(/usr/local/lib)

# Specify the source files
set(SOURCES main.cpp Derivatives.cpp Utilities.cpp)

# Add the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link the libraries
target_link_libraries(${PROJECT_NAME} OpenXLSX boost_system boost_filesystem)
