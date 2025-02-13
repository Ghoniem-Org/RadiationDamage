cmake_minimum_required(VERSION 3.10)
project(ode_solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find SUNDIALS
find_package(SUNDIALS REQUIRED)

# Add the executable
# add_executable(main_spatial main_spatial.cpp)
add_executable(main_spatial_from_curve main_spatial_from_curve.cpp)

# Link SUNDIALS libraries
# target_link_libraries(main_spatial SUNDIALS::cvode SUNDIALS::nvecserial) 
target_link_libraries(main_spatial_from_curve SUNDIALS::cvode SUNDIALS::nvecserial) 

