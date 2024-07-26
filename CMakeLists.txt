cmake_minimum_required(VERSION 3.10)
project(ode_solver)

# Find SUNDIALS
find_package(SUNDIALS REQUIRED)

# Add the executable
add_executable(main main.cpp)

# Link SUNDIALS libraries
target_link_libraries(main SUNDIALS::cvode SUNDIALS::nvecserial)