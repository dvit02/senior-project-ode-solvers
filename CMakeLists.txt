cmake_minimum_required(VERSION 3.16)
project(Senior_Project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----------------------------
# Library: the framework
# ----------------------------
add_library(ode_solvers
        src/Newton1D.cpp
        src/ScalarFunction.cpp
        src/Euler1D.cpp
        src/RK4_1D.cpp
        src/ExponentialDecayODE.cpp
        src/Euler.cpp
)

target_include_directories(ode_solvers
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ----------------------------
# Example 1: your current main (convergence demo)
# ----------------------------
add_executable(convergence_example
        examples/main.cpp
)
target_link_libraries(convergence_example PRIVATE ode_solvers)

# ----------------------------
# Example 2: oscillator system demo
# ----------------------------
add_executable(oscillator_example
        examples/oscillator_example.cpp
)
target_link_libraries(oscillator_example PRIVATE ode_solvers)

# ----------------------------
# Example 3: Newton quadratic test
# ----------------------------
add_executable(test_quadratic_example
        examples/TestQuadraticEq.cpp
)
target_link_libraries(test_quadratic_example PRIVATE ode_solvers)
