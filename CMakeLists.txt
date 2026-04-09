cmake_minimum_required(VERSION 3.16)
project(Senior_Project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----------------------------
# Library: the framework
# ----------------------------
add_library(ode_solvers
        src/ExponentialDecayODE.cpp
        src/Euler.cpp
        src/RK4.cpp
        src/HarmonicOscillator.cpp
        src/VanDerPol.cpp
        src/Integrator.cpp
        src/VelocityVerlet.cpp
        src/DormandPrince.cpp
        src/DoublePendulum.cpp
)

target_include_directories(ode_solvers
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ----------------------------
# Example 1: convergence demo
# ----------------------------
add_executable(convergence_example
        examples/convergence_example.cpp
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
# Example 4: Van der Pol demo
# ----------------------------
add_executable(vanderpol_example
        examples/vanderpol_example.cpp
)
target_link_libraries(vanderpol_example PRIVATE ode_solvers)

# ----------------------------
# Example 5: Velocity Verlet oscillator demo
# ----------------------------
add_executable(verlet_oscillator_example
        examples/verlet_oscillator_example.cpp
)
target_link_libraries(verlet_oscillator_example PRIVATE ode_solvers)

# ----------------------------
# Example 6: Dormand-Prince RK45 vs adaptive RK4
# ----------------------------
add_executable(dp45_example
        examples/dp45_example.cpp
)
target_link_libraries(dp45_example PRIVATE ode_solvers)


add_executable(benchmark_example examples/benchmark_example.cpp)
target_link_libraries(benchmark_example ode_solvers)

add_executable(double_pendulum_example
        examples/double_pendulum_example.cpp
)
target_link_libraries(double_pendulum_example PRIVATE ode_solvers)

# ----------------------------
# Example 8: Double pendulum tolerance sweep for accuracy ranking
# ----------------------------
add_executable(double_pendulum_sweep
        examples/double_pendulum_sweep.cpp
)
target_link_libraries(double_pendulum_sweep PRIVATE ode_solvers)


# ----------------------------
# Unit tests
# ----------------------------
enable_testing()

add_executable(fixed_steppers_unit_test
        tests/test_fixed_steppers.cpp
)
target_link_libraries(fixed_steppers_unit_test PRIVATE ode_solvers)

add_executable(dp45_unit_test
        tests/test_dormand_prince.cpp
)
target_link_libraries(dp45_unit_test PRIVATE ode_solvers)

add_executable(velocity_verlet_unit_test
        tests/test_velocity_verlet.cpp
)
target_link_libraries(velocity_verlet_unit_test PRIVATE ode_solvers)

add_executable(test_double_pendulum tests/test_double_pendulum.cpp
        src/DoublePendulum.cpp
        src/Integrator.cpp
        src/RK4.cpp)

target_include_directories(test_double_pendulum PRIVATE include)


add_test(NAME fixed_steppers_unit_test COMMAND fixed_steppers_unit_test)
add_test(NAME dp45_unit_test COMMAND dp45_unit_test)
add_test(NAME velocity_verlet_unit_test COMMAND velocity_verlet_unit_test)
add_test(NAME test_double_pendulum COMMAND test_double_pendulum)





