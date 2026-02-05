cmake_minimum_required(VERSION 3.26)
project(Senior_Project)

set(CMAKE_CXX_STANDARD 17)

add_executable(Senior_Project main.cpp
        Newton1D.cpp
        Newton1D.h
        ScalarFunction.cpp
        ScalarFunction.h
        TestQuadraticEq.cpp
        TestQuadraticEq.h
        ODE1D.cpp
        ODE1D.h
        Euler1D.cpp
        Euler1D.h
        ExponentialDecayODE.cpp
        ExponentialDecayODE.h
        Integrator1D.cpp
        Integrator1D.h
        Solution1D.cpp
        Solution1D.h
)
