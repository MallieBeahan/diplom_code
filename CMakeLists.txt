cmake_minimum_required(VERSION 3.17)
project(diplom_code)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(diplom_code
        include/GAS_CONSTANTS.h
        include/GLOBAL_VARS.H
        include/HEADERS.h
        include/SYSTEM_PARAMS.h
        Help_functions.cpp
        Main.cpp
        Speed_generation.cpp
        Mathemathical_modelling.cpp
        Verlet_scheme.cpp
        Energy_calc.cpp)
