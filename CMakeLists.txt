cmake_minimum_required(VERSION 3.17)
project(diplom_code)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(diplom_code
        Gas_constants.h
        Global_vars.h
        Help_functions.h
        main.cpp
        Molecule.h
        Params.h
        Speed_generation.h
        Start_conditions.h)
