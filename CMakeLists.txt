cmake_minimum_required(VERSION 3.29)
set(PROJECT EMPPIS_LABS)
project(${PROJECT})

set(CMAKE_CXX_STANDARD 23)

add_executable(${PROJECT}
        main.cpp
        Plot.h
        Bocchiemark.hpp
        Random.hpp
        lab1.hpp
        lab2.hpp
        lab3.hpp
        lab4.hpp
        lab5.hpp
        lab6.hpp
        lab7.hpp
        lab8.hpp
)
