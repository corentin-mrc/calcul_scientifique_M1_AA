cmake_minimum_required(VERSION 3.0.0)

project("calculsci" CXX)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wextra -pedantic-errors")

add_executable(main main.cpp noeud.cpp triangle.cpp maillage.cpp boite_a_outils.cpp)
