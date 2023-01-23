cmake_minimum_required(VERSION 3.0.0)

project("calculsci" CXX)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra")

add_executable(main main.cpp maillage.cpp noeud.cpp triangle.cpp vApprox.cpp)