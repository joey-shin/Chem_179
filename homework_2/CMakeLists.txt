cmake_minimum_required(VERSION 3.27)
project(homework_2)

set(CMAKE_CXX_STANDARD 17)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(homework_2 main.cpp)
target_link_libraries(homework_2 ${ARMADILLO_LIBRARIES})
