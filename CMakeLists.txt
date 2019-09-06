###############################################
# Preliminaries  
###############################################

cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(odeintWrapper CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Boost REQUIRED)

add_executable(main_tests boost_odeint_wrapper.cpp)

target_link_libraries(main_tests Boost::boost)
target_compile_features(main_tests PUBLIC cxx_std_17)
