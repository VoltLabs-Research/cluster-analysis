cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(cluster-analysis VERSION 1.0.0 LANGUAGES CXX)

if(NOT TARGET coretoolkit::coretoolkit)
    find_package(coretoolkit REQUIRED)
endif()

volt_add_plugin(cluster-analysis)
