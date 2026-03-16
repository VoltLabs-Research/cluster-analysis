#
# Copyright (C) Rodolfo Herrera Hernandez. All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(cluster-analysis VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)

find_package(coretoolkit REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

set(CLUSTER_ANALYSIS_BOOST_TARGET "")
if(TARGET Boost::headers)
    set(CLUSTER_ANALYSIS_BOOST_TARGET Boost::headers)
elseif(TARGET Boost::boost)
    set(CLUSTER_ANALYSIS_BOOST_TARGET Boost::boost)
elseif(TARGET boost::headers)
    set(CLUSTER_ANALYSIS_BOOST_TARGET boost::headers)
elseif(TARGET boost::boost)
    set(CLUSTER_ANALYSIS_BOOST_TARGET boost::boost)
else()
    message(FATAL_ERROR "ClusterAnalysis requires a Boost headers target")
endif()

# Collect all library sources (exclude main.cpp)
file(GLOB LIB_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/*.cpp)
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

# Create static library with algorithm implementation
add_library(${PROJECT_NAME}_lib STATIC ${LIB_SOURCES})
target_include_directories(${PROJECT_NAME}_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME}_lib PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC 
    ${CLUSTER_ANALYSIS_BOOST_TARGET}
    coretoolkit::coretoolkit 
    spdlog::spdlog
    nlohmann_json::nlohmann_json
)
set_target_properties(${PROJECT_NAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Create executable
add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME} PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${PROJECT_NAME}_lib 
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(TARGETS ${PROJECT_NAME}_lib DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
