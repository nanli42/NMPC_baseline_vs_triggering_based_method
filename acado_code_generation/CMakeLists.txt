#
# Project settings
#

# Minimum required version of cmake
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

# Project name and programming languages used
PROJECT( code_generation_ACADO_for_f110_nmpc CXX )

# CMake module(s) path
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR} )

#
# Prerequisites
#
FIND_PACKAGE( ACADO REQUIRED )

#
# Include directories
#
INCLUDE_DIRECTORIES( . ${ACADO_INCLUDE_DIRS} )

FILE(GLOB SOURCES_N15 "src/formulation_N15.cpp")
FILE(GLOB SOURCES_N30 "src/formulation_N30.cpp")

#
# Build an executable
#
ADD_EXECUTABLE( main_N15 ${SOURCES_N15} )
ADD_EXECUTABLE( main_N30 ${SOURCES_N30} )

# Linux and OS X
TARGET_LINK_LIBRARIES( main_N15 ${ACADO_SHARED_LIBRARIES} )
TARGET_LINK_LIBRARIES( main_N30 ${ACADO_SHARED_LIBRARIES} )
