# ================================================================================
# 
# This file is part of the MagBottleC project.
# MagBottleC : RootDirectory
# 
# File: CMakeLists.txt
# 
# Original Author(s) of this File:
# 	Tom Gade, 2017
# 
# Author(s) of Significant Updates/Modifications to the File:
# 	
# 
# This file is heavily modeled after the U of M's MinVR Project CMakeLists.txt. 
# -----------------------------------------------------------------------------------
 
project(MagBottleC)

# CMakeLists files in this project can refer to the root source directory of the 
# project as ${MINVR_SOURCE_DIR} and to the root binary directory of the project 
# as ${MINVR_BINARY_DIR}.

cmake_minimum_required (VERSION 2.8.2) 

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${${PROJECT_NAME}_SOURCE_DIR}/build/install" CACHE PATH "default install path" FORCE )
endif()

# Let's follow the convention of putting all of the libs in "build/lib" and all
# of the exes in "build/bin".  To distinguish between debugging, release, and
# other builds, we'll add a postfix to the name of the lib or exe that we generate.

if(NOT CMAKE_DEBUG_POSTFIX)
  set(CMAKE_DEBUG_POSTFIX "d")
endif()
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_RELWITHDEBINFO_POSTFIX "rd")
set(CMAKE_MINSIZEREL_POSTFIX "s")

if (CMAKE_BUILD_TYPE MATCHES "Release")
  set(CMAKE_BUILD_POSTFIX "${CMAKE_RELEASE_POSTFIX}")
elseif (CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
  set(CMAKE_BUILD_POSTFIX "${CMAKE_MINSIZEREL_POSTFIX}")
elseif (CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
  set(CMAKE_BUILD_POSTFIX "${CMAKE_RELWITHDEBINFO_POSTFIX}")
elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CMAKE_BUILD_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
else()
  set(CMAKE_BUILD_POSTFIX "")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
add_definitions(-DINSTALLPATH="${CMAKE_INSTALL_PREFIX}")

make_directory(${CMAKE_CURRENT_SOURCE_DIR}/../lib)
make_directory(${CMAKE_CURRENT_SOURCE_DIR}/build/bin)
make_directory(${CMAKE_CURRENT_SOURCE_DIR}/build/lib)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
foreach (CONF ${CMAKE_CONFIGURATION_TYPES})
  string (TOUPPER ${CONF} CONF)
  set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONF} ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
  set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONF} ${CMAKE_CURRENT_SOURCE_DIR}/build/lib)
  set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONF} ${CMAKE_CURRENT_SOURCE_DIR}/build/lib)
endforeach(CONF CMAKE_CONFIGURATION_TYPES)


if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Linux specific code
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC") 
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # Apple-specific code
  set(CMAKE_CXX_FLAGS "-DOSX")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if (WIN32)
  # Windows-specific
endif (WIN32)


enable_testing()

set(mbc_base_dir ${CMAKE_CURRENT_SOURCE_DIR})
set(mbc_src_dir ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(mbc_inc_dir ${CMAKE_CURRENT_SOURCE_DIR}/include)

include_directories(
  ${mbc_inc_dir}
)

add_subdirectory(src)

# Organize folder structure for Xcode, Visual Studio, etc.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

