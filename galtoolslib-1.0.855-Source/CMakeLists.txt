cmake_minimum_required(VERSION 2.8.8)

PROJECT(galtoolslib CXX Fortran C)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

set(galtoolslib_VERSION_MAJOR "1")
set(galtoolslib_VERSION_MINOR "0")
set(galtoolslib_VERSION_PATCH "855")
set(galtoolslib_VERSION "${galtoolslib_VERSION_MAJOR}.${galtoolslib_VERSION_MINOR}.${galtoolslib_VERSION_PATCH}")

#Enable rpaths in install directories
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#Make position independent code the default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#Find the required packages
find_package(CCfits REQUIRED)
find_package(HEALPix REQUIRED)
find_package(WCS)
find_package(CLHEP REQUIRED)
find_package(GSL REQUIRED)
find_package(Xerces REQUIRED)


SET(ENABLE_OPENCL FALSE CACHE BOOL "Enable OpenCL functionality")
if (ENABLE_OPENCL)
   find_package(OpenCL)
endif()

if (CCFITS_FOUND)
   set(HAVE_CCFITS ON)
endif()

if (HEALPIX_FOUND)
   set(HAVE_HEALPIX ON)
endif()

if (OpenCL_FOUND)
   #The OpenCL version check is not reliable.  Implementing my own.
   include(CheckSymbolExists)
   set(CMAKE_REQUIRED_INCLUDES ${OpenCL_INCLUDE_DIRS})
   foreach(VERSION "2_0" "1_2")
      if (APPLE)
         CHECK_SYMBOL_EXISTS("CL_VERSION_${VERSION}" "OpenCL/opencl.h" CL_VERSION_${VERSION})
      else()
         CHECK_SYMBOL_EXISTS("CL_VERSION_${VERSION}" "CL/opencl.h" CL_VERSION_${VERSION})
      endif()
      if (CL_VERSION_${VERSION})
         string(REPLACE "_" "" VERSION "${VERSION}")
         set(OPENCL_VERSION "${VERSION}0")
         add_definitions("-DCL_HPP_TARGET_OPENCL_VERSION=${OPENCL_VERSION}")
         add_definitions("-DCL_HPP_MINIMUM_OPENCL_VERSION=120")
         break()
      endif()
   endforeach()
   if (NOT OPENCL_VERSION)
      message(WARNING "Could not detect openCL version, defaulting to 1.2")
      set(OPENCL_VERSION "120")
      add_definitions("-DCL_HPP_TARGET_OPENCL_VERSION=${OPENCL_VERSION}")
      add_definitions("-DCL_HPP_MINIMUM_OPENCL_VERSION=120")
   endif()
   set(HAVE_OPENCL ON)
   add_definitions("-DHAVE_OPENCL")
   include_directories(${OpenCL_INCLUDE_DIRS})
else ()
   set(OpenCL_INCLUDE_DIRS "")
   set(OpenCL_LIBRARIES "")
endif ()

if (WCS_FOUND)
   set(HAVE_WCS ON)
   include_directories(${WCS_INCLUDE_DIR})
   add_definitions(-DHAVE_WCS)
else ()
   set(WCS_LIBRARY "")
   set(WCS_INCLUDE_DIR "")
endif()

#Turn testing on
enable_testing()
  
#Set the default build type.  User can over-ride in cache or on command line.
if (NOT CMAKE_BUILD_TYPE)
   set(DEFAULT_BUILD_TYPE RelWithDebInfo)
   message(STATUS "CMAKE_BUILD_TYPE not set, defaulting to ${DEFAULT_BUILD_TYPE}.\n Set in cache or on command line to override.")
   set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE})
endif (NOT CMAKE_BUILD_TYPE)

#Check for debug build type
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   add_definitions( "-DDEBUG" )
   set(DEBUG ON)
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

#Set the compiler flags, print a message to the user about the over-ruling of the Cache variables
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
if (CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo|Release")
   message(STATUS "The configure script will now override the cached CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER} variable to a hand crafted value for most compilers.\nCMAKE_CXX_FLAGS can be used to add compiler flags.  To manually define the build flags a CMAKE_BUILD_TYPE other than RelWithDebInfo (the default) and Release should be used.")
endif ()

if (CMAKE_CXX_COMPILER MATCHES ".*g[+][+].*" OR CMAKE_COMPILER_IS_GNUCXX)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -pipe -O3 -Wall")
   set(CMAKE_CXX_FLAGS_RELEASE "-pipe -O3 -Wall")
   set(CMAKE_CXX_FLAGS_DEBUG "-pipe -O0 -Wall -g")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (CMAKE_CXX_COMPILER MATCHES ".*ic(c|pc).*")
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -pipe -O3 -no-ipo -Wall -gcc")
   set(CMAKE_CXX_FLAGS_RELEASE "-pipe -O3 -no-ipo -Wall -gcc")
   set(CMAKE_CXX_FLAGS_DEBUG "-pipe -O0 -g -Wall -gcc")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# This is a check for OSX version: versions 13 and higher need -stdlib=libc++ appended to the CXXFLAGS
  IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
     EXEC_PROGRAM(uname ARGS -v OUTPUT_VARIABLE DARWIN_VERSION)
     STRING(REGEX MATCH "[0-9]+" DARWIN_VERSION ${DARWIN_VERSION})
     MESSAGE(STATUS "DARWIN_VERSION=${DARWIN_VERSION}")
     IF (DARWIN_VERSION GREATER 12)
     	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -stdlib=libc++")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -stdlib=libc++")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -stdlib=libc++")
     ENDIF (DARWIN_VERSION GREATER 12)
  ENDIF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#elseif (CMAKE_CXX_COMPILER MATCHES ".*open[cC][cC].*")
#   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -pipe -O3 -Wall -apo -ffast-math -mso")
#   set(CMAKE_CXX_FLAGS_RELEASE "-pipe -O3 -Wall -apo -ffast-math -mso")
#elseif (CMAKE_CXX_COMPILER MATCHES ".*sun.*")
#   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -fast")
#   set(CMAKE_CXX_FLAGS_RELEASE "-fast")
endif ()

if (CMAKE_Fortran_COMPILER MATCHES ".*(gfortran|g95).*" OR CMAKE_COMPILER_IS_GNUFortran)
   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -pipe -O3 -Wall")
   set(CMAKE_Fortran_FLAGS_RELEASE "-pipe -O3 -Wall")
   set(CMAKE_Fortran_FLAGS_DEBUG "-pipe -O0 -g -Wall")
elseif (CMAKE_Fortran_COMPILER MATCHES ".*ifort.*")
   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O3 -no-ipo -warn -gcc -cxxlib")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -no-ipo -warn -gcc -cxxlib")
   set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -warn -gcc -cxxlib")
#elseif (CMAKE_Fortran_COMPILER MATCHES ".*openf9.*")
#   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -pipe -O3 -Wall -apo -ffast-math -mso")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-pipe -O3 -Wall -apo -ffast-math -mso")
#elseif (CMAKE_Fortran_COMPILER MATCHES ".*sun.*")
#   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -fast")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
endif ()

#Use openmp if available and not specifically disabled
set(ENABLE_OPENMP ON CACHE BOOL "Multi-threaded parallel version using OpenMP. Requires a compiler with OpenMP support.")
if (ENABLE_OPENMP)
   find_package(OpenMP)
   if (OPENMP_FOUND)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
      set(HAVE_OPENMP ON)
      add_definitions(-DHAVE_OPENMP)
   endif (OPENMP_FOUND)
endif (ENABLE_OPENMP)

set(galtoolslib_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
set(galtoolslib_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "Using CXX flags ${galtoolslib_CXX_FLAGS} for compiler ${CMAKE_CXX_COMPILER}") 
message(STATUS "Using Fortran flags ${galtoolslib_Fortran_FLAGS} for compiler ${CMAKE_Fortran_COMPILER}") 

# add name mangling scheme for fortran compilers
enable_language(C)
include(FortranCInterface)
FortranCInterface_HEADER(fortranMangling.h MACRO_NAMESPACE "FC_")
#This creates a "fortranMangling.h" header that defines mangling macros
# FC_GLOBAL(), FC_GLOBAL_(), FC_MODULE(), and FC_MODULE_().
include_directories(${CMAKE_BINARY_DIR}/src)

# Make sure the mangling is working
FortranCInterface_VERIFY(CXX)

#Set the compression flag if needed
set(ENABLE_COMPRESSION ON CACHE BOOL "Enable fits file compression.")
if (ENABLE_COMPRESSION)
   add_definitions(-DENABLE_COMPRESSION)
endif()

#Enable user to set the installation directories for binaries and libraries
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_DATA_DIR share CACHE PATH "Installation directory for data files")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

#Test for native regular expressions
include(CheckForRegex)
check_cxx_native_regex_works(native_regex_found)
set(USE_NATIVE_REGEX ${native_regex_found} CACHE BOOL "Override native regex lookup")
add_definitions("-DUSE_NATIVE_REGEX=${USE_NATIVE_REGEX}")

#Check for size of size_t and define to unsigned int if needed
#NOT WORKING
#include(CheckTypeSize)
#list(APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS .c)
#CHECK_TYPE_SIZE("size_t" SIZEOF_SIZE_T)
#if (SIZEOF_SIZE_T EQUAL "0")
#   #I have no idea if this actually works
#   #Should only be useful on macs that have multiple definitions of size_t with universal binaries
#   set(size_t ${SIZEOF_SIZE_T_CODE})
#elseif(NOT SIZEOF_SIZE_T)
#   #Define it to unsigned int
#   set(size_t "unsigned int")
#endif()

#Create the config.h
configure_file ("${PROJECT_SOURCE_DIR}/cmake_config.h.in" "${PROJECT_BINARY_DIR}/config.h" )

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

#Add the include directories from HEALPix and CCfits
include_directories(${CCfits_INCLUDE_DIRS} ${HEALPix_INCLUDE_DIRS} ${XERCESC_INCLUDE} ${CLHEP_INCLUDE_DIRS} ${GSL_INCLUDE_DIR})

#Default to include current directory
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#Add the install lib dir to rpaths
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

#Add the subdirectories
add_subdirectory(utils)
add_subdirectory(galstruct)
add_subdirectory(libskymap)
add_subdirectory(nuclei)
add_subdirectory(processes)
add_subdirectory(random)
add_subdirectory(cppUnitTests)
add_subdirectory(manualTests)

#install the cl.hpp file 
if (HAVE_OPENCL)
   install( FILES ${PROJECT_SOURCE_DIR}/CL/cl.hpp DESTINATION ${INSTALL_INCLUDE_DIR}/CL COMPONENT dev )
   install( FILES ${PROJECT_SOURCE_DIR}/CL/cl2.hpp DESTINATION ${INSTALL_INCLUDE_DIR}/CL COMPONENT dev )
endif ()

#export the targets
export(TARGETS skymap utils nuclei galstruct processes random FILE "${PROJECT_BINARY_DIR}/galtoolslibLibraryDepends.cmake")

#Export the package so it is available in the CMake cache without manual interventions by the user
export(PACKAGE galtoolslib)

# Create the galtoolslibConfig.cmake and galtoolslibConfigVersion files
include(CMakePackageConfigHelpers)
configure_package_config_file(galtoolslibConfig.cmake.in ${PROJECT_BINARY_DIR}/galtoolslibConfig.cmake 
   INSTALL_DESTINATION ${INSTALL_CMAKE_DIR} 
   PATH_VARS INSTALL_INCLUDE_DIR INSTALL_CMAKE_DIR INSTALL_DATA_DIR INSTALL_BIN_DIR)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/galtoolslibConfigVersion.cmake 
   VERSION ${galtoolslib_VERSION}
   COMPATIBILITY SameMajorVersion )

#Install the config files and the export file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/galtoolslibConfig.cmake 
   ${CMAKE_CURRENT_BINARY_DIR}/galtoolslibConfigVersion.cmake 
   DESTINATION ${INSTALL_CMAKE_DIR} )

install(EXPORT galtoolslibLibraryDepends DESTINATION ${INSTALL_CMAKE_DIR} COMPONENT dev)

#Install the cmake macros in cmake/Modules directory
file(GLOB macros cmake/Modules/*.cmake)
install(FILES ${macros} DESTINATION ${INSTALL_CMAKE_DIR})

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${galtoolslib_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${galtoolslib_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${galtoolslib_VERSION_PATCH}")
set (CPACK_SOURCE_GENERATOR TBZ2)
#List of files to ignore.  These are actually a list of regular expressions.
#set (CPACK_SOURCE_IGNORE_FILES "")
include (CPack)
