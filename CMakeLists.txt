cmake_minimum_required(VERSION 2.8.8)

PROJECT(GALPROP CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

set(GALPROP_VERSION_MAJOR "56")
set(GALPROP_VERSION_MINOR "0")
set(GALPROP_VERSION_PATCH "2870")
set(GALPROP_VERSION "${GALPROP_VERSION_MAJOR}.${GALPROP_VERSION_MINOR}.${GALPROP_VERSION_PATCH}")

#Enable rpaths in install directories
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#Find the required packages
set(CMAKE_FIND_ROOT_PATH /home/stephan/galprop/galtoolslib-1.0.855-Source/build)
find_package(galtoolslib REQUIRED ONLY_CMAKE_FIND_ROOT_PATH)
message(STATUS "galtoolslib_lib (${galtoolslib_DIR})")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${galtoolslib_CMAKE_DIR})
message(STATUS "galtoolslib (${galtoolslib_CMAKE_DIR})")

#Make position independent code the default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#Propagate the build type from galtoolslib unless it is specifically set
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "CMAKE_BUILD_TYPE not set, defaulting to build type from galtoolslib (${galtoolslib_BUILD_TYPE})")
   set(CMAKE_BUILD_TYPE ${galtoolslib_BUILD_TYPE})
endif (NOT CMAKE_BUILD_TYPE)

#Check for debug build type
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   add_definitions( "-DDEBUG" )
   set(DEBUG ON)
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

#Set the compiler flags to the same ones as galtoolslib unless the build_type or compiler differs
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
if (CMAKE_BUILD_TYPE STREQUAL galtoolslib_BUILD_TYPE AND CMAKE_CXX_COMPILER STREQUAL galtoolslib_CXX AND CMAKE_Fortran_COMPILER STREQUAL galtoolslib_Fortran)
   message(STATUS "Using compiler flags from galtoolslib.  If that is not desired, please use a CMAKE_BUILD_TYPE different from the one used in galtoolslib.")
   set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER} ${galtoolslib_CXX_FLAGS})
   set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER} ${galtoolslib_Fortran_FLAGS})
else()
   message(STATUS "To ensure compatibility with galtoolslib please use the same compiler and build type")
endif ()

#if (CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo|Release")
#   message(STATUS "The configure script will now override the cached CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER} variable to a hand crafted value for most compilers.\nCMAKE_Fortran_FLAGS can be used to add compiler flags.  To manually define the build flags a CMAKE_BUILD_TYPE other than RelWithDebInfo and Release should be used.")
#endif ()

#if (CMAKE_Fortran_COMPILER MATCHES ".*(gfortran|g95).*" OR CMAKE_COMPILER_IS_GNUFortran)
#   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -pipe -O3 -Wall -mtune=native")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-pipe -O3 -Wall -mtune=native")
#elseif (CMAKE_Fortran_COMPILER MATCHES ".*ifort.*")
#   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -pipe -O3 -xhost -no-ipo -parallel -warn -gcc -cxxlib")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-pipe -O3 -xhost -no-ipo -parallel -warn -gcc -cxxlib")
#elseif (CMAKE_Fortran_COMPILER MATCHES ".*openf9.*")
#   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -pipe -O3 -Wall -apo -ffast-math -mso")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-pipe -O3 -Wall -apo -ffast-math -mso")
#elseif (CMAKE_Fortran_COMPILER MATCHES ".*sun.*")
#   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -fast")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
#endif ()

#Use openmp if available and not specifically disabled
set(ENABLE_OPENMP ON CACHE BOOL "Multi-threaded parallel version using OpenMP. Requires a compiler with OpenMP support.")
if (ENABLE_OPENMP)
   find_package(OpenMP)
   if (OPENMP_FOUND)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
      set(HAVE_OPENMP ON)
   endif (OPENMP_FOUND)
else()
   if (CMAKE_BUILD_TYPE STREQUAL galtoolslib_BUILD_TYPE)
      message(STATUS "OpenMP may still be enabled for GALPROP if it was enabled for galtoolslib")
   endif()
endif ()

# add name mangling scheme for fortran compilers
enable_language(C)
include(FortranCInterface)
FortranCInterface_HEADER(fortranMangling.h MACRO_NAMESPACE "FC_")
#This creates a "fortranMangling.h" header that defines mangling macros
# FC_GLOBAL(), FC_GLOBAL_(), FC_MODULE(), and FC_MODULE_().
include_directories(${CMAKE_BINARY_DIR}/src)

# Make sure the mangling is working
FortranCInterface_VERIFY(CXX)

#Default paths for fits data and galdef files
set(FITSDATA_PATH "" CACHE PATH "Default path to fits data files")
set(GALDEF_PATH "" CACHE PATH "Default path to GALDEF files")

if (FITSDATA_PATH)
   set (HAVE_FITSDATA 1)
   add_definitions("-DFITSDATA_PATH=\"${FITSDATA_PATH}\"")
endif()
if (GALDEF_PATH)
   set (HAVE_GALDEF 1)
   add_definitions("-DGALDEF_PATH=\"${GALDEF_PATH}\"")
endif()
add_definitions("-DGALTOOLSLIB_DATA_PATH=\"${galtoolslib_DATA_DIR}\"")

set(GALPROP_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "Using CXX flags ${GALPROP_CXX_FLAGS} for compiler ${CMAKE_CXX_COMPILER}") 
set(GALPROP_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "Using Fortran flags ${GALPROP_Fortran_FLAGS} for compiler ${CMAKE_Fortran_COMPILER}") 

#Add the definitions from galtoolslib
add_definitions(${galtoolslib_DEFINITIONS})

#Enable user to set the installation directories for binaries and libraries
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_DATA_DIR share CACHE PATH "Installation directory for data files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

#Add the path to data in definitions
if (IS_ABSOLUTE ${INSTALL_DATA_DIR})
   set(DATA_PATH ${INSTALL_DATA_DIR})
else()
   set(DATA_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_DATA_DIR})
endif()
add_definitions("-DDATA_PATH=\"${DATA_PATH}\"")

#Add the install lib dir to rpaths
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

#Create the config.h
configure_file ("${PROJECT_SOURCE_DIR}/cmake_config.h.in" "${PROJECT_BINARY_DIR}/config.h" )

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

#Add the include directories from HEALPix and CCfits
include_directories(${galtoolslib_INCLUDE_DIRS})

#Default to include current directory
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#Add the install lib dir to rpaths
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

#Add the subdirectories
add_subdirectory(source)
#add_subdirectory(tests)

#Export the targets
export(TARGETS galprop galpropapp FILE "${PROJECT_BINARY_DIR}/GALPROPLibraryDepends.cmake")

#Export the package so it is available in the CMake cache without manual interventions by the user
export(PACKAGE GALPROP)

# Create the GALPROPConfig.cmake and GALPROPConfigVersion files
include(CMakePackageConfigHelpers)
configure_package_config_file(GALPROPConfig.cmake.in ${PROJECT_BINARY_DIR}/GALPROPConfig.cmake 
   INSTALL_DESTINATION ${INSTALL_CMAKE_DIR} 
   PATH_VARS INSTALL_INCLUDE_DIR INSTALL_CMAKE_DIR INSTALL_BIN_DIR)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/GALPROPConfigVersion.cmake 
   VERSION ${GALPROP_VERSION}
   COMPATIBILITY SameMajorVersion )

#Install the config files and the export file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GALPROPConfig.cmake 
   ${CMAKE_CURRENT_BINARY_DIR}/GALPROPConfigVersion.cmake 
   DESTINATION ${INSTALL_CMAKE_DIR} )

install(EXPORT GALPROPLibraryDepends DESTINATION ${INSTALL_CMAKE_DIR} COMPONENT dev)

#Install the cmake macros in cmake/Modules directory
file(GLOB macros cmake/Modules/*.cmake)
install(FILES ${macros} DESTINATION ${INSTALL_CMAKE_DIR})

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${GALPROP_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${GALPROP_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${GALPROP_VERSION_PATCH}")
set (CPACK_SOURCE_GENERATOR TBZ2)
#List of files to ignore.  These are actually a list of regular expressions.
#The documentation for this is completely wrong, there is no escape and no /s
set (CPACK_SOURCE_IGNORE_FILES "/galplotdef;/galdef_5[2345];FITS/;\\\\.swp$;galprop_explanatory_supplement;test_source_SNR_event.cc;source_SNR_event_vec.cc;stokes.cc;plot_luminosity.cc;ppfrag.cc;protri.cc;notes.dat;modulate.cc;ProductionMatrices.h;ProductionMatrices.cc;Bremsstrahlung.cc;Bremsstrahlung.h;B_field.pro;cr_luminosity.cc;tests/;tools/")
include (CPack)

