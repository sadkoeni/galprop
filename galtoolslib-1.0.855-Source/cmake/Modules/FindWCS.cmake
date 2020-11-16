# Locate the astro library and include files from xephem package
# This module defines
#  WCS_LIBRARY       - The fitsio library
#  WCS_INCLUDE_DIR   - The include directory
#  WCS_FOUND         - True if found

# The search can be assisted by using the variables WCS_PREFIX, 
# WCS_INC, and WCS_LIB.  Note that those are absolute directories.
# Those can also be set as environment variables which then override the cached ones.
# If WCS_INC and WCS_LIB are in cache, it is not enough to specify
# WCS_PREFIX to override the cache, both WCS_INC and WCS_LIB must be set.

set(WCS_PREFIX /usr/local CACHE PATH "Prefix where cfitsio is installed")
set(WCS_INC "" CACHE PATH "Path to cfitsio header files [WCS_PREFIX/include]")
set(WCS_LIB "" CACHE PATH "Path to cfitsio library [WCS_PREFIX/lib]")

include(helpfulMacros)
setFromEnvVariable(WCS_PREFIX)
setFromEnvVariable(WCS_INC)
setFromEnvVariable(WCS_LIB)

if (NOT WCS_INC)
  set(WCS_INC ${WCS_PREFIX}/include)
endif(NOT WCS_INC)

if (NOT WCS_LIB)
  set(WCS_LIB ${WCS_PREFIX}/lib)
endif(NOT WCS_LIB)

find_path(WCS_INCLUDE_DIR wcslib/wcs.h PATHS ${WCS_INC} NO_DEFAULT_PATH)
find_path(WCS_INCLUDE_DIR wcslib/wcs.h PATHS )
find_library(WCS_LIBRARY wcs PATHS ${WCS_LIB} NO_DEFAULT_PATH)
find_library(WCS_LIBRARY wcs PATHS ${WCS_LIB})

if (NOT WCS_LIBRARY OR NOT WCS_INCLUDE_DIR)
  message(STATUS "WCS not found.  Use WCS_PREFIX or WCS_INC and WCS_LIB to specify its location")
endif (NOT WCS_LIBRARY OR NOT WCS_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set WCS_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(WCS  DEFAULT_MSG  WCS_LIBRARY  WCS_INCLUDE_DIR)
