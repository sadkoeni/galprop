# Locate the astro library and include files from xephem package
# This module defines
#  ASTRO_LIBRARY       - The fitsio library
#  ASTRO_INCLUDE_DIR   - The include directory
#  ASTRO_FOUND         - True if found

# The search can be assisted by using the variables ASTRO_PREFIX, 
# ASTRO_INC, and ASTRO_LIB.  Note that those are absolute directories.
# Those can also be set as environment variables which then override the cached ones.
# If ASTRO_INC and ASTRO_LIB are in cache, it is not enough to specify
# ASTRO_PREFIX to override the cache, both ASTRO_INC and ASTRO_LIB must be set.

set(ASTRO_PREFIX /usr/local CACHE PATH "Prefix where cfitsio is installed")
set(ASTRO_INC "" CACHE PATH "Path to cfitsio header files [ASTRO_PREFIX/include]")
set(ASTRO_LIB "" CACHE PATH "Path to cfitsio library [ASTRO_PREFIX/lib]")

include(helpfulMacros)
setFromEnvVariable(ASTRO_PREFIX)
setFromEnvVariable(ASTRO_INC)
setFromEnvVariable(ASTRO_LIB)

if (NOT ASTRO_INC)
  set(ASTRO_INC ${ASTRO_PREFIX}/include)
endif(NOT ASTRO_INC)

if (NOT ASTRO_LIB)
  set(ASTRO_LIB ${ASTRO_PREFIX}/lib)
endif(NOT ASTRO_LIB)

find_path(ASTRO_INCLUDE_DIR astro.h PATHS ${ASTRO_INC} NO_DEFAULT_PATH)
find_path(ASTRO_INCLUDE_DIR astro.h PATHS )
find_library(ASTRO_LIBRARY astro PATHS ${ASTRO_LIB} NO_DEFAULT_PATH)
find_library(ASTRO_LIBRARY astro PATHS ${ASTRO_LIB})

if (NOT ASTRO_LIBRARY OR NOT ASTRO_INCLUDE_DIR)
  message(STATUS "ASTRO not found.  Use ASTRO_PREFIX or ASTRO_INC and ASTRO_LIB to specify its location")
endif (NOT ASTRO_LIBRARY OR NOT ASTRO_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set ASTRO_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ASTRO  DEFAULT_MSG  ASTRO_LIBRARY  ASTRO_INCLUDE_DIR)
