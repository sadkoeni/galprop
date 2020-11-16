# Locate the cfitsio library and include files
# This module defines
#  Cfitsio_LIBRARY       - The fitsio library
#  Cfitsio_INCLUDE_DIR   - The include directory
#  CFITSIO_FOUND         - True if found

# The search can be assisted by using the variables Cfitsio_PREFIX, 
# Cfitsio_INC, and Cfitsio_LIB.  Note that those are absolute directories.
# Those can also be set as environment variables which then override the cached ones.
# If Cfitsio_INC and Cfitsio_LIB are in cache, it is not enough to specify
# Cfitsio_PREFIX to override the cache, both Cfitsio_INC and Cfitsio_LIB must be set.

set(Cfitsio_PREFIX /usr/local CACHE PATH "Prefix where cfitsio is installed")
set(Cfitsio_INC "" CACHE PATH "Path to cfitsio header files [Cfitsio_PREFIX/include]")
set(Cfitsio_LIB "" CACHE PATH "Path to cfitsio library [Cfitsio_PREFIX/lib]")

include(helpfulMacros)
setFromEnvVariable(Cfitsio_PREFIX)
setFromEnvVariable(Cfitsio_INC)
setFromEnvVariable(Cfitsio_LIB)

if (NOT Cfitsio_INC)
  set(Cfitsio_INC ${Cfitsio_PREFIX}/include)
endif(NOT Cfitsio_INC)

if (NOT Cfitsio_LIB)
  set(Cfitsio_LIB ${Cfitsio_PREFIX}/lib)
endif(NOT Cfitsio_LIB)

find_path(Cfitsio_INCLUDE_DIR fitsio.h PATHS ${Cfitsio_INC} NO_DEFAULT_PATH)
find_path(Cfitsio_INCLUDE_DIR fitsio.h PATHS /usr/include /usr/local/include /usr/include/libcfitsio0)
find_library(Cfitsio_LIBRARY cfitsio PATHS ${Cfitsio_LIB} NO_DEFAULT_PATH)
find_library(Cfitsio_LIBRARY cfitsio PATHS ${Cfitsio_LIB})

if (NOT Cfitsio_LIBRARY OR NOT Cfitsio_INCLUDE_DIR)
  message(STATUS "Cfitsio not found.  Use Cfitsio_PREFIX or Cfitsio_INC and Cfitsio_LIB to specify its location")
endif (NOT Cfitsio_LIBRARY OR NOT Cfitsio_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set Cfitsio_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cfitsio  DEFAULT_MSG  Cfitsio_LIBRARY  Cfitsio_INCLUDE_DIR)
