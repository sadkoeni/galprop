# Locate the CCfits library and include files
# This module requires cfitsio and defines
#  CCfits_LIBRARIES
#  CCfits_INCLUDE_DIRS
#  CCfits_FOUND
# if it finds  cfitsio. It also defines
#  CCfits_INCLUDE_DIR
#  CCfits_LIBRARY
# which contains the ccfits include only

# To assist in the search use the variables CCfits_PREFIX, CCfits_INC, and CCfits_LIB.
# Those can also be set as environment variables which take presedence over the cached ones.
# Note that if CCfits_INC and CCfits_LIB are set in cache, then it is not enough
# to specify only CCfits_PREFIX as an environment variable

find_package(Cfitsio REQUIRED)

if (CFITSIO_FOUND)
  set(CCfits_PREFIX /usr/local CACHE PATH "Prefix where CCfits is installed")
  set(CCfits_INC "" CACHE PATH "CCfits include directory [CCfits_PREFIX/include]")
  set(CCfits_LIB "" CACHE PATH "CCfits library directory [CCfits_PREFIX/lib]")

  include(helpfulMacros)
  setFromEnvVariable(CCfits_PREFIX)
  setFromEnvVariable(CCfits_INC)
  setFromEnvVariable(CCfits_LIB)

  if (NOT CCfits_INC)
    set(CCfits_INC ${CCfits_PREFIX}/include)
  endif(NOT CCfits_INC)
  if (NOT CCfits_LIB)
    set(CCfits_LIB ${CCfits_PREFIX}/lib)
  endif(NOT CCfits_LIB)
  
  find_path(CCfits_INCLUDE_DIR CCfits/CCfits PATHS ${CCfits_INC} NO_DEFAULT_PATH)
  find_path(CCfits_INCLUDE_DIR CCfits/CCfits)
  find_library(CCfits_LIBRARY CCfits PATHS ${CCfits_LIB} NO_DEFAULT_PATH)
  find_library(CCfits_LIBRARY CCfits)

  if (CCfits_INCLUDE_DIR)
     set(CCfits_INCLUDE_DIRS ${CCfits_INCLUDE_DIR} ${Cfitsio_INCLUDE_DIR})
  endif (CCfits_INCLUDE_DIR)
  if (CCfits_LIBRARY)
    set(CCfits_LIBRARIES ${CCfits_LIBRARY} ${Cfitsio_LIBRARY})
  endif (CCfits_LIBRARY)
endif (CFITSIO_FOUND)

if (NOT CCfits_INCLUDE_DIRS OR NOT CCfits_LIBRARIES)
  message(STATUS "CCfits not found.  Use CCfits_PREFIX or CCfits_INC and CCfits_LIB to specify its location")
endif ()

# handle the QUIETLY and REQUIRED arguments and set CCfits_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CCfits  DEFAULT_MSG  CCfits_LIBRARIES CCfits_INCLUDE_DIRS)

mark_as_advanced(CCfits_INCLUDE_DIR CCfits_LIBRARY)
