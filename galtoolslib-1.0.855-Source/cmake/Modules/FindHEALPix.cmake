# - find HEALPix

# HEALPix_INCLUDE_DIR       - Where to find HEALPix header files
# HEALPix_INCLUDE_DIRS      - Where to find HEALPix header files and Cfitsio header files
# HEALPix_LIBRARIES         - HEALPix libraries
# HEALPix_FOUND             - True if HEALPix is found
# HAVE_LSCONSTANTS_H        - True if lsconstants.h is found
# HAVE_ALM_HEALPix_TOOLS_H  - True if alm_healpix_tools.h is found

# HEALPix_PREFIX can be used to point to the healpix distribution, either in cache or as an ENVIRONMENT variable that will override the cache.
# The environment variables HEALPix_INCLUDE_PATH and HEALPix_LIBRARY_PATH can also be set to point to the HEALPix include directory and library directory
# Note that HEALPix_INCLUDE_PATH and HEALPix_LIBRARY_PATH take presedence over the HEALPix_PREFIX

find_package(Cfitsio REQUIRED)

set(HEALPix_PREFIX "/usr/local" CACHE PATH "Location where Healpix is installed")

if(DEFINED ENV{HEALPix_PREFIX})
   set(HEALPix_PREFIX $ENV{HEALPix_PREFIX})
endif()

if(DEFINED ENV{HEALPix_INCLUDE_PATH})
   find_path(HEALPix_INCLUDE_DIR "healpix_base.h" PATHS $ENV{HEALPix_INCLUDE_PATH} NO_DEFAULT_PATH)
endif()
find_path(HEALPix_INCLUDE_DIR "healpix_base.h" PATHS ${HEALPix_PREFIX}/include NO_DEFAULT_PATH)
find_path(HEALPix_INCLUDE_DIR "healpix_base.h")

if(DEFINED ENV{HEALPix_LIBRARY_PATH})
   find_library(HEALPix_healpix_cxx NAMES "healpix_cxx" PATHS $ENV{HEALPix_LIBRARY_PATH} NO_DEFAULT_PATH)
   find_library(HEALPix_fftpack NAMES "fftpack" PATHS $ENV{HEALPix_LIBRARY_PATH} NO_DEFAULT_PATH)
   find_library(HEALPix_cxxsupport NAMES "cxxsupport" PATHS $ENV{HEALPix_LIBRARY_PATH} NO_DEFAULT_PATH)
   find_library(HEALPix_psht NAMES "psht" PATHS $ENV{HEALPix_LIBRARY_PATH} NO_DEFAULT_PATH)
   find_library(HEALPix_c_utils NAMES "c_utils" PATHS $ENV{HEALPix_LIBRARY_PATH} NO_DEFAULT_PATH)
   find_library(HEALPix_sharp NAMES "sharp" PATHS $ENV{HEALPix_LIBRARY_PATH} NO_DEFAULT_PATH)
endif()
find_library(HEALPix_healpix_cxx NAMES "healpix_cxx" PATHS ${HEALPix_PREFIX}/lib NO_DEFAULT_PATH)
find_library(HEALPix_fftpack NAMES "fftpack" PATHS ${HEALPix_PREFIX}/lib NO_DEFAULT_PATH)
find_library(HEALPix_cxxsupport NAMES "cxxsupport" PATHS ${HEALPix_PREFIX}/lib NO_DEFAULT_PATH)
find_library(HEALPix_psht NAMES "psht" PATHS ${HEALPix_PREFIX}/lib NO_DEFAULT_PATH)
find_library(HEALPix_c_utils NAMES "c_utils" PATHS ${HEALPix_PREFIX}/lib NO_DEFAULT_PATH)
find_library(HEALPix_sharp NAMES "sharp" PATHS ${HEALPix_PREFIX}/lib NO_DEFAULT_PATH)
find_library(HEALPix_healpix_cxx NAMES "healpix_cxx")
find_library(HEALPix_fftpack NAMES "fftpack")
find_library(HEALPix_cxxsupport NAMES "cxxsupport")
find_library(HEALPix_psht NAMES "psht")
find_library(HEALPix_c_utils NAMES "c_utils")
find_library(HEALPix_sharp NAMES "sharp")

#Different Versions of HEALPix have different libraries.
set(HEALPix_LIBRARIES ${HEALPix_healpix_cxx})
if (HEALPix_healpix_cxxsupport)
   set(HEALPix_LIBRARIES ${HEALPix_healpix_cxx} ${HEALPix_fftpack} ${HEALPix_cxxsupport})
endif()
if (HEALPix_psht)
   set(HEALPix_LIBRARIES ${HEALPix_healpix_cxx} ${HEALPix_psht} ${HEALPix_c_utils} ${HEALPix_fftpack} ${HEALPix_cxxsupport})
endif()
if (HEALPix_sharp)
   set(HEALPix_LIBRARIES ${HEALPix_healpix_cxx} ${HEALPix_sharp} ${HEALPix_c_utils} ${HEALPix_fftpack} ${HEALPix_cxxsupport})
endif()

#Check for a few header files and add the Cfitsio include directories
if(HEALPix_INCLUDE_DIR)
   set(HEALPix_INCLUDE_DIRS ${HEALPix_INCLUDE_DIR} ${Cfitsio_INCLUDE_DIR})
   include(CheckIncludeFileCXX)

   set(CMAKE_REQUIRED_INCLUDES ${HEALPix_INCLUDE_DIR})

   CHECK_INCLUDE_FILE_CXX("alm_healpix_tools.h" HAVE_ALM_HEALPIX_TOOLS_H)
   CHECK_INCLUDE_FILE_CXX("lsconstants.h" HAVE_LSCONSTANTS_H)
   CHECK_INCLUDE_FILE_CXX("alm_powspec_tools.h" HAVE_ALM_POWSPEC_TOOLS_H)

   #Version check, version 3 has this file
   CHECK_INCLUDE_FILE_CXX("alloc_utils.h" HAVE_ALLOC_UTILS_H)
   if (HAVE_ALLOC_UTILS_H)
      add_definitions(-DHEALPIX_V3)
      #Add a check for std::complex
      include(CheckCXXSourceCompiles)
      set(CMAKE_REQUIRED_INCLUDES ${HEALPix_INCLUDE_DIRS})
      check_cxx_source_compiles("#include <xcomplex.h>\n int main() { dcomplex var; var.imag(); return 0;}" HEALPIX_STD_COMPLEX)
   else()
      add_definitions(-DHEALPIX_V2)
   endif()
endif()


#Now everything should be ready
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set HEALPix_FOUND to TRUE
# if all listed variablesare TRUE
set(errorMsg "Could NOT find HEALPix.  Use the HEALPix_PREFIX cache variable or HEALPix_INCLUDE_PATH and HEALPix_LIBRARY_PATH environment variables to point to the library.")
find_package_handle_standard_args(HEALPix ${errorMsg} HEALPix_LIBRARIES HEALPix_INCLUDE_DIRS)
