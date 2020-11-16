
#Add utils as an include directory
include_directories(${PROJECT_SOURCE_DIR}/utils)

#Set up the library
add_library(skymap Coordinate.cc HealpixBaseExtended.cc Region.cc SkySelection.cc Skymap.cc BaseSky.cc BaseSkyFitsIO.cc Skymap.h Coordinate.h HealpixBaseExtended.h Region.h SkySelection.h ArraySlice.h SparseSkymap.h BaseSky.h FullSky.h SparseSky.h BaseSkyFitsIO.h)

set_target_properties(skymap PROPERTIES PUBLIC_HEADER "Skymap.h;SparseSkymap.h;Coordinate.h;HealpixBaseExtended.h;Region.h;SkySelection.h;ArraySlice.h;BaseSky.h;FullSky.h;SparseSky.h;BaseSkyFitsIO.h")

target_link_libraries(skymap utils ${HEALPix_LIBRARIES} ${CCfits_LIBRARIES} ${WCS_LIBRARY})

install(TARGETS skymap EXPORT galtoolslibLibraryDepends 
   RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin 
   LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib 
   ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT stlib
   PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)
