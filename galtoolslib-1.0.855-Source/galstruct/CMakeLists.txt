
#Add utils as an include directory
include_directories(${PROJECT_SOURCE_DIR}/utils)

#Set up the library
add_library(galstruct radialprofiles.cpp spiralarms.cpp planeprofiles.cpp cylindricalprofiles.cpp radialprofiles.h planeprofiles.h cylindricalprofiles.h spiralarms.h)

set_target_properties(galstruct PROPERTIES PUBLIC_HEADER "radialprofiles.h;planeprofiles.h;spiralarms.h;cylindricalprofiles.h")

target_link_libraries(galstruct utils ${GSL_LIBRARIES} ${XERCES_LIBRARIES} ${CLHEP_LIBRARIES} ${OpenCL_LIBRARIES})

install(TARGETS galstruct EXPORT galtoolslibLibraryDepends 
   RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin 
   LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib 
   ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT stlib
   PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)
