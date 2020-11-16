
#Add utils as an include directory
include_directories(${PROJECT_SOURCE_DIR}/utils)

#Set up the library
add_library(geometry )

set_target_properties(geometry PROPERTIES PUBLIC_HEADER "")

target_link_libraries(geometry utils)

install(TARGETS geometry EXPORT galtoolslibLibraryDepends 
   RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin 
   LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib 
   ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT stlib
   PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)
