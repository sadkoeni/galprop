include_directories(${PROJECT_SOURCE_DIR}/utils)

#Add the path to data in definitions
if (IS_ABSOLUTE ${INSTALL_DATA_DIR})
   set(DATA_PATH ${INSTALL_DATA_DIR})
else()
   set(DATA_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_DATA_DIR})
endif()
add_definitions("-DDATA_PATH=\"${DATA_PATH}\"")

add_library(processes Processes_Interface.cc pp_meson.f antiproton.f synchrotron.f bremss_spec.f crn6.f nucleon_cs.cc e_loss_compton.f energy_losses.cc kamae.cc inter.f ap_sergey.f Processes_Interface.h constants.h)

target_link_libraries(processes utils)

set_target_properties(processes PROPERTIES PUBLIC_HEADER "Processes_Interface.h;constants.h")

install(TARGETS processes EXPORT galtoolslibLibraryDepends 
RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin 
LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib 
ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT stlib
PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)

install(FILES barpol.dat ap-table.dat
   DESTINATION ${INSTALL_DATA_DIR})
