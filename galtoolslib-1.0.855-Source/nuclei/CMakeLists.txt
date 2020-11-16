include_directories(${PROJECT_SOURCE_DIR}/utils)

#Add the path to data in definitions
if (IS_ABSOLUTE ${INSTALL_DATA_DIR})
   set(DATA_PATH ${INSTALL_DATA_DIR})
else()
   set(DATA_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_DATA_DIR})
endif()
add_definitions("-DDATA_PATH=\"${DATA_PATH}\"")

add_library(nuclei Nuclei_Interface.cc nuc_package.cc sigma_boron_dec_heinbach_simon.cc He_to_H_CS.cc WNEWTR_FUNC_aws.f YIELDX_011000_imos.f Kcapture_cs.cc Nuclei_Interface.h constants.h)

target_link_libraries(nuclei utils)

set_target_properties(nuclei PROPERTIES PUBLIC_HEADER "Nuclei_Interface.h")

install(TARGETS nuclei EXPORT galtoolslibLibraryDepends 
RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin 
LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib 
ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT stlib
PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)

install(FILES eval_iso_cs.dat isotope_cs.dat p_cs_fits.dat nucdata.dat WNEWTR_082693.CDR.dat
   DESTINATION ${INSTALL_DATA_DIR})
