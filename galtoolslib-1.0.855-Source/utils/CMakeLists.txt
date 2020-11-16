include_directories(${PROJECT_SOURCE_DIR}/utils)

add_library(utils ErrorLogger.cc Timer.cc Parameters.cc StatusIndicator.cc Interpolation.cc ReaderErrorReporter.cc Reader.cc Variables.cc md5.cc Parameters.h PhysicalConstants.h Units.h ErrorLogger.h Singleton.h Timer.h StatusIndicator.h Interpolation.h ReaderErrorReporter.h Reader.h ShadowPtr.h ShadowPtr_fwd.h Variables.h Registry.h md5.h Malloc.cc Malloc.h)

set_target_properties(utils PROPERTIES PUBLIC_HEADER "PhysicalConstants.h;Units.h;ErrorLogger.h;Singleton.h;Timer.h;Parameters.h;StatusIndicator.h;Interpolation.h;Reader.h;ReaderErrorReporter.h;ShadowPtr.h;ShadowPtr_fwd.h;Variables.h;Registry.h;md5.h;Malloc.h")

target_link_libraries(utils ${XERCESC_LIBRARY} ${CLHEP_LIBRARIES})

install(TARGETS utils EXPORT galtoolslibLibraryDepends 
RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin 
LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib 
ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT stlib
PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)
