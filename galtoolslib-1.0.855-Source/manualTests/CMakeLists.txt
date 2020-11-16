
include_directories(${PROJECT_SOURCE_DIR}/utils ${PROJECT_SOURCE_DIR}/galstruct ${PROJECT_SOURCE_DIR}/libskymap)

add_executable(testStatusIndicator testStatusIndicator.cc)
target_link_libraries(testStatusIndicator utils)
  
add_executable(armTests armTests.cpp)
target_link_libraries(armTests utils galstruct)

add_executable(testBaseSkyFitsIO testBaseSkyFitsIO.cc)
target_link_libraries(testBaseSkyFitsIO skymap)

add_executable(convertWCStoFullSky convertWCStoFullSky.cc)
target_link_libraries(convertWCStoFullSky skymap)

add_executable(convertToFullSky convertToFullSky.cc)
target_link_libraries(convertToFullSky skymap)

add_executable(testInterpolationAndConversion testInterpolationAndConversion.cc)
target_link_libraries(testInterpolationAndConversion skymap)
