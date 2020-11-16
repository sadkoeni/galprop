
find_package(CppUnit)

include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/galstruct)
include_directories(${PROJECT_SOURCE_DIR}/libskymap)

#Unit tests if found
if (CPPUNIT_FOUND)
  file(GLOB cppUnitTest_cc "*.cc")
  file(GLOB cppUnitTest_h "*.h")
  
  include_directories(${CPPUNIT_INCLUDE_DIR})
  
  add_executable(cppUnitTest ${cppUnitTest_cc} ${cppUnitTest_h})
  target_link_libraries(cppUnitTest utils)
  target_link_libraries(cppUnitTest galstruct)
  target_link_libraries(cppUnitTest skymap)
  target_link_libraries(cppUnitTest ${CPPUNIT_LIBRARY})
  
  add_test(unitTest cppUnitTest)
endif (CPPUNIT_FOUND)

