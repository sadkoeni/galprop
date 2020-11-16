include_directories(${PROJECT_SOURCE_DIR}/random)

add_library(random Random.cpp Random.hpp RandomAlgorithm.hpp RandomCanonical.hpp RandomEngine.hpp RandomMixer.hpp RandomNumber.hpp RandomPower2.hpp RandomSeed.hpp RandomSelect.hpp RandomType.hpp NormalDistribution.hpp LeadingZeros.hpp ExponentialProb.hpp ExponentialDistribution.hpp ExactPower.hpp ExactExponential.hpp Config.h )

set_target_properties(random PROPERTIES PUBLIC_HEADER "Random.hpp;RandomAlgorithm.hpp;RandomCanonical.hpp;RandomEngine.hpp;RandomMixer.hpp;RandomNumber.hpp;RandomPower2.hpp;RandomSeed.hpp;RandomSelect.hpp;RandomType.hpp;NormalDistribution.hpp;LeadingZeros.hpp;ExponentialProb.hpp;ExponentialDistribution.hpp;ExactPower.hpp;ExactExponential.hpp;Config.h")

install(TARGETS random EXPORT galtoolslibLibraryDepends 
RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin 
LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib 
ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT stlib
PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)
