#Set up the library
add_library(galprop B_field_model.cc Configure.cc D_pp.cc D_xx.cc Distribution.cc 
   GalacticRadiationField.cc RadiationField.cc Galaxy.cc Galdef.cc Galprop.cc 
   Particle.cc  
   blattnig.cc cfactor.f create_SNR.cc create_galaxy.cc
   create_gcr.cc fill_transport_arrays.cc decayed_cross_sections.cc 
   e_KN_loss.cc electrons_normalize.cc 
   eprop.cc gauss.cc 
   gen_DM_source.cc gen_IC_emiss.cc gen_bremss_emiss.cc 
   gen_ionization_rate.cc 
   gen_isrf_energy_density.cc gen_knock_on_electron_source.cc gen_pi0_decay_emiss.cc 
   gen_secondary_antiproton_source.cc 
   gen_secondary_positron_source.cc gen_secondary_proton_source.cc 
   gen_secondary_pair_source.cc gen_secondary_source.cc gen_synch_emiss.cc 
   gen_tertiary_antiproton_source.cc integ.cc 
   ionization_bethe.cc isrf_energy_density.cc kinematic.cc 
   nH.cc nuclei_normalize.cc print_BC.cc 
   propagate_particles.cc propel.cc propel_diagnostics.cc read_gcr.cc 
   read_isrf.cc read_gas_maps.cc source_SNR_event.cc 
   source_distribution.cc store_DM_emiss.cc store_DM_skymap.cc 
   store_IC_skymap.cc store_IC_skymap_comp.cc store_bremss_H2R_skymap.cc 
   store_bremss_HIR_skymap.cc store_bremss_emiss.cc store_bremss_ionized_skymap.cc 
   store_bremss_skymap.cc store_bremss_HII_skymap.cc store_pi0_decay_HII_skymap.cc 
   store_gcr.cc store_gcr_full.cc store_ionization_rate.cc store_mapcube_skymap.cc 
   store_pi0_decay_H2R_skymap.cc store_pi0_decay_HIR_skymap.cc store_IC_emiss.cc 
   store_pi0_decay_emiss.cc store_pi0_decay_skymap.cc store_skymap.cc 
   store_synch_skymap.cc test_Distribution.cc test_Particle.cc 
   test_float_accuracy.cc test_isotope_cs.cc test_nH.cc 
   tridag.cc tridag_ext.cc tridag_sym.cc tridag_sym_ext.cc  
   B_field_3D_model.cc synchrotron_emissivity.cc 
   synchrotron_emissivity_B_field.cc synchrotron_emissivity_aws.cc test_sync_package.cc 
   store_synch_emiss.cc galprop_internal.cc gen_luminosity.cc store_luminosity.cc 
   store_gcr_source_functions.cc DistributionFunction.cc GasFunction.cc 
   GasEmissFunction.cc los_integration.cc AvXCO.cc los_intersection.cc 
   gen_skymaps.cc kappa_free_free.cc SourceClass.cc SourceClass_Compatibility.cc
   frag.f SpectralDistribution.cc PowerLawSD.cc SpatialDistribution.cc
   GalstructXMLDistribution.cc SourceClass_UniformSpectra.cc
   SmoothPowerLawSD.cc PropelOperatorSplitting.cc
   Configure.h Distribution.h Galaxy.h Galdef.h Galprop.h Particle.h Spectrum.h 
   galprop_internal.h 
   galprop_classes.h GalacticRadiationField.h RadiationField.h 
   synchrotron_emissivity_aws.h synchrotron_emissivity_B_field.h 
   synchrotron_emissivity.h B_field_3D_model.h DistributionFunction.h 
   los_integration.h integ.h ArmData.h StellarGeometryDistributions.h  
   los_intersection.h kappa_free_free.h SourceClass.h SourceClass_Compatibility.h
   SpectralDistribution.h PowerLawSD.h SpatialDistribution.h GalstructXMLDistribution.h
   SourceClass_UniformSpectra.h SmoothPowerLawSD.h PropelBase.h PropelOperatorSplitting.h
   )

set_target_properties(galprop PROPERTIES PUBLIC_HEADER "Configure.h;Distribution.h;Galaxy.h;Galdef.h;Galprop.h;Particle.h;Spectrum.h;galprop_internal.h;galprop_classes.h;GalacticRadiationField.h;RadiationField.h;synchrotron_emissivity_aws.h;synchrotron_emissivity_B_field.h;synchrotron_emissivity.h;B_field_3D_model.h;DistributionFunction.h;los_integration.h;integ.h;ArmData.h;StellarGeometryDistributions.h;los_intersection.h;kappa_free_free.h;SourceClass.h;SourceClass_Compatibility.h;SpectralDistribution.h;PowerLawSD.h;SmoothPowerLawSD.h;SpatialDistribution.h;GalstructXMLDistribution.h;SourceClass_UniformSpectra.h")

target_link_libraries(galprop processes galstruct ${galtoolslib_LIBRARIES} ${GSL_LIBRARIES} ${CLHEP_LIBRARIES})

add_executable(galpropapp galpropMain.cc)
target_link_libraries(galpropapp galprop)

add_executable(changeGasMaps changeGasMaps.cc)
target_link_libraries(changeGasMaps galprop)


SET_TARGET_PROPERTIES(galpropapp PROPERTIES OUTPUT_NAME galprop)

install(TARGETS changeGasMaps galprop galpropapp EXPORT GALPROPLibraryDepends 
   RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin 
   LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib 
   ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT stlib 
   PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev)

#Install the data files
install(FILES j0zero500.dat 
   DESTINATION ${INSTALL_DATA_DIR})

#Run the post install script to create a link to the average X_CO program
install(SCRIPT postInstallScript.cmake)
