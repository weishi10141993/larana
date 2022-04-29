cet_enable_asserts()

cet_make_library(LIBRARY_NAME HitParticleAssociationsTool INTERFACE
  SOURCE IHitParticleAssociations.h
  LIBRARIES INTERFACE
  lardataobj::AnalysisBase
  canvas::canvas
)

cet_write_plugin_builder(lar::HitParticleAssociationsTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larana::HitParticleAssociationsTool
)

include(lar::HitParticleAssociationsTool)

cet_build_plugin(DirectHitParticleAssns lar::HitParticleAssociationsTool
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(IndirectHitParticleAssns lar::HitParticleAssociationsTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_headers()
install_source()
install_fhicl()
