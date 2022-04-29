add_subdirectory(AssociationsTools)

cet_build_plugin(MCParticleHitMatching art::EDProducer
  LIBRARIES PRIVATE
  larana::HitParticleAssociationsTool
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art_plugin_support::toolMaker
  fhiclcpp::fhiclcpp
)

cet_build_plugin(MCParticleShowerMatching art::EDProducer
  LIBRARIES PRIVATE
  lardataobj::AnalysisBase
  nusimdata::SimulationBase
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(MCParticleTrackMatching art::EDProducer
  LIBRARIES PRIVATE
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(MCTruthT0Matching art::EDProducer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::AssociationUtil
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(PhotonCounterT0Matching art::EDProducer
  LIBRARIES PRIVATE
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
