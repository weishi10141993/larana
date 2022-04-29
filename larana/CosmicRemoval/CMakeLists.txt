add_subdirectory(TrackContainment)

cet_make_library(SOURCE
  BeamFlashTrackMatchTaggerAlg.cxx
  HitTagAssociatorAlg.cxx
  LIBRARIES
  PUBLIC
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  larcorealg::headers
  nusimdata::SimulationBase
  PRIVATE
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larcore::Geometry_Geometry_service
  lardata::LArPropertiesService
  lardataalg::DetectorInfo
  larcorealg::Geometry
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Physics
  ROOT::Tree
)

cet_build_plugin(BeamFlashTrackMatchTagger art::EDProducer
  LIBRARIES PRIVATE
  larana::CosmicRemoval
  larana::OpticalDetector_OpDigiProperties_service
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larcore::Geometry_Geometry_service
  lardata::LArPropertiesService
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
  fhiclcpp::fhiclcpp
)

cet_build_plugin(CRHitRemovalByPCA art::EDProducer
  LIBRARIES PRIVATE
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
)

cet_build_plugin(CRHitRemoval art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  messagefacility::MF_MessageLogger
)

cet_build_plugin(CosmicClusterTagger art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
)

cet_build_plugin(CosmicPCAxisTagger art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg_Cluster3DAlgs
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  ROOT::Physics
)

cet_build_plugin(CosmicPFParticleTagger art::EDProducer
  LIBRARIES PRIVATE
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  ROOT::Physics
  
)

cet_build_plugin(CosmicRemovalAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(CosmicTrackTagger art::EDProducer
  LIBRARIES PRIVATE
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  ROOT::Core
  canvas::canvas
)

cet_build_plugin(TrackContainmentAnalyzer art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::CosmicRemoval_TrackContainment
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  art::Framework_Principal
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(TrackContainmentTagger art::EDProducer
  LIBRARIES PRIVATE
  larana::CosmicRemoval_TrackContainment
  lardata::AssociationUtil
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  art::Framework_Principal
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(TrackPFParticleMatch art::EDProducer
  LIBRARIES PRIVATE
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
)

install_headers()
install_fhicl()
install_source()
