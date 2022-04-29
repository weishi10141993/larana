cet_make_library(SOURCE TrackContainmentAlg.cxx
  LIBRARIES
  PUBLIC
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  PRIVATE
  larcorealg::Geometry
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

install_headers()
install_source()
