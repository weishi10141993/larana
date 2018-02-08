include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( 

  LIB_LIBRARIES
  lardataobj_AnalysisBase
  larreco_Calorimetry
  lardata_Utilities
  lardataobj_RecoBase
  larcorealg_Geometry
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Persistency_Provenance
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  cetlib cetlib_except
  ${CLHEP}
  ${ROOT_BASIC_LIB_LIST}
  ${ROOT_GEOM}
  ${Boost_SYSTEM_LIBRARY}
  MODULE_LIBRARIES 
  nusimdata_SimulationBase
  larana_ParticleIdentification
  lardataobj_RecoBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larreco_Calorimetry
  lardataobj_AnalysisBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  cetlib cetlib_except
  )

install_headers()
install_fhicl()
install_source()

add_subdirectory(Algorithm)
