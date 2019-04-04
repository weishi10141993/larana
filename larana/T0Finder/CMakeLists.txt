add_subdirectory(AssociationsTools)

art_make(MODULE_LIBRARIES lardataobj_RawData
                          lardataobj_AnalysisBase
                          lardataobj_RecoBase
                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities
                          larsim_MCCheater_BackTrackerService_service
                          larsim_MCCheater_ParticleInventoryService_service
                          nusimdata_SimulationBase
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                          art_Persistency_Common
                          art_Persistency_Provenance
                          art_Utilities
                          canvas
                          ${MF_MESSAGELOGGER}
                          
                          ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
install_scripts()
