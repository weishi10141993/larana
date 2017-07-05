
# OpFlashFinder_module.cc is not built - it seems specific to uboone
art_make( 
          
          LIB_LIBRARIES 
                       lardataobj_RawData
                       lardataobj_RecoBase
                       larcorealg_Geometry
                       lardata_Utilities
                       ${MF_MESSAGELOGGER}
                       ${Boost_SYSTEM_LIBRARY}
                       ${FHICLCPP}
                       cetlib cetlib_except
                       ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_source()
