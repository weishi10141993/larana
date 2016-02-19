
# OpFlashFinder_module.cc is not built - it seems specific to uboone
art_make( 
          
          LIB_LIBRARIES 
                       lardata_RawData
                       lardata_RecoBase
                       larcore_Geometry
                       lardata_Utilities
                       ${MF_MESSAGELOGGER}
                       ${Boost_SYSTEM_LIBRARY}
                       ${FHICLCPP}
                       ${CETLIB}
                       ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_source()
