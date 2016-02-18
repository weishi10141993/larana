
# OpFlashFinder_module.cc is not built - it seems specific to uboone
art_make( BASENAME_ONLY
          LIBRARY_NAME OpHitFinder
          LIB_LIBRARIES 
                       RawData
                       RecoBase
                       Geometry
                       Utilities
                       ${MF_MESSAGELOGGER}
                       ${Boost_SYSTEM_LIBRARY}
                       ${FHICLCPP}
                       ${CETLIB}
                       ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_source()
