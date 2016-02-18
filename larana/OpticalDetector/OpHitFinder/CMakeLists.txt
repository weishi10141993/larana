
# OpFlashFinder_module.cc is not built - it seems specific to uboone
art_make( 
          
          LIB_LIBRARIES 
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_source()
