# basic source code CMakeLists.txt

art_make( 
	LIB_LIBRARIES
		Geometry
		RecoBase
                AnalysisBase
                ${FHICLCPP}
                ${ROOT_BASIC_LIB_LIST}
                ${CETLIB}
                ${BOOST_SYSTEM} #why the hell do I need this?
)



install_headers()
install_source()
