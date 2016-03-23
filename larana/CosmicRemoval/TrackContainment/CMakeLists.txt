# basic source code CMakeLists.txt

art_make( 
	LIB_LIBRARIES
		larcore_Geometry
		lardata_RecoBase
                lardata_AnalysisBase
                ${FHICLCPP}
                ${ROOT_BASIC_LIB_LIST}
                ${CETLIB}
                ${BOOST_SYSTEM} #why the hell do I need this?
)



install_headers()
install_source()
