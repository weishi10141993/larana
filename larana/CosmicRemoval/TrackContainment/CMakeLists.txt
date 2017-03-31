# basic source code CMakeLists.txt

art_make( 
	LIB_LIBRARIES
		larcore_Geometry
		lardataobj_RecoBase
                lardataobj_AnalysisBase
                ${FHICLCPP}
                ${ROOT_BASIC_LIB_LIST}
                cetlib cetlib_except
                ${BOOST_SYSTEM} #why the hell do I need this?
)



install_headers()
install_source()
