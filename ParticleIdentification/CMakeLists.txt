
art_make( BASENAME_ONLY
	  LIBRARY_NAME PID
	  LIB_LIBRARIES
                        AnalysisBase
                        Utilities
                        ${MF_MESSAGELOGGER}
		        ${MF_UTILITIES}
		        ${FHICLCPP}
		        ${CETLIB}
		        ${CLHEP}
		        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
                        RecoBase
                        Geometry
                        Geometry_service
                        AnalysisAlg
                        AnalysisBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
         )

install_headers()
install_fhicl()
install_source()
