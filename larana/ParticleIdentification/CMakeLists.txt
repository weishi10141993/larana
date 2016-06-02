include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( 
	  
	  LIB_LIBRARIES
			nusimdata_SimulationBase
                        lardataobj_AnalysisBase
			lardata_AnalysisAlg
                        lardata_Utilities
                        larcore_Geometry
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
		        ${MF_UTILITIES}
		        ${FHICLCPP}
		        ${CETLIB}
		        ${CLHEP}
		        ${ROOT_BASIC_LIB_LIST}
			${ROOT_TMVA}
                        ${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES 
			nusimdata_SimulationBase
	  		larana_ParticleIdentification
                        lardataobj_RecoBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_AnalysisAlg
                        lardataobj_AnalysisBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)
