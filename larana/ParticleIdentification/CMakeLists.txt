include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( 
	  
	  LIB_LIBRARIES
			nusimdata_SimulationBase
                        lardataobj_AnalysisBase
			larreco_Calorimetry
                        lardata_Utilities
                        lardataobj_RecoBase
                        larcore_Geometry
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
		        ${MF_UTILITIES}
		        ${FHICLCPP}
		        cetlib cetlib_except
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
                        larreco_Calorimetry
                        lardataobj_AnalysisBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			cetlib cetlib_except
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)
