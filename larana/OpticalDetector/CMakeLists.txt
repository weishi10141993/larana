
# larana_OpticalDetector_OpFlashFinder_module.cc is not built - it seems specific to uboone

art_make_library( LIBRARY_NAME FlashHypothesis
		  	     SOURCE FlashHypothesis.cxx
		                           FlashHypothesisCalculator.cxx
		                           FlashUtilities.cxx
			     LIBRARIES
			                             ${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME SimPhotonCounter
		  	     SOURCE SimPhotonCounter.cxx
			     LIBRARIES
			                ${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME FlashHypothesisComparison
		  	     SOURCE FlashHypothesisComparison.cxx
			     LIBRARIES
					FlashHypothesis
					SimPhotonCounter
					MCSTReco
			                ${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME OpFlashAnaAlg
		  	     SOURCE OpFlashAnaAlg.cxx
			     LIBRARIES
					RecoBase
			                ${ROOT_BASIC_LIB_LIST}
)


art_make( 
          EXCLUDE FlashHypothesis.cxx
			 FlashHypothesisCalculator.cxx
		         FlashUtilities.cxx
		         FlashHypothesisComparison.cxx
			 SimPhotonCounter.cxx
			 OpFlashAnaAlg.cxx
          LIBRARY_NAME OpticalDetector
          LIB_LIBRARIES 
			RecoBase
	  		FlashHypothesis
			SimPhotonCounter
	  		FlashHypothesisComparison
			OpFlashAnaAlg
			MCBase
                        Geometry
                        larcore_Geometry_Geometry_service
                        larsim_PhotonPropagation_PhotonVisibilityService_service
		        larana_OpticalDetector_OpDigiProperties_service
		        Utilities
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES larana_OpticalDetector_OpDigiProperties_service
	                OpticalDetector
                        larsim_PhotonPropagation_PhotonVisibilityService_service
			larsim_MCCheater_BackTracker_service
                        Filters
                        RecoObjects
                        RecoBase
			MCBase
                        RecoAlg
                        Geometry
                        larcore_Geometry_Geometry_service
                        Simulation
                        larsim_TriggerAlgo_TriggerAlgoMicroBoone_service
                        Utilities
                        RawData
                        AnalysisAlg
                        AnalysisBase
                        SeedService_service # artextensions
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
       SERVICE_LIBRARIES
                        Geometry
                        larcore_Geometry_Geometry_service
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
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
install_fw( LIST toyWaveform.txt
            SUBDIRNAME OpticalDetector )
