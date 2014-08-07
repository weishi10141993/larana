# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

cet_test(AlgoThreshold_test USE_BOOST_UNIT
			    LIBRARIES OpticalDetector
					     ${FHICLCPP}
)

cet_test(OpFlashAlg_test USE_BOOST_UNIT
			 LIBRARIES OpticalDetector
)

#cet_test(standalone_test)
