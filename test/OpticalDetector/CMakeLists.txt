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

#cet_test(standalone_test)
