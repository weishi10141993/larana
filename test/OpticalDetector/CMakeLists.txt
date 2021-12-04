# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

find_package(Boost REQUIRED)

cet_test(OpFlashAlg_test USE_BOOST_UNIT
                         LIBRARIES larana_OpticalDetector
)
