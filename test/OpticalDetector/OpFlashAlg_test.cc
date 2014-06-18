#define BOOST_TEST_MODULE ( OpFlashAlg_test )
#include "boost/test/auto_unit_test.hpp"

#include "OpticalDetector/OpFlashAlg.h"
#include "OpticalDetectorData/OpticalTypes.h"

const float HitThreshold = 3;
const int Channel = 1;
const uint32_t TimeSlice = 1000;
const unsigned short Frame = 1;
const optdata::TimeSlice_t TimeSlicesPerFrame = 102400;
const double opdigi_SampleFreq = 64;
const double TrigTimeAbs = 0.8;
const double SPESize = 20;

const double tolerance = 1e-6;

BOOST_AUTO_TEST_SUITE(OpFlashAlg_test)

BOOST_AUTO_TEST_CASE(ConstructHit_checkDefaultPulse)
{
  pmtana::pulse_param pulse;
  std::vector<recob::OpHit> HitVector;

  opdet::ConstructHit(HitThreshold,
		      Channel,
		      TimeSlice,
		      Frame,
		      &pulse,
		      TimeSlicesPerFrame,
		      opdigi_SampleFreq,
		      TrigTimeAbs,
		      SPESize,
		      HitVector);
  
  BOOST_CHECK_EQUAL(HitVector.size(),0ul);

}

BOOST_AUTO_TEST_CASE(ConstructHit_checkPulseBelowThreshold)
{
  pmtana::pulse_param pulse;
  pulse.t_start=3; pulse.t_max = 6; pulse.t_end = 9;
  pulse.area = 4; pulse.peak = 2;

  std::vector<recob::OpHit> HitVector;

  opdet::ConstructHit(HitThreshold,
		      Channel,
		      TimeSlice,
		      Frame,
		      &pulse,
		      TimeSlicesPerFrame,
		      opdigi_SampleFreq,
		      TrigTimeAbs,
		      SPESize,
		      HitVector);
  
  BOOST_CHECK_EQUAL(HitVector.size(),0ul);
}

BOOST_AUTO_TEST_CASE(ConstructHit_checkPulseAbovThreshold)
{



  pmtana::pulse_param pulse;
  pulse.t_start=3; pulse.t_max = 6; pulse.t_end = 9;
  pulse.area = 40; pulse.peak = 20;

  std::vector<recob::OpHit> HitVector;

  opdet::ConstructHit(HitThreshold,
		      Channel,
		      TimeSlice,
		      Frame,
		      &pulse,
		      TimeSlicesPerFrame,
		      opdigi_SampleFreq,
		      TrigTimeAbs,
		      SPESize,
		      HitVector);
  
  double peak_time_abs = (pulse.t_max + TimeSlice + Frame*TimeSlicesPerFrame)/opdigi_SampleFreq;

  BOOST_CHECK_EQUAL(HitVector.size(),1ul);
  BOOST_CHECK_EQUAL(HitVector[0].OpChannel(),Channel);
  BOOST_CHECK_CLOSE(HitVector[0].PeakTimeAbs(),peak_time_abs,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].PeakTime(),peak_time_abs-TrigTimeAbs,tolerance);
  BOOST_CHECK_EQUAL(HitVector[0].Frame(),Frame);
  BOOST_CHECK_CLOSE(HitVector[0].Width(),pulse.t_end-pulse.t_start,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].Area(),pulse.area,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].Amplitude(),pulse.peak,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].PE(),pulse.peak/SPESize,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].FastToTotal(),0.,tolerance);

}

BOOST_AUTO_TEST_SUITE_END()
