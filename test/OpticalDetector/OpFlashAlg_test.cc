#define BOOST_TEST_MODULE ( OpFlashAlg_test )
#include "boost/test/auto_unit_test.hpp"

#include "OpticalDetector/OpFlashAlg.h"
#include "OpticalDetectorData/OpticalTypes.h"

const float HitThreshold = 3;
const float FlashThreshold = 50;
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

BOOST_AUTO_TEST_CASE(checkGetAccumIndex)
{

  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,0,1,0),10ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,0,2,0),5ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,0,3,0),3ul);

  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,5,1,0),15ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,5,2,0),7ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,5,3,0),5ul);

  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,5,1,0.5),15ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,5,2,1),8ul);

}

BOOST_AUTO_TEST_CASE(FillAccumulator_checkBelowThreshold){

  const size_t vector_size = 1;
  const double PE_base = 10;

  unsigned int AccumIndex = 0;
  unsigned int HitIndex = 0;
  double PE = 10;
  std::vector<double> Binned(vector_size,PE_base);
  std::vector< std::vector<int> > Contributors(vector_size);
  std::vector<int> FlashesInAccumulator;

  opdet::FillAccumulator(AccumIndex,HitIndex,PE,FlashThreshold,
			 Binned,Contributors,FlashesInAccumulator);

  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).size() , 1);
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(0) , HitIndex );
  BOOST_CHECK_EQUAL( Binned.at(AccumIndex) , PE+PE_base );
  BOOST_CHECK_EQUAL( FlashesInAccumulator.size() , 0);

}

BOOST_AUTO_TEST_CASE(FillAccumulator_checkAboveThreshold){

  const size_t vector_size = 1;
  const double PE_base = 10;

  unsigned int AccumIndex = 0;
  unsigned int HitIndex = 0;
  double PE = 40;
  std::vector<double> Binned(vector_size,PE_base);
  std::vector< std::vector<int> > Contributors(vector_size);
  std::vector<int> FlashesInAccumulator;

  opdet::FillAccumulator(AccumIndex,HitIndex,PE,FlashThreshold,
			 Binned,Contributors,FlashesInAccumulator);

  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).size() , 1);
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(0) , HitIndex );
  BOOST_CHECK_EQUAL( Binned.at(AccumIndex) , PE+PE_base );
  BOOST_CHECK_EQUAL( FlashesInAccumulator.size() , 1);

}

BOOST_AUTO_TEST_CASE(FillAccumulator_checkMultipleHits){

  const size_t vector_size = 1;
  const double PE_base = 10;

  unsigned int AccumIndex = 0;
  unsigned int HitIndex = 0;
  double PE = 25;
  std::vector<double> Binned(vector_size,PE_base);
  std::vector< std::vector<int> > Contributors(vector_size);
  std::vector<int> FlashesInAccumulator;

  opdet::FillAccumulator(AccumIndex,HitIndex,PE,FlashThreshold,
			 Binned,Contributors,FlashesInAccumulator);

  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).size() , 1);
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(0) , HitIndex );
  BOOST_CHECK_EQUAL( Binned.at(AccumIndex) , PE+PE_base );
  BOOST_CHECK_EQUAL( FlashesInAccumulator.size() , 0);


  unsigned int HitIndex2 = 1;
  opdet::FillAccumulator(AccumIndex,HitIndex2,PE,FlashThreshold,
			 Binned,Contributors,FlashesInAccumulator);

  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).size() , 2);
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(0) , HitIndex );
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(1) , HitIndex2 );
  BOOST_CHECK_EQUAL( Binned.at(AccumIndex) , PE*2+PE_base );
  BOOST_CHECK_EQUAL( FlashesInAccumulator.size() , 1);

}

BOOST_AUTO_TEST_CASE(FillFlashesBySizeMap_checkNoFlash)
{
  const size_t vector_size = 1;
  const double PE_vals = 10;

  std::vector<double> Binned(vector_size,PE_vals);
  std::vector<int> FlashesInAccumulator;

  

}

BOOST_AUTO_TEST_SUITE_END()
