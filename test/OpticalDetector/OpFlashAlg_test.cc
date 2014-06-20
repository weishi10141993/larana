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
  const size_t vector_size = 10;
  const double PE_vals = 10;

  std::vector<double> BinnedPE(vector_size,PE_vals);
  std::vector<int> FlashesInAccumulator;
  std::map<double, std::map<int,std::vector<int> >, std::greater<double> > FlashesBySize;

  opdet::FillFlashesBySizeMap(FlashesInAccumulator,BinnedPE,1,FlashesBySize);

  BOOST_CHECK_EQUAL( FlashesBySize.size() , 0 );

}

BOOST_AUTO_TEST_CASE(FillFlashesBySizeMap_checkOneFlash)
{
  const size_t vector_size = 10;
  const double PE_vals = 10;

  std::vector<double> BinnedPE(vector_size,PE_vals);
  BinnedPE[2] = 50;

  std::vector<int> FlashesInAccumulator(1,2);

  std::map<double, std::map<int,std::vector<int> >, std::greater<double> > FlashesBySize;

  opdet::FillFlashesBySizeMap(FlashesInAccumulator,BinnedPE,1,FlashesBySize);

  BOOST_CHECK_EQUAL( FlashesBySize.size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize.count(50) , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize.count(10) , 0 );
  BOOST_CHECK_EQUAL( FlashesBySize[50].size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[50].count(1) , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1].size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1][0] , 2 );

}

BOOST_AUTO_TEST_CASE(FillFlashesBySizeMap_checkTwoFlashes)
{
  const size_t vector_size = 10;
  const double PE_vals = 10;

  std::vector<double> BinnedPE(vector_size,PE_vals);
  BinnedPE[2] = 50;
  BinnedPE[8] = 50;

  std::vector<int> FlashesInAccumulator { 2, 8 };

  std::map<double, std::map<int,std::vector<int> >, std::greater<double> > FlashesBySize;

  opdet::FillFlashesBySizeMap(FlashesInAccumulator,BinnedPE,1,FlashesBySize);

  BOOST_CHECK_EQUAL( FlashesBySize.size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize.count(50) , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize.count(10) , 0 );

  BOOST_CHECK_EQUAL( FlashesBySize[50].size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[50].count(1) , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1].size() , 2 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1][0] , 2 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1][1] , 8 );

}

BOOST_AUTO_TEST_CASE(FillFlashesBySizeMap_checkTwoAccumulators)
{
  const size_t vector_size = 10;
  const double PE_vals = 10;

  std::vector<double> BinnedPE1(vector_size,PE_vals);
  BinnedPE1[2] = 50;
  BinnedPE1[5] = 60;
  std::vector<double> BinnedPE2(vector_size,PE_vals);
  BinnedPE2[8] = 50;

  std::vector<int> FlashesInAccumulator1 { 2, 5 };
  std::vector<int> FlashesInAccumulator2(1,8);

  std::map<double, std::map<int,std::vector<int> >, std::greater<double> > FlashesBySize;

  opdet::FillFlashesBySizeMap(FlashesInAccumulator1,BinnedPE1,1,FlashesBySize);
  opdet::FillFlashesBySizeMap(FlashesInAccumulator2,BinnedPE2,2,FlashesBySize);

  BOOST_CHECK_EQUAL( FlashesBySize.size() , 2 );

  BOOST_CHECK_EQUAL( FlashesBySize.count(50) , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize.count(60) , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize.count(10) , 0 );

  auto map_begin = FlashesBySize.begin();
  BOOST_CHECK_EQUAL( map_begin->first , 60 );

  auto map_last = FlashesBySize.end();
  map_last--;
  BOOST_CHECK_EQUAL( map_last->first , 50 );

  BOOST_CHECK_EQUAL( FlashesBySize[50].size() , 2 );
  BOOST_CHECK_EQUAL( FlashesBySize[50].count(1) , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[50].count(2) , 1 );

  BOOST_CHECK_EQUAL( FlashesBySize[50][1].size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1][0] , 2 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][2].size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][2][0] , 8 );

  BOOST_CHECK_EQUAL( FlashesBySize[60].size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[60].count(1) , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[60].count(2) , 0 );
  BOOST_CHECK_EQUAL( FlashesBySize[60][1].size() , 1 );
  BOOST_CHECK_EQUAL( FlashesBySize[60][1][0] , 5 );
}

BOOST_AUTO_TEST_CASE(FillHitsThisFlash_EmptyContributors)
{

  size_t NHits_prev = 0;
  size_t NHits = 10;

  const size_t vector_size = 1;
  std::vector< std::vector<int> > Contributors(vector_size);
  int Bin = 0;

  std::vector<int> HitClaimedByFlash(NHits,-1);
  std::vector<int> HitsThisFlash;

  opdet::FillHitsThisFlash(Contributors,
			   Bin,
			   NHits_prev,
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 0);
}

BOOST_AUTO_TEST_CASE(FillHitsThisFlash_NoPrevClaimedHits)
{

  size_t NHits_prev = 0;
  size_t NHits = 10;

  const size_t vector_size = 1;
  std::vector< std::vector<int> > Contributors(vector_size);
  int Bin = 0;
  Contributors[0].push_back(1); Contributors[0].push_back(3);

  std::vector<int> HitClaimedByFlash(NHits,-1);
  std::vector<int> HitsThisFlash;

  opdet::FillHitsThisFlash(Contributors,
			   Bin,
			   NHits_prev,
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 2 );
  BOOST_CHECK_EQUAL( HitsThisFlash[0] , 1 );
  BOOST_CHECK_EQUAL( HitsThisFlash[1] , 3 );
}

BOOST_AUTO_TEST_CASE(FillHitsThisFlash_PrevClaimedHits)
{

  size_t NHits_prev = 0;
  size_t NHits = 10;

  const size_t vector_size = 1;
  std::vector< std::vector<int> > Contributors(vector_size);
  int Bin = 0;
  Contributors[0].push_back(1); Contributors[0].push_back(3);

  std::vector<int> HitClaimedByFlash(NHits,-1);
  HitClaimedByFlash[1] = 0;
  HitClaimedByFlash[2] = 1;
  std::vector<int> HitsThisFlash;

  opdet::FillHitsThisFlash(Contributors,
			   Bin,
			   NHits_prev,
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 1 );
  BOOST_CHECK_EQUAL( HitsThisFlash[0] , 3 );
}

BOOST_AUTO_TEST_CASE(FillHitsThisFlash_PrevHitsOffset)
{

  size_t NHits_prev = 10;
  size_t NHits = 10;

  const size_t vector_size = 1;
  std::vector< std::vector<int> > Contributors(vector_size);
  int Bin = 0;
  Contributors[0].push_back(11); Contributors[0].push_back(13);

  std::vector<int> HitClaimedByFlash(NHits,-1);
  HitClaimedByFlash[1] = 0;
  HitClaimedByFlash[2] = 1;
  std::vector<int> HitsThisFlash;

  opdet::FillHitsThisFlash(Contributors,
			   Bin,
			   NHits_prev,
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 1 );
  BOOST_CHECK_EQUAL( HitsThisFlash[0] , 13 );
}

BOOST_AUTO_TEST_CASE(FillHitsThisFlash_MultipleContributorVectors)
{

  size_t NHits_prev = 0;
  size_t NHits = 10;

  const size_t vector_size = 2;
  std::vector< std::vector<int> > Contributors(vector_size);
  Contributors[0].push_back(1); Contributors[0].push_back(3);
  Contributors[1].push_back(5); Contributors[1].push_back(6);

  std::vector<int> HitClaimedByFlash(NHits,-1);
  HitClaimedByFlash[1] = 0;
  HitClaimedByFlash[2] = 1;
  std::vector<int> HitsThisFlash;

  opdet::FillHitsThisFlash(Contributors,
			   1,
			   NHits_prev,
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 2 );
  BOOST_CHECK_EQUAL( HitsThisFlash[0] , 5 );
  BOOST_CHECK_EQUAL( HitsThisFlash[1] , 6 );
}

BOOST_AUTO_TEST_CASE(ClaimHits_NoHitsThisFlash)
{
  size_t NHits_prev = 0;
  size_t NHits = 1;

  //need this part to make a hit...
  double hit_pe = 10;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++)
    HitVector.emplace_back(0,0,0,0,0,0,0,hit_pe,0);

  std::vector<int> HitsThisFlash;
  std::vector< std::vector<int> > HitsPerFlash;

  std::vector<int> HitClaimedByFlash(NHits,-1);

  opdet::ClaimHits(HitVector,
		   HitsThisFlash,
		   FlashThreshold,
		   HitsPerFlash,
		   NHits_prev,
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], -1);
}

BOOST_AUTO_TEST_CASE(ClaimHits_BelowFlashThreshold)
{
  size_t NHits_prev = 0;
  size_t NHits = 1;

  //need this part to make a hit...
  double hit_pe = 10;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++)
    HitVector.emplace_back(0,0,0,0,0,0,0,hit_pe,0);

  std::vector<int> HitsThisFlash(1,0);
  std::vector< std::vector<int> > HitsPerFlash;

  std::vector<int> HitClaimedByFlash(NHits,-1);

  opdet::ClaimHits(HitVector,
		   HitsThisFlash,
		   FlashThreshold,
		   HitsPerFlash,
		   NHits_prev,
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], -1);
}

BOOST_AUTO_TEST_CASE(ClaimHits_AboveFlashThreshold)
{
  size_t NHits_prev = 0;
  size_t NHits = 1;

  //need this part to make a hit...
  double hit_pe = 100;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++)
    HitVector.emplace_back(0,0,0,0,0,0,0,hit_pe,0);

  std::vector<int> HitsThisFlash(1,0);
  std::vector< std::vector<int> > HitsPerFlash;

  std::vector<int> HitClaimedByFlash(NHits,-1);

  opdet::ClaimHits(HitVector,
		   HitsThisFlash,
		   FlashThreshold,
		   HitsPerFlash,
		   NHits_prev,
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 1);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], 0);
}

BOOST_AUTO_TEST_CASE(ClaimHits_OneHitThisFlash)
{
  size_t NHits_prev = 0;
  size_t NHits = 2;

  //need this part to make a hit...
  double hit_pe = 100;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++)
    HitVector.emplace_back(0,0,0,0,0,0,0,hit_pe,0);

  std::vector<int> HitsThisFlash(1,0);
  std::vector< std::vector<int> > HitsPerFlash;

  std::vector<int> HitClaimedByFlash(NHits,-1);

  opdet::ClaimHits(HitVector,
		   HitsThisFlash,
		   FlashThreshold,
		   HitsPerFlash,
		   NHits_prev,
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 1);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[1], -1);
}

BOOST_AUTO_TEST_CASE(ClaimHits_TwoHitsThisFlash_BelowThreshold)
{
  size_t NHits_prev = 0;
  size_t NHits = 2;

  //need this part to make a hit...
  double hit_pe = 10;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++)
    HitVector.emplace_back(0,0,0,0,0,0,0,hit_pe,0);

  std::vector<int> HitsThisFlash;
  HitsThisFlash.push_back(0); HitsThisFlash.push_back(1);
  std::vector< std::vector<int> > HitsPerFlash;

  std::vector<int> HitClaimedByFlash(NHits,-1);

  opdet::ClaimHits(HitVector,
		   HitsThisFlash,
		   FlashThreshold,
		   HitsPerFlash,
		   NHits_prev,
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], -1);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[1], -1);
}

BOOST_AUTO_TEST_CASE(ClaimHits_TwoHitsThisFlash_AboveThreshold)
{
  size_t NHits_prev = 0;
  size_t NHits = 2;

  //need this part to make a hit...
  double hit_pe = 30;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++)
    HitVector.emplace_back(0,0,0,0,0,0,0,hit_pe,0);

  std::vector<int> HitsThisFlash;
  HitsThisFlash.push_back(0); HitsThisFlash.push_back(1);
  std::vector< std::vector<int> > HitsPerFlash;

  std::vector<int> HitClaimedByFlash(NHits,-1);

  opdet::ClaimHits(HitVector,
		   HitsThisFlash,
		   FlashThreshold,
		   HitsPerFlash,
		   NHits_prev,
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 1);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][0], 0);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][1], 1);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[1], 0);
}

BOOST_AUTO_TEST_CASE(ClaimHits_WithPrevHits)
{
  size_t NHits_prev = 10;
  size_t NHits = 2;

  //need this part to make a hit...
  double hit_pe = 30;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits+NHits_prev; i++)
    HitVector.emplace_back(0,0,0,0,0,0,0,hit_pe,0);

  std::vector<int> HitsThisFlash;
  HitsThisFlash.push_back(10); HitsThisFlash.push_back(11);
  std::vector< std::vector<int> > HitsPerFlash;

  std::vector<int> HitClaimedByFlash(NHits,-1);

  opdet::ClaimHits(HitVector,
		   HitsThisFlash,
		   FlashThreshold,
		   HitsPerFlash,
		   NHits_prev,
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 1);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][0], 10);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][1], 11);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[1], 0);
}

BOOST_AUTO_TEST_SUITE_END()
