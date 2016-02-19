#define BOOST_TEST_MODULE ( OpFlashAlg_test )
#include "boost/test/auto_unit_test.hpp"

#include "larana/OpticalDetector/OpFlashAlg.h"
#include "lardata/OpticalDetectorData/OpticalTypes.h"

const float HitThreshold = 3;
const float FlashThreshold = 50;
const double WidthTolerance = 0.5;
const int Channel = 1;
const uint32_t TimeSlice = 1000;
const unsigned short Frame = 1;
const optdata::TimeSlice_t TimeSlicesPerFrame = 102400;
const double opdigi_SampleFreq = 64;
const double TrigTimeAbs = 0.8;
const double SPESize = 20;

const double tolerance = 1e-6;

BOOST_AUTO_TEST_SUITE(OpFlashAlg_test)

/*
BOOST_AUTO_TEST_CASE(ConstructHit_checkDefaultPulse)
{
  pmtana::pulse_param pulse;
  std::vector<recob::OpHit> HitVector;
  art::Service<detinfo::DetectorClocksService> ts;
  //::detinfo::DetectorClocksService ts;

  opdet::ConstructHit(HitThreshold,
		      Channel,
		      TimeSlice,
		      Frame,
		      &pulse,
		      *ts,
		      SPESize,
		      HitVector);
  
  BOOST_CHECK_EQUAL(HitVector.size(),0ul);

}
*/
/*
BOOST_AUTO_TEST_CASE(ConstructHit_checkPulseBelowThreshold)
{
  pmtana::pulse_param pulse;
  pulse.t_start=3; pulse.t_max = 6; pulse.t_end = 9;
  pulse.area = 4; pulse.peak = 2;

  std::vector<recob::OpHit> HitVector;
  art::Service<detinfo::DetectorClocksService> ts;

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
*/
 /*
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
  double pulse_end = (pulse.t_end + TimeSlice + Frame*TimeSlicesPerFrame)/opdigi_SampleFreq;
  double pulse_begin = (pulse.t_start + TimeSlice + Frame*TimeSlicesPerFrame)/opdigi_SampleFreq;

  BOOST_CHECK_EQUAL(HitVector.size(),1ul);
  BOOST_CHECK_EQUAL(HitVector[0].OpChannel(),Channel);
  BOOST_CHECK_CLOSE(HitVector[0].PeakTimeAbs(),peak_time_abs,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].PeakTime(),peak_time_abs-TrigTimeAbs,tolerance);
  BOOST_CHECK_EQUAL(HitVector[0].Frame(),Frame);
  BOOST_CHECK_CLOSE(HitVector[0].Width(),pulse_end-pulse_begin,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].Area(),pulse.area,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].Amplitude(),pulse.peak,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].PE(),pulse.peak/SPESize,tolerance);
  BOOST_CHECK_CLOSE(HitVector[0].FastToTotal(),0.,tolerance);

}
 */
BOOST_AUTO_TEST_CASE(checkGetAccumIndex)
{

  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,0,1,0),10ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,0,2,0),5ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,0,3,0),3ul);

  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,-5,1,0),15ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,-5,2,0),7ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,-5,3,0),5ul);

  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,-5,1,0.5),15ul);
  BOOST_CHECK_EQUAL(opdet::GetAccumIndex(10.0,-5,2,1),8ul);

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

  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).size() , 1U);
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(0) , (int)HitIndex );
  BOOST_CHECK_EQUAL( Binned.at(AccumIndex) , PE+PE_base );
  BOOST_CHECK_EQUAL( FlashesInAccumulator.size() , 0U);

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

  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).size() , 1U);
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(0) , (int)HitIndex );
  BOOST_CHECK_EQUAL( Binned.at(AccumIndex) , PE+PE_base );
  BOOST_CHECK_EQUAL( FlashesInAccumulator.size() , 1U);

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

  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).size() , 1U);
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(0) , (int)HitIndex );
  BOOST_CHECK_EQUAL( Binned.at(AccumIndex) , PE+PE_base );
  BOOST_CHECK_EQUAL( FlashesInAccumulator.size() , 0U);


  unsigned int HitIndex2 = 1;
  opdet::FillAccumulator(AccumIndex,HitIndex2,PE,FlashThreshold,
			 Binned,Contributors,FlashesInAccumulator);

  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).size() , 2U);
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(0) , (int)HitIndex );
  BOOST_CHECK_EQUAL( Contributors.at(AccumIndex).at(1) , (int)HitIndex2 );
  BOOST_CHECK_EQUAL( Binned.at(AccumIndex) , PE*2+PE_base );
  BOOST_CHECK_EQUAL( FlashesInAccumulator.size() , 1U);

}

BOOST_AUTO_TEST_CASE(FillFlashesBySizeMap_checkNoFlash)
{
  const size_t vector_size = 10;
  const double PE_vals = 10;

  std::vector<double> BinnedPE(vector_size,PE_vals);
  std::vector<int> FlashesInAccumulator;
  std::map<double, std::map<int,std::vector<int> >, std::greater<double> > FlashesBySize;

  opdet::FillFlashesBySizeMap(FlashesInAccumulator,BinnedPE,1,FlashesBySize);

  BOOST_CHECK_EQUAL( FlashesBySize.size() , 0U );

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

  BOOST_CHECK_EQUAL( FlashesBySize.size() , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize.count(50) , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize.count(10) , 0U );
  BOOST_CHECK_EQUAL( FlashesBySize[50].size() , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[50].count(1) , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1].size() , 1U );
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

  BOOST_CHECK_EQUAL( FlashesBySize.size() , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize.count(50) , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize.count(10) , 0U );

  BOOST_CHECK_EQUAL( FlashesBySize[50].size() , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[50].count(1) , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1].size() , 2U );
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

  BOOST_CHECK_EQUAL( FlashesBySize.size() , 2U );

  BOOST_CHECK_EQUAL( FlashesBySize.count(50) , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize.count(60) , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize.count(10) , 0U );

  auto map_begin = FlashesBySize.begin();
  BOOST_CHECK_EQUAL( map_begin->first , 60 );

  auto map_last = FlashesBySize.end();
  map_last--;
  BOOST_CHECK_EQUAL( map_last->first , 50 );

  BOOST_CHECK_EQUAL( FlashesBySize[50].size() , 2U );
  BOOST_CHECK_EQUAL( FlashesBySize[50].count(1) , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[50].count(2) , 1U );

  BOOST_CHECK_EQUAL( FlashesBySize[50][1].size() , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[50][1][0] , 2 );
  BOOST_CHECK_EQUAL( FlashesBySize[50][2].size() , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[50][2][0] , 8 );

  BOOST_CHECK_EQUAL( FlashesBySize[60].size() , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[60].count(1) , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[60].count(2) , 0U );
  BOOST_CHECK_EQUAL( FlashesBySize[60][1].size() , 1U );
  BOOST_CHECK_EQUAL( FlashesBySize[60][1][0] , 5 );
}

BOOST_AUTO_TEST_CASE(FillHitsThisFlash_EmptyContributors)
{

  size_t NHits = 10;

  const size_t vector_size = 1;
  std::vector< std::vector<int> > Contributors(vector_size);
  int Bin = 0;

  std::vector<int> HitClaimedByFlash(NHits,-1);
  std::vector<int> HitsThisFlash;

  opdet::FillHitsThisFlash(Contributors,
			   Bin,
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 0U);
}

BOOST_AUTO_TEST_CASE(FillHitsThisFlash_NoPrevClaimedHits)
{

  size_t NHits = 10;

  const size_t vector_size = 1;
  std::vector< std::vector<int> > Contributors(vector_size);
  int Bin = 0;
  Contributors[0].push_back(1); Contributors[0].push_back(3);

  std::vector<int> HitClaimedByFlash(NHits,-1);
  std::vector<int> HitsThisFlash;

  opdet::FillHitsThisFlash(Contributors,
			   Bin,
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 2U );
  BOOST_CHECK_EQUAL( HitsThisFlash[0] , 1 );
  BOOST_CHECK_EQUAL( HitsThisFlash[1] , 3 );
}

BOOST_AUTO_TEST_CASE(FillHitsThisFlash_PrevClaimedHits)
{

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
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 1U );
  BOOST_CHECK_EQUAL( HitsThisFlash[0] , 3 );
}
/*
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

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 1U );
  BOOST_CHECK_EQUAL( HitsThisFlash[0] , 13 );
}
*/
BOOST_AUTO_TEST_CASE(FillHitsThisFlash_MultipleContributorVectors)
{

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
			   HitClaimedByFlash,
			   HitsThisFlash);

  BOOST_CHECK_EQUAL( HitsThisFlash.size() , 2U );
  BOOST_CHECK_EQUAL( HitsThisFlash[0] , 5 );
  BOOST_CHECK_EQUAL( HitsThisFlash[1] , 6 );
}

BOOST_AUTO_TEST_CASE(ClaimHits_NoHitsThisFlash)
{
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
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 0U);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], -1);
}

BOOST_AUTO_TEST_CASE(ClaimHits_BelowFlashThreshold)
{
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
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 0U);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], -1);
}

BOOST_AUTO_TEST_CASE(ClaimHits_AboveFlashThreshold)
{
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
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 1U);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], 0);
}

BOOST_AUTO_TEST_CASE(ClaimHits_OneHitThisFlash)
{
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
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 1U);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[1], -1);
}

BOOST_AUTO_TEST_CASE(ClaimHits_TwoHitsThisFlash_BelowThreshold)
{
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
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 0U);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], -1);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[1], -1);
}

BOOST_AUTO_TEST_CASE(ClaimHits_TwoHitsThisFlash_AboveThreshold)
{
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
		   HitClaimedByFlash);

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 1U);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][0], 0);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][1], 1);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[1], 0);
}
/*
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

  BOOST_CHECK_EQUAL( HitsPerFlash.size(), 1U);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][0], 10);
  BOOST_CHECK_EQUAL( HitsPerFlash[0][1], 11);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[0], 0);
  BOOST_CHECK_EQUAL( HitClaimedByFlash[1], 0);
}
*/
BOOST_AUTO_TEST_CASE(FindSeedHit_AllUsed)
{

  size_t NHits = 5;

  std::vector<int> HitsThisRefinedFlash;
  double PEAccumulated=0, FlashMaxTime=0, FlashMinTime=0;
  
  double hit_pe = 20;
  double peak_time = 0;
  double width = 10;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++){
    HitVector.emplace_back(0,peak_time,0,0,width,0,0,hit_pe,0);
    hit_pe+=10;
    peak_time+=1;
  }

  std::map<double, std::vector<int>, std::greater<double> > HitsBySize;
  for(size_t i=0; i<NHits; i++)
    HitsBySize[HitVector.at(i).PE()].push_back(i);

  std::vector<bool> HitsUsed(NHits,true);

  opdet::FindSeedHit(HitsBySize,
		     HitsUsed,
		     HitVector,
		     HitsThisRefinedFlash,
		     PEAccumulated,
		     FlashMaxTime,
		     FlashMinTime);

  BOOST_CHECK_EQUAL( HitsThisRefinedFlash.size() , 0U );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , (int)NHits);
  BOOST_CHECK_EQUAL( PEAccumulated , 0 );
  BOOST_CHECK_EQUAL( FlashMaxTime , 0 );
  BOOST_CHECK_EQUAL( FlashMinTime , 0 );

}

BOOST_AUTO_TEST_CASE(FindSeedHit_NoneUsed)
{

  size_t NHits = 5;

  std::vector<int> HitsThisRefinedFlash;
  double PEAccumulated=0, FlashMaxTime=0, FlashMinTime=0;
  
  double hit_pe = 20;
  double peak_time = 0;
  double width = 10;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++){
    HitVector.emplace_back(0,peak_time,0,0,width,0,0,hit_pe,0);
    hit_pe+=10;
    peak_time+=1;
  }

  std::map<double, std::vector<int>, std::greater<double> > HitsBySize;
  for(size_t i=0; i<NHits; i++)
    HitsBySize[HitVector.at(i).PE()].push_back(i);

  std::vector<bool> HitsUsed(NHits,false);

  opdet::FindSeedHit(HitsBySize,
		     HitsUsed,
		     HitVector,
		     HitsThisRefinedFlash,
		     PEAccumulated,
		     FlashMaxTime,
		     FlashMinTime);

  BOOST_CHECK_EQUAL( HitsThisRefinedFlash.size() , 1U );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 1);
  BOOST_CHECK_EQUAL( HitsUsed.at(NHits-1) , true);
  BOOST_CHECK_EQUAL( PEAccumulated , 60 );
  BOOST_CHECK_EQUAL( FlashMaxTime , 9. );
  BOOST_CHECK_EQUAL( FlashMinTime , -1. );

}

BOOST_AUTO_TEST_CASE(FindSeedHit_FirstUsed)
{

  size_t NHits = 5;

  std::vector<int> HitsThisRefinedFlash;
  double PEAccumulated=0, FlashMaxTime=0, FlashMinTime=0;
  
  double hit_pe = 20;
  double peak_time = 0;
  double width = 10;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++){
    HitVector.emplace_back(0,peak_time,0,0,width,0,0,hit_pe,0);
    hit_pe+=10;
    peak_time+=1;
  }

  std::map<double, std::vector<int>, std::greater<double> > HitsBySize;
  for(size_t i=0; i<NHits; i++)
    HitsBySize[HitVector.at(i).PE()].push_back(i);

  std::vector<bool> HitsUsed(NHits,false);
  HitsUsed[NHits-1] = true;

  opdet::FindSeedHit(HitsBySize,
		     HitsUsed,
		     HitVector,
		     HitsThisRefinedFlash,
		     PEAccumulated,
		     FlashMaxTime,
		     FlashMinTime);

  BOOST_CHECK_EQUAL( HitsThisRefinedFlash.size() , 1U );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 2);
  BOOST_CHECK_EQUAL( HitsUsed.at(NHits-1) , true);
  BOOST_CHECK_EQUAL( HitsUsed.at(NHits-2) , true);
  BOOST_CHECK_EQUAL( PEAccumulated , 50 );
  BOOST_CHECK_EQUAL( FlashMaxTime , 8. );
  BOOST_CHECK_EQUAL( FlashMinTime , -2. );

}

BOOST_AUTO_TEST_CASE(AddHitToFlash_UsedHit)
{

  size_t NHits = 5;

  std::vector<int> HitsThisRefinedFlash;
  double PEAccumulated=0, FlashMaxTime=0, FlashMinTime=0;
  
  double hit_pe = 20;
  double peak_time = 0;
  double width = 10;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++){
    HitVector.emplace_back(0,peak_time,0,0,width,0,0,hit_pe,0);
    hit_pe+=10;
    peak_time+=1;
  }

  std::map<double, std::vector<int>, std::greater<double> > HitsBySize;
  for(size_t i=0; i<NHits; i++)
    HitsBySize[HitVector.at(i).PE()].push_back(i);
  std::vector<bool> HitsUsed(NHits,false);

  opdet::FindSeedHit(HitsBySize,
		     HitsUsed,
		     HitVector,
		     HitsThisRefinedFlash,
		     PEAccumulated,
		     FlashMaxTime,
		     FlashMinTime);

  int HitID = 4;
  opdet::AddHitToFlash( HitID,
			HitsUsed,
			HitVector.at(HitID),
			WidthTolerance,
			HitsThisRefinedFlash,
			PEAccumulated,
			FlashMaxTime,
			FlashMinTime);

  BOOST_CHECK_EQUAL( HitsThisRefinedFlash.size() , 1U );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 1);
  BOOST_CHECK_EQUAL( HitsUsed.at(NHits-1) , true);
  BOOST_CHECK_EQUAL( PEAccumulated , 60 );
  BOOST_CHECK_EQUAL( FlashMaxTime , 9. );
  BOOST_CHECK_EQUAL( FlashMinTime , -1. );

}

BOOST_AUTO_TEST_CASE(AddHitToFlash_NewHit)
{

  size_t NHits = 5;

  std::vector<int> HitsThisRefinedFlash;
  double PEAccumulated=0, FlashMaxTime=0, FlashMinTime=0;
  
  double hit_pe = 20;
  double peak_time = 0;
  double width = 10;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++){
    HitVector.emplace_back(0,peak_time,0,0,width,0,0,hit_pe,0);
    hit_pe+=10;
    peak_time+=1;
  }

  std::map<double, std::vector<int>, std::greater<double> > HitsBySize;
  for(size_t i=0; i<NHits; i++)
    HitsBySize[HitVector.at(i).PE()].push_back(i);
  std::vector<bool> HitsUsed(NHits,false);

  opdet::FindSeedHit(HitsBySize,
		     HitsUsed,
		     HitVector,
		     HitsThisRefinedFlash,
		     PEAccumulated,
		     FlashMaxTime,
		     FlashMinTime);

  int HitID = 3;
  opdet::AddHitToFlash( HitID,
			HitsUsed,
			HitVector.at(HitID),
			WidthTolerance,
			HitsThisRefinedFlash,
			PEAccumulated,
			FlashMaxTime,
			FlashMinTime);

  BOOST_CHECK_EQUAL( HitsThisRefinedFlash.size() , 2U );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 2);
  BOOST_CHECK_EQUAL( HitsUsed.at(NHits-1) , true);
  BOOST_CHECK_EQUAL( HitsUsed.at(NHits-2) , true);
  BOOST_CHECK_EQUAL( PEAccumulated , 110 );
  BOOST_CHECK_EQUAL( FlashMaxTime , 9. );
  BOOST_CHECK_EQUAL( FlashMinTime , -2. );

}

BOOST_AUTO_TEST_CASE(AddHitToFlash_OutsideWidth)
{

  size_t NHits = 5;

  std::vector<int> HitsThisRefinedFlash;
  double PEAccumulated=0, FlashMaxTime=0, FlashMinTime=0;
  
  double hit_pe = 20;
  double peak_time = 0;
  double width = 0.1;
  std::vector<recob::OpHit> HitVector;
  for(size_t i=0; i<NHits; i++){
    HitVector.emplace_back(0,peak_time,0,0,width,0,0,hit_pe,0);
    hit_pe+=10;
    peak_time+=1;
  }

  std::map<double, std::vector<int>, std::greater<double> > HitsBySize;
  for(size_t i=0; i<NHits; i++)
    HitsBySize[HitVector.at(i).PE()].push_back(i);
  std::vector<bool> HitsUsed(NHits,false);

  opdet::FindSeedHit(HitsBySize,
		     HitsUsed,
		     HitVector,
		     HitsThisRefinedFlash,
		     PEAccumulated,
		     FlashMaxTime,
		     FlashMinTime);

  int HitID = 3;
  opdet::AddHitToFlash( HitID,
			HitsUsed,
			HitVector.at(HitID),
			WidthTolerance,
			HitsThisRefinedFlash,
			PEAccumulated,
			FlashMaxTime,
			FlashMinTime);

  BOOST_CHECK_EQUAL( HitsThisRefinedFlash.size() , 1U );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 1);
  BOOST_CHECK_EQUAL( HitsUsed.at(NHits-1) , true);
  BOOST_CHECK_EQUAL( HitsUsed.at(NHits-2) , false);
  BOOST_CHECK_EQUAL( PEAccumulated , 60 );
  BOOST_CHECK_EQUAL( FlashMaxTime , 4.05 );
  BOOST_CHECK_EQUAL( FlashMinTime , 3.95 );

}

BOOST_AUTO_TEST_CASE(CheckAndStoreFlash_AboveThreshold)
{
  std::vector< std::vector<int> > RefinedHitsPerFlash;

  std::vector<int> HitsThisRefinedFlash{0,1,2};
  std::vector<bool> HitsUsed{true,true,true,false,false};
  double PEAccumulated = 60;

  opdet::CheckAndStoreFlash( RefinedHitsPerFlash,
			     HitsThisRefinedFlash,
			     PEAccumulated,
			     FlashThreshold,
			     HitsUsed );

  BOOST_CHECK_EQUAL( RefinedHitsPerFlash.size(), 1U );
  BOOST_CHECK_EQUAL( RefinedHitsPerFlash[0].size(), 3U );
  BOOST_CHECK_EQUAL( RefinedHitsPerFlash[0][0], 0 );
  BOOST_CHECK_EQUAL( RefinedHitsPerFlash[0][1], 1 );
  BOOST_CHECK_EQUAL( RefinedHitsPerFlash[0][2], 2 );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 3);

}

BOOST_AUTO_TEST_CASE(CheckAndStoreFlash_AboveThreshold_OneHit)
{
  std::vector< std::vector<int> > RefinedHitsPerFlash;

  std::vector<int> HitsThisRefinedFlash{0};
  std::vector<bool> HitsUsed{true,false,false,false,false};
  double PEAccumulated = 60;

  opdet::CheckAndStoreFlash( RefinedHitsPerFlash,
			     HitsThisRefinedFlash,
			     PEAccumulated,
			     FlashThreshold,
			     HitsUsed );

  BOOST_CHECK_EQUAL( RefinedHitsPerFlash.size(), 1U );
  BOOST_CHECK_EQUAL( RefinedHitsPerFlash[0].size(), 1U );
  BOOST_CHECK_EQUAL( RefinedHitsPerFlash[0][0], 0 );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 1);

}

BOOST_AUTO_TEST_CASE(CheckAndStoreFlash_BelowThreshold_OneHit)
{
  std::vector< std::vector<int> > RefinedHitsPerFlash;

  std::vector<int> HitsThisRefinedFlash{0};
  std::vector<bool> HitsUsed{true,false,false,false,false};
  double PEAccumulated = 30;

  opdet::CheckAndStoreFlash( RefinedHitsPerFlash,
			     HitsThisRefinedFlash,
			     PEAccumulated,
			     FlashThreshold,
			     HitsUsed );

  BOOST_CHECK_EQUAL( RefinedHitsPerFlash.size(), 0U );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 1);
  BOOST_CHECK_EQUAL( HitsUsed[0] , true);

}

BOOST_AUTO_TEST_CASE(CheckAndStoreFlash_BelowThreshold_MultipleHits)
{
  std::vector< std::vector<int> > RefinedHitsPerFlash;

  std::vector<int> HitsThisRefinedFlash{0,1,2};
  std::vector<bool> HitsUsed{true,true,true,false,false};
  double PEAccumulated = 30;

  opdet::CheckAndStoreFlash( RefinedHitsPerFlash,
			     HitsThisRefinedFlash,
			     PEAccumulated,
			     FlashThreshold,
			     HitsUsed );

  BOOST_CHECK_EQUAL( RefinedHitsPerFlash.size(), 0U );
  BOOST_CHECK_EQUAL( std::count(HitsUsed.begin(),HitsUsed.end(),true) , 1);
  BOOST_CHECK_EQUAL( HitsUsed[0] , true);
  BOOST_CHECK_EQUAL( HitsUsed[1] , false);
  BOOST_CHECK_EQUAL( HitsUsed[2] , false);

}

BOOST_AUTO_TEST_CASE(AddHitContribution_AddFirstHit)
{
    double MaxTime = -1e9, MinTime = 1e9;
    double TotalPE=0, AveTime=0, AveAbsTime=0, FastToTotal=0;

    size_t NOpChannels = 5;
    std::vector<double> PEs(NOpChannels,0);    

    double hit_pe = 20;
    double peak_time = 1;
    double width = 5;
    int op_channel = 0;
    recob::OpHit currentHit(op_channel,peak_time,0,0,width,0,0,hit_pe,0);

    opdet::AddHitContribution( currentHit,
			       MaxTime,
			       MinTime,
			       AveTime,
			       FastToTotal,
			       AveAbsTime,
			       TotalPE,
			       PEs);

    BOOST_CHECK_EQUAL( MaxTime , peak_time );
    BOOST_CHECK_EQUAL( MinTime , peak_time );
    BOOST_CHECK_EQUAL( AveTime , peak_time*hit_pe );
    BOOST_CHECK_EQUAL( AveAbsTime , 0 );
    BOOST_CHECK_EQUAL( FastToTotal , 0 );
    BOOST_CHECK_EQUAL( TotalPE , hit_pe );
    BOOST_CHECK_EQUAL( PEs.at(op_channel) , hit_pe );
    BOOST_CHECK_EQUAL( PEs.at(4) , 0 );
}

BOOST_AUTO_TEST_CASE(AddHitContribution_AddSecondHit)
{
    double MaxTime = 1, MinTime = 1;
    double TotalPE=20, AveTime=20, AveAbsTime=0, FastToTotal=0;

    size_t NOpChannels = 5;
    std::vector<double> PEs(NOpChannels,0);
    PEs.at(0)+=20;

    double hit_pe = 30;
    double peak_time = 5;
    double width = 5;
    int op_channel = 2;
    recob::OpHit currentHit(op_channel,peak_time,0,0,width,0,0,hit_pe,0);

    opdet::AddHitContribution( currentHit,
			       MaxTime,
			       MinTime,
			       AveTime,
			       FastToTotal,
			       AveAbsTime,
			       TotalPE,
			       PEs);

    BOOST_CHECK_EQUAL( MaxTime , peak_time );
    BOOST_CHECK_EQUAL( MinTime , 1 );
    BOOST_CHECK_EQUAL( AveTime , peak_time*hit_pe+20 );
    BOOST_CHECK_EQUAL( AveAbsTime , 0 );
    BOOST_CHECK_EQUAL( FastToTotal , 0 );
    BOOST_CHECK_EQUAL( TotalPE , hit_pe+20 );
    BOOST_CHECK_EQUAL( PEs.at(op_channel) , hit_pe );
    BOOST_CHECK_EQUAL( PEs.at(0) , 20 );
    BOOST_CHECK_EQUAL( PEs.at(4) , 0 );
}


BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_BackwardsTime)
{

  double iPE = 100; double iTime = 0; double iWidth=0.5;
  double jPE = 100; double jTime = -1; double jWidth=0.5;
  
  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
						jPE, jTime, jWidth);
  
  BOOST_CHECK_CLOSE( result, 1e6, tolerance);
}

BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_EqualFlashes)
{

  double iPE = 100; double iTime = 0; double iWidth=0.5;
  double jPE = 100; double jTime = 0; double jWidth=0.5;
  
  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
						jPE, jTime, jWidth);
  
  BOOST_CHECK_CLOSE( result, 0, tolerance);
}

BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_LateFlash)
{

  double iPE = 100; double iTime = 0; double iWidth=0.5;
  double jPE = 10; double jTime = 1.6; double jWidth=0.5;
  
  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
						jPE, jTime, jWidth);
  
  double good_result = (jPE - std::exp(-1)*iPE)/(std::sqrt(std::exp(-1)*iPE));

  BOOST_CHECK_CLOSE( result, good_result, tolerance);
}

BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_VeryLateFlash)
{

  double iPE = 100; double iTime = 0; double iWidth=0.5;
  double jPE = 10; double jTime = 16; double jWidth=0.5;
  
  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
						jPE, jTime, jWidth);
  
  double good_result = (jPE - std::exp(-10)*iPE)/(std::sqrt(std::exp(-10)*iPE));

  BOOST_CHECK_CLOSE( result, good_result, tolerance);
}

BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_UnequalWidths)
{

  double iPE = 100; double iTime = 0; double iWidth=1;
  double jPE = 10; double jTime = 1.6; double jWidth=0.5;
  
  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
						jPE, jTime, jWidth);
  
  double good_result = (jPE - std::exp(-1)*iPE*0.5)/(std::sqrt(std::exp(-1)*iPE*0.5));

  BOOST_CHECK_CLOSE( result, good_result, tolerance);
}

BOOST_AUTO_TEST_CASE(MarkFlashesForRemoval_NoFlashes)
{
  size_t NFlashes=5;
  size_t BeginFlash=5; //that is, no new flashes to look at

  std::vector<recob::OpFlash> FlashVector(NFlashes);
  std::vector<bool> MarkedForRemoval(NFlashes-BeginFlash,false);

  opdet::MarkFlashesForRemoval(FlashVector,
			       BeginFlash,
			       MarkedForRemoval);

  BOOST_CHECK_EQUAL( MarkedForRemoval.size() , 0U );

}

BOOST_AUTO_TEST_CASE(MarkFlashesForRemoval_OneFlash)
{
  size_t NFlashes=1;
  size_t BeginFlash=0;

  std::vector<double> PEs(30,0);
  PEs.at(0) = 100;
  std::vector<double> WireCenters(3,0);
  std::vector<double> WireWidths(3,0);

  std::vector<recob::OpFlash> FlashVector;
  FlashVector.emplace_back(0, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  std::vector<bool> MarkedForRemoval(NFlashes-BeginFlash,false);

  opdet::MarkFlashesForRemoval(FlashVector,
			       BeginFlash,
			       MarkedForRemoval);

  BOOST_CHECK_EQUAL( MarkedForRemoval.size() , 1U );
  BOOST_CHECK_EQUAL( MarkedForRemoval[0] , false );
  BOOST_CHECK_EQUAL( FlashVector.size() , 1U );
  BOOST_CHECK_EQUAL( FlashVector[0].Time() , 0 );
  BOOST_CHECK_EQUAL( FlashVector[0].TimeWidth() , 0.5 );

}

BOOST_AUTO_TEST_CASE(MarkFlashesForRemoval_TwoIndieFlashes)
{
  size_t NFlashes=2;
  size_t BeginFlash=0;

  std::vector<double> PEs(30,0);
  PEs.at(0) = 100;
  std::vector<double> WireCenters(3,0);
  std::vector<double> WireWidths(3,0);

  std::vector<recob::OpFlash> FlashVector;
  FlashVector.emplace_back(0, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(1e6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  std::vector<bool> MarkedForRemoval(NFlashes-BeginFlash,false);

  opdet::MarkFlashesForRemoval(FlashVector,
			       BeginFlash,
			       MarkedForRemoval);

  BOOST_CHECK_EQUAL( MarkedForRemoval.size() , 2U );
  BOOST_CHECK_EQUAL( MarkedForRemoval[0] , false );
  BOOST_CHECK_EQUAL( MarkedForRemoval[1] , false );
  BOOST_CHECK_EQUAL( FlashVector.size() , 2U );
  BOOST_CHECK_EQUAL( FlashVector[0].Time() , 0 );
  BOOST_CHECK_EQUAL( FlashVector[1].Time() , 1e6 );

}

BOOST_AUTO_TEST_CASE(MarkFlashesForRemoval_RemoveOneFlash)
{
  size_t NFlashes=3;
  size_t BeginFlash=0;

  std::vector<double> PEs(30,0);
  PEs.at(0) = 100;
  std::vector<double> PEs_Small(30,0);
  PEs_Small.at(0) = 5;
  std::vector<double> WireCenters(3,0);
  std::vector<double> WireWidths(3,0);

  std::vector<recob::OpFlash> FlashVector;
  FlashVector.emplace_back(0, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(1.6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs_Small, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(1e6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  std::vector<bool> MarkedForRemoval(NFlashes-BeginFlash,false);

  opdet::MarkFlashesForRemoval(FlashVector,
			       BeginFlash,
			       MarkedForRemoval);

  BOOST_CHECK_EQUAL( MarkedForRemoval.size() , 3U );
  BOOST_CHECK_EQUAL( MarkedForRemoval[0] , false );
  BOOST_CHECK_EQUAL( MarkedForRemoval[1] , true );
  BOOST_CHECK_EQUAL( MarkedForRemoval[2] , false );
  BOOST_CHECK_EQUAL( FlashVector.size() , 3U );
  BOOST_CHECK_EQUAL( FlashVector[0].Time() , 0 );
  BOOST_CHECK_EQUAL( FlashVector[0].TotalPE() , 100 );
  BOOST_CHECK_EQUAL( FlashVector[1].Time() , 1.6 );
  BOOST_CHECK_EQUAL( FlashVector[1].TotalPE() , 5 );
  BOOST_CHECK_EQUAL( FlashVector[2].Time() , 1e6 );
  BOOST_CHECK_EQUAL( FlashVector[2].TotalPE() , 100 );

}

BOOST_AUTO_TEST_CASE(MarkFlashesForRemoval_IgnoreFirstFlash)
{
  size_t NFlashes=4;
  size_t BeginFlash=1;

  std::vector<double> PEs(30,0);
  PEs.at(0) = 100;
  std::vector<double> PEs_Small(30,0);
  PEs_Small.at(0) = 5;
  std::vector<double> WireCenters(3,0);
  std::vector<double> WireWidths(3,0);

  std::vector<recob::OpFlash> FlashVector;
  FlashVector.emplace_back(-1e6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(0, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(1.6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs_Small, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(1e6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  std::vector<bool> MarkedForRemoval(NFlashes-BeginFlash,false);

  opdet::MarkFlashesForRemoval(FlashVector,
			       BeginFlash,
			       MarkedForRemoval);

  BOOST_CHECK_EQUAL( MarkedForRemoval.size() , 3U );
  BOOST_CHECK_EQUAL( MarkedForRemoval[0] , false );
  BOOST_CHECK_EQUAL( MarkedForRemoval[1] , true );
  BOOST_CHECK_EQUAL( MarkedForRemoval[2] , false );
  BOOST_CHECK_EQUAL( FlashVector.size() , 4U );
  BOOST_CHECK_EQUAL( FlashVector[0].Time() , -1e6 );
  BOOST_CHECK_EQUAL( FlashVector[0].TotalPE() , 100 );
  BOOST_CHECK_EQUAL( FlashVector[1].Time() , 0 );
  BOOST_CHECK_EQUAL( FlashVector[1].TotalPE() , 100 );
  BOOST_CHECK_EQUAL( FlashVector[2].Time() , 1.6 );
  BOOST_CHECK_EQUAL( FlashVector[2].TotalPE() , 5 );
  BOOST_CHECK_EQUAL( FlashVector[3].Time() , 1e6 );
  BOOST_CHECK_EQUAL( FlashVector[3].TotalPE() , 100 );

}

BOOST_AUTO_TEST_CASE(RemoveFlashesFromVectors_IgnoreFirstFlash)
{
  size_t NFlashes=4;
  size_t BeginFlash=1;

  std::vector<double> PEs(30,0);
  PEs.at(0) = 100;
  std::vector<double> PEs_Small(30,0);
  PEs_Small.at(0) = 5;
  std::vector<double> WireCenters(3,0);
  std::vector<double> WireWidths(3,0);

  std::vector<recob::OpFlash> FlashVector;
  FlashVector.emplace_back(-1e6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(0, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(1.6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs_Small, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  FlashVector.emplace_back(1e6, //time
			   0.5, //TimeWidth,
			   0, //AveAbsTime,
			   0, //Frame,
			   PEs, 
			   0, //InBeamFrame,
			   0, //OnBeamTime,
			   0, //FastToTotal,
			   0, //meany, 
			   0, //widthy, 
			   0, //meanz, 
			   0, //widthz, 
			   WireCenters, 
			   WireWidths);
  std::vector<bool> MarkedForRemoval{false, true, false};
  std::vector< std::vector<int> > RefinedHitsPerFlash(NFlashes-BeginFlash);
  RefinedHitsPerFlash[0].push_back(0);
  RefinedHitsPerFlash[1].push_back(1);
  RefinedHitsPerFlash[2].push_back(2);

  opdet::RemoveFlashesFromVectors(MarkedForRemoval,
				  FlashVector,
				  BeginFlash,
				  RefinedHitsPerFlash);

  BOOST_CHECK_EQUAL( FlashVector.size() , 3U );
  BOOST_CHECK_EQUAL( FlashVector[0].Time() , -1e6 );
  BOOST_CHECK_EQUAL( FlashVector[0].TotalPE() , 100 );
  BOOST_CHECK_EQUAL( FlashVector[1].Time() , 0 );
  BOOST_CHECK_EQUAL( FlashVector[1].TotalPE() , 100 );
  BOOST_CHECK_EQUAL( FlashVector[2].Time() , 1e6 );
  BOOST_CHECK_EQUAL( FlashVector[2].TotalPE() , 100 );

  BOOST_CHECK_EQUAL( RefinedHitsPerFlash.size() , 2U );
  BOOST_CHECK_EQUAL( RefinedHitsPerFlash[0][0] , 0 );
  BOOST_CHECK_EQUAL( RefinedHitsPerFlash[1][0] , 2 );

}

BOOST_AUTO_TEST_CASE(RemoveFlashesFromVectors_NoFlashes)
{
  size_t NFlashes=5;
  size_t BeginFlash=5; //that is, no new flashes to look at

  std::vector<recob::OpFlash> FlashVector(NFlashes);
  std::vector<bool> MarkedForRemoval(NFlashes-BeginFlash,false);
  std::vector< std::vector<int> > RefinedHitsPerFlash(NFlashes-BeginFlash);

  opdet::RemoveFlashesFromVectors(MarkedForRemoval,
				  FlashVector,
				  BeginFlash,
				  RefinedHitsPerFlash);

  BOOST_CHECK_EQUAL( FlashVector.size() , NFlashes );

}


BOOST_AUTO_TEST_SUITE_END()
