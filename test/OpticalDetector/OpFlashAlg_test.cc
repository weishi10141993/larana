#define BOOST_TEST_MODULE ( OpFlashAlg_test )
#include "boost/test/unit_test.hpp"

#include "larana/OpticalDetector/OpFlashAlg.h"

constexpr float FlashThreshold = 50;
constexpr double WidthTolerance = 0.5;

auto const tolerance = 1e-4% boost::test_tools::tolerance();

BOOST_AUTO_TEST_SUITE(OpFlashAlg_test)

BOOST_AUTO_TEST_CASE(checkGetAccumIndex)
{

  BOOST_TEST(opdet::GetAccumIndex(10.0,0,1,0) == 10ul);
  BOOST_TEST(opdet::GetAccumIndex(10.0,0,2,0) == 5ul);
  BOOST_TEST(opdet::GetAccumIndex(10.0,0,3,0) == 3ul);

  BOOST_TEST(opdet::GetAccumIndex(10.0,-5,1,0) == 15ul);
  BOOST_TEST(opdet::GetAccumIndex(10.0,-5,2,0) == 7ul);
  BOOST_TEST(opdet::GetAccumIndex(10.0,-5,3,0) == 5ul);

  BOOST_TEST(opdet::GetAccumIndex(10.0,-5,1,0.5) == 15ul);
  BOOST_TEST(opdet::GetAccumIndex(10.0,-5,2,1) == 8ul);

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

  BOOST_TEST( Contributors.at(AccumIndex).size() == 1U);
  BOOST_TEST( Contributors.at(AccumIndex).at(0) == (int)HitIndex );
  BOOST_TEST( Binned.at(AccumIndex) == PE+PE_base );
  BOOST_TEST( FlashesInAccumulator.size() == 0U);

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

  BOOST_TEST( Contributors.at(AccumIndex).size() == 1U);
  BOOST_TEST( Contributors.at(AccumIndex).at(0) == (int)HitIndex );
  BOOST_TEST( Binned.at(AccumIndex) == PE+PE_base );
  BOOST_TEST( FlashesInAccumulator.size() == 1U);

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

  BOOST_TEST( Contributors.at(AccumIndex).size() == 1U);
  BOOST_TEST( Contributors.at(AccumIndex).at(0) == (int)HitIndex );
  BOOST_TEST( Binned.at(AccumIndex) == PE+PE_base );
  BOOST_TEST( FlashesInAccumulator.size() == 0U);


  unsigned int HitIndex2 = 1;
  opdet::FillAccumulator(AccumIndex,HitIndex2,PE,FlashThreshold,
                         Binned,Contributors,FlashesInAccumulator);

  BOOST_TEST( Contributors.at(AccumIndex).size() == 2U);
  BOOST_TEST( Contributors.at(AccumIndex).at(0) == (int)HitIndex );
  BOOST_TEST( Contributors.at(AccumIndex).at(1) == (int)HitIndex2 );
  BOOST_TEST( Binned.at(AccumIndex) == PE*2+PE_base );
  BOOST_TEST( FlashesInAccumulator.size() == 1U);

}

BOOST_AUTO_TEST_CASE(FillFlashesBySizeMap_checkNoFlash)
{
  const size_t vector_size = 10;
  const double PE_vals = 10;

  std::vector<double> BinnedPE(vector_size,PE_vals);
  std::vector<int> FlashesInAccumulator;
  std::map<double, std::map<int,std::vector<int> >, std::greater<double> > FlashesBySize;

  opdet::FillFlashesBySizeMap(FlashesInAccumulator,BinnedPE,1,FlashesBySize);

  BOOST_TEST( FlashesBySize.size() == 0U );

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

  BOOST_TEST( FlashesBySize.size() == 1U );
  BOOST_TEST( FlashesBySize.count(50) == 1U );
  BOOST_TEST( FlashesBySize.count(10) == 0U );
  BOOST_TEST( FlashesBySize[50].size() == 1U );
  BOOST_TEST( FlashesBySize[50].count(1) == 1U );
  BOOST_TEST( FlashesBySize[50][1].size() == 1U );
  BOOST_TEST( FlashesBySize[50][1][0] == 2 );

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

  BOOST_TEST( FlashesBySize.size() == 1U );
  BOOST_TEST( FlashesBySize.count(50) == 1U );
  BOOST_TEST( FlashesBySize.count(10) == 0U );

  BOOST_TEST( FlashesBySize[50].size() == 1U );
  BOOST_TEST( FlashesBySize[50].count(1) == 1U );
  BOOST_TEST( FlashesBySize[50][1].size() == 2U );
  BOOST_TEST( FlashesBySize[50][1][0] == 2 );
  BOOST_TEST( FlashesBySize[50][1][1] == 8 );

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

  BOOST_TEST( FlashesBySize.size() == 2U );

  BOOST_TEST( FlashesBySize.count(50) == 1U );
  BOOST_TEST( FlashesBySize.count(60) == 1U );
  BOOST_TEST( FlashesBySize.count(10) == 0U );

  auto map_begin = FlashesBySize.begin();
  BOOST_TEST( map_begin->first == 60 );

  auto map_last = FlashesBySize.end();
  map_last--;
  BOOST_TEST( map_last->first == 50 );

  BOOST_TEST( FlashesBySize[50].size() == 2U );
  BOOST_TEST( FlashesBySize[50].count(1) == 1U );
  BOOST_TEST( FlashesBySize[50].count(2) == 1U );

  BOOST_TEST( FlashesBySize[50][1].size() == 1U );
  BOOST_TEST( FlashesBySize[50][1][0] == 2 );
  BOOST_TEST( FlashesBySize[50][2].size() == 1U );
  BOOST_TEST( FlashesBySize[50][2][0] == 8 );

  BOOST_TEST( FlashesBySize[60].size() == 1U );
  BOOST_TEST( FlashesBySize[60].count(1) == 1U );
  BOOST_TEST( FlashesBySize[60].count(2) == 0U );
  BOOST_TEST( FlashesBySize[60][1].size() == 1U );
  BOOST_TEST( FlashesBySize[60][1][0] == 5 );
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

  BOOST_TEST( HitsThisFlash.size() == 0U);
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

  BOOST_TEST( HitsThisFlash.size() == 2U );
  BOOST_TEST( HitsThisFlash[0] == 1 );
  BOOST_TEST( HitsThisFlash[1] == 3 );
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

  BOOST_TEST( HitsThisFlash.size() == 1U );
  BOOST_TEST( HitsThisFlash[0] == 3 );
}

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

  BOOST_TEST( HitsThisFlash.size() == 2U );
  BOOST_TEST( HitsThisFlash[0] == 5 );
  BOOST_TEST( HitsThisFlash[1] == 6 );
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

  BOOST_TEST( HitsPerFlash.size() == 0U);
  BOOST_TEST( HitClaimedByFlash[0] == -1);
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

  BOOST_TEST( HitsPerFlash.size() == 0U);
  BOOST_TEST( HitClaimedByFlash[0] == -1);
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

  BOOST_TEST( HitsPerFlash.size() == 1U);
  BOOST_TEST( HitsPerFlash[0][0] == 0);
  BOOST_TEST( HitClaimedByFlash[0] == 0);
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

  BOOST_TEST( HitsPerFlash.size() == 1U);
  BOOST_TEST( HitsPerFlash[0][0] == 0);
  BOOST_TEST( HitClaimedByFlash[0] == 0);
  BOOST_TEST( HitClaimedByFlash[1] == -1);
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

  BOOST_TEST( HitsPerFlash.size() == 0U);
  BOOST_TEST( HitClaimedByFlash[0] == -1);
  BOOST_TEST( HitClaimedByFlash[1] == -1);
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

  BOOST_TEST( HitsPerFlash.size() == 1U);
  BOOST_TEST( HitsPerFlash[0][0] == 0);
  BOOST_TEST( HitsPerFlash[0][1] == 1);
  BOOST_TEST( HitClaimedByFlash[0] == 0);
  BOOST_TEST( HitClaimedByFlash[1] == 0);
}

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

  BOOST_TEST( HitsThisRefinedFlash.size() == 0U );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == (int)NHits);
  BOOST_TEST( PEAccumulated == 0 );
  BOOST_TEST( FlashMaxTime == 0 );
  BOOST_TEST( FlashMinTime == 0 );

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

  BOOST_TEST( HitsThisRefinedFlash.size() == 1U );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 1);
  BOOST_TEST( HitsUsed.at(NHits-1) == true);
  BOOST_TEST( PEAccumulated == 60 );
  BOOST_TEST( FlashMaxTime == 9. );
  BOOST_TEST( FlashMinTime == -1. );

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

  BOOST_TEST( HitsThisRefinedFlash.size() == 1U );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 2);
  BOOST_TEST( HitsUsed.at(NHits-1) == true);
  BOOST_TEST( HitsUsed.at(NHits-2) == true);
  BOOST_TEST( PEAccumulated == 50 );
  BOOST_TEST( FlashMaxTime == 8. );
  BOOST_TEST( FlashMinTime == -2. );

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

  BOOST_TEST( HitsThisRefinedFlash.size() == 1U );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 1);
  BOOST_TEST( HitsUsed.at(NHits-1) == true);
  BOOST_TEST( PEAccumulated == 60 );
  BOOST_TEST( FlashMaxTime == 9. );
  BOOST_TEST( FlashMinTime == -1. );

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

  BOOST_TEST( HitsThisRefinedFlash.size() == 2U );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 2);
  BOOST_TEST( HitsUsed.at(NHits-1) == true);
  BOOST_TEST( HitsUsed.at(NHits-2) == true);
  BOOST_TEST( PEAccumulated == 110 );
  BOOST_TEST( FlashMaxTime == 9. );
  BOOST_TEST( FlashMinTime == -2. );

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

  BOOST_TEST( HitsThisRefinedFlash.size() == 1U );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 1);
  BOOST_TEST( HitsUsed.at(NHits-1) == true);
  BOOST_TEST( HitsUsed.at(NHits-2) == false);
  BOOST_TEST( PEAccumulated == 60 );
  BOOST_TEST( FlashMaxTime == 4.05 );
  BOOST_TEST( FlashMinTime == 3.95 );

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

  BOOST_TEST( RefinedHitsPerFlash.size() == 1U );
  BOOST_TEST( RefinedHitsPerFlash[0].size() == 3U );
  BOOST_TEST( RefinedHitsPerFlash[0][0] == 0 );
  BOOST_TEST( RefinedHitsPerFlash[0][1] == 1 );
  BOOST_TEST( RefinedHitsPerFlash[0][2] == 2 );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 3);

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

  BOOST_TEST( RefinedHitsPerFlash.size() == 1U );
  BOOST_TEST( RefinedHitsPerFlash[0].size() == 1U );
  BOOST_TEST( RefinedHitsPerFlash[0][0] == 0 );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 1);

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

  BOOST_TEST( RefinedHitsPerFlash.size() == 0U );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 1);
  BOOST_TEST( HitsUsed[0] == true);

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

  BOOST_TEST( RefinedHitsPerFlash.size() == 0U );
  BOOST_TEST( std::count(HitsUsed.begin(),HitsUsed.end(),true) == 1);
  BOOST_TEST( HitsUsed[0] == true);
  BOOST_TEST( HitsUsed[1] == false);
  BOOST_TEST( HitsUsed[2] == false);

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

    BOOST_TEST( MaxTime == peak_time );
    BOOST_TEST( MinTime == peak_time );
    BOOST_TEST( AveTime == peak_time*hit_pe );
    BOOST_TEST( AveAbsTime == 0 );
    BOOST_TEST( FastToTotal == 0 );
    BOOST_TEST( TotalPE == hit_pe );
    BOOST_TEST( PEs.at(op_channel) == hit_pe );
    BOOST_TEST( PEs.at(4) == 0 );
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

    BOOST_TEST( MaxTime == peak_time );
    BOOST_TEST( MinTime == 1 );
    BOOST_TEST( AveTime == peak_time*hit_pe+20 );
    BOOST_TEST( AveAbsTime == 0 );
    BOOST_TEST( FastToTotal == 0 );
    BOOST_TEST( TotalPE == hit_pe+20 );
    BOOST_TEST( PEs.at(op_channel) == hit_pe );
    BOOST_TEST( PEs.at(0) == 20 );
    BOOST_TEST( PEs.at(4) == 0 );
}


BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_BackwardsTime)
{

  double iPE = 100; double iTime = 0; double iWidth=0.5;
  double jPE = 100; double jTime = -1; double jWidth=0.5;

  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
                                                jPE, jTime, jWidth);

  BOOST_TEST( result == 1e6, tolerance);
}

BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_EqualFlashes)
{

  double iPE = 100; double iTime = 0; double iWidth=0.5;
  double jPE = 100; double jTime = 0; double jWidth=0.5;

  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
                                                jPE, jTime, jWidth);

  BOOST_TEST( result == 0, tolerance);
}

BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_LateFlash)
{

  double iPE = 100; double iTime = 0; double iWidth=0.5;
  double jPE = 10; double jTime = 1.6; double jWidth=0.5;

  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
                                                jPE, jTime, jWidth);

  double good_result = (jPE - std::exp(-1)*iPE)/(std::sqrt(std::exp(-1)*iPE));

  BOOST_TEST( result == good_result, tolerance);
}

BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_VeryLateFlash)
{

  double iPE = 100; double iTime = 0; double iWidth=0.5;
  double jPE = 10; double jTime = 16; double jWidth=0.5;

  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
                                                jPE, jTime, jWidth);

  double good_result = (jPE - std::exp(-10)*iPE)/(std::sqrt(std::exp(-10)*iPE));

  BOOST_TEST( result == good_result, tolerance);
}

BOOST_AUTO_TEST_CASE(GetLikelihoodLateLight_UnequalWidths)
{

  double iPE = 100; double iTime = 0; double iWidth=1;
  double jPE = 10; double jTime = 1.6; double jWidth=0.5;

  double result = opdet::GetLikelihoodLateLight(iPE, iTime, iWidth,
                                                jPE, jTime, jWidth);

  double good_result = (jPE - std::exp(-1)*iPE*0.5)/(std::sqrt(std::exp(-1)*iPE*0.5));

  BOOST_TEST( result == good_result, tolerance);
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

  BOOST_TEST( MarkedForRemoval.size() == 0U );

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

  BOOST_TEST( MarkedForRemoval.size() == 1U );
  BOOST_TEST( MarkedForRemoval[0] == false );
  BOOST_TEST( FlashVector.size() == 1U );
  BOOST_TEST( FlashVector[0].Time() == 0 );
  BOOST_TEST( FlashVector[0].TimeWidth() == 0.5 );

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

  BOOST_TEST( MarkedForRemoval.size() == 2U );
  BOOST_TEST( MarkedForRemoval[0] == false );
  BOOST_TEST( MarkedForRemoval[1] == false );
  BOOST_TEST( FlashVector.size() == 2U );
  BOOST_TEST( FlashVector[0].Time() == 0 );
  BOOST_TEST( FlashVector[1].Time() == 1e6 );

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

  BOOST_TEST( MarkedForRemoval.size() == 3U );
  BOOST_TEST( MarkedForRemoval[0] == false );
  BOOST_TEST( MarkedForRemoval[1] == true );
  BOOST_TEST( MarkedForRemoval[2] == false );
  BOOST_TEST( FlashVector.size() == 3U );
  BOOST_TEST( FlashVector[0].Time() == 0 );
  BOOST_TEST( FlashVector[0].TotalPE() == 100 );
  BOOST_TEST( FlashVector[1].Time() == 1.6 );
  BOOST_TEST( FlashVector[1].TotalPE() == 5 );
  BOOST_TEST( FlashVector[2].Time() == 1e6 );
  BOOST_TEST( FlashVector[2].TotalPE() == 100 );

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

  BOOST_TEST( MarkedForRemoval.size() == 3U );
  BOOST_TEST( MarkedForRemoval[0] == false );
  BOOST_TEST( MarkedForRemoval[1] == true );
  BOOST_TEST( MarkedForRemoval[2] == false );
  BOOST_TEST( FlashVector.size() == 4U );
  BOOST_TEST( FlashVector[0].Time() == -1e6 );
  BOOST_TEST( FlashVector[0].TotalPE() == 100 );
  BOOST_TEST( FlashVector[1].Time() == 0 );
  BOOST_TEST( FlashVector[1].TotalPE() == 100 );
  BOOST_TEST( FlashVector[2].Time() == 1.6 );
  BOOST_TEST( FlashVector[2].TotalPE() == 5 );
  BOOST_TEST( FlashVector[3].Time() == 1e6 );
  BOOST_TEST( FlashVector[3].TotalPE() == 100 );

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

  BOOST_TEST( FlashVector.size() == 3U );
  BOOST_TEST( FlashVector[0].Time() == -1e6 );
  BOOST_TEST( FlashVector[0].TotalPE() == 100 );
  BOOST_TEST( FlashVector[1].Time() == 0 );
  BOOST_TEST( FlashVector[1].TotalPE() == 100 );
  BOOST_TEST( FlashVector[2].Time() == 1e6 );
  BOOST_TEST( FlashVector[2].TotalPE() == 100 );

  BOOST_TEST( RefinedHitsPerFlash.size() == 2U );
  BOOST_TEST( RefinedHitsPerFlash[0][0] == 0 );
  BOOST_TEST( RefinedHitsPerFlash[1][0] == 2 );

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

  BOOST_TEST( FlashVector.size() == NFlashes );

}


BOOST_AUTO_TEST_SUITE_END()
