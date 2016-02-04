#define BOOST_TEST_MODULE ( AlgoThreshold_test )
#include "boost/test/auto_unit_test.hpp"

#include "fhiclcpp/ParameterSet.h"
#include "OpticalDetector/OpHitFinder/AlgoThreshold.h"


struct AlgoThresholdFixture{

  AlgoThresholdFixture() : myAlgoThreshold() {};
  pmtana::AlgoThreshold myAlgoThreshold;

};

double const tolerance = 1e-6;

BOOST_FIXTURE_TEST_SUITE(AlgoThreshold_test, AlgoThresholdFixture)

BOOST_AUTO_TEST_CASE(checkZeroVector)
{ 

  std::vector<short> wf(20,0);
  std::vector<double> ped_mean(20,0);
  std::vector<double> ped_sigma(20,0.1);
  myAlgoThreshold.Reconstruct(wf,ped_mean,ped_sigma);

  BOOST_CHECK_EQUAL(myAlgoThreshold.GetNPulse(),0ul);

  //BOOST_FAIL("This test fails"); 

}

BOOST_AUTO_TEST_CASE(checkNPulse)
{
  std::vector<short> wf(20,0);
  std::vector<double> ped_mean(20,0);
  std::vector<double> ped_sigma(20,0.1);
  wf[10]=10;

  myAlgoThreshold.Reconstruct(wf,ped_mean,ped_sigma);
  BOOST_CHECK_EQUAL(myAlgoThreshold.GetNPulse(),1ul);
  //BOOST_CHECK_NE(myAlgoThreshold.GetPulse(0),(void*)0);
  //BOOST_CHECK_EQUAL(myAlgoThreshold.GetPulse(1),(void*)0);

}
BOOST_AUTO_TEST_CASE(checkSquarePulse)
{ 

  std::vector<short> wf(20,0);
  std::vector<double> ped_mean(20,0);
  std::vector<double> ped_sigma(20,0.1);
  //

  double area = 0;
  for(size_t iter=0; iter<wf.size(); iter++){
    if(iter>5 && iter<15) {
      wf[iter] += 10;
      area += wf[iter];
    }
  }

  myAlgoThreshold.Reconstruct(wf,ped_mean,ped_sigma);

  BOOST_CHECK_EQUAL(myAlgoThreshold.GetNPulse(),1ul);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_start,5,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_end,15,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_max,6,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).area,area,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).peak,10.0,tolerance);
}

BOOST_AUTO_TEST_CASE(checkTrianglePulse)
{ 

  std::vector<short> wf(20,0);
  std::vector<double> ped_mean(20,0);
  std::vector<double> ped_sigma(20,0.1);
  //

  double area = 0;
  for(size_t iter=0; iter<wf.size(); iter++){
    if(iter<=10)
      wf[iter] += iter;
    else if(iter>10)
      wf[iter] += 20-iter;

    if(wf[iter]>=3) 
      area += wf[iter];

  }

  myAlgoThreshold.Reconstruct(wf,ped_mean,ped_sigma);
  BOOST_CHECK_EQUAL(myAlgoThreshold.GetNPulse(),1ul);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_start,1,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_end,19,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_max,10,tolerance);
  //BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).area,area,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).peak,10.0,tolerance);
}

BOOST_AUTO_TEST_CASE(checkNonZeroPed)
{ 

  double ped = 2;
  std::vector<short> wf(20,(short)ped);
  std::vector<double> ped_mean(20,ped);
  std::vector<double> ped_sigma(20,0.1);

  double area = 0;
  for(size_t iter=0; iter<wf.size(); iter++){
    if(iter<=10)
      wf[iter] += iter;
    else if(iter>10)
      wf[iter] += 20-iter;

    if(wf[iter]>=3+ped) 
      area += wf[iter] - ped;

  }

  myAlgoThreshold.Reconstruct(wf,ped_mean,ped_sigma);
  BOOST_CHECK_EQUAL(myAlgoThreshold.GetNPulse(),1ul);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_start,0,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_end,19,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_max,10,tolerance);
  //BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).area,area,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).peak,10.0,tolerance);
}

BOOST_AUTO_TEST_CASE(checkPulseOffEnd)
{ 

  std::vector<short> wf(20,0);
  std::vector<double> ped_mean(20,0);
  std::vector<double> ped_sigma(20,0.1);
  

  wf[18] = 5; wf[19] = 10;
  double area = 15;

  myAlgoThreshold.Reconstruct(wf,ped_mean,ped_sigma);
  BOOST_CHECK_EQUAL(myAlgoThreshold.GetNPulse(),1ul);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_start,17,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_end,19,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_max,19,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).area,area,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).peak,10.0,tolerance);
}

BOOST_AUTO_TEST_CASE(checkPulseOffFront)
{ 

  std::vector<short> wf(20,0);
  std::vector<double> ped_mean(20,0);
  std::vector<double> ped_sigma(20,0.1);
  

  wf[0] = 10; wf[1] = 5;
  double area = 15;

  myAlgoThreshold.Reconstruct(wf,ped_mean,ped_sigma);
  BOOST_CHECK_EQUAL(myAlgoThreshold.GetNPulse(),1ul);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_start,0,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_end,2,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_max,0,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).area,area,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).peak,10.0,tolerance);
}

BOOST_AUTO_TEST_CASE(checkDoublePulse)
{ 

  std::vector<short> wf(20,0);
  std::vector<double> ped_mean(20,0);
  std::vector<double> ped_sigma(20,0.1);
  

  wf[4] = 5; wf[5] = 10; wf[6] = 5;
  wf[14] = 5; wf[15] = 10; wf[16] = 5;
  double area = 20;

  myAlgoThreshold.Reconstruct(wf,ped_mean,ped_sigma);
  BOOST_CHECK_EQUAL(myAlgoThreshold.GetNPulse(),2ul);

  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_start,3,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_end,7,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).t_max,5,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).area,area,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(0).peak,10.0,tolerance);

  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(1).t_start,13,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(1).t_end,17,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(1).t_max,15,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(1).area,area,tolerance);
  BOOST_CHECK_CLOSE(myAlgoThreshold.GetPulse(1).peak,10.0,tolerance);

}

BOOST_AUTO_TEST_SUITE_END()
