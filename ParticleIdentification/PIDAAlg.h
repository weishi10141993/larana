#ifndef PIDAALG_H
#define PIDAALG_H
/*!
 * Title:   PIDA Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on ideas/code from Bruce Baller
 *
 * Description: Algorithm that calculates the PIDA from a calorimetry object
 * Input:       anab::Calorimetry
 * Output:      PIDA information 
*/
#include "fhiclcpp/ParameterSet.h"
#include "AnalysisBase/Calorimetry.h"

namespace util{
  class NormalDistribution;
}

class util::NormalDistribution{

 public:
  NormalDistribution() {}
  NormalDistribution(float,float);

  float getValue(float);

 private:
  float fStepSize;
  float fMaxSigma;
  std::vector<float> fValues;
  

};

namespace pid{
  class PIDAAlg;
}

class pid::PIDAAlg{
 public:
 PIDAAlg(fhicl::ParameterSet const& p):
  fPIDA_BOGUS(-9999)
    { this->reconfigure(p); }

  void reconfigure(fhicl::ParameterSet const& p);

  void RunPIDAAlg(std::vector<double> const&, std::vector<double> const&);
  void RunPIDAAlg(anab::Calorimetry const&);
  void RunPIDAAlg(anab::Calorimetry const&, float&, float&);

  float getPIDAMean();
  float getPIDASigma();

  void setExponentConstant(float const& ex) { fExponentConstant = ex; }

 private:

  const float fPIDA_BOGUS;

  float fExponentConstant;
  float fMinResRange;
  float fMaxResRange;

  std::vector<float> fpida_values;
  std::vector<float> fpida_errors;
  float fpida_mean;
  float fpida_sigma;

  void calculatePIDAMean();
  void calculatePIDASigma();

  void ClearInternalData();

  util::NormalDistribution fnormalDist;

};

#endif
