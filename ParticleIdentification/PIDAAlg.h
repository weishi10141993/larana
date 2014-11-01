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

namespace pid{
  class PIDAAlg;
}

class pid::PIDAAlg{
 public:
 PIDAAlg(fhicl::ParameterSet const& p):
  fPIDA_BOGUS(-9999) 
    { this->reconfigure(p); }

  void reconfigure(fhicl::ParameterSet const& p);

  void RunPIDAAlg(anab::Calorimetry const&);
  void RunPIDAAlg(anab::Calorimetry const&, float&, float&);
  float getPIDAMean();
  float getPIDASigma();

 private:

  const float fPIDA_BOGUS;

  float fExponentConstant;
  float fMaxResRange;

  std::vector<float> fpida_values;
  float fpida_mean;
  float fpida_sigma;

  void calculatePIDAMean();
  void calculatePIDASigma();

};

#endif
