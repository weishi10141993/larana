/*!
 * Title:   PIDA Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on ideas/code from Bruce Baller
 *
 * Description: Algorithm that calculates the PIDA from a calorimetry object
 * Input:       anab::Calorimetry
 * Output:      PIDA information 
*/

#include <iostream>
#include <cmath>

#include "PIDAAlg.h"

void pid::PIDAAlg::reconfigure(fhicl::ParameterSet const& p){
  fExponentConstant = p.get<float>("ExponentConstant",0.42);
  fMinResRange      = p.get<float>("MinResRange",0);
  fMaxResRange      = p.get<float>("MaxResRange",30);

  fnormalDist = util::NormalDistribution(p.get<float>("KDEMaxSigma",3),p.get<float>("KDEStepSize",0.01));

  ClearInternalData();
}

void pid::PIDAAlg::ClearInternalData(){
  fpida_mean = fPIDA_BOGUS;
  fpida_sigma = fPIDA_BOGUS;
  fpida_values.clear();
  fpida_errors.clear();
}

float pid::PIDAAlg::getPIDAMean(){
  if(fpida_mean==fPIDA_BOGUS)
    calculatePIDAMean();

  return fpida_mean;
}

float pid::PIDAAlg::getPIDASigma(){
  if(fpida_sigma==fPIDA_BOGUS)
    calculatePIDASigma();

  return fpida_sigma;
}

void pid::PIDAAlg::RunPIDAAlg(anab::Calorimetry const& calo){
  std::vector<double> const& resRange = calo.ResidualRange();
  std::vector<double> const& dEdx     = calo.dEdx();
  RunPIDAAlg(resRange,dEdx);
}

void pid::PIDAAlg::RunPIDAAlg(anab::Calorimetry const& calo,
			      float& mean,
			      float& sigma)
{
  RunPIDAAlg(calo);
  mean = getPIDAMean();
  sigma = getPIDASigma();
}

void pid::PIDAAlg::RunPIDAAlg(std::vector<double> const& resRange,
			      std::vector<double> const& dEdx){

  ClearInternalData();

  fpida_values.reserve( resRange.size() );
  fpida_errors.reserve( resRange.size() );

  for(size_t i_r=0; i_r<resRange.size(); i_r++){
    if(resRange[i_r]>fMaxResRange || resRange[i_r]<fMinResRange) continue;
    fpida_values.push_back(dEdx[i_r]*std::pow(resRange[i_r],fExponentConstant));
    fpida_errors.push_back(0);
  }

}

void pid::PIDAAlg::calculatePIDAMean(){

  if(fpida_values.size()==0)
    throw "pid::PIDAAlg --- PIDA Values not filled!";

  fpida_mean = 0;
  for(auto const& val : fpida_values)
    fpida_mean += val;
  fpida_mean /= fpida_values.size();

}

void pid::PIDAAlg::calculatePIDASigma(){

  if(fpida_mean==fPIDA_BOGUS)
    calculatePIDAMean();

  fpida_sigma = 0;
  for(auto const& val : fpida_values)
    fpida_sigma += (fpida_mean-val)*(fpida_mean-val);

  fpida_sigma = std::sqrt(fpida_sigma)/fpida_values.size();
}

util::NormalDistribution::NormalDistribution(float max_sigma, float step_size){

  if(step_size==0)
    throw "util::NormalDistribution --- Cannot have zero step size!";

  const size_t vector_size = (size_t)(max_sigma / step_size);
  fValues.resize(vector_size);

  const float AMPLITUDE = 1. / std::sqrt(2*M_PI);

  for(size_t i_step=0; i_step<vector_size; i_step++){
    float diff = i_step*step_size - 1.;
    fValues[i_step] = AMPLITUDE * std::exp(-0.5*diff*diff);
  }

  fStepSize = step_size;
  fMaxSigma = fStepSize * vector_size;
  
}

float util::NormalDistribution::getValue(float x){

  x = std::abs(x);
  if(x > fMaxSigma) return 0;

  size_t bin_low = x / fStepSize;
  float remainder = (x - (bin_low*fStepSize)) / fStepSize;
  
  return fValues[bin_low]*(1-remainder) + remainder*fValues[bin_low+1];

}
