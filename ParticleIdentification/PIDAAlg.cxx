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
  fpida_mean = fPIDA_BOGUS;
  fpida_sigma = fPIDA_BOGUS;
  fpida_values.clear();
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

  fpida_values.clear();
  fpida_values.reserve( resRange.size() );

  for(size_t i_r=0; i_r<resRange.size(); i_r++){
    if(resRange[i_r]>fMaxResRange || resRange[i_r]<fMinResRange) continue;
    fpida_values.push_back(dEdx[i_r]*std::pow(resRange[i_r],fExponentConstant));
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

