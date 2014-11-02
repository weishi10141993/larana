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
  fMaxPIDAValue     = p.get<float>("MaxPIDAValue",50);

  fKDEEvalMaxSigma      = p.get<float>("KDEEvalMaxSigma",3);
  fKDEEvalStepSize      = p.get<float>("KDEEvalStepSize",0.01);
  fKDEDefaultBandwidth  = p.get<float>("KDEDefaultBandwidth",-1);

  fnormalDist = util::NormalDistribution(fKDEEvalMaxSigma,fKDEEvalStepSize);

  ClearInternalData();
}

void pid::PIDAAlg::ClearInternalData(){
  fpida_mean = fPIDA_BOGUS;
  fpida_sigma = fPIDA_BOGUS;
  fpida_values.clear();
  fpida_errors.clear();

  fpida_kde_mp = fPIDA_BOGUS;
  fpida_kde_fwhm = fPIDA_BOGUS;
  fkde_distribution.clear();
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

float pid::PIDAAlg::getPIDAKDEMostProbable(){
  if(fpida_kde_mp==fPIDA_BOGUS)
    calculatePIDAKDEMostProbable();

  return fpida_kde_mp;
}

float pid::PIDAAlg::getPIDAKDEFullWidthHalfMax(){
  if(fpida_kde_fwhm==fPIDA_BOGUS)
    calculatePIDAKDEFullWidthHalfMax();

  return fpida_kde_fwhm;
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
    float val = dEdx[i_r]*std::pow(resRange[i_r],fExponentConstant);
    if(val > fMaxPIDAValue)
      fpida_values.push_back(val);
    //fpida_errors.push_back(0);
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

void pid::PIDAAlg::createKDE(){

  if(fpida_values.size()==0)
    throw "pid::PIDAAlg --- PIDA Values not filled!";

  if( fkde_distribution.size()!=0 ) return;

  if(fpida_errors.size()==0){
    if(fKDEDefaultBandwidth<=0) {
      calculatePIDASigma();
      fpida_errors = std::vector<float>(fpida_values.size(),fpida_sigma);
    }
    else
      fpida_errors = std::vector<float>(fpida_values.size(),fKDEDefaultBandwidth);
  }

  const auto min_pida_iterator = std::min_element(fpida_values.begin(),fpida_values.end());
  const size_t min_pida_location = std::distance(fpida_values.begin(),min_pida_iterator);
  const float min_pida_value = fpida_values[min_pida_location] - fKDEEvalMaxSigma*fpida_errors[min_pida_location];

  const auto max_pida_iterator = std::max_element(fpida_values.begin(),fpida_values.end());
  const size_t max_pida_location = std::distance(fpida_values.begin(),max_pida_iterator);
  const float max_pida_value = fpida_values[max_pida_location] + fKDEEvalMaxSigma*fpida_errors[max_pida_location];

  const size_t kde_dist_size = (size_t)( (max_pida_value - min_pida_value)/fKDEEvalStepSize ) + 1;
  fkde_distribution.resize(kde_dist_size);
  float kde_max=0;
  for(size_t i_step=0; i_step<kde_dist_size; i_step++){
    float pida_val = min_pida_value + i_step*fKDEEvalStepSize;
    fkde_distribution[i_step]=0;

    for(size_t i_pida=0; i_pida<fpida_values.size(); i_pida++)
      fkde_distribution[i_step] += fnormalDist.getValue( (fpida_values[i_pida]-pida_val)/fpida_errors[i_pida] );

    if(fkde_distribution[i_step]>kde_max){
      kde_max = fkde_distribution[i_step];
      fpida_kde_mp = pida_val;
    }
  }
  
}

void pid::PIDAAlg::calculatePIDAKDEMostProbable(){
  if(fkde_distribution.size()==0) createKDE();
}

void pid::PIDAAlg::calculatePIDAKDEFullWidthHalfMax(){
  if(fkde_distribution.size()==0) createKDE();

  float half_max = 0.5*fpida_kde_mp;
  const auto max_kde_iterator = std::max_element(fkde_distribution.begin(), fkde_distribution.end());

  float low_width=0;
  for(std::vector<float>::iterator iter=max_kde_iterator;
      iter!=fkde_distribution.begin();
      iter--)
    {
      if(*iter < half_max) break;
      low_width += fKDEEvalStepSize;
    }

  float high_width=0;
  for(std::vector<float>::iterator iter=max_kde_iterator;
      iter!=fkde_distribution.end();
      iter++)
    {
      if(*iter < half_max) break;
      high_width += fKDEEvalStepSize;
    }

  fpida_kde_fwhm = low_width+high_width;  

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
