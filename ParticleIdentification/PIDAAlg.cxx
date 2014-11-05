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

  fPIDAHistNbins = p.get<unsigned int>("PIDAHistNbins",100);
  fPIDAHistMin   = p.get<float>("PIDAHistMin",0.0);
  fPIDAHistMax   = p.get<float>("PIDAHistMin",30.0);

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

void pid::PIDAAlg::SetPIDATree(TTree *tree, TH1F* hist_vals, TH1F* hist_kde){
  fPIDATree = tree;
  hPIDAvalues = hist_vals;
  hPIDAKDE = hist_kde;

  hPIDAvalues->SetNameTitle("hPIDAvalues","PIDA Distribution");
  hPIDAvalues->SetBins(fPIDAHistNbins,fPIDAHistMin,fPIDAHistMax);
  hPIDAKDE->SetNameTitle("hPIDAKDE","PIDA KDE-smoothed Distribution");
  hPIDAKDE->SetBins(fPIDAHistNbins,fPIDAHistMin,fPIDAHistMax);

  fPIDATree->Branch("pida",&fPIDAProperties,fPIDAProperties.leaf_structure.c_str());
  fPIDATree->Branch("hpida_vals","TH1F",hPIDAvalues);
  fPIDATree->Branch("hpida_kde","TH1F",hPIDAKDE);
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
    if(val < fMaxPIDAValue)
      fpida_values.push_back(val);
    //fpida_errors.push_back(0);
  }

  if(fpida_values.size()==0)
    fpida_values.push_back(-99);

}

void pid::PIDAAlg::FillPIDATree(unsigned int run, 
				unsigned int event, 
				unsigned int calo_index, 
				anab::Calorimetry const& calo){
  RunPIDAAlg(calo);
  createKDE();
  FillPIDAProperties(run,event,calo_index,calo.PlaneID().Plane);
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
      float bandwidth = fpida_sigma*1.06*std::pow((float)(fpida_values.size()),-0.2);
      fpida_errors = std::vector<float>(fpida_values.size(),bandwidth);
    }
    else
      fpida_errors = std::vector<float>(fpida_values.size(),fKDEDefaultBandwidth);
  }

  const auto min_pida_iterator = std::min_element(fpida_values.begin(),fpida_values.end());
  const size_t min_pida_location = std::distance(fpida_values.begin(),min_pida_iterator);
  fkde_dist_min = fpida_values[min_pida_location] - fKDEEvalMaxSigma*fpida_errors[min_pida_location];

  const auto max_pida_iterator = std::max_element(fpida_values.begin(),fpida_values.end());
  const size_t max_pida_location = std::distance(fpida_values.begin(),max_pida_iterator);
  fkde_dist_max = fpida_values[max_pida_location] + fKDEEvalMaxSigma*fpida_errors[max_pida_location];

  const size_t kde_dist_size = (size_t)( (fkde_dist_max - fkde_dist_min)/fKDEEvalStepSize ) + 1;
  fkde_distribution.resize(kde_dist_size);
  float kde_max=0;
  for(size_t i_step=0; i_step<kde_dist_size; i_step++){
    float pida_val = fkde_dist_min + i_step*fKDEEvalStepSize;
    fkde_distribution[i_step]=0;

    for(size_t i_pida=0; i_pida<fpida_values.size(); i_pida++)
      fkde_distribution[i_step] += fnormalDist.getValue((fpida_values[i_pida]-pida_val)/fpida_errors[i_pida])/fpida_errors[i_pida];

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

void pid::PIDAAlg::FillPIDAProperties(unsigned int run,
				      unsigned int event,
				      unsigned int calo_index,
				      unsigned int planeid){
  fPIDAProperties.run = run;
  fPIDAProperties.event = event;
  fPIDAProperties.calo_index = calo_index;
  fPIDAProperties.planeid = planeid;

  calculatePIDASigma();
  calculatePIDAKDEFullWidthHalfMax();
  fPIDAProperties.pida_mean = fpida_mean;
  fPIDAProperties.pida_sigma = fpida_sigma;
  fPIDAProperties.pida_kde_mp = fpida_kde_mp;
  fPIDAProperties.pida_kde_fwhm = fpida_kde_fwhm;

  hPIDAvalues->Reset();
  hPIDAKDE->Reset();
  for(auto const& val: fpida_values)
    hPIDAvalues->Fill(val);
  for(size_t i_step=0; i_step<fkde_distribution.size(); i_step++)
    hPIDAKDE->AddBinContent(hPIDAKDE->FindBin(i_step*fKDEEvalStepSize+fkde_dist_min),
			    fkde_distribution[i_step]);
  
  fPIDATree->Fill();
}

void pid::PIDAAlg::PrintPIDAValues(){
  for(size_t i_pida=0; i_pida<fpida_values.size(); i_pida++)
    std::cout << "\tPIDA --- " << i_pida << "\t" << fpida_values[i_pida] << std::endl;
}

util::NormalDistribution::NormalDistribution(float max_sigma, float step_size){

  if(step_size==0)
    throw "util::NormalDistribution --- Cannot have zero step size!";

  const size_t vector_size = (size_t)(max_sigma / step_size);
  fValues.resize(vector_size);

  const float AMPLITUDE = 1. / std::sqrt(2*M_PI);

  float integral=0;
  for(size_t i_step=0; i_step<vector_size; i_step++){
    float diff = i_step*step_size;
    fValues[i_step] = AMPLITUDE * std::exp(-0.5*diff*diff);
    integral+= fValues[i_step];
  }

  for(size_t i_step=0; i_step<vector_size; i_step++)
    fValues[i_step] /= (integral*2);

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
