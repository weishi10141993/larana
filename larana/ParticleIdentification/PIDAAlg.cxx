/*!
 * Title:   PIDA Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on ideas/code from Bruce Baller
 *
 * Description: Algorithm that calculates the PIDA from a calorimetry object
 * Input:       anab::Calorimetry
 * Output:      PIDA information 
*/

#include <iostream>
#include <sstream>
#include <cmath>

#include "PIDAAlg.h"

void pid::PIDAAlg::reconfigure(fhicl::ParameterSet const& p){
  fExponentConstant = p.get<float>("ExponentConstant",0.42);
  fMinResRange      = p.get<float>("MinResRange",0);
  fMaxResRange      = p.get<float>("MaxResRange",30);
  fMaxPIDAValue     = p.get<float>("MaxPIDAValue",50);

  fKDEEvalMaxSigma      = p.get<float>("KDEEvalMaxSigma",3);
  fKDEEvalStepSize      = p.get<float>("KDEEvalStepSize",0.01);
  fKDEBandwidths        = p.get< std::vector<float> >("KDEBandwidths");

  fnormalDist = util::NormalDistribution(fKDEEvalMaxSigma,fKDEEvalStepSize);

  fPIDAHistNbins = p.get<unsigned int>("PIDAHistNbins",100);
  fPIDAHistMin   = p.get<float>("PIDAHistMin",0.0);
  fPIDAHistMax   = p.get<float>("PIDAHistMax",50.0);

  ClearInternalData();
}

void pid::PIDAAlg::ClearInternalData(){
  fpida_mean = fPIDA_BOGUS;
  fpida_sigma = fPIDA_BOGUS;
  fpida_integral_dedx = fPIDA_BOGUS;
  fpida_integral_pida = fPIDA_BOGUS;
  fpida_values.clear();
  fpida_errors.clear();

  fpida_kde_mp      = std::vector<float>(fKDEBandwidths.size(),fPIDA_BOGUS);
  fpida_kde_fwhm    = std::vector<float>(fKDEBandwidths.size(),fPIDA_BOGUS);
  fpida_kde_b       = std::vector<float>(fKDEBandwidths.size(),fPIDA_BOGUS);
  fkde_distribution = std::vector< std::vector<float> >(fKDEBandwidths.size());
  fkde_dist_min     = std::vector<float>(fKDEBandwidths.size(),fPIDA_BOGUS);
  fkde_dist_max     = std::vector<float>(fKDEBandwidths.size(),fPIDA_BOGUS);
}

void pid::PIDAAlg::SetPIDATree(TTree *tree, TH1F* hist_vals, std::vector<TH1F*> hist_kde){

  if(hist_kde.size()>MAX_BANDWIDTHS)
    throw "Error: input histograms larger than max allowed bandwidths.";
  if(hist_kde.size()!=fKDEBandwidths.size())
    throw "Error: input histograms do not have same size as bandwidths.";

  fPIDATree = tree;

  hPIDAvalues = hist_vals;
  hPIDAvalues->SetNameTitle("hPIDAvalues","PIDA Distribution");
  hPIDAvalues->SetBins(fPIDAHistNbins,fPIDAHistMin,fPIDAHistMax);

  for(size_t i_hist=0; i_hist<hist_kde.size(); i_hist++){
    hPIDAKDE[i_hist] = hist_kde[i_hist];

    std::stringstream hname,htitle;
    hname << "hPIDAKDE_" << i_hist;
    htitle << "PIDA KDE-smoothed Distribution, Bandwidth=" << fKDEBandwidths.at(i_hist);

    hPIDAKDE[i_hist]->SetNameTitle(hname.str().c_str(),htitle.str().c_str());
    hPIDAKDE[i_hist]->SetBins(fPIDAHistNbins,fPIDAHistMin,fPIDAHistMax);
  }

  fPIDATree->Branch("pida",&fPIDAProperties,fPIDAProperties.leaf_structure.c_str());
  fPIDATree->Branch("hpida_vals","TH1F",hPIDAvalues);
  fPIDATree->Branch("n_bandwidths",&(fPIDAProperties.n_bandwidths),"n_bandwidths/i");
  fPIDATree->Branch("kde_bandwidth",fPIDAProperties.kde_bandwidth,"kde_bandwidth[n_bandwidths]/F");
  fPIDATree->Branch("kde_mp",fPIDAProperties.kde_mp,"kde_mp[n_bandwidths]/F");
  fPIDATree->Branch("kde_fwhm",fPIDAProperties.kde_fwhm,"kde_fwhm[n_bandwidths]/F");
  for(size_t i_hist=0; i_hist<hist_kde.size(); i_hist++){  
    std::stringstream bname;
    bname << "hpida_kde_" << i_hist;
    fPIDATree->Branch(bname.str().c_str(),"TH1F",hPIDAKDE[i_hist]);
  }

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

float pid::PIDAAlg::getPIDAKDEMostProbable(const size_t i_b){
  if(fpida_kde_mp[i_b]==fPIDA_BOGUS)
    createKDE(i_b);

  return fpida_kde_mp[i_b];
}

float pid::PIDAAlg::getPIDAKDEFullWidthHalfMax(const size_t i_b){
  if(fpida_kde_fwhm[i_b]==fPIDA_BOGUS)
    createKDE(i_b);
  
  return fpida_kde_fwhm[i_b];
}

const std::vector<float>& pid::PIDAAlg::getPIDAValues(){
  return fpida_values;
}

const std::vector<float>& pid::PIDAAlg::getPIDAErrors(){
  return fpida_errors;
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

  std::map<double,double> range_dEdx_map;

  for(size_t i_r=0; i_r<resRange.size(); i_r++){
    if(resRange[i_r]>fMaxResRange || resRange[i_r]<fMinResRange) continue;

    range_dEdx_map[ resRange[i_r] ] = dEdx[i_r];
    
    float val = dEdx[i_r]*std::pow(resRange[i_r],fExponentConstant);
    if(val < fMaxPIDAValue){
      fpida_values.push_back(val);
      //fpida_errors.push_back(0);
    }
    
  }

  calculatePIDAIntegral(range_dEdx_map);

  if(fpida_values.size()==0)
    fpida_values.push_back(-99);

}

void pid::PIDAAlg::FillPIDATree(unsigned int run, 
				unsigned int event, 
				unsigned int calo_index, 
				anab::Calorimetry const& calo){
  RunPIDAAlg(calo);
  FillPIDAProperties(run,event,calo_index,calo);
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

void pid::PIDAAlg::calculatePIDAIntegral(std::map<double,double> const& range_dEdx_map){
  
  if(range_dEdx_map.size()<2) return;

  fpida_integral_dedx = 0;
  
  for( std::map<double,double>::const_iterator map_iter = range_dEdx_map.begin();
       map_iter != std::prev(range_dEdx_map.end());
       map_iter++)
    {
      double range_width = std::next(map_iter)->first - map_iter->first;
      fpida_integral_dedx +=  range_width*( std::next(map_iter)->second + 0.5*(map_iter->second-std::next(map_iter)->second));
    }

  fpida_integral_pida = fpida_integral_dedx * (1-fExponentConstant) * 
    std::pow( (std::prev(range_dEdx_map.end())->first - range_dEdx_map.begin()->first),(fExponentConstant-1) );
}

void pid::PIDAAlg::createKDE(const size_t i_b){

  if(fpida_values.size()==0)
    throw "pid::PIDAAlg --- PIDA Values not filled!";

  //if( fkde_distribution[i_b].size()!=0 ) return;

  if(fKDEBandwidths[i_b]<=0) {
    calculatePIDASigma();
    float bandwidth = fpida_sigma*1.06*std::pow((float)(fpida_values.size()),-0.2);
    fpida_errors = std::vector<float>(fpida_values.size(),bandwidth);
    fpida_kde_b[i_b] = bandwidth;
  }
  else{
    fpida_errors = std::vector<float>(fpida_values.size(),fKDEBandwidths[i_b]);
    fpida_kde_b[i_b] = fKDEBandwidths[i_b];
  }

  const auto min_pida_iterator = std::min_element(fpida_values.begin(),fpida_values.end());
  const size_t min_pida_location = std::distance(fpida_values.begin(),min_pida_iterator);
  fkde_dist_min[i_b] = fpida_values[min_pida_location] - fKDEEvalMaxSigma*fpida_errors[min_pida_location];

  const auto max_pida_iterator = std::max_element(fpida_values.begin(),fpida_values.end());
  const size_t max_pida_location = std::distance(fpida_values.begin(),max_pida_iterator);
  fkde_dist_max[i_b] = fpida_values[max_pida_location] + fKDEEvalMaxSigma*fpida_errors[max_pida_location];

  //make the kde distribution, and get the max value
  const size_t kde_dist_size = (size_t)( (fkde_dist_max[i_b] - fkde_dist_min[i_b])/fKDEEvalStepSize ) + 1;
  fkde_distribution[i_b].resize(kde_dist_size);
  float kde_max=0;
  size_t step_max=0;
  for(size_t i_step=0; i_step<kde_dist_size; i_step++){
    float pida_val = fkde_dist_min[i_b] + i_step*fKDEEvalStepSize;
    fkde_distribution[i_b][i_step]=0;

    for(size_t i_pida=0; i_pida<fpida_values.size(); i_pida++)
      fkde_distribution[i_b][i_step] += fnormalDist.getValue((fpida_values[i_pida]-pida_val)/fpida_errors[i_pida])/fpida_errors[i_pida];

    if(fkde_distribution[i_b][i_step]>kde_max){
      kde_max = fkde_distribution[i_b][i_step];
      step_max = i_step;
      fpida_kde_mp[i_b] = pida_val;
    }
  }

  //now get fwhm
  float half_max = 0.5*fpida_kde_mp[i_b];
  float low_width=0;
  for(size_t i_step=step_max; i_step>0; i_step--){
    if(fkde_distribution[i_b][i_step] < half_max) break;
    low_width += fKDEEvalStepSize;
  }
  float high_width=0;
  for(size_t i_step=step_max; i_step<kde_dist_size; i_step++){
    if(fkde_distribution[i_b][i_step] < half_max) break;
    high_width += fKDEEvalStepSize;
  }
  fpida_kde_fwhm[i_b] = low_width+high_width;  

}

void pid::PIDAAlg::createKDEs(){
  for(size_t i_b=0; i_b < fKDEBandwidths.size(); i_b++)
    createKDE(i_b);
}

void pid::PIDAAlg::FillPIDAProperties(unsigned int run,
				      unsigned int event,
				      unsigned int calo_index,
				      anab::Calorimetry const& calo){

  fPIDAProperties.run = run;
  fPIDAProperties.event = event;
  fPIDAProperties.calo_index = calo_index;
  fPIDAProperties.planeid = calo.PlaneID().Plane;
  fPIDAProperties.trk_range = calo.Range();
  fPIDAProperties.calo_KE = calo.KineticEnergy();

  fPIDAProperties.n_bandwidths = fKDEBandwidths.size();
  for(size_t i_b=0; i_b<fPIDAProperties.n_bandwidths; i_b++)
  
  calculatePIDASigma();
  fPIDAProperties.n_pid_pts = fpida_values.size();
  fPIDAProperties.mean = fpida_mean;
  fPIDAProperties.sigma = fpida_sigma;

  fPIDAProperties.integral_dedx = fpida_integral_dedx;
  fPIDAProperties.integral_pida = fpida_integral_pida;
  
  createKDEs();
  for(size_t i_b=0; i_b<fPIDAProperties.n_bandwidths; i_b++){
    fPIDAProperties.kde_mp[i_b]        = fpida_kde_mp[i_b];
    fPIDAProperties.kde_fwhm[i_b]      = fpida_kde_fwhm[i_b];
    fPIDAProperties.kde_bandwidth[i_b] = fpida_kde_b[i_b];
  }

  hPIDAvalues->Reset();
  for(auto const& val: fpida_values)
    hPIDAvalues->Fill(val);

  for(size_t i_b=0; i_b<fPIDAProperties.n_bandwidths; i_b++){
    hPIDAKDE[i_b]->Reset();
    for(size_t i_step=0; i_step<fkde_distribution[i_b].size(); i_step++)
      hPIDAKDE[i_b]->AddBinContent(hPIDAKDE[i_b]->FindBin(i_step*fKDEEvalStepSize+fkde_dist_min[i_b]),
				   fkde_distribution[i_b][i_step]);
  }

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
