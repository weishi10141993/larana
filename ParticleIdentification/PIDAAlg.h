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

#include "TTree.h"
#include "TH1F.h"

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

const unsigned int MAX_BANDWIDTHS=100;

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

  const std::vector<float>& getPIDAValues();
  const std::vector<float>& getPIDAErrors();

  size_t getNKDEBandwidths() { return fKDEBandwidths.size(); }
  float getKDEBandwidth(const size_t i_b) { return fKDEBandwidths.at(i_b); }
  float getPIDAKDEMostProbable(const size_t);
  float getPIDAKDEFullWidthHalfMax(const size_t);

  void PrintPIDAValues();

  void setExponentConstant(float const& ex) { fExponentConstant = ex; }

  void SetPIDATree(TTree*,TH1F*,std::vector<TH1F*>);
  void FillPIDATree(unsigned int, unsigned int, unsigned int, anab::Calorimetry const&);

 private:

  const float fPIDA_BOGUS;

  float fExponentConstant;
  float fMinResRange;
  float fMaxResRange;
  float fMaxPIDAValue;
  float fKDEEvalMaxSigma;
  float fKDEEvalStepSize;
  std::vector<float> fKDEBandwidths;

  std::vector<float> fpida_values;
  std::vector<float> fpida_errors;
  float fpida_mean;
  float fpida_sigma;
  float fpida_integral_dedx;
  float fpida_integral_pida;

  void calculatePIDAMean();
  void calculatePIDASigma();
  void calculatePIDAIntegral(std::map<double,double> const&);

  void ClearInternalData();

  void createKDEs();
  void createKDE(const size_t);
  void calculatePIDAKDEMostProbable();
  void calculatePIDAKDEFullWidthHalfMax();
  std::vector<float> fpida_kde_mp;
  std::vector<float> fpida_kde_fwhm;
  std::vector<float> fpida_kde_b;

  //this is only for making a histogram later ...
  std::vector< std::vector<float> > fkde_distribution;
  std::vector<float> fkde_dist_min;
  std::vector<float> fkde_dist_max;
  
  util::NormalDistribution fnormalDist;

  TTree*         fPIDATree;
  TH1F*          hPIDAvalues;
  TH1F*          hPIDAKDE[MAX_BANDWIDTHS];
  unsigned int   fPIDAHistNbins;
  float          fPIDAHistMin;
  float          fPIDAHistMax;  
  typedef struct PIDAProperties{
    unsigned int run;
    unsigned int event;
    unsigned int calo_index;
    unsigned int planeid;
    float        trk_range;
    float        calo_KE;

    unsigned int n_pid_pts;
    float mean;
    float sigma;
    float integral_dedx;
    float integral_pida;

    unsigned int n_bandwidths;
    float kde_bandwidth[MAX_BANDWIDTHS];
    float kde_mp[MAX_BANDWIDTHS];
    float kde_fwhm[MAX_BANDWIDTHS];

    std::string leaf_structure;
    PIDAProperties():
    leaf_structure("run/i:event/i:calo_index/i:planeid/i:trk_range/F:calo_KE/F:n_pid_pts/i:mean/F:sigma/F:integral_dedx/F:integral_pida/F"){}

  } PIDAProperties_t;
  PIDAProperties_t fPIDAProperties;
  void FillPIDAProperties(unsigned int, unsigned int, unsigned int, anab::Calorimetry const&);

};

#endif
