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

  float getPIDAKDEMostProbable();
  float getPIDAKDEFullWidthHalfMax();

  void PrintPIDAValues();

  void setExponentConstant(float const& ex) { fExponentConstant = ex; }

  void SetPIDATree(TTree*,TH1F*,TH1F*);
  void FillPIDATree(unsigned int, unsigned int, unsigned int, anab::Calorimetry const&);

 private:

  const float fPIDA_BOGUS;

  float fExponentConstant;
  float fMinResRange;
  float fMaxResRange;
  float fMaxPIDAValue;
  float fKDEEvalMaxSigma;
  float fKDEEvalStepSize;
  float fKDEDefaultBandwidth;

  std::vector<float> fpida_values;
  std::vector<float> fpida_errors;
  float fpida_mean;
  float fpida_sigma;

  void calculatePIDAMean();
  void calculatePIDASigma();

  void ClearInternalData();

  void createKDE();
  void calculatePIDAKDEMostProbable();
  void calculatePIDAKDEFullWidthHalfMax();
  std::vector<float> fkde_distribution;
  float fpida_kde_mp;
  float fpida_kde_fwhm;

  float fkde_dist_min;
  float fkde_dist_max;
  
  util::NormalDistribution fnormalDist;


  TTree*             fPIDATree;
  TH1F*              hPIDAvalues;
  TH1F*              hPIDAKDE;
  unsigned int       fPIDAHistNbins;
  float              fPIDAHistMin;
  float              fPIDAHistMax;  
  typedef struct PIDAProperties{
    unsigned int run;
    unsigned int event;
    unsigned int calo_index;
    unsigned int planeid;

    unsigned int n_pid_pts;
    float pida_mean;
    float pida_sigma;
    float pida_kde_mp;
    float pida_kde_fwhm;

    std::string leaf_structure;
    PIDAProperties():
    leaf_structure("run/i:event/i:calo_index/i:planeid/i:n_pid_pts/i:pida_mean/F:pida_sigma/F:pida_kde_mp/F:pida_kde_fwhm/F"){}

  } PIDAProperties_t;
  PIDAProperties_t fPIDAProperties;
  void FillPIDAProperties(unsigned int, unsigned int, unsigned int, unsigned int);

};

#endif
