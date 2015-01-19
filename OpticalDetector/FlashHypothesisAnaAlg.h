#ifndef FLASHHYPOTHESISANAALG_H
#define FLASHHYPOTHESISANAALG_H

/*!
 * Title:   FlashHypothesisAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Alg that compares the flash hypotheses with truth photons and stores the 
 * results in a TTree.
 * 
*/

#include "fhiclcpp/ParameterSet.h"

#include "FlashHypothesisAlg.h"
#include "SimPhotonCounterAlg.h"

#include "TTree.h"

namespace opdet{

  class FlashHypothesisAnaAlg{

  public:
    FlashHypothesisAnaAlg(fhicl::ParameterSet const& p):
      fSPCAlg(p.get<fhicl::ParameterSet>("SimPhotonCounterAlgParams")),
      fCounterIndex(p.get<fhicl::ParameterSet>("SimPhotonCounterIndex")){}


    void SetOutputTree(TTree*,
		       TH1F*, TH1F*,
		       TH1F*, TH1F*,
		       TH1F*, TH1F*,
		       TH1F*, TH1F*,
		       TH1F*, TH1F*,
		       geo::Geometry const&);
		       
    void FillAnaTree(std::vector<TVector3> const& trajVector, 
		     std::vector<float> const& dEdxVector,
		     sim::SimPhotonsCollection const&,
		     geo::Geometry const& geom,
		     opdet::OpDigiProperties const& opdigip,
		     phot::PhotonVisibilityService const& pvs,
		     util::LArProperties const& larp,
		     float XOffset=0);
		     
    
  private:

    unsigned int        fCounterIndex;

    FlashHypothesisAlg  fFHAlg;
    SimPhotonCounterAlg fSPCAlg;

    TTree *fTree;

    typedef struct {
      unsigned int run;
      unsigned int event;

      float comp_tot_prompt;
      float comp_tot_late;
    } FlashHypothesisAnaAlg_t;
    
    FlashHypothesisAnaAlg_t fStruct;

    TH1F* hhyp_prompt;
    TH1F* hhyp_prompt_err;
    TH1F* hhyp_late;
    TH1F* hhyp_late_err;
    TH1F* hhyp_total;
    TH1F* hhyp_total_err;

    TH1F* hsim_prompt;
    TH1F* hsim_late;

    TH1F* hcomp_prompt;
    TH1F* hcomp_late;

    /*
    std::vector<float> fhyp_prompt;
    std::vector<float> fhyp_prompt_err;
    std::vector<float> fhyp_late;
    std::vector<float> fhyp_late_err;
    std::vector<float> fhyp_total;
    std::vector<float> fhyp_total_err;

    std::vector<float> fsim_prompt;
    std::vector<float> fsim_late;

    std::vector<float> fcomp_prompt;
    std::vector<float> fcomp_late;
    */
    void FillFlashHypothesis(FlashHypothesisCollection const&);
    void FillSimPhotonCounter(std::vector<float> const&, std::vector<float> const&);
    void FillComparisonInfo(std::vector<float> const&, float const&,
			    std::vector<float> const&, float const&);
    
  };
  
}


#endif
