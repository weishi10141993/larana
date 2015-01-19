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


    void SetOutputTree(TTree*);

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

    void FillFlashHypothesis(FlashHypothesisCollection const&);
    void FillSimPhotonCounter(std::vector<float> const&, std::vector<float> const&);
    void FillComparisonInfo(std::vector<float> const&, float const&,
			    std::vector<float> const&, float const&);
    
  };
  
}


#endif
