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
#include "Geometry/Geometry.h"
#include "Geometry/OpDetGeo.h"

#include "FlashHypothesis.h"
#include "FlashHypothesisCreator.h"
#include "SimPhotonCounterAlg.h"
#include "FlashHypothesisComparison.h"

#include "TTree.h"

namespace opdet{

  class FlashHypothesisAnaAlg{

  public:
  FlashHypothesisAnaAlg(fhicl::ParameterSet const& p):
    fCounterIndex(p.get<unsigned int>("SimPhotonCounterIndex",0)),
      fXOffset(p.get<float>("HypothesisXOffset",0.0)),
      fSPCAlg(p.get<fhicl::ParameterSet>("SimPhotonCounterAlgParams")) {}
    
    
    void SetOutputObjects(TTree*,
			  TH1F*,TH1F*,TH1F*,
			  TH1F*,TH1F*,TH1F*,
			  geo::Geometry const&);

    void FillOpDetPositions(geo::Geometry const&);
    
    void RunComparison(const unsigned int run,
		       const unsigned int event,
		       std::vector<TVector3> const& trajVector, 
		       std::vector<float> const& dEdxVector,
		       sim::SimPhotonsCollection const&,
		       geo::Geometry const& geom,
		       opdet::OpDigiProperties const& opdigip,
		       phot::PhotonVisibilityService const& pvs,
		       util::LArProperties const& larp);
    
    
  private:

    unsigned int        fCounterIndex;
    float               fXOffset;

    FlashHypothesisCreator    fFHCreator;
    SimPhotonCounterAlg       fSPCAlg;
    FlashHypothesisComparison fFHCompare;

    std::vector<float> fOpDetPositions_Y;
    std::vector<float> fOpDetPositions_Z;
    
  };
  
}


#endif
