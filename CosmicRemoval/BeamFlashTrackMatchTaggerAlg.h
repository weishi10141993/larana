#ifndef BEAMFLASHTRACKMATCHTAGGERALG_H
#define BEAMFLASHTRACKMATCHTAGGERALG_H
/*!
 * Title:   Beam Flash<-->Track Match Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on code from Ben Jones
 *
 * Description: Algorithm that compares all tracks to the flash during the 
 *              beam gate, and determines if that track is consistent with
 *              having produced that flash.
 * Input:       recob::OpFlash, recob::Track
 * Output:      anab::CosmicTag (and Assn<anab::CosmicTag,recob::Track>) 
*/
#include <iostream>

#include "fhiclcpp/ParameterSet.h"

#include "RecoBase/OpFlash.h"
#include "RecoBase/Track.h"
#include "AnalysisBase/CosmicTag.h"

#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "Utilities/LArProperties.h"
#include "OpticalDetector/OpDigiProperties.h"

#include "TTree.h"

namespace cosmic{
  class BeamFlashTrackMatchTaggerAlg;
}

class cosmic::BeamFlashTrackMatchTaggerAlg{
 public:
  BeamFlashTrackMatchTaggerAlg(fhicl::ParameterSet const& p);
  void reconfigure(fhicl::ParameterSet const& p);

  //how to run the algorithm
  void RunCompatibilityCheck(std::vector<recob::OpFlash> const&,
			     std::vector<recob::Track> const&,
			     std::vector<anab::CosmicTag>&,
			     std::vector<size_t>&,
			     geo::Geometry const&,
			     phot::PhotonVisibilityService const&,
			     util::LArProperties const&,
			     opdet::OpDigiProperties const&);

  void SetHypothesisComparisonTree(TTree*);

  void RunHypothesisComparison(std::vector<recob::OpFlash> const&,
			       std::vector<recob::Track> const&,
			       TTree &,
			       geo::Geometry const&,
			       phot::PhotonVisibilityService const&,
			       util::LArProperties const&,
			       opdet::OpDigiProperties const&);
  
 private:

  const anab::CosmicTagID_t COSMIC_TYPE_FLASHMATCH;
  const anab::CosmicTagID_t COSMIC_TYPE_OUTSIDEDRIFT;
  const bool DEBUG_FLAG;

  float fMinTrackLength;
  float fMIPdQdx;
  float fSingleChannelCut;
  float fCumulativeChannelThreshold;
  unsigned int fCumulativeChannelCut;
  float fIntegralCut;

  bool fMakeOutsideDriftTags;
  bool fNormalizeHypothesisToFlash;

  TTree*             cTree;
  std::vector<float> cOpDetVector;
  std::vector<float> cTrackHypVector;
  typedef struct cFlashProperties{
  cFlashProperties(float pe,float y,float sy,float z,float sz):
    totalPE(pe),y(y),sigmay(sy),z(z),sigmaz(sz) {}
    float totalPE;
    float y;
    float sigmay;
    float z;
    float sigmaz;
  } cFlashProperties_t;
  cFlashProperties_t cFlash_p;
  cFlashProperties_t cTrkHyp_p;

  typedef struct cEventProperties{
    int run;
    int event;
    float chi2;
  } cEventProperties_t;

  typedef enum CompatibilityResultType{
    kCompatible = 0,
    kSingleChannelCut,
    kCumulativeChannelCut,
    kIntegralCut
  } CompatibilityResultType;

  //core functions
  std::vector<float> GetMIPHypotheses(recob::Track const& track, 
				      geo::Geometry const& geom,
				      phot::PhotonVisibilityService const& pvs,
				      util::LArProperties const&,
				      opdet::OpDigiProperties const&,
				      float XOffset=0);

  CompatibilityResultType CheckCompatibility(std::vector<float> const& lightHypothesis, 
					     const recob::OpFlash* flashPointer);

  bool InDriftWindow(double, double, geo::Geometry const&);

  //debugging functions
  void PrintTrackProperties(recob::Track const&, std::ostream* output=&std::cout);
  void PrintFlashProperties(recob::OpFlash const&, std::ostream* output=&std::cout);
  void PrintHypothesisFlashComparison(std::vector<float> const&,
				      const recob::OpFlash*,
				      CompatibilityResultType,
				      std::ostream* output=&std::cout);

};

#endif
