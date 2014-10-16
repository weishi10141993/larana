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
#include "TObject.h"


namespace cosmic{
  class BeamFlashTrackMatchTaggerAlg;
  //class TestProperties;
}


class TestProperties : public TObject{
 public:
  TestProperties(){}
  virtual ~TestProperties(){}
  unsigned int run;
  unsigned int event;
  ClassDef(TestProperties,1);
};

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

  void RunHypothesisComparison(unsigned int const,
			       unsigned int const,
			       std::vector<recob::OpFlash> const&,
			       std::vector<recob::Track> const&,
			       geo::Geometry const&,
			       phot::PhotonVisibilityService const&,
			       util::LArProperties const&,
			       opdet::OpDigiProperties const&);
  
  /*
  class FlashProperties : public TObject{

  public:
    FlashProperties(){}
    virtual ~FlashProperties(){}
  FlashProperties(unsigned int i, float pe,float y,float sy,float z,float sz, std::vector<float> opdet):
    index(i), totalPE(pe),y(y),sigmay(sy),z(z),sigmaz(sz),opdetVector(opdet) {}
    unsigned int index;
    float totalPE;
    float y;
    float sigmay;
    float z;
    float sigmaz;
    std::vector<float> opdetVector;
    ClassDef(FlashProperties,1);
  };

  class ComparisonProperties : public TObject{
  public:
    ComparisonProperties(){}
    virtual ~ComparisonProperties(){}
  ComparisonProperties(FlashProperties const& f, FlashProperties const& t, float c):
    flash(f), trkhyp(t),chi2(c) {}
    FlashProperties flash;
    FlashProperties trkhyp;
    float chi2;
    ClassDef(ComparisonProperties,1);
  };

  class EventProperties : public TObject{
  public:
    EventProperties(){}
    virtual ~EventProperties(){}
    unsigned int run;
    unsigned int event;
    std::vector<ComparisonProperties> comps;
    ClassDef(EventProperties,1);
  };
  */
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
  
  //EventProperties *cEvent_p;
  TestProperties *cTest_p;

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
  /*
  void FillFlashProperties(std::vector<float> const& opdetVector,
			   FlashProperties & flash_p,
			   geo::Geometry const& geom);
  */
  //debugging functions
  void PrintTrackProperties(recob::Track const&, std::ostream* output=&std::cout);
  void PrintFlashProperties(recob::OpFlash const&, std::ostream* output=&std::cout);
  void PrintHypothesisFlashComparison(std::vector<float> const&,
				      const recob::OpFlash*,
				      CompatibilityResultType,
				      std::ostream* output=&std::cout);

};

#endif
/*
#ifdef _CINT_
#pragma link C++ class FlashProperties+
#pragma link C++ class ComparisonProperties+
#pragma link C++ class EventProperties+
#endif
*/
