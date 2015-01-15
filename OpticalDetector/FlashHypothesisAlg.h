#ifndef FLASHHYPOTHESISALG_H
#define FLASHHYPOTHESISALG_H

/*!
 * Title:   FlashHypothesis Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Algorithm that produces a flash hypothesis for a trajectory.
 * Input:       Trajectory (std::vector<TVector3> objects)
 * Output:      FlashHypotheses
*/

#include <iostream>
#include <numeric>

#include "fhiclcpp/ParameterSet.h"

#include "RecoBase/Track.h"
#include "SimulationBase/MCParticle.h"

#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "Utilities/LArProperties.h"
#include "OpticalDetector/OpDigiProperties.h"

#include "TVector3.h"

#include "FlashHypothesis.h"
#include "FlashHypothesisCalculator.h"

namespace opdet{


class FlashHypothesisAlg{

 public:
  FlashHypothesisAlg() {}

  FlashHypothesisCollection GetFlashHypothesisCollection(recob::Track const& track, 
							 std::vector<float> dEdxVector,
							 geo::Geometry const& geom,
							 phot::PhotonVisibilityService const& pvs,
							 util::LArProperties const& larp,
							 opdet::OpDigiProperties const& opdigip,
							 float XOffset=0);
  
  FlashHypothesisCollection GetFlashHypothesisCollection(std::vector<TVector3> const& trajVector, 
							 std::vector<float> dEdxVector,
							 geo::Geometry const& geom,
							 phot::PhotonVisibilityService const& pvs,
							 util::LArProperties const& larp,
							 opdet::OpDigiProperties const& opdigip,
							 float XOffset=0);
  
  FlashHypothesisCollection GetFlashHypothesisCollection(TVector3 const& pt1, TVector3 const& pt2, 
							 float dEdx,
							 geo::Geometry const& geom,
							 phot::PhotonVisibilityService const& pvs,
							 util::LArProperties const& larp,
							 opdet::OpDigiProperties const& opdigip,
							 float XOffset=0);

 private: 
  void CreateFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
					float dEdx,
					geo::Geometry const& geom,
					phot::PhotonVisibilityService const& pvs,
					util::LArProperties const& larp,
					opdet::OpDigiProperties const& opdigip,
					float XOffset,
					FlashHypothesis &prompt_flash_hyp,
					FlashHypothesis &late_flash_hyp);
  
  void AddFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
				     float dEdx,
				     geo::Geometry const& geom,
				     phot::PhotonVisibilityService const& pvs,
				     util::LArProperties const& larp,
				     opdet::OpDigiProperties const& opdigip,
				     float XOffset,
				     FlashHypothesis &prompt_flash_hyp,
				     FlashHypothesis &late_flash_hyp);
  
  FlashHypothesisCalculator _calc;
  
};

 class FlashHypothesisCollection{
  
 public:

 FlashHypothesisCollection(FlashHypothesis const& prompt,
			   FlashHypothesis const& late):
  _prompt_hypothesis(prompt), _late_hypothesis(late)
  {
    _total_hypothesis = _prompt_hypothesis + _late_hypothesis;
  }
  
  void Initialize(size_t s)
  {
    _prompt_hypothesis = FlashHypothesis(s);
    _late_hypothesis = FlashHypothesis(s);
    _total_hypothesis = _prompt_hypothesis + _late_hypothesis;
  }

  void Normalize(float const& totalPE, util::LArProperties const& larp);
  
  FlashHypothesis const& GetPromptHypothesis() { return _prompt_hypothesis; }
  FlashHypothesis const& GetLateHypothesis() { return _late_hypothesis; }
  FlashHypothesis const& GetTotalHypothesis() { return _total_hypothesis; }

  float GetPromptPEs() const { return _prompt_hypothesis.GetTotalPEs(); }
  float GetPromptPEsError() const { return _prompt_hypothesis.GetTotalPEsError(); }
  float GetLatePEs() const { return _late_hypothesis.GetTotalPEs(); }
  float GetLatePEsError() const { return _late_hypothesis.GetTotalPEsError(); }
  float GetTotalPEs() const { return _total_hypothesis.GetTotalPEs(); }
  float GetTotalPEsError() const { return _total_hypothesis.GetTotalPEsError(); }
  
 private:
  FlashHypothesis _prompt_hypothesis;
  FlashHypothesis _late_hypothesis;
  FlashHypothesis _total_hypothesis;

  FlashHypothesisAlg _alg;
};


}

#endif
