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
    FlashHypothesisCollection CreateFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
							       float dEdx,
							       geo::Geometry const& geom,
							       phot::PhotonVisibilityService const& pvs,
							       util::LArProperties const& larp,
							       opdet::OpDigiProperties const& opdigip,
							       float XOffset);
    
    void AddFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
				       float dEdx,
				       geo::Geometry const& geom,
				       phot::PhotonVisibilityService const& pvs,
				       util::LArProperties const& larp,
				       opdet::OpDigiProperties const& opdigip,
				       float XOffset,
				       FlashHypothesisCollection &hyp_collection);
    
    FlashHypothesisCalculator _calc;
    
  };
  
}

#endif
