#ifndef FLASHHYPOTHESISCREATOR_H
#define FLASHHYPOTHESISCREATOR_H

/*!
 * Title:   FlashHypothesis Creator Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: class that produces a flash hypothesis for a trajectory.
 * Input:       Trajectory (std::vector<TVector3> objects)
 * Output:      FlashHypotheses
*/

#include <iostream>
#include <numeric>

#include "fhiclcpp/ParameterSet.h"

#include "RecoBase/Track.h"
#include "MCBase/MCTrack.h"

#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "Utilities/LArProperties.h"
#include "OpticalDetector/OpDigiProperties.h"

#include "TVector3.h"

#include "FlashHypothesis.h"
#include "FlashHypothesisCalculator.h"

namespace opdet{
  
  class FlashHypothesisCreator{
    
  public:
    FlashHypothesisCreator() {}
    
    FlashHypothesisCollection GetFlashHypothesisCollection(recob::Track const& track, 
							   std::vector<float> const& dEdxVector,
							   geo::Geometry const& geom,
							   phot::PhotonVisibilityService const& pvs,
							   util::LArProperties const& larp,
							   opdet::OpDigiProperties const& opdigip,
							   float XOffset=0);
    
    FlashHypothesisCollection GetFlashHypothesisCollection(sim::MCTrack const& mctrack, 
							   std::vector<float> const& dEdxVector,
							   geo::Geometry const& geom,
							   phot::PhotonVisibilityService const& pvs,
							   util::LArProperties const& larp,
							   opdet::OpDigiProperties const& opdigip,
							   float XOffset=0);

    FlashHypothesisCollection GetFlashHypothesisCollection(std::vector<TVector3> const& trajVector, 
							   std::vector<float> const& dEdxVector,
							   geo::Geometry const& geom,
							   phot::PhotonVisibilityService const& pvs,
							   util::LArProperties const& larp,
							   opdet::OpDigiProperties const& opdigip,
							   float XOffset=0);
    
    FlashHypothesisCollection GetFlashHypothesisCollection(TVector3 const& pt1, TVector3 const& pt2, 
							   float const& dEdx,
							   geo::Geometry const& geom,
							   phot::PhotonVisibilityService const& pvs,
							   util::LArProperties const& larp,
							   opdet::OpDigiProperties const& opdigip,
							   float XOffset=0);
    
  private: 
    FlashHypothesisCollection CreateFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
							       float const& dEdx,
							       geo::Geometry const& geom,
							       phot::PhotonVisibilityService const& pvs,
							       util::LArProperties const& larp,
							       opdet::OpDigiProperties const& opdigip,
							       float XOffset);
    
    FlashHypothesisCalculator _calc;
    
  };
  
}

#endif
