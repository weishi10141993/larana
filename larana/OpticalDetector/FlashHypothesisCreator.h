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

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "lardata/DetectorInfo/LArProperties.h"
#include "larcorealg/CoreUtils/ProviderPack.h"
#include "larana/OpticalDetector/OpDigiProperties.h"

#include "TVector3.h"

#include "FlashHypothesis.h"
#include "FlashHypothesisCalculator.h"

namespace opdet{
  
  class FlashHypothesisCreator{
    
  public:
    
    /// Set of service providers used in the common(est) interface
    using Providers_t = lar::ProviderPack<geo::GeometryCore, detinfo::LArProperties>;
    
    FlashHypothesisCreator() {}
    
    FlashHypothesisCollection GetFlashHypothesisCollection(recob::Track const& track, 
							   std::vector<float> const& dEdxVector,
							   Providers_t providers,
							   phot::PhotonVisibilityService const& pvs,
							   opdet::OpDigiProperties const& opdigip,
							   float XOffset=0);
    
    FlashHypothesisCollection GetFlashHypothesisCollection(sim::MCTrack const& mctrack, 
							   std::vector<float> const& dEdxVector,
							   Providers_t providers,
							   phot::PhotonVisibilityService const& pvs,
							   opdet::OpDigiProperties const& opdigip,
							   float XOffset=0);

    FlashHypothesisCollection GetFlashHypothesisCollection(std::vector<TVector3> const& trajVector, 
							   std::vector<float> const& dEdxVector,
							   Providers_t providers,
							   phot::PhotonVisibilityService const& pvs,
							   opdet::OpDigiProperties const& opdigip,
							   float XOffset=0);
    
    FlashHypothesisCollection GetFlashHypothesisCollection(TVector3 const& pt1, TVector3 const& pt2, 
							   float const& dEdx,
							   Providers_t providers,
							   phot::PhotonVisibilityService const& pvs,
							   opdet::OpDigiProperties const& opdigip,
							   float XOffset=0);
    
  private: 
    FlashHypothesisCollection CreateFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
							       float const& dEdx,
							       Providers_t providers,
							       phot::PhotonVisibilityService const& pvs,
							       opdet::OpDigiProperties const& opdigip,
							       float XOffset);
    
    FlashHypothesisCalculator _calc;
    
  };
  
}

#endif
