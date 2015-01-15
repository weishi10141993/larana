/*!
 * Title:   FlashHypothesis Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Algorithm that produces a flash hypothesis for a trajectory.
 * Input:       Trajectory (std::vector<TVector3> objects)
 * Output:      FlashHypotheses
*/

#include "FlashHypothesisAlg.h"

opdet::FlashHypothesisCollection 
opdet::FlashHypothesisAlg::GetFlashHypothesisCollection(recob::Track const& track, 
							std::vector<float> dEdxVector,
							geo::Geometry const& geom,
							phot::PhotonVisibilityService const& pvs,
							util::LArProperties const& larp,
							opdet::OpDigiProperties const& opdigip,
							float XOffset)
{
  if(track.NumberTrajectoryPoints() != dEdxVector.size())
    throw "ERROR in FlashHypothesisAlg: dEdx vector size not same as track size.";

  FlashHypothesisCollection fhc(geom.NOpDet());
  for(size_t pt=1; pt<track.NumberTrajectoryPoints(); pt++)
    fhc = fhc + _alg.CreateFlashHypothesesFromSegment(track.LocationAtPoint(pt-1),
						      track.LocationAtPoint(pt),
						      dEdxVector[pt],
						      geom,pvs,larp,opdigip,XOffset);
}

opdet::FlashHypothesisCollection 
opdet::FlashHypothesisCollection::FlashHypothesisCollection(std::vector<TVector3> const& trajVector, 
							    std::vector<float> dEdxVector,
							    geo::Geometry const& geom,
							    phot::PhotonVisibilityService const& pvs,
							    util::LArProperties const& larp,
							    opdet::OpDigiProperties const& opdigip,
							    float XOffset)
{
  if(trajVector.size() != dEdxVector.size())
    throw "ERROR in FlashHypothesisCollection: dEdx vector size not same as track size.";

  FlashHypothesisCollection fhc(geom.NOpDet());
  for(size_t pt=1; pt<trajVector.size(); pt++)
    fhc = fhc + _alg.CreateFlashHypothesesFromSegment(trajVector[pt-1],
						      trajVector[pt],
						      dEdxVector[pt],
						      geom,pvs,larp,opdigip,XOffset);

}

opdet::FlashHypothesisCollection 
opdet::FlashHypothesisCollection::GetFlashHypothesisCollection(TVector3 const& pt1, TVector3 const& pt2, 
							       float dEdx,
							       geo::Geometry const& geom,
							       phot::PhotonVisibilityService const& pvs,
							       util::LArProperties const& larp,
							       opdet::OpDigiProperties const& opdigip,
							       float XOffset)
{
  return _alg.CreateFlashHypothesesFromSegment(pt1,pt2,dEdx,geom,pvs,larp,opdigip,XOffset);
}

opdet::FlashHypothesisCollection
opdet::FlashHypothesisAlg::CreateFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
							    float dEdx,
							    geo::Geometry const& geom,
							    phot::PhotonVisibilityService const& pvs,
							    util::LArProperties const& larp,
							    opdet::OpDigiProperties const& opdigip,
							    float XOffset)
{
  FlashHypothesis prompt_hyp = FlashHypothesis(geom.NOpDet());
  
  std::vector<double> xyz_segment(_calc.SegmentMidpoint(pt1,pt2,XOffset));
  
  //get the visibility vector
  const std::vector<float>* PointVisibility = pvs.GetAllVisibilities(xyz_segment[0]);
  
  //check vector size, as it may be zero if given a y/z outside some range
  if(PointVisibility->size()!=geom.NOpDet()) return;
  
  //klugey ... right now, set a qe_vector that gives constant qe across all opdets
  std::vector<float> qe_vector(geom.NOpDet(),opdigip.QE());
  _calc.FillFlashHypotheses(larp.ScintYield()*larp.ScientYieldRatio(),
			    dEdx,
			    pt1,pt2,
			    qe_vector,
			    *PointVisibility,
			    prompt_hyp);
  
  FlashHypothesisCollection fhc;
  fhc.SetPromptHypAndPromptFrac(prompt_hyp,larp.ScientYieldRatio());
  return fhc;
}
