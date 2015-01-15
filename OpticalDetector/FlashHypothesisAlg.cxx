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

  _fhc.Initialize(geom.NOpDet());
  for(size_t pt=1; pt<track.NumberTrajectoryPoints(); pt++)
    _alg.AddFlashHypothesesFromSegment(track.LocationAtPoint(pt-1),
				       track.LocationAtPoint(pt),
				       dEdxVector[pt],
				       geom,pvs,larp,opdigip,XOffset,
				       _fhc._prompt_hypothesis,
				       _fhc._late_hypothesis);

  _total_hypothesis = _prompt_hypothesis+_late_hypothesis;
  
}

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

  Initialize(geom.NOpDet());
  for(size_t pt=1; pt<trajVector.size(); pt++)
    _alg.AddFlashHypothesesFromSegment(trajVector[pt-1],
				       trajVector[pt],
				       dEdxVector[pt],
				       geom,pvs,larp,opdigip,XOffset,
				       _prompt_hypothesis,
				       _late_hypothesis);

  _total_hypothesis = _prompt_hypothesis+_late_hypothesis;
  
}

opdet::FlashHypothesisCollection::FlashHypothesisCollection(TVector3 const& pt1, TVector3 const& pt2, 
							    float dEdx,
							    geo::Geometry const& geom,
							    phot::PhotonVisibilityService const& pvs,
							    util::LArProperties const& larp,
							    opdet::OpDigiProperties const& opdigip,
							    float XOffset)
{
  _alg.CreateFlashHypothesesFromSegment(pt1,
					pt2,
					dEdx,
					geom,pvs,larp,opdigip,XOffset,
					_prompt_hypothesis,
					_late_hypothesis);
  _total_hypothesis = _prompt_hypothesis+_late_hypothesis;
}

void opdet::FlashHypothesisCollection::Normalize(float const& totalPE_target, util::LArProperties const& larp)
{
  if( GetTotalPEs() < std::numeric_limits<float>::epsilon() ) return;

  const float promptPE_target = totalPE_target * larp.ScintYieldRatio();
  const float latePE_target = totalPE_target * (1.-larp.ScintYieldRatio());

  _prompt_hypothesis.Normalize(promptPE_target);
  _late_hypothesis.Normalize(latePE_target);  
  _total_hypothesis = _prompt_hypothesis + _late_hypothesis;

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

void opdet::FlashHypothesisAlg::AddFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
							      float dEdx,
							      geo::Geometry const& geom,
							      phot::PhotonVisibilityService const& pvs,
							      util::LArProperties const& larp,
							      opdet::OpDigiProperties const& opdigip,
							      float XOffset,
							      FlashHypothesisCollection &hyp_collection)
{
  FlashHypothesisCollection fhc =
    CreateFlashHypothesesFromSegment(pt1,pt2,dEdx,geom,pvs,larp,opdigip,XOffset);
  hyp_collection = hyp_collection + fhc;
}
