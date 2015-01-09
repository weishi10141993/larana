/*!
 * Title:   FlashHypothesis Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Algorithm that produces a flash hypothesis for a trajectory.
 * Input:       Trajectory (std::vector<TVector3> objects)
 * Output:      FlashHypotheses
*/

#include "FlashHypothesisAlg.h"

opdet::FlashHypothesisCollection::FlashHypothesisCollection(recob::Track const& track, 
							    std::vector<float> dEdxVector,
							    geo::Geometry const& geom,
							    phot::PhotonVisibilityService const& pvs,
							    util::LArProperties const& larp,
							    opdet::OpDigiProperties const& opdigip,
							    float XOffset)
{
  if(track.NumberTrajectoryPoints() != dEdxVector.size())
    throw "ERROR in FlashHypothesisCollection: dEdx vector size not same as track size.";

  Initialize(geom.NOpDet());
  for(size_t pt=1; pt<track.NumberTrajectoryPoints(); pt++)
    _alg.AddFlashHypothesesFromSegment(track.LocationAtPoint(pt-1),
				       track.LocationAtPoint(pt),
				       dEdxVector[pt],
				       geom,pvs,larp,opdigip,XOffset,
				       _prompt_hypothesis,
				       _late_hypothesis);

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

void opdet::FlashHypothesisAlg::CreateFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
								 float dEdx,
								 geo::Geometry const& geom,
								 phot::PhotonVisibilityService const& pvs,
								 util::LArProperties const& larp,
								 opdet::OpDigiProperties const& opdigip,
								 float XOffset,
								 FlashHypothesis &prompt_flash_hyp,
								 FlashHypothesis &late_flash_hyp)
{
  prompt_flash_hyp = FlashHypothesis(geom.NOpDet());
  late_flash_hyp = FlashHypothesis(geom.NOpDet());

  double xyz_segment[3];
  xyz_segment[0] = 0.5*(pt2.x()+pt1.x()) + XOffset;
  xyz_segment[1] = 0.5*(pt2.y()+pt1.y());
  xyz_segment[2] = 0.5*(pt2.z()+pt1.z());
  
  //get the visibility vector
  const std::vector<float>* PointVisibility = pvs.GetAllVisibilities(xyz_segment);
  
  //check vector size, as it may be zero if given a y/z outside some range
  if(PointVisibility->size()!=geom.NOpDet()) return;
  
  //get the amount of light
  const float PromptLightAmount = larp.ScintYield()*larp.ScintYieldRatio()*opdigip.QE()*dEdx*(pt2-pt1).Mag();
  const float LateLightAmount = larp.ScintYield()*(1.-larp.ScintYieldRatio())*opdigip.QE()*dEdx*(pt2-pt1).Mag();
  
  for(size_t opdet_i=0; opdet_i<geom.NOpDet(); opdet_i++){
    float vis = PointVisibility->at(opdet_i);
    prompt_flash_hyp.SetHypothesisAndError(opdet_i,vis*PromptLightAmount);
    late_flash_hyp.SetHypothesisAndError(opdet_i,vis*LateLightAmount);

    //WK, IIT Workshop --- Ignoring saturation limit for now.
    //
    //apply saturation limit
    //if(lightHypothesis[opdet_i]>fOpDetSaturation){
    //totalHypothesisPE -= (lightHypothesis[opdet_i]-fOpDetSaturation);
    //lightHypothesis[opdet_i] = fOpDetSaturation;
    //}
  }

}

void opdet::FlashHypothesisAlg::AddFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
							      float dEdx,
							      geo::Geometry const& geom,
							      phot::PhotonVisibilityService const& pvs,
							      util::LArProperties const& larp,
							      opdet::OpDigiProperties const& opdigip,
							      float XOffset,
							      FlashHypothesis &prompt_flash_hyp,
							      FlashHypothesis &late_flash_hyp)
{
  FlashHypothesis this_segment_prompt(geom.NOpDet());
  FlashHypothesis this_segment_late(geom.NOpDet());
  
  CreateFlashHypothesesFromSegment(pt1,pt2,dEdx,
				   geom,pvs,larp,opdigip,XOffset,
				   this_segment_prompt,this_segment_late);

  prompt_flash_hyp = prompt_flash_hyp + this_segment_prompt;
  late_flash_hyp = late_flash_hyp + this_segment_late;
  
}
