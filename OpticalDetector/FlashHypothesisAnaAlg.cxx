/*!
 * Title:   FlashHypothesisAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Alg that compares the flash hypotheses with truth photons and stores the 
 * results in a TTree.
 * 
*/

#include "FlashHypothesisAnaAlg.h"

void opdet::FlashHypothesisAnaAlg::SetOutputObjects(TTree *tree,
						    TH1F* h_h_p, TH1F* h_s_p, TH1F* h_c_p,
						    TH1F* h_h_l, TH1F* h_s_l, TH1F* h_c_l,
						    TH1F* h_h_t, TH1F* h_s_t, TH1F* h_c_t,
						    geo::Geometry const& geom)
{
  fFHCompare.SetOutputObjects(tree,
			      h_h_p,h_s_p,h_c_p,
			      h_h_l,h_s_l,h_c_l,
			      h_h_t,h_s_t,h_c_t,
			      geom.NOpDet());
}

void opdet::FlashHypothesisAnaAlg::FillOpDetPositions(geo::Geometry const& geom)
{

  fOpDetPositions_Y.resize(geom.NOpDet());
  fOpDetPositions_Z.resize(geom.NOpDet());

  double xyz[3];
  for(size_t i_opdet=0; i_opdet<geom.NOpDet(); i_opdet++){
    geom.Cryostat(0).OpDet(i_opdet).GetCenter(xyz);
    fOpDetPositions_Y[i_opdet] = (float)xyz[1];
    fOpDetPositions_Z[i_opdet] = (float)xyz[2];
  }

}

void opdet::FlashHypothesisAnaAlg::RunComparison(const unsigned int run,
						 const unsigned int event,
						 std::vector<sim::MCTrack> const& mctrackVec,
						 std::vector<sim::SimPhotons> const& simPhotonsVec,
						 geo::Geometry const& geom,
						 opdet::OpDigiProperties const& opdigip,
						 phot::PhotonVisibilityService const& pvs,
						 util::LArProperties const& larp)
{
  FlashHypothesisCollection fhc(geom.NOpDet());
  for(auto const& mctrack : mctrackVec){
    if(mctrack.size()==0) continue;
    std::vector<float> dEdxVector(mctrack.size()-1,fdEdx);
    fhc = fhc + fFHCreator.GetFlashHypothesisCollection(mctrack, 
							dEdxVector,
							geom,
							pvs,
							larp,
							opdigip,
							fXOffset);
  }
  
  fSPCAlg.InitializeCounters(geom,opdigip);
  fSPCAlg.AddSimPhotonsVector(simPhotonsVec);

  fFHCompare.RunComparison(run,event,
			   fhc,fSPCAlg.GetSimPhotonCounter(fCounterIndex),
			   fOpDetPositions_Y,fOpDetPositions_Z);
}
