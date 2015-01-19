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

void opdet::FlashHypothesisAnaAlg::SetOutputTree(TTree* tree,
						 TH1F* hh_p, TH1F* hh_pe,
						 TH1F* hh_l, TH1F* hh_le,
						 TH1F* hh_t, TH1F* hh_te,
						 TH1F* hs_p, TH1F* hs_l,
						 TH1F* hc_p, TH1F* hc_l,
						 geo::Geometry const& geom)
{
  fTree = tree;
  int n_opdets = geom.NOpDet();
  hhyp_prompt = hh_p; hhyp_prompt_err = hh_pe
  hh_p->SetNameTitle(
}

void opdet::FlashHypothesisAnaAlg::FillAnaTree(std::vector<TVector3> const& trajVector, 
					       std::vector<float> const& dEdxVector,
					       sim::SimPhotonsCollection const& simPhotonCol,
					       geo::Geometry const& geom,
					       opdet::OpDigiProperties const& opdigip,
					       phot::PhotonVisibilityService const& pvs,
					       util::LArProperties const& larp,
					       float XOffset=0)
{
  FlashHypothesisCollection fhc = fFHAlg.GetFlashHypothesisCollection(trajVector, 
								      dEdxVector,
								      geom,
								      pvs,
								      larp,
								      opdigip,
								      XOffset=0);
  FillFlashHypothesis(fhc);

  fSPCAlg.InitializeCounters(geom,opdigip);
  fSPCAlg.AddSimPhotonCollection(simPhotonCol);
  FillSimPhotonCounter(fSPCAlg.PromptPhotonVector(fCounterIndex),fSPCAlg.LatePhotonVector(fCounterIndex));

  std::vector<float> resultVec_PromptCompare(geom.NOpDet());
  std::vector<float> resultVec_LateCompare(geom.NOpDet());
  float result_PromptCompare = fhc.GetPromptHypothesis().CompareByError(fSPCAlg.PromptPhotonVector(fCounterIndex),
									resultVec_PromptCompare);
  float result_PromptCompare = fhc.GetLateHypothesis().CompareByError(fSPCAlg.LatePhotonVector(fCounterIndex),
								      resultVec_LateCompare);
  FillComparisonInfo(resultVec_PromptCompare,result_PromptCompare,
		     resultVec_LateCompare,result_LateCompare);
}

void opdet::FlashHypothesisAnaAlg::FillFlashHypothesis(FlashHypothesisCollection const& fhc)
{
}

void opdet::FlashHypothesisAnaAlg::FillSimPhotonCounter(std::vector<float> const& promptVec,
							std::vector<float> const& lateVec)
{
}

void opdet::FlashHypothesisAnaAlg::FillComparisonInfo(std::vector<float> const& promptVec, float const& promptTotal
						      std::vector<float> const& lateVec, float const& lateTotal)
{
}
