#include "RecoBase/OpHit.h"
#include "RecoBase/OpFlash.h"

#include "OpticalDetector/OpFlashAnaAlg.h"

void opdet::OpFlashAnaAlg::SetOpFlashTree(TTree* tree,bool makeOpHitHist)
{
  fOpFlashTree = tree;
  
  fOpFlashDataPtr.reset(new recob::OpFlash);  
  recob::OpFlash* flashptr = fOpFlashDataPtr.get();
  fOpFlashTree->Branch(fOpFlashTree->GetName(),
		       "recob::OpFlash",
		       &flashptr);

  fMakeOpHitHist = makeOpHitHist;
  if(fMakeOpHitHist){
    fOpFlashHitHistPtr.reset(new TH1D);
    TH1D* hithistptr = fOpFlashHitHistPtr.get();
    hithistptr->SetName("hOpHits");
    hithistptr->SetTitle("PEs per Optical Detector; OpDet Number; PE");
    fOpFlashTree->Branch("flash_OpHitHist","TH1D",&hithistptr);
  }

}

void opdet::OpFlashAnaAlg::SetOpHitTree(TTree* tree)
{
  fOpHitTree = tree;
  
  fOpHitDataPtr.reset(new recob::OpHit);  
  recob::OpHit* hitptr = fOpHitDataPtr.get();
  fOpHitTree->Branch(fOpHitTree->GetName(),
		       "recob::OpHit",
		       &hitptr);

}

void opdet::OpFlashAnaAlg::FillOpFlashes(const std::vector<recob::OpFlash>& flashVector)
{
  if(fOpFlashTree) FillOpFlashTree(flashVector);
}

void opdet::OpFlashAnaAlg::FillOpHits(const std::vector<recob::OpHit>& hitVector)
{
  if(fOpHitTree) FillOpHitTree(hitVector);
}

void opdet::OpFlashAnaAlg::FillOpFlashTree(const std::vector<recob::OpFlash>& flashVector)
{
  for(auto const& flash : flashVector){
    *fOpFlashDataPtr = flash;

    if(fMakeOpHitHist){
      TH1D* hithistptr = fOpFlashHitHistPtr.get();
      int n_opdets = flash.WireCenters().size();
      hithistptr->SetBins(n_opdets,-0.5,(float)n_opdets-0.5);
      for(int i=0; i<n_opdets; i++)
	hithistptr->SetBinContent(i+1,flash.PE(i));
    }

    fOpFlashTree->Fill();
  }
}

void opdet::OpFlashAnaAlg::FillOpHitTree(const std::vector<recob::OpHit>& hitVector)
{
  for(auto const& hit : hitVector){
    *fOpHitDataPtr = hit;
    fOpHitTree->Fill();
  }
}
