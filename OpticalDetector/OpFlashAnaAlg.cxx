#include "OpticalDetector/OpFlashAnaAlg.h"

void opdet::OpFlashAnaAlg::SetOpFlashTree(TTree* tree)
{

  fOpFlashDataPtr.reset(new recob::OpFlash);
  
  recob::OpFlash* flashptr = fOpFlashDataPtr.get();

  fOpFlashTree->Branch(fOpFlashTree->GetName(),
		       "recob::OpFlash",
		       &flashptr);

}

void opdet::OpFlashAnaAlg::SetFlashTimeHist(TH1F* hist)
{
  fFlashTimeHist = hist;
}

void opdet::OpFlashAnaAlg::FillOpFlashes(const std::vector<recob::OpFlash>& flashVector)
{

  if(!fOpFlashTree) return;
  
  for(auto const& flash : flashVector){
    *fOpFlashDataPtr = flash;
    fOpFlashTree->Fill();
  }

}
