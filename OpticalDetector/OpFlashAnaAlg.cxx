#include "OpticalDetector/OpFlashAnaAlg.h"
#include <iostream>

void opdet::OpFlashAnaAlg::SetOpFlashTree(TTree* tree,bool makeOpDetPEHist)
{
  fOpFlashTree = tree;
  fMakeOpDetPEHist = makeOpDetPEHist;
  
  fOpFlashTree->Branch(fOpFlashTree->GetName(),
		       &fOpFlashAnaStruct,
		       fOpFlashAnaStruct.LeafList.c_str());
  if(fMakeOpDetPEHist){
    fOpFlashAnaStruct.FlashOpDetPEHist->SetName("hOpDetPEs");
    fOpFlashAnaStruct.FlashOpDetPEHist->SetTitle("PEs per Optical Detector; OpDet Number; PE");
    fOpFlashTree->Branch("flash_OpDetPEHist","TH1D",&(fOpFlashAnaStruct.FlashOpDetPEHist));
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

    std::cout << "flash y = " << flash.YCenter() << "\tflash total pe=" << flash.TotalPE() << std::endl;

    fOpFlashAnaStruct.FlashTime = flash.Time();
    fOpFlashAnaStruct.FlashTimeWidth = flash.TimeWidth();
    fOpFlashAnaStruct.FlashAbsTime = flash.AbsTime();
    fOpFlashAnaStruct.FlashFrame = flash.Frame();
    fOpFlashAnaStruct.FlashY = flash.YCenter();
    fOpFlashAnaStruct.FlashYWidth = flash.YWidth();
    fOpFlashAnaStruct.FlashZ = flash.ZCenter();
    fOpFlashAnaStruct.FlashZWidth = flash.ZWidth();
    fOpFlashAnaStruct.FlashInBeamFrame = flash.InBeamFrame();
    fOpFlashAnaStruct.FlashOnBeamTime = flash.OnBeamTime();
    fOpFlashAnaStruct.FlashFastToTotal = flash.FastToTotal();
    fOpFlashAnaStruct.FlashTotalPE = flash.TotalPE();

    if(fMakeOpDetPEHist){
      int n_opdets = flash.WireCenters().size();
      fOpFlashAnaStruct.FlashOpDetPEHist->SetBins(n_opdets,-0.5,(float)n_opdets-0.5);
      for(int i=0; i<n_opdets; i++)
	fOpFlashAnaStruct.FlashOpDetPEHist->SetBinContent(i+1,flash.PE(i));
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
