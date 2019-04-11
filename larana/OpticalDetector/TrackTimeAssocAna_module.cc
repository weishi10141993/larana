// This analyzer writes out a TTree containing the properties of
// each reconstructed flash
//

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// C++ includes
#include <cstring>
#include <sstream>
#include "math.h"

// LArSoft includes
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "lardata/Utilities/AssociationUtil.h"


namespace opdet {
 
  class TrackTimeAssocAna : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    TrackTimeAssocAna(const fhicl::ParameterSet&);

    // The analyzer routine, called once per event. 
    void analyze (const art::Event&); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fMatchModuleLabel;

    TTree * fMatchTree;
    
    // Match variables
    Int_t   fEventID;
    Int_t   fMatchID;
    Int_t   fFlashID;
    Int_t   fTrackID;
    Float_t fChi2;

    // Flash variables
    Float_t fFFlashTime; 
    Float_t fFAbsTime;
    Float_t fFTimeWidth;
    bool    fFInBeamFrame;
    int     fFOnBeamTime;
    Float_t fFTotalPE;
    Float_t fFYCenter;
    Float_t fFYWidth;
    Float_t fFZCenter;
    Float_t fFZWidth;
    Float_t fFUCenter;
    Float_t fFUWidth;
    Float_t fFVCenter;
    Float_t fFVWidth;
    Float_t fFWCenter;
    Float_t fFWWidth;
    Float_t fFFastToTotal;
    
    // Track variables
    Float_t fTEnd1X;
    Float_t fTEnd2X;
    Float_t fTEnd1Y;
    Float_t fTEnd2Y;
    Float_t fTEnd1Z;
    Float_t fTEnd2Z;
    Float_t fTLength;
 

  };

} 

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  TrackTimeAssocAna::TrackTimeAssocAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {


    // Indicate that the Input Module comes from .fcl
    fMatchModuleLabel = pset.get<std::string>("MatchModuleLabel");

    art::ServiceHandle<art::TFileService const> tfs;

    fMatchTree = tfs->make<TTree>("MatchTree","MatchTree");
    
    fMatchTree->Branch("EventID", &fEventID, "EventID/I");
    fMatchTree->Branch("FlashID", &fFlashID, "FlashID/I");
    fMatchTree->Branch("TrackID", &fTrackID, "TrackID/I");
    fMatchTree->Branch("Chi2",    &fChi2,    "Chi2/F");
    

    // Flash variables
    fMatchTree->Branch("FFlashTime",    &fFFlashTime,    "FFlashTime/F");
    fMatchTree->Branch("FAbsTime",      &fFAbsTime,      "FAbsTime/F");
    fMatchTree->Branch("FTimeWidth",    &fFTimeWidth,    "FTimeWidth/F");
    fMatchTree->Branch("FInBeamFrame",  &fFInBeamFrame,  "FInBeamFrame/B");
    fMatchTree->Branch("FOnBeamTime",   &fFOnBeamTime,   "FOnBeamTime/I");
    fMatchTree->Branch("FTotalPE",      &fFTotalPE,      "FTotalPE/F");
    fMatchTree->Branch("FYCenter",      &fFYCenter,      "FYCenter/F");
    fMatchTree->Branch("FYWidth",       &fFYWidth,       "FYWidth/F");
    fMatchTree->Branch("FZCenter",      &fFZCenter,      "FZCenter/F");
    fMatchTree->Branch("FZWidth",       &fFZWidth,       "FZWidth/F");
    fMatchTree->Branch("FUCenter",      &fFUCenter,      "FUCenter/F");
    fMatchTree->Branch("FUWidth",       &fFUWidth,       "FUWidth/F");
    fMatchTree->Branch("FVCenter",      &fFVCenter,      "FVCenter/F");
    fMatchTree->Branch("FVWidth",       &fFVWidth,       "FVWidth/F");
    fMatchTree->Branch("FWCenter",      &fFWCenter,      "FWCenter/F");
    fMatchTree->Branch("FWWidth",       &fFWWidth,       "FWWidth/F"); 
    fMatchTree->Branch("FFastToTotal",  &fFFastToTotal,  "FFastToTotal/F"); 

    fMatchTree->Branch("TEnd1X",        &fTEnd1X,        "TEnd1X/F"); 
    fMatchTree->Branch("TEnd1Y",        &fTEnd1Y,        "TEnd1Y/F"); 
    fMatchTree->Branch("TEnd1Z",        &fTEnd1Z,        "TEnd1Z/F"); 
    fMatchTree->Branch("TEnd2X",        &fTEnd2X,        "TEnd2X/F"); 
    fMatchTree->Branch("TEnd2Y",        &fTEnd2Y,        "TEnd2Y/F"); 
    fMatchTree->Branch("TEnd2Z",        &fTEnd2Z,        "TEnd2Z/F"); 
    fMatchTree->Branch("TLength",       &fTLength,       "TLength/F"); 

  }

  //-----------------------------------------------------------------------
  void TrackTimeAssocAna::analyze(const art::Event& evt) 
  {
    
    // Get flashes from event
    art::Handle< std::vector< anab::FlashMatch > > MatchHandle;
    evt.getByLabel(fMatchModuleLabel, MatchHandle);
    
    art::PtrVector<anab::FlashMatch> MatchVec;
    for(unsigned int i = 0; i < MatchHandle->size(); ++i){
      art::Ptr<anab::FlashMatch> match(MatchHandle, i);
      MatchVec.push_back(match);
    }

    art::FindManyP<recob::OpFlash> FlashesFMH(MatchHandle, evt, fMatchModuleLabel);
    
    art::FindManyP<recob::Track>   TracksFMH(MatchHandle,  evt, fMatchModuleLabel);



    for(size_t iMatch = 0; iMatch!=MatchVec.size(); ++iMatch)
      {
	fEventID  = evt.id().event();
	fMatchID  = iMatch;
	fTrackID  = MatchVec.at(iMatch)->SubjectID();
	fFlashID  = MatchVec.at(iMatch)->FlashID();
	fChi2     = MatchVec.at(iMatch)->Chi2();
	
	std::vector<art::Ptr<recob::Track> >   Tracks  = TracksFMH.at(iMatch);
	std::vector<art::Ptr<recob::OpFlash> > Flashes = FlashesFMH.at(iMatch);

	if(Tracks.size()!=1)
	  {
	    mf::LogError("TrackTimeAssocAna")<<"ERROR! Match to " << Tracks.size()<<" tracks - should be one!";
	    continue;
	  }

	if(Flashes.size()!=1)
	  {
	    mf::LogError("TrackTimeAssocAna")<<"ERROR! Match to " << Flashes.size()<<" flashes - should be one!";
	    continue;
	  }
	
	// Fill flash variables
	fFFlashTime   = Flashes.at(0)->Time(); 
        fFAbsTime     = Flashes.at(0)->AbsTime();
	fFTimeWidth   = Flashes.at(0)->TimeWidth();
	fFInBeamFrame = Flashes.at(0)->InBeamFrame();
	fFOnBeamTime  = Flashes.at(0)->OnBeamTime();
        fFTotalPE     = Flashes.at(0)->TotalPE();
        fFYCenter     = Flashes.at(0)->YCenter();
	fFYWidth      = Flashes.at(0)->YWidth();
	fFZCenter     = Flashes.at(0)->ZCenter();
	fFZWidth      = Flashes.at(0)->ZWidth();
	fFUCenter     = Flashes.at(0)->WireCenters()[0];
        fFUWidth      = Flashes.at(0)->WireWidths()[0];
	fFVCenter     = Flashes.at(0)->WireCenters()[1];
	fFVWidth      = Flashes.at(0)->WireWidths()[1];
	fFWCenter     = Flashes.at(0)->WireCenters()[2];
	fFWWidth      = Flashes.at(0)->WireWidths()[2];
        fFFastToTotal = Flashes.at(0)->FastToTotal();
	
	
	// Fill track variables
	fTEnd1X = Tracks.at(0)->Vertex().X();
	fTEnd1Y = Tracks.at(0)->Vertex().Y();
	fTEnd1Z = Tracks.at(0)->Vertex().Z();
	fTEnd2X = Tracks.at(0)->End().X();
        fTEnd2Y = Tracks.at(0)->End().Y();
	fTEnd2Z = Tracks.at(0)->End().Z();

        fTLength = Tracks.at(0)->Length();	

	fMatchTree->Fill();
      } 

    
  }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(TrackTimeAssocAna)
}
