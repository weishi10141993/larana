// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree containing the properties of
// each reconstructed flash
//

#ifndef OpFlashAna_H
#define OpFlashAna_H 1

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RawData/OpDetPulse.h"
#include "RecoBase/OpFlash.h"
#include "RecoBase/OpHit.h"
#include "DetectorInfoServices/DetectorClocksService.h"
//#include "OpticalDetector/OpDigiProperties.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"

namespace opdet {
 
  class OpFlashAna : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    OpFlashAna(const fhicl::ParameterSet&);
    virtual ~OpFlashAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The analyzer routine, called once per event. 
    void analyze (const art::Event&); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fOpFlashModuleLabel;              // Input tag for OpFlash collection
    std::string fOpHitModuleLabel;              // Input tag for OpHit collection
    float fSampleFreq;                     // in MHz
    float fTimeBegin;                      // in us
    float fTimeEnd;                        // in us
    
    float fYMin, fYMax, fZMin, fZMax;

    int PosHistYRes, PosHistZRes;

    bool fMakeFlashTimeHist;
    bool fMakeFlashPosHist;
    bool fMakePerFlashHists;

    bool fMakePerFlashTree;
    bool fMakeFlashBreakdownTree;
    bool fMakePerOpHitTree;
    bool fMakeFlashHitMatchTree;
      
    TTree * fPerFlashTree;
    TTree * fPerOpHitTree;
    TTree * fFlashBreakdownTree;
    TTree * fFlashHitMatchTree;
   
    Int_t fEventID;
    Int_t fFlashID;
    Int_t fHitID;
    Float_t fFlashTime; 
    Float_t fAbsTime;
    bool  fInBeamFrame;
    int   fOnBeamTime;
    Float_t fTotalPE;
    Int_t   fFlashFrame;
    
    Float_t fNPe;
    Float_t fYCenter;
    Float_t fYWidth;
    Float_t fZCenter;
    Float_t fZWidth;

    Int_t   fOpChannel;
    Float_t fPeakTimeAbs;
    Float_t fPeakTime;
    Int_t   fFrame;
    Float_t fWidth;
    Float_t fArea;
    Float_t fAmplitude;
    Float_t fPE;
    Float_t fFastToTotal;
  };

} 

#endif // OpFlashAna_H

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpFlashAna::OpFlashAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {


    // Indicate that the Input Module comes from .fcl
    fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel   = pset.get<std::string>("OpHitModuleLabel");


//    art::ServiceHandle<OpDigiProperties> odp;
//    fTimeBegin  = odp->TimeBegin();
//    fTimeEnd    = odp->TimeEnd();
//    fSampleFreq = odp->SampleFreq();

    const detinfo::DetectorClocks* clks = lar::providerFrom<detinfo::DetectorClocksService>();
    fTimeBegin  = clks->OpticalClock().Time();
    fTimeEnd    = clks->OpticalClock().FramePeriod();
    fSampleFreq = clks->OpticalClock().Frequency();
    
   
    fYMin =  pset.get<float>("YMin");
    fYMax =  pset.get<float>("YMax");
    fZMin =  pset.get<float>("ZMin");
    fZMax =  pset.get<float>("ZMax");
    
    fMakeFlashTimeHist = pset.get<bool>("MakeFlashTimeHist");
    fMakeFlashPosHist  = pset.get<bool>("MakeFlashPosHist");
    fMakePerFlashHists = pset.get<bool>("MakePerFlashHists");

    fMakePerFlashTree       =  pset.get<bool>("MakePerFlashTree");
    fMakeFlashBreakdownTree =  pset.get<bool>("MakeFlashBreakdownTree");
    fMakePerOpHitTree       =  pset.get<bool>("MakePerOpHitTree");
    fMakeFlashHitMatchTree  =  pset.get<bool>("MakeFlashHitMatchTree");

    
    PosHistYRes=100;
    PosHistZRes=100;

    art::ServiceHandle<art::TFileService> tfs;

    if(fMakeFlashBreakdownTree)
      {
	fFlashBreakdownTree = tfs->make<TTree>("FlashBreakdownTree","FlashBreakdownTree");
	fFlashBreakdownTree->Branch("EventID",   &fEventID,    "EventID/I");
	fFlashBreakdownTree->Branch("FlashID",   &fFlashID,    "FlashID/I");
	fFlashBreakdownTree->Branch("OpChannel", &fOpChannel,  "OpChannel/I");
	fFlashBreakdownTree->Branch("FlashTime", &fFlashTime,  "FlashTime/F");
	fFlashBreakdownTree->Branch("NPe",       &fNPe,        "NPe/F");
	fFlashBreakdownTree->Branch("AbsTime",   &fAbsTime,    "AbsTime/F");
	fFlashBreakdownTree->Branch("FlashFrame",&fFlashFrame, "FlashFrame/I");
	fFlashBreakdownTree->Branch("InBeamFrame",&fInBeamFrame, "InBeamFrame/B");
	fFlashBreakdownTree->Branch("OnBeamTime",&fOnBeamTime, "OnBeamTime/I");
	fFlashBreakdownTree->Branch("YCenter",   &fYCenter,    "YCenter/F");
	fFlashBreakdownTree->Branch("ZCenter",   &fZCenter,    "ZCenter/F");
	fFlashBreakdownTree->Branch("YWidth",    &fYWidth,     "YWidth/F");
	fFlashBreakdownTree->Branch("ZWidth",    &fZWidth,     "ZWidth/F");
	fFlashBreakdownTree->Branch("TotalPE",   &fTotalPE,    "TotalPE/F");
      }

    if(fMakePerOpHitTree)
      {
	fPerOpHitTree = tfs->make<TTree>("PerOpHitTree","PerOpHitTree");
	fPerOpHitTree->Branch("EventID",     &fEventID,      "EventID/I");
	fPerOpHitTree->Branch("HitID",       &fHitID,        "HitID/I");
	fPerOpHitTree->Branch("OpChannel",   &fOpChannel,    "OpChannel/I");
	fPerOpHitTree->Branch("PeakTimeAbs", &fPeakTimeAbs,  "PeakTimeAbs/F");
	fPerOpHitTree->Branch("PeakTime",    &fPeakTime,     "PeakTime/F");
	fPerOpHitTree->Branch("Frame",       &fFrame,        "Frame/I");
	fPerOpHitTree->Branch("Width",       &fWidth,        "Width/F");
	fPerOpHitTree->Branch("Area",        &fArea,         "Area/F");
	fPerOpHitTree->Branch("Amplitude",   &fAmplitude,    "Amplitude/F");
	fPerOpHitTree->Branch("PE",          &fPE,           "PE/F");
	fPerOpHitTree->Branch("FastToTotal", &fFastToTotal,  "FastToTotal/F");
      }
    if(fMakePerFlashTree)
      {
	fPerFlashTree = tfs->make<TTree>("PerFlashTree","PerFlashTree");
	fPerFlashTree->Branch("EventID",   &fEventID,    "EventID/I");
	fPerFlashTree->Branch("FlashID",   &fFlashID,    "FlashID/I");
	fPerFlashTree->Branch("YCenter",   &fYCenter,    "YCenter/F");
	fPerFlashTree->Branch("ZCenter",   &fZCenter,    "ZCenter/F");
	fPerFlashTree->Branch("YWidth",    &fYWidth,     "YWidth/F");
	fPerFlashTree->Branch("ZWidth",    &fZWidth,     "ZWidth/F");
	fPerFlashTree->Branch("FlashTime", &fFlashTime,  "FlashTime/F");
	fPerFlashTree->Branch("AbsTime",   &fAbsTime,    "AbsTime/F");
	fPerFlashTree->Branch("FlashFrame",&fFlashFrame, "FlashFrame/I");
	fPerFlashTree->Branch("InBeamFrame",&fInBeamFrame, "InBeamFrame/B");
	fPerFlashTree->Branch("OnBeamTime",&fOnBeamTime, "OnBeamTime/I");
	fPerFlashTree->Branch("TotalPE",   &fTotalPE,    "TotalPE/F");
      }
    if(fMakeFlashHitMatchTree)
      {
	fFlashHitMatchTree = tfs->make<TTree>("FlashHitMatchTree","FlashHitMatchTree");
	fFlashHitMatchTree->Branch("EventID",       &fEventID,    "EventID/I");
	fFlashHitMatchTree->Branch("FlashID",       &fFlashID,    "FlashID/I");
	fFlashHitMatchTree->Branch("HitID",         &fHitID,      "HitID/I");
	fFlashHitMatchTree->Branch("OpChannel",     &fOpChannel,  "OpChannel/I");
	fFlashHitMatchTree->Branch("HitPeakTimeAbs",&fPeakTimeAbs,"HitPeakTimeAbs/F");
	fFlashHitMatchTree->Branch("HitPeakTime",   &fPeakTime,   "HitPeakTime/F");
	fFlashHitMatchTree->Branch("HitPE",         &fPE,         "HitPE/F");
	fFlashHitMatchTree->Branch("FlashPE",       &fTotalPE,    "FlashPE/F");
	fFlashHitMatchTree->Branch("FlashTimeAbs",  &fAbsTime,    "FlashTimeAbs/F");
	fFlashHitMatchTree->Branch("FlashTime",     &fFlashTime,  "FlashTime/F");
	fFlashHitMatchTree->Branch("HitFrame",      &fFrame,      "HitFrame/I");
	fFlashHitMatchTree->Branch("FlashFrame",    &fFlashFrame, "FlashFrame/I");
      }
    
    fFlashID=0;
  }

  //-----------------------------------------------------------------------
  // Destructor
  OpFlashAna::~OpFlashAna() 
  {}
   
  //-----------------------------------------------------------------------
  void OpFlashAna::beginJob()
  {
  }
   

  //-----------------------------------------------------------------------
  void OpFlashAna::analyze(const art::Event& evt) 
  {
    
    // Get flashes from event
    art::Handle< std::vector< recob::OpFlash > > FlashHandle;
    evt.getByLabel(fOpFlashModuleLabel, FlashHandle);

    // Get assosciations between flashes and hits
    art::FindManyP<recob::OpHit> Assns(FlashHandle, evt, fOpFlashModuleLabel);

    // Create string for histogram name
    char HistName[50];
    
    fFlashID=0;


    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
    art::ServiceHandle<art::TFileService> tfs;
    
    std::vector<TH1D*> FlashHist;
    
    fEventID=evt.id().event();

    sprintf(HistName, "Event %d Flash Times", evt.id().event());
    TH1D * FlashTimes = nullptr;
    if(fMakeFlashTimeHist)
      {
	FlashTimes = tfs->make<TH1D>(HistName, ";t (ns);", 
					    int((fTimeEnd - fTimeBegin) * fSampleFreq), 
					    fTimeBegin * 1000., 
					    fTimeEnd * 1000.);
      }

    TH2D * FlashPositions = nullptr;
    if(fMakeFlashPosHist)
      {
	sprintf(HistName, "Event %d All Flashes YZ", evt.id().event());
	
	FlashPositions = tfs->make<TH2D>(HistName, ";y ;z ", 
						PosHistYRes, fYMin, fYMax,
						PosHistZRes, fZMin, fZMax);
      }
    
    art::ServiceHandle<geo::Geometry> geom;
    unsigned int NOpChannels = geom->NOpChannels();


    // For every OpFlash in the vector
    for(unsigned int i = 0; i < FlashHandle->size(); ++i)
      {

	// Get OpFlash
	art::Ptr< recob::OpFlash > TheFlashPtr(FlashHandle, i);
	recob::OpFlash TheFlash = *TheFlashPtr;


	fFlashTime = TheFlash.Time();
	fFlashID = i; //++;
	
	TH2D * ThisFlashPosition = nullptr;
	if(fMakePerFlashHists)
	  {
	    sprintf(HistName, "Event %d t = %f", evt.id().event(), fFlashTime);
	    FlashHist.push_back ( tfs->make<TH1D>(HistName, ";OpChannel;PE", 
						  NOpChannels, 0, NOpChannels));
	
	    sprintf(HistName, "Event %d Flash %f YZ", evt.id().event(), fFlashTime);
	    
	    ThisFlashPosition = tfs->make<TH2D>(HistName, ";y ;z ", 
						PosHistYRes, fYMin, fYMax,
						PosHistZRes, fZMin, fZMax);
	  }
	fYCenter     = TheFlash.YCenter();
	fZCenter     = TheFlash.ZCenter();
	fYWidth      = TheFlash.YWidth();
	fZWidth      = TheFlash.ZWidth();
	fInBeamFrame = TheFlash.InBeamFrame();
	fOnBeamTime  = TheFlash.OnBeamTime();
	fAbsTime     = TheFlash.AbsTime();
	fFlashFrame  = TheFlash.Frame();
	fTotalPE     = TheFlash.TotalPE();
	
	
	for(unsigned int j=0; j!=NOpChannels; ++j)
	  {
	    if(fMakePerFlashHists) FlashHist.at(FlashHist.size()-1)->Fill(j, TheFlash.PE(j));
	    fNPe = TheFlash.PE(j);
	    fOpChannel=j;
	    
	    if((fMakeFlashBreakdownTree)&&(fNPe>0)) fFlashBreakdownTree->Fill();
	    
	  }

	for(int y=0; y!=PosHistYRes; ++y)
	  for(int z=0; z!=PosHistZRes; ++z)
	    {
	      float ThisY = fYMin + (fYMax-fYMin)*float(y)/PosHistYRes + 0.0001;
	      float ThisZ = fZMin + (fZMax-fZMin)*float(z)/PosHistZRes + 0.0001;
	      if (fMakePerFlashHists) ThisFlashPosition->Fill(ThisY, ThisZ, fTotalPE * exp(-pow((ThisY-fYCenter)/fYWidth,2)/2.-pow((ThisZ-fZCenter)/fZWidth,2)/2.));
	      if (fMakeFlashPosHist) FlashPositions->Fill(ThisY, ThisZ, fTotalPE * exp(-pow((ThisY-fYCenter)/fYWidth,2)-pow((ThisZ-fZCenter)/fZWidth,2)));
	      
	    }
	
	
	

	if(fMakeFlashTimeHist) FlashTimes->Fill(fFlashTime, fTotalPE);
	
	if(fMakePerFlashTree)  fPerFlashTree->Fill();

	if(fMakeFlashHitMatchTree) {
	  // Extract the assosciations
	  fHitID = 0;
	  std::vector< art::Ptr<recob::OpHit> > matchedHits = Assns.at(i);
	  for (auto ophit : matchedHits) {
	    fOpChannel   = ophit->OpChannel();
	    fPeakTimeAbs = ophit->PeakTimeAbs();
	    fPeakTime    = ophit->PeakTime();
	    fFrame       = ophit->Frame();
	    fWidth       = ophit->Width();
	    fArea        = ophit->Area();
	    fAmplitude   = ophit->Amplitude();
	    fPE          = ophit->PE();
	    fFastToTotal = ophit->FastToTotal();
	    fFlashHitMatchTree->Fill();
	    fHitID++;
	  }
	}
      }
    
    if(fMakePerOpHitTree)
      {
	art::Handle< std::vector< recob::OpHit > > OpHitHandle;
	evt.getByLabel(fOpHitModuleLabel, OpHitHandle);
	
	for(size_t i=0; i!=OpHitHandle->size(); ++i)
	  {
	    fOpChannel   = OpHitHandle->at(i).OpChannel();
	    fPeakTimeAbs = OpHitHandle->at(i).PeakTimeAbs();
	    fPeakTime    = OpHitHandle->at(i).PeakTime();
	    fFrame       = OpHitHandle->at(i).Frame();
	    fWidth       = OpHitHandle->at(i).Width();
	    fArea        = OpHitHandle->at(i).Area();
	    fAmplitude   = OpHitHandle->at(i).Amplitude();
	    fPE          = OpHitHandle->at(i).PE();
	    fFastToTotal = OpHitHandle->at(i).FastToTotal();
	    fHitID       = i;
	    fPerOpHitTree->Fill();
	  }

      }

  }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(OpFlashAna)
}

