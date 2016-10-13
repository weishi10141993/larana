// This analyzer writes out a TTree containing the properties of
// each reconstructed flash
//
#ifndef FlashPurityCheckAna_H
#define FlashPurityCheckAna_H 1

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
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larana/OpticalDetector/OpDigiProperties.h"

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
#include "larreco/Deprecated/BezierTrack.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"


#include "nusimdata/SimulationBase/MCTruth.h"

namespace opdet {
 
  class FlashPurityCheckAna : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    FlashPurityCheckAna(const fhicl::ParameterSet&);
    virtual ~FlashPurityCheckAna();

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
    std::string fTrackModuleLabel;                // Input tag for OpFlash collection
    std::string fMatchModuleLabel;                // Input tag for OpFlash collection
    std::string fGenieGenModuleLabel;              // Input tag for OpFlash collection

    
    TTree * fPerEventTree;

    Float_t fEventID;

    Float_t VertexX, VertexY, VertexZ;
    Int_t fNVtxTracks;
    Int_t fNVtxTracksRejected;
    Int_t fNVtxTracks20cm;
    Int_t fNVtxTracksRejected20cm;
    Int_t fNNonVtxTracks;
    Int_t fNNonVtxTracksRejected;
    Int_t fNNonVtxTracks20cm;
    Int_t fNNonVtxTracksRejected20cm;

    
  };

} 

#endif // FlashPurityCheckAna_H

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  FlashPurityCheckAna::FlashPurityCheckAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {


    // Indicate that the Input Module comes from .fcl
    fOpFlashModuleLabel  = pset.get<std::string>("OpFlashModuleLabel");
    fTrackModuleLabel    = pset.get<std::string>("TrackModuleLabel");
    fMatchModuleLabel    = pset.get<std::string>("MatchModuleLabel");
    fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");


    art::ServiceHandle<art::TFileService> tfs;

    
    
    fPerEventTree = tfs->make<TTree>("PerEventTree","PerEventTree");
    fPerEventTree->Branch("NVtxTracks",             &fNVtxTracks,             "NVtxTracks/I");
    fPerEventTree->Branch("NVtxTracksRejected",     &fNVtxTracksRejected,     "NVtxTracksRejected/I");
    fPerEventTree->Branch("NVtxTracks20cm",         &fNVtxTracks20cm,         "NVtxTracks20cm/I");
    fPerEventTree->Branch("NVtxTracksRejected20cm", &fNVtxTracksRejected20cm, "NVtxTracksRejected20cm/I");
    fPerEventTree->Branch("NNonVtxTracks",             &fNNonVtxTracks,             "NNonVtxTracks/I");
    fPerEventTree->Branch("NNonVtxTracksRejected",     &fNNonVtxTracksRejected,     "NNonVtxTracksRejected/I");
    fPerEventTree->Branch("NNonVtxTracks20cm",         &fNNonVtxTracks20cm,         "NNonVtxTracks20cm/I");
    fPerEventTree->Branch("NNonVtxTracksRejected20cm", &fNNonVtxTracksRejected20cm, "NNonVtxTracksRejected20cm/I");


    
  }

  //-----------------------------------------------------------------------
  // Destructor
  FlashPurityCheckAna::~FlashPurityCheckAna() 
  {}
   
  //-----------------------------------------------------------------------
  void FlashPurityCheckAna::beginJob()
  {
  }
   

  //-----------------------------------------------------------------------
  void FlashPurityCheckAna::analyze(const art::Event& evt) 
  {
    
    // Get flashes from event
    art::Handle< std::vector< recob::OpFlash > > FlashHandle;
    evt.getByLabel(fOpFlashModuleLabel, FlashHandle);

    // Get matches
    art::Handle< std::vector< anab::FlashMatch > > MatchHandle;
    evt.getByLabel(fMatchModuleLabel, MatchHandle);


    // Get tracks

    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);

    std::vector<art::Ptr<recob::Track> >  Tracks;
    for(unsigned int i=0; i < trackh->size(); ++i)
      {
	art::Ptr<recob::Track> track(trackh,i);
	Tracks.push_back(track);
      }
    std::vector<trkf::BezierTrack*> BTracks;
    BTracks.clear();
    for(size_t i=0; i!=Tracks.size(); i++)
      BTracks.push_back(new trkf::BezierTrack(*Tracks.at(i)));
    std::cout<<"N Tracks : " << BTracks.size()<<std::endl;
    



    // Get MC Truth
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);
    
    if(mclist.size()==0)
      std::cout<<"confused! MC list is zero length!"<<std::endl;
    else
      {
	art::Ptr<simb::MCTruth> mctruth = mclist[0];
      
    
	//    art::FindManyP<recob::OpFlash> FlashesFMH(MatchHandle, evt, fMatchModuleLabel);
	art::FindManyP<anab::FlashMatch>   MatchFMH(trackh,  evt, fMatchModuleLabel);
	
	std::cout<<"No of FMH entries : " << MatchFMH.size()<<std::endl;
	std::vector<bool> Rejected(Tracks.size(), false);
	
	for(size_t i=0; i!=Tracks.size(); ++i)
	  {
	    std::cout<<"FMH at " << i << " is " << MatchFMH.at(i).size()<<std::endl; 
	    for(size_t j=0; j!=MatchFMH.at(i).size(); ++j)
	      {
		if(!MatchFMH.at(i).at(j)->InBeam())
		  Rejected[i]=true;
	      }
	  }
	
	
	
	fEventID=evt.id().event();
	
	
	fNVtxTracks20cm =0 ;
	fNVtxTracksRejected20cm =0;
	fNVtxTracks =0 ;
	fNVtxTracksRejected =0;
	fNNonVtxTracks =0;
	fNNonVtxTracksRejected =0;
	fNNonVtxTracks20cm =0;
	fNNonVtxTracksRejected20cm =0;
	
	
	//	if (mctruth->NeutrinoSet())
	//	  {
	    VertexX = mctruth->GetNeutrino().Nu().Vx();
	    VertexY = mctruth->GetNeutrino().Nu().Vy();
	    VertexZ = mctruth->GetNeutrino().Nu().Vz();
	    
	    TVector3 Vertex = TVector3(VertexX,VertexY, VertexZ);
	    //std::cout<<"Vertex is at " << vtxx_truth<<", " << vtxy_truth<<", " << vtxz_truth<<std::endl;
	    double s=0, d=0;
	    for(size_t i=0; i!=BTracks.size(); ++i)
	      {
		BTracks.at(i)->GetClosestApproach(Vertex, s, d);
		if(d<2) 
		  {
		    fNVtxTracks++;
		    if(Rejected[i]) fNVtxTracksRejected++;
		    if(BTracks.at(i)->Length()>20)
		      {
			fNVtxTracks20cm++;
			if(Rejected[i]) fNVtxTracksRejected20cm++;
		      }
		  }
		else
		  {
		    fNNonVtxTracks++;
		    if(Rejected[i]) fNNonVtxTracksRejected++;
		    if(BTracks.at(i)->Length()>20)
		      {
			fNNonVtxTracks20cm++;
			if(Rejected[i]) fNNonVtxTracksRejected20cm++;
		      }
		  }
	      }
	    
      }
    fPerEventTree->Fill();
	//  }
  }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(FlashPurityCheckAna)
}

