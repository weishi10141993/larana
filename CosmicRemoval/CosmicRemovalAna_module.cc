/////////////////////////////////////////////////////////////
// Cosmic Removal Module Ana
//
// Module Designed to loop over tracks / clusters / hits that
// have the cosmic tag association and remove or ignore the 
// hits and check to see what is left over
//
// Yale Workshop (Cosmics Removal Group)
//
/////////////////////////////////////////////////////////////

#ifndef COSMICSREMOVALANA_H
#define COSMICSREMOVALANA_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Principal/View.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Core/FindMany.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft Includes
#include "Geometry/Geometry.h"
#include "SimulationBase/MCTruth.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoObjects/BezierTrack.h"
#include "AnalysisBase/CosmicTag.h"

// ROOT Includes
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"


namespace microboone {
   
  class CosmicRemovalAna : public art::EDAnalyzer {

  public:
          
    explicit CosmicRemovalAna(fhicl::ParameterSet const& pset); 
    virtual ~CosmicRemovalAna();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

    
    private:
    std::string fGenieGenModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fHitsModuleLabel;
    std::string fClusterModuleLabel; 
    std::string fTrackModuleLabel;
    std::string fCosmicTagAssocLabel;
    };//<---End     
}

// =====================================================
// fhicl::ParameterSet
// =====================================================
microboone::CosmicRemovalAna::CosmicRemovalAna(fhicl::ParameterSet const& pset):
EDAnalyzer(pset),
fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      ),
fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")	      ),
fCosmicTagAssocLabel      (pset.get< std::string >("CosmicTagAssocLabel")     )

{
}



// =====================================================
// Deconstructor
// =====================================================
microboone::CosmicRemovalAna::~CosmicRemovalAna()
{
}

// =====================================================
// BeginJob
// =====================================================

void microboone::CosmicRemovalAna::beginJob(){
art::ServiceHandle<art::TFileService> tfs;


}


// =====================================================
// Event Loop
// =====================================================
void microboone::CosmicRemovalAna::analyze(const art::Event& evt)
{

// #################################################
// ### Picking up track information on the event ###
// #################################################
art::Handle< std::vector<recob::Track> > trackh; //<---Track Handle
evt.getByLabel(fTrackModuleLabel, trackh); 

std::vector<art::Ptr<recob::Track> > tracklist;//<---Ptr Vector
art::fill_ptr_vector(tracklist,trackh);//<---Fill the vector


unsigned int trklist = trackh->size();

// ####################################################
// ### Getting the collection of cosmic tag objects ###
// ####################################################
art::FindManyP<anab::CosmicTag> cosmictag( tracklist ,evt,fCosmicTagAssocLabel);

//std::cout<<"cosmictag.size() "<<cosmictag.size()<<std::endl;
//std::cout<<"trklist = "<<trklist<<std::endl;


// ###########################
// ### Looping over tracks ###
// ###########################
for(unsigned int trk = 0; trk < trklist; trk++)
	{
	//std::cout<<"Running Over Tracks"<<std::endl;
	//auto it = trackh.begin();
	
	//art::FindManyP<anab::CosmicTag> costag(trk, evt, fCosmicTagAssocLabel);
	  
	  //std::cerr << "what is this? " << cosmictag.at(trk)->cosmicScore << std::endl;
  
  
  

	}//<---end trk loop
	
	
// ###################################################
// ### Picking up cluster information on the event ###
// ###################################################	
art::Handle< std::vector<recob::Cluster> > clusterh; //<---Cluster Handle
evt.getByLabel(fClusterModuleLabel, clusterh); 

std::vector<art::Ptr<recob::Cluster> > clusterlist;//<---Ptr Vector
art::fill_ptr_vector(clusterlist,clusterh);//<---Fill the vector


}


namespace microboone{

DEFINE_ART_MODULE(CosmicRemovalAna)
}

#endif
