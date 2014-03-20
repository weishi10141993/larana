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
#include "MCCheater/BackTracker.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoObjects/BezierTrack.h"
#include "AnalysisBase/CosmicTag.h"
#include "AnalysisBase/FlashMatch.h"


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
#include <iterator>


namespace microboone {
   
  class CosmicRemovalAna : public art::EDAnalyzer {

  public:
          
    explicit CosmicRemovalAna(fhicl::ParameterSet const& pset); 
    virtual ~CosmicRemovalAna();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

    
    private:
    unsigned int nCosmicTags;
    
    std::string fGenieGenModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fHitsModuleLabel;
    std::string fClusterModuleLabel; 
    std::string fTrackModuleLabel;
    std::vector <std::string> fCosmicTagAssocLabel;
    std::vector <float> fCosmicScoreThresholds;
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
fCosmicTagAssocLabel      (pset.get<std::vector< std::string > >("CosmicTagAssocLabel") ),
fCosmicScoreThresholds    (pset.get<std::vector<float> > ("CosmicScoreThresholds") )
{
}


/////
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


nCosmicTags = fCosmicTagAssocLabel.size();

}


// =====================================================
// Event Loop
// =====================================================
void microboone::CosmicRemovalAna::analyze(const art::Event& evt)
{


/*
// #################################################
// ### Picking up track information on the event ###
// #################################################
art::Handle< std::vector<recob::Track> > trackh; //<---Track Handle
evt.getByLabel(fTrackModuleLabel, trackh); 

std::vector<art::Ptr<recob::Track> > tracklist;//<---Ptr Vector
art::fill_ptr_vector(tracklist,trackh);//<---Fill the vector

// ###################################################
// ### Picking up cluster information on the event ###
// ###################################################	
art::Handle< std::vector<recob::Cluster> > clusterh; //<---Cluster Handle
evt.getByLabel(fClusterModuleLabel, clusterh); 

std::vector<art::Ptr<recob::Cluster> > clusterlist;//<---Ptr Vector
art::fill_ptr_vector(clusterlist,clusterh);//<---Fill the vector


// ### Determining the number of tracks in the event ###
unsigned int trklist = trackh->size();
// ### Determining the number of clusters in the event ###
//unsigned int clulist = clusterh->size();


// #################################################
// ### Looping over the collection of CosmicTags ###
// #################################################
 for(unsigned int nCT = 0; nCT < nCosmicTags; nCT++){
   
   // ### Getting current cosmic tag associations ###
   art::FindManyP<anab::CosmicTag> cosmicTrackTag( tracklist, evt, fCosmicTagAssocLabel[nCT]);
   
   // ============================================================================================================
   // ============================================== LOOKING AT TRACKS ===========================================
   // ============================================================================================================
   // ### Getting hits associatied with tracks ###	
   art::FindManyP<recob::Hit> TrkHit( tracklist, evt, fHitsModuleLabel);
   
   art::Ptr<anab::CosmicTag> currentTag;
   
   std::vector< art::Ptr<recob::Hit> > associatedTPCHits;
   // ###########################
   // ### Looping over tracks ###
   // ###########################
   for(unsigned int trk = 0; trk < trklist; trk++)
     {
       // ############################################################
       // ### Start assuming the track is not associated w/ cosmic ###
       // ############################################################
       bool TrkMatchToCosmicTag = false;
       
       //if(cosmicTrackTag.at(nCT).size()>1) 
	 //mf::LogInfo("CosmicRemovalAna") << "Warning : more than one cosmic tag per track in module " << fCosmicTagAssocLabel[nCT] << ". Confused, but just taking the first one anyway.";
       
       if(cosmicTrackTag.at(nCT).size()==0) continue;
       
      // float CosmicScore = cosmicTrackTag.at(nCT).at(0).CosmicScore();
      // if(CosmicScore > fCosmicScoreThresholds[nCT]) TrkMatchToCosmicTag = true;
       
       //currentTag = cosmictag.at(trk).at(0);
       currentTag = cosmicTrackTag.at(trk).at(0);
       float Score = currentTag->CosmicScore();
       std::cerr << "Score = " << Score << std::endl;
	
	
	// ############################################
	// ###   If Track is matched to a CosmicTag ###
	// ### for the TPC then find associated hits###
	// ############################################
	if(TrkMatchToCosmicTag)
	  {
	    associatedTPCHits = TrkHit.at(trk);
	    
	    
	    
	  }//<---TrkMatchToCosmicTag
	
	
     }//<---end trk loop
   
	

// ============================================================================================================
// ============================================ LOOKING AT CLUSTERS ===========================================
// ============================================================================================================

// ####################################################
// ### Getting the hits associated with the cluster ###
// ####################################################
art::FindManyP<recob::Hit> clusterHit(clusterlist, evt, fHitsModuleLabel);

 art::FindManyP<anab::CosmicTag> cosmicClusterTag( clusterlist, evt, fCosmicTagAssocLabel[nCT]);

std::vector< art::Ptr<recob::Hit> > associatedClusterHits;

for(unsigned int clu = 0; clu < clulist; clu++)
	{
	// ##############################################################
	// ### Start assuming the cluster is not associated w/ cosmic ###
	// ##############################################################
	bool CluMatchToCosmicTag = false;
	
	//if(cosmicClusterTag.at(nCT).size()>1) 
	  //mf::LogInfo("CosmicRemovalAna") << "Warning : more than one cosmic tag per cluster in module " << fCosmicTagAssocLabel[nCT] << ". Confused, but just taking the first one anyway.";
	
       if(cosmicClusterTag.at(nCT).size()==0) continue;
       
       //float CosmicScore = cosmicClusterTag.at(nCT).at(0).CosmicScore();
       
       currentTag = cosmicClusterTag.at(clu).at(0);
       float CosmicScore = currentTag->CosmicScore();
       if(CosmicScore > fCosmicScoreThresholds[nCT]) CluMatchToCosmicTag = true;	
	
	
	// #############################################
	// ### If Cluster is matched to a Flash then ###
	// ###    find the hits associated with it   ###
	// #############################################
	if(CluMatchToCosmicTag)
		{
		associatedClusterHits = clusterHit.at(clu);
		
		
		
		}
	
	}//<---End clu loop
	
		
}//<---End nCT (number of cosmic tag) loop

*/

// ===========================================================================================================
// ============================================ LOOKING AT MCTRUTH ===========================================
// ===========================================================================================================

// ######################################
// ### Picking up BackTracker Service ###
// ######################################
art::ServiceHandle<cheat::BackTracker> bt;

// ###################################################################
// ### Defining a std::map of trackIDEs and vector of hit indicies ###
// ###################################################################
std::map<int, std::vector<int> > trkIDEsHitIndex;

// ##################################################################
// ### Grabbing ALL HITS in the event to monitor the backtracking ###
// ##################################################################
art::Handle< std::vector<recob::Hit> > hitListHandle;
std::vector<art::Ptr<recob::Hit> > hitlist;

if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
	{art::fill_ptr_vector(hitlist, hitListHandle);}


/*std::set<int> backtracker;

backtracker = bt->GetSetOfTrackIDs(hitlist);

for(std::set<int>::iterator it=backtracker.begin(); it!=backtracker.end(); ++it)
{
	std::cout<<"set element : " << *it<<std::endl;
}*/


// ###########################################
// ### Looping over hits to get TrackIDE's ###
// ###########################################
int nhits = hitlist.size();

int counter = 0;

std::cout<<"nhits = "<<nhits<<std::endl;
std::cout<<std::endl;
int ntimethrough = 0;
for ( auto const& itr : hitlist )
	{
	
	
	std::vector<cheat::TrackIDE> eveIDs    = bt->HitToEveID(itr);
	std::vector<cheat::TrackIDE> trackides = bt->HitToTrackID(itr);
	

	
	//bt->ChannelToTrackID(trkIDE, channel, startTime, endTime);
	
	if(eveIDs.size() == 0)
		{
		//std::cout<<"Something is wrong, trkIDs.size() = "<<trackides.size()<<std::endl;
		continue;
		}
	// ############################
	// ### Loop over eveIDs's ###
	// ############################
	//for (size_t e = 0; e<trackides.size(); e++)
	for (size_t e = 0; e<eveIDs.size(); e++)
		{
		
		//art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(trackides[e].trackID);
		art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(eveIDs[e].trackID);
		
		int origin = mctruth->Origin();
		
		// Origin == 1 Not Cosmic
		// Origin == 2 Cosmic
		// ########################################################
		// ### If the Origin of this track is a Cosmic (i.e. 2) ###
		// ### then fill a map keyed by MCParticleID and the 
		// ### the second entry is a vector of hit indicies 
		// ### associated with that particle
		// ### Origin == 1 Not Cosmic
		// ### Origin == 2 Cosmic
		// ########################################################
		if(origin == 2)
			{
			//std::cout<<std::endl;
			//std::cout<<"### This is a cosmic ###"<<std::endl;
			//std::cout<<"Trackides[e].trackID = "<<trackides[e].trackID<<std::endl;
			//std::cout<<"eveIDs[e].trackID = "<<eveIDs[e].trackID<<std::endl;
			//std::cout<<"counter = "<<counter<<std::endl;
			
			trkIDEsHitIndex[eveIDs[e].trackID].push_back(counter);
			
			//std::cout<<"Current hit = "<<std::distance(hitlist.begin(), itr)<<std::endl;
			
			}//<---
		ntimethrough++;
		}//<---End e loop

	counter++;
	}//<--End nhits loop

std::cout<<"Made it through the trackides loop "<<ntimethrough<<" times"<<std::endl;

}


namespace microboone{

DEFINE_ART_MODULE(CosmicRemovalAna)
}

#endif
