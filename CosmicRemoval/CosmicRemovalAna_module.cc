/////////////////////////////////////////////////////////////
// Cosmic Removal Module Ana
//
// Module Designed to loop over tracks / clusters / hits that
// have the cosmic tag association and removeor ignore the 
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
#include "art/Framework/Core/FindOne.h"
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
#include <sstream>


#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TStopwatch.h"
#include <iterator>


namespace microboone {
   
  class CosmicRemovalAna : public art::EDAnalyzer {

  public:
          
    explicit CosmicRemovalAna(fhicl::ParameterSet const& pset); 
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

    enum hit_origin_t {
      hit_origin_Unknown = -1,
      hit_origin_Cosmic = 0,
      hit_origin_NonCosmic = 1
    };
    
    private:

    unsigned int nCosmicTags;
    
    TTree *tTagTree;
    TTree *tEventTree;

    std::string fHitsModuleLabel;
    std::string fClusterModuleLabel; 
    std::string fTrackModuleLabel;
    std::vector <std::string> fCosmicTagAssocLabel;
    std::vector <float> fCosmicScoreThresholds;

    void InitEventTree(int run_number, int event_number);
    
    void FillMCInfo( std::vector< art::Ptr<recob::Hit> > const& hitlist,
		     std::vector<hit_origin_t> & hitOrigins);
    
    void FillTrackInfo( size_t const& hit_iter,
			hit_origin_t const& origin,
			float const& charge,
			art::FindManyP<recob::Track> const& tracks_per_hit,
			art::Event const& evt,
			std::vector<bool> & hitsAccounted_per_tag);
    
    void FillClusterInfo( size_t const& hit_iter,
			  hit_origin_t const& origin,
			  float const& charge,
			  art::FindManyP<recob::Cluster> const& clusters_per_hit,
			  art::Event const& evt,
			  std::vector<bool> & hitsAccounted_per_tag);
    
    std::vector<float> cTaggedCharge_Cosmic;
    std::vector<float> cTaggedCharge_NonCosmic;
    std::vector<int>   cTaggedHits_Cosmic;
    std::vector<int>   cTaggedHits_NonCosmic;
    std::vector<std::string> *cTagAlgorithmNames;

    typedef struct {
      int eventNumber;
      int tagType;
      float x0;
      float x1;
      float y0;
      float y1;
      float z0;
      float z1;
      int nHits;
      int nGoodHits;
      int origin;
      float score;
      int pdg;
      float energy;
    } cTagProperties_t;
    cTagProperties_t cTagVals;
    
  };//<---End class
}

    typedef struct {
      int runNumber;
      int eventNumber;
      
      int nHitsTotal_Unknown;
      int nHitsTotal_Cosmic;
      int nHitsTotal_NonCosmic;
      
      float qTotal_Unknown;
      float qTotal_Cosmic;
      float qTotal_NonCosmic;
      
      int nHitsTrack;
      int nHitsTrack_Cosmic;
      int nHitsTrack_NonCosmic;
      
      float qTrack;
      float qTrack_Cosmic;
      float qTrack_NonCosmic;
      
      int nHitsCluster;
      int nHitsCluster_Cosmic;
      int nHitsCluster_NonCosmic;
      
      float qCluster;
      float qCluster_Cosmic;
      float qCluster_NonCosmic;
      
    } cEventProperties_t;
    cEventProperties_t cEventVals;


// =====================================================
// fhicl::ParameterSet
// =====================================================
microboone::CosmicRemovalAna::CosmicRemovalAna(fhicl::ParameterSet const& pset):
  EDAnalyzer(pset),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")	      ),
  fCosmicTagAssocLabel      (pset.get<std::vector< std::string > >("CosmicTagAssocLabel") ),
  fCosmicScoreThresholds    (pset.get<std::vector<float> > ("CosmicScoreThresholds") )
{
}


// =====================================================
// BeginJob
// =====================================================

void microboone::CosmicRemovalAna::beginJob()
{
 
  //static cEventProperties_t cEventVals;
 
  nCosmicTags = fCosmicTagAssocLabel.size();
  cTaggedCharge_Cosmic.resize(nCosmicTags);
  cTaggedCharge_NonCosmic.resize(nCosmicTags);
  cTaggedHits_Cosmic.resize(nCosmicTags);
  cTaggedHits_NonCosmic.resize(nCosmicTags);
  
  art::ServiceHandle<art::TFileService> tfs;
  tEventTree = (TTree*)tfs->make<TTree>("CosmicEventTree","CosmicEventTree");
  
  //for(size_t i=0; i!=nCosmicTags; ++i)
  //cTagAlgorithmNames->push_back(fCosmicTagAssocLabel.at(i));
  
  //tTagTree->Branch("tag", &cTagVals, "eventNumber/I:tagType/I:x0/F:x1/F:y0/F:y1/F:z0/F,z1/F:nHits/I:nGoodHits/I:score/F:origin/I:pdg/I:energy/F");
  
  tEventTree->Branch("event", &cEventVals, "runNumber/I:eventNumber/I:nHitsTotal_Unknown/I:nHitsTotal_Cosmic/I:nHitsTotal_NonCosmic/I:qTotal_Unknown/F:qTotal_Cosmic/F:qTotal_NonCosmic/F:nHitsTrack/I:nHitsTrack_Cosmic/I:nHitsTrack_NonCosmic/I:qTrack/F:qTrack_Cosmic/F:qTrack_NonCosmic/F:nHitsCluster/I:nHitsCluster_Cosmic/I:nHitsCluster_NonCosmic/I:qCluster/F:qCluster_Cosmic/F:qCluster_NonCosmic/F");
  tEventTree->Branch("TaggedCharge_Cosmic",&cTaggedCharge_Cosmic);
  tEventTree->Branch("TaggedCharge_NonCosmic",&cTaggedCharge_NonCosmic);
  tEventTree->Branch("TaggedHits_Cosmic",&cTaggedHits_Cosmic);
  tEventTree->Branch("TaggedHits_NonCosmic",&cTaggedHits_NonCosmic);
  //tEventTree->Branch("TagNames",&cTagAlgorithmNames);

}


// =====================================================
// Event Loop
// =====================================================
void microboone::CosmicRemovalAna::analyze(const art::Event& evt)
{

  std::vector< std::pair<std::string,double> > times;
  TStopwatch sw;

  InitEventTree(evt.run(), evt.event());

  sw.Stop();
  times.push_back( std::make_pair<std::string,double>("initTree",sw.RealTime()) );
  sw.Start();

  // ##################################################################
  // ### Grabbing ALL HITS in the event to monitor the backtracking ###
  // ##################################################################
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);

  std::vector<art::Ptr<recob::Hit> > hitlist;  
  art::fill_ptr_vector(hitlist, hitListHandle);

  sw.Stop();
  times.push_back( std::make_pair<std::string,double>("readHits",sw.RealTime()) );
  sw.Start();

  std::vector<hit_origin_t> hitOrigins; hitOrigins.resize(hitlist.size());
  FillMCInfo(hitlist, hitOrigins);  
 
  sw.Stop();
  times.push_back( std::make_pair<std::string,double>("FillMCInfo",sw.RealTime()) );
  sw.Start();
  
  art::FindManyP<recob::Track> tracks_per_hit(hitlist, evt, fTrackModuleLabel);
  art::FindManyP<recob::Cluster> clusters_per_hit(hitlist, evt, fClusterModuleLabel);

  sw.Stop();
  times.push_back( std::make_pair<std::string,double>("GrabAssociations",sw.RealTime()) );
  sw.Start();

  std::vector<bool> false_vector(fCosmicTagAssocLabel.size(),false);
  std::vector< std::vector<bool> > hitsAccounted(hitlist.size(),false_vector);

  for(size_t hit_iter=0; hit_iter<hitlist.size(); hit_iter++){

    float charge = hitlist.at(hit_iter)->Charge();
    hit_origin_t origin = hitOrigins.at(hit_iter);

    FillTrackInfo(hit_iter,
		  origin,
		  charge,
		  tracks_per_hit,
		  evt,
		  hitsAccounted.at(hit_iter));

    FillClusterInfo(hit_iter,
		    origin,
		    charge,
		    clusters_per_hit,
		    evt,
		    hitsAccounted.at(hit_iter));

  }//end loop over all the hits

  sw.Stop();
  times.push_back( std::make_pair<std::string,double>("GetRecoInfo",sw.RealTime()) );
  sw.Start();
  

  tEventTree->Fill();

  sw.Stop();
  times.push_back( std::make_pair<std::string,double>("FillTree",sw.RealTime()) );
  sw.Start();


  for(auto const& time_check : times)
    std::cout << time_check.first << ":" << time_check.second << std::endl;
}

//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::InitEventTree(int run_number, int event_number){

  cEventVals.runNumber = run_number;//evt.run();
  cEventVals.eventNumber = event_number;//evt.event();
  
  cEventVals.nHitsTotal_Unknown = 0;
  cEventVals.nHitsTotal_Cosmic = 0;
  cEventVals.nHitsTotal_NonCosmic = 0;

  cEventVals.qTotal_Unknown = 0;
  cEventVals.qTotal_Cosmic = 0;
  cEventVals.qTotal_NonCosmic = 0;

  cEventVals.nHitsTrack = 0;
  cEventVals.nHitsTrack_Cosmic = 0;
  cEventVals.nHitsTrack_NonCosmic = 0;
  
  cEventVals.qTrack = 0;
  cEventVals.qTrack_Cosmic = 0;
  cEventVals.qTrack_NonCosmic = 0;

  cEventVals.nHitsCluster = 0;
  cEventVals.nHitsCluster_Cosmic = 0;
  cEventVals.nHitsCluster_NonCosmic = 0;
  
  cEventVals.qCluster = 0;
  cEventVals.qCluster_Cosmic = 0;
  cEventVals.qCluster_NonCosmic = 0;
  
  for(size_t iter=0; iter<fCosmicTagAssocLabel.size(); iter++){
    cTaggedHits_Cosmic.at(iter) = 0;
    cTaggedCharge_Cosmic.at(iter) = 0;
    cTaggedHits_NonCosmic.at(iter) = 0;
    cTaggedCharge_NonCosmic.at(iter) = 0;
  }

}

//-------------------------------------------------------------------------------------------------------------------
// take in a list of hits, and determine the origin for those hits (and fill in the tree info)
void microboone::CosmicRemovalAna::FillMCInfo( std::vector< art::Ptr<recob::Hit> > const& hitlist,
					       std::vector<hit_origin_t> & hitOrigins){

  art::ServiceHandle<cheat::BackTracker> bt;
  
  for(size_t itr=0; itr<hitlist.size(); itr++){
    
    art::Ptr<recob::Hit> const& hitptr = hitlist.at(itr);
    std::vector<cheat::TrackIDE> eveIDs = bt->HitToEveID(hitptr);
    
    if(eveIDs.size()==0){
      hitOrigins.at(itr) = hit_origin_Unknown;
      cEventVals.nHitsTotal_Unknown++;
      cEventVals.qTotal_Unknown += hitptr->Charge();
      continue;
    }
    
    float cosmic_energy=0;
    float non_cosmic_energy=0;
    for(auto const& id : eveIDs){
      
      int origin = (bt->TrackIDToMCTruth(id.trackID))->Origin();  
      if(origin == simb::kBeamNeutrino)
	non_cosmic_energy += id.energy;
      else
	cosmic_energy += id.energy;
      
    }
    
    if(non_cosmic_energy > cosmic_energy){
      hitOrigins.at(itr) = hit_origin_NonCosmic;
      cEventVals.nHitsTotal_NonCosmic++;
      cEventVals.qTotal_NonCosmic += hitptr->Charge();
    }
    else{
      hitOrigins.at(itr) = hit_origin_Cosmic;
      cEventVals.nHitsTotal_Cosmic++;
      cEventVals.qTotal_Cosmic += hitptr->Charge();
    }
  }

}//end FillMCInfo


//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillTrackInfo(size_t const& hit_iter,
						 hit_origin_t const& origin,
						 float const& charge,
						 art::FindManyP<recob::Track> const& tracks_per_hit,
						 art::Event const& evt,
						 std::vector<bool> & hitsAccounted_per_tag){

  std::vector< art::Ptr<recob::Track> > const& tracks_this_hit(tracks_per_hit.at(hit_iter));
  if(tracks_this_hit.size()==0) return;
  
  cEventVals.nHitsTrack++;
  cEventVals.qTrack += charge;
  
  if(origin==hit_origin_Cosmic){
    cEventVals.nHitsTrack_Cosmic++;
    cEventVals.qTrack_Cosmic += charge;
  }
  else if(origin==hit_origin_NonCosmic){
    cEventVals.nHitsTrack_NonCosmic++;
    cEventVals.qTrack_NonCosmic += charge;
  }
  
  for(unsigned int nCT = 0; nCT < fCosmicTagAssocLabel.size(); nCT++){//<---This loops over the vector of cosmicTags in stored in the event
    
    if(hitsAccounted_per_tag.at(nCT)) continue;
    
        try{ 
	  art::FindManyP<anab::CosmicTag> cosmic_tags_per_track(tracks_this_hit,evt,fCosmicTagAssocLabel.at(nCT));
	  //for(auto const& cosmic_tags_this_track : cosmic_tags_per_track){ 
	  for(size_t tag_iter=0; tag_iter<cosmic_tags_per_track.size(); tag_iter++){
	    if(cosmic_tags_per_track.at(tag_iter).size()==0) continue;
	    
	    art::Ptr<anab::CosmicTag> const& currentTag(cosmic_tags_per_track.at(tag_iter).at(0));
	    if( currentTag->CosmicScore() > fCosmicScoreThresholds.at(nCT) ){
	      hitsAccounted_per_tag.at(nCT) = true;
	      if(origin==hit_origin_Cosmic){
		cTaggedHits_Cosmic.at(nCT)++;
		cTaggedCharge_Cosmic.at(nCT) += charge;
	      }
	      else if(origin==hit_origin_NonCosmic){
		cTaggedHits_NonCosmic.at(nCT)++;
		cTaggedCharge_NonCosmic.at(nCT) += charge;
	      }
	    }
	    
	  }
	}
	catch (...){ return; }

    
  }
  
}//end FillTrackInfo


//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillClusterInfo(size_t const& hit_iter,
						   hit_origin_t const& origin,
						   float const& charge,
						   art::FindManyP<recob::Cluster> const& clusters_per_hit,
						   art::Event const& evt,
						   std::vector<bool> & hitsAccounted_per_tag){

  std::vector< art::Ptr<recob::Cluster> > const& clusters_this_hit(clusters_per_hit.at(hit_iter));
  if(clusters_this_hit.size()==0) return;
  
  cEventVals.nHitsCluster++;
  cEventVals.qCluster += charge;
  
  if(origin==hit_origin_Cosmic){
    cEventVals.nHitsCluster_Cosmic++;
    cEventVals.qCluster_Cosmic += charge;
  }
  else if(origin==hit_origin_NonCosmic){
    cEventVals.nHitsCluster_NonCosmic++;
    cEventVals.qCluster_NonCosmic += charge;
  }
  
  for(unsigned int nCT = 0; nCT < fCosmicTagAssocLabel.size(); nCT++){//<---This loops over the vector of cosmicTags in stored in the event
    
    if(hitsAccounted_per_tag.at(nCT)) continue;
    
    try{ 
      art::FindManyP<anab::CosmicTag> cosmic_tags_per_cluster(clusters_this_hit,evt,fCosmicTagAssocLabel.at(nCT));
      for(size_t tag_iter=0; tag_iter<cosmic_tags_per_cluster.size(); tag_iter++){
	if(cosmic_tags_per_cluster.at(tag_iter).size()==0) continue;
	
	art::Ptr<anab::CosmicTag> const& currentTag(cosmic_tags_per_cluster.at(tag_iter).at(0));
	if( currentTag->CosmicScore() > fCosmicScoreThresholds.at(nCT) ){
	  hitsAccounted_per_tag.at(nCT) = true;
	  if(origin==hit_origin_Cosmic){
	    cTaggedHits_Cosmic.at(nCT)++;
	    cTaggedCharge_Cosmic.at(nCT) += charge;
	  }
	  else if(origin==hit_origin_NonCosmic){
	    cTaggedHits_NonCosmic.at(nCT)++;
	    cTaggedCharge_NonCosmic.at(nCT) += charge;
	  }
	}
	
      }
    }
    catch(...) { return; }

    
  }
  
}//end FillClusterInfo


namespace microboone{
  
  DEFINE_ART_MODULE(CosmicRemovalAna)  
}








#endif
