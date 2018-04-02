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
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"


// ROOT Includes
#include "TTree.h"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>
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
    
  //  TTree *tTagTree;
    TTree *tEventTree;
    
    std::string fHitsModuleLabel;
    std::string fMCModuleLabel;
    std::string fMCHitsModuleLabel;
    std::string fClusterModuleLabel; 
    std::string fTrackModuleLabel;
    float       fHitCompareCut;
    std::vector <std::string> fCosmicTagAssocLabel;
    std::vector <float> fCosmicScoreThresholds;

    void InitEventTree(int run_number, int event_number);
    
    void FillMCInfo( std::vector<recob::Hit> const& hitlist,
		     std::vector<hit_origin_t> & hitOrigins,
		     std::vector<sim::MCHitCollection> const& mchitCollectionVector,
		     std::map<int,const simb::MCTruth* > const& trackIDToTruthMap);

    void FillTrackInfo(size_t const& hit_iter,
		       hit_origin_t const& origin,
		       float const& charge,
		       std::vector<size_t> const& track_indices_this_hit,
		       std::vector< std::vector< const anab::CosmicTag* > > const& tags_per_cluster,
		       std::vector<bool> & hitsAccounted_per_tag,
		       std::vector<bool> & hitsAllTags);
    
    void FillClusterInfo(size_t const& hit_iter,
			 hit_origin_t const& origin,
			 float const& charge,
			 std::vector<size_t> const& cluster_indices_this_hit,
			 std::vector< std::vector< const anab::CosmicTag* > > const& tags_per_cluster,
			 std::vector<bool> & hitsAccounted_per_tag,
			 std::vector<bool> & hitsAllTags);
    
    void FillAllTagsInfo(recob::Hit const& hit,
			 hit_origin_t const& origin);

    std::vector<float> cTaggedCharge_Cosmic;
    std::vector<float> cTaggedCharge_NonCosmic;
    std::vector<int>   cTaggedHits_Cosmic;
    std::vector<int>   cTaggedHits_NonCosmic;
  //  std::vector<std::string> *cTagAlgorithmNames;

  /*
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
  */
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
      
      float TotalTaggedCharge_Cosmic;
      float TotalTaggedCharge_NonCosmic;
      int   TotalTaggedHits_Cosmic;
      int   TotalTaggedHits_NonCosmic;


    } cEventProperties_t;
    cEventProperties_t cEventVals;

// =====================================================
// fhicl::ParameterSet
// =====================================================
microboone::CosmicRemovalAna::CosmicRemovalAna(fhicl::ParameterSet const& pset):
  EDAnalyzer(pset),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fMCModuleLabel            (pset.get< std::string >("MCModuleLabel")         ),
  fMCHitsModuleLabel        (pset.get< std::string >("MCHitsModuleLabel")         ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")	      ),
  fHitCompareCut            (pset.get< float >("HitCompareCut")	      ),
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
    
  tEventTree->Branch("event", &cEventVals, "runNumber/I:eventNumber/I:nHitsTotal_Unknown/I:nHitsTotal_Cosmic/I:nHitsTotal_NonCosmic/I:qTotal_Unknown/F:qTotal_Cosmic/F:qTotal_NonCosmic/F:nHitsTrack/I:nHitsTrack_Cosmic/I:nHitsTrack_NonCosmic/I:qTrack/F:qTrack_Cosmic/F:qTrack_NonCosmic/F:nHitsCluster/I:nHitsCluster_Cosmic/I:nHitsCluster_NonCosmic/I:qCluster/F:qCluster_Cosmic/F:qCluster_NonCosmic/F:TotalTaggedCharge_Cosmic/F:TotalTaggedCharge_NonCosmic/F:TotalTaggedHits_Cosmic/I:TotalTaggedHits_NonCosmic/I");
  tEventTree->Branch("TaggedCharge_Cosmic",&cTaggedCharge_Cosmic);
  tEventTree->Branch("TaggedCharge_NonCosmic",&cTaggedCharge_NonCosmic);
  tEventTree->Branch("TaggedHits_Cosmic",&cTaggedHits_Cosmic);
  tEventTree->Branch("TaggedHits_NonCosmic",&cTaggedHits_NonCosmic);

}


// =====================================================
// Event Loop
// =====================================================
void microboone::CosmicRemovalAna::analyze(const art::Event& evt)
{

  InitEventTree(evt.run(), evt.event());

  // ##################################################################
  // ### Grabbing ALL HITS in the event to monitor the backtracking ###
  // ##################################################################
  art::Handle< std::vector<raw::RawDigit> > rawdigitHandle;
  evt.getByLabel("daq",rawdigitHandle);

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);
  std::vector<recob::Hit> const& hitVector(*hitListHandle);  

  std::vector<hit_origin_t> hitOrigins(hitVector.size());

  //get mcHitCollection
  art::Handle< std::vector<sim::MCHitCollection> > mchitListHandle;
  evt.getByLabel(fMCHitsModuleLabel,mchitListHandle);
  std::vector<sim::MCHitCollection> const& mchitcolVector(*mchitListHandle);

  //get mcparticles out of the event
  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
  evt.getByLabel(fMCModuleLabel,mcParticleHandle);
  std::vector<simb::MCParticle> const& mcParticleVector(*mcParticleHandle);

  //get associations of mc particles to mc truth
  art::Handle< art::Assns<simb::MCParticle,simb::MCTruth> > assnMCParticleTruthHandle;
  evt.getByLabel(fMCModuleLabel,assnMCParticleTruthHandle);

  std::vector< const simb::MCTruth* > 
    particle_to_truth = util::GetAssociatedVectorOneP(assnMCParticleTruthHandle,
						      mcParticleHandle);

  //make trackId to MCTruth map
  std::map<int,const simb::MCTruth* > trackIDToTruthMap;
  for(size_t p_iter=0; p_iter<mcParticleVector.size(); p_iter++)
    trackIDToTruthMap[ mcParticleVector[p_iter].TrackId() ] = particle_to_truth[p_iter];

  FillMCInfo(hitVector, hitOrigins, mchitcolVector, trackIDToTruthMap);  
 
  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);
  std::vector<recob::Track> const& trackVector(*trackListHandle);  

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);
  std::vector<recob::Cluster> const& clusterVector(*clusterListHandle);  

  art::Handle< art::Assns<recob::Hit,recob::Track> > assnHitTrackHandle;
  evt.getByLabel(fTrackModuleLabel,assnHitTrackHandle);
  std::vector< std::vector<size_t> > 
    track_indices_per_hit = util::GetAssociatedVectorManyI(assnHitTrackHandle,
							   hitListHandle);

  art::Handle< art::Assns<recob::Hit,recob::Cluster> > assnHitClusterHandle;
  evt.getByLabel(fClusterModuleLabel,assnHitClusterHandle);
  std::vector< std::vector<size_t> > 
    cluster_indices_per_hit = util::GetAssociatedVectorManyI(assnHitClusterHandle,
							     hitListHandle);

  std::vector< art::Handle< std::vector<anab::CosmicTag> > > cosmicTagHandlesVector(fCosmicTagAssocLabel.size());
  std::vector< art::Handle< art::Assns<recob::Track,anab::CosmicTag> > > assnTrackTagHandlesVector(fCosmicTagAssocLabel.size());
  std::vector< std::vector< const anab::CosmicTag* > > tags_per_track(trackVector.size(), std::vector<const anab::CosmicTag*>(fCosmicTagAssocLabel.size()));
  std::vector< art::Handle< art::Assns<recob::Cluster,anab::CosmicTag> > > assnClusterTagHandlesVector(fCosmicTagAssocLabel.size());
  std::vector< std::vector< const anab::CosmicTag* > > tags_per_cluster(clusterVector.size(), std::vector<const anab::CosmicTag*>(fCosmicTagAssocLabel.size()));

  for(size_t label_i=0; label_i<fCosmicTagAssocLabel.size(); label_i++){
    try{ evt.getByLabel(fCosmicTagAssocLabel[label_i],cosmicTagHandlesVector[label_i]); }
    catch(...){ continue; }
    try{ 
      evt.getByLabel(fCosmicTagAssocLabel[label_i],assnTrackTagHandlesVector[label_i]); 
      for(auto const& pair : *assnTrackTagHandlesVector[label_i])
	tags_per_track.at(pair.first.key())[label_i] = &(*(pair.second));
    }
    catch(...){}
    try{ 
      evt.getByLabel(fCosmicTagAssocLabel[label_i],assnClusterTagHandlesVector[label_i]); 
      for(auto const& pair : *assnClusterTagHandlesVector[label_i])
	tags_per_cluster.at(pair.first.key())[label_i] = &(*(pair.second));
    }
    catch(...){}
  }

  std::vector< std::vector<bool> > hitsAccounted(hitVector.size(),std::vector<bool>(fCosmicTagAssocLabel.size(),false));
  std::vector<bool> hitsAllTags(hitVector.size(),false);

  for(size_t hit_iter=0; hit_iter<hitVector.size(); hit_iter++){

    float charge = hitVector[hit_iter].Integral();
    hit_origin_t origin = hitOrigins[hit_iter];

    if(track_indices_per_hit[hit_iter].size()!=0)
      FillTrackInfo(hit_iter,
		    origin,
		    charge,
		    track_indices_per_hit[hit_iter],
		    tags_per_track,
		    hitsAccounted[hit_iter],
		    hitsAllTags);

    if(cluster_indices_per_hit[hit_iter].size()!=0)
      FillClusterInfo(hit_iter,
		      origin,
		      charge,
		      cluster_indices_per_hit[hit_iter],
		      tags_per_cluster,
		      hitsAccounted[hit_iter],
		      hitsAllTags);
    
    if(hitsAllTags[hit_iter]) FillAllTagsInfo(hitVector[hit_iter],origin);

  }//end loop over all the hits

  tEventTree->Fill();
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
  
  cEventVals.TotalTaggedCharge_Cosmic = 0.;
  cEventVals.TotalTaggedCharge_NonCosmic = 0.;
  cEventVals.TotalTaggedHits_Cosmic = 0;
  cEventVals.TotalTaggedHits_NonCosmic = 0;

  for(size_t iter=0; iter<fCosmicTagAssocLabel.size(); iter++){
    cTaggedHits_Cosmic.at(iter) = 0;
    cTaggedCharge_Cosmic.at(iter) = 0;
    cTaggedHits_NonCosmic.at(iter) = 0;
    cTaggedCharge_NonCosmic.at(iter) = 0;
  }

}

//-------------------------------------------------------------------------------------------------------------------
// take in a list of hits, and determine the origin for those hits (and fill in the tree info)
void microboone::CosmicRemovalAna::FillMCInfo( std::vector<recob::Hit> const& hitlist,
					       std::vector<hit_origin_t> & hitOrigins,
					       std::vector<sim::MCHitCollection> const& mchitCollectionVector,
					       std::map<int,const simb::MCTruth* > const& trackIdToTruthMap){

  auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();

  for(size_t itr=0; itr<hitlist.size(); itr++){

    recob::Hit const& this_hit = hitlist[itr];
    
    std::vector<int> trackIDs;
    std::vector<double> energy;    

    for( auto const& mchit : mchitCollectionVector[this_hit.Channel()] ){
      if( std::abs(ts->TPCTDC2Tick(mchit.PeakTime()) - this_hit.PeakTime()) < fHitCompareCut){
	trackIDs.push_back(mchit.PartTrackId());
	energy.push_back(mchit.PartEnergy());
      }
    }
    
    if(trackIDs.size()==0){
      hitOrigins[itr] = hit_origin_Unknown;
      cEventVals.nHitsTotal_Unknown++;
      cEventVals.qTotal_Unknown += this_hit.Integral();
      continue;
    }
    
    float cosmic_energy=0;
    float non_cosmic_energy=0;

    for(size_t iter=0; iter<trackIDs.size(); iter++){ 
      auto map_element = trackIdToTruthMap.find(std::abs(trackIDs[iter]));
      if(map_element==trackIdToTruthMap.end()) continue;
      int origin = map_element->second->Origin();  
      if(origin == simb::kBeamNeutrino)
	non_cosmic_energy += energy[iter];
      else
	cosmic_energy += energy[iter];
    }
    
    if(non_cosmic_energy > cosmic_energy){
      hitOrigins[itr] = hit_origin_NonCosmic;
      cEventVals.nHitsTotal_NonCosmic++;
      cEventVals.qTotal_NonCosmic += this_hit.Integral();
    }
    else{
      hitOrigins[itr] = hit_origin_Cosmic;
      cEventVals.nHitsTotal_Cosmic++;
      cEventVals.qTotal_Cosmic += this_hit.Integral();
    }
  }

}//end FillMCInfo

//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillTrackInfo(size_t const& hit_iter,
						 hit_origin_t const& origin,
						 float const& charge,
						 std::vector<size_t> const& track_indices_this_hit,
						 std::vector< std::vector< const anab::CosmicTag* > > const& tags_per_track,
						 std::vector<bool> & hitsAccounted_per_tag,
						 std::vector<bool> & hitsAllTags){

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
    if(hitsAccounted_per_tag[nCT]) continue;

    for(auto const& track_index : track_indices_this_hit){
      if(!tags_per_track[track_index][nCT]) continue;
      const anab::CosmicTag* currentTag(tags_per_track[track_index][nCT]);
      if( currentTag->CosmicScore() > fCosmicScoreThresholds[nCT] ){

	hitsAccounted_per_tag[nCT] = true;
	hitsAllTags[hit_iter] = true;
	if(origin==hit_origin_Cosmic){
	  cTaggedHits_Cosmic[nCT]++;
	  cTaggedCharge_Cosmic[nCT] += charge;
	}
	else if(origin==hit_origin_NonCosmic){
	  cTaggedHits_NonCosmic[nCT]++;
	  cTaggedCharge_NonCosmic[nCT] += charge;
	}
      }
    }
    
  }
  
}//end FillTrackInfo

//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillClusterInfo(size_t const& hit_iter,
						   hit_origin_t const& origin,
						   float const& charge,
						   std::vector<size_t> const& cluster_indices_this_hit,
						   std::vector< std::vector< const anab::CosmicTag* > > const& tags_per_cluster,
						   std::vector<bool> & hitsAccounted_per_tag,
						   std::vector<bool> & hitsAllTags){

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
    if(hitsAccounted_per_tag[nCT]) continue;

    for(auto const& cluster_index : cluster_indices_this_hit){
      if(!tags_per_cluster[cluster_index][nCT]) continue;
      const anab::CosmicTag* currentTag(tags_per_cluster[cluster_index][nCT]);
      if( currentTag->CosmicScore() > fCosmicScoreThresholds[nCT] ){

	hitsAccounted_per_tag[nCT] = true;
	hitsAllTags[hit_iter] = true;
	if(origin==hit_origin_Cosmic){
	  cTaggedHits_Cosmic[nCT]++;
	  cTaggedCharge_Cosmic[nCT] += charge;
	}
	else if(origin==hit_origin_NonCosmic){
	  cTaggedHits_NonCosmic[nCT]++;
	  cTaggedCharge_NonCosmic[nCT] += charge;
	}
      }
    }
    

  }
  
}//end FillClusterInfo


//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillAllTagsInfo(recob::Hit const& hit,
						   hit_origin_t const& origin){
  if(origin==hit_origin_Cosmic){
    cEventVals.TotalTaggedCharge_Cosmic += hit.Integral();
    cEventVals.TotalTaggedHits_Cosmic++;
  }
  else if(origin==hit_origin_NonCosmic){
    cEventVals.TotalTaggedCharge_NonCosmic += hit.Integral();
    cEventVals.TotalTaggedHits_NonCosmic++;
  }

}

namespace microboone{
  
  DEFINE_ART_MODULE(CosmicRemovalAna)  
}








#endif
