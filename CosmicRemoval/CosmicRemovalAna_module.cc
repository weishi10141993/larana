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
#include "Utilities/TimeService.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "MCBase/MCHitCollection.h"
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
    
    TH1F *th0;
    TH1F *th1;
    TH1F *th2;
    TH1F *th3;
    TH1F *th4;
    TH1F *th5;
    TH1F *th6;
    TH1F *th7;
    TH1F *th12;
    TH1F *th17;
    TH1F *th18;
    
    std::string fHitsModuleLabel;
    std::string fMCModuleLabel;
    std::string fMCHitsModuleLabel;
    std::string fClusterModuleLabel; 
    std::string fTrackModuleLabel;
    std::vector <std::string> fCosmicTagAssocLabel;
    std::vector <float> fCosmicScoreThresholds;

    void InitEventTree(int run_number, int event_number);
    
    void FillMCInfo( std::vector< art::Ptr<recob::Hit> > const& hitlist,
		     std::vector<hit_origin_t> & hitOrigins,
		     double & time);
    
    void FillMCInfoNew( std::vector< art::Ptr<recob::Hit> > const& hitlist,
			std::vector<hit_origin_t> & hitOrigins,
			double & time,
			std::vector<sim::MCHitCollection> const& mchitCollectionVector,
			std::map<int,const simb::MCTruth* > const& trackIDToTruthMap);

    void FillTrackInfo( size_t const& hit_iter,
			hit_origin_t const& origin,
			float const& charge,
			art::FindManyP<recob::Track> const& tracks_per_hit,
			art::Event const& evt,
			std::vector<bool> & hitsAccounted_per_tag,
			double & time);
    
    void FillTrackInfoNew(size_t const& hit_iter,
			    hit_origin_t const& origin,
			    float const& charge,
			    std::vector< art::Ptr<recob::Track> > const& tracks_this_hit,
			    std::vector< std::vector< art::Ptr<anab::CosmicTag> > > const& tags_per_cluster,
			    std::vector<bool> & hitsAccounted_per_tag,
			    double & time);

    void FillClusterInfo( size_t const& hit_iter,
			  hit_origin_t const& origin,
			  float const& charge,
			  art::FindManyP<recob::Cluster> const& clusters_per_hit,
			  art::Event const& evt,
			  std::vector<bool> & hitsAccounted_per_tag,
			  double & time);
    
    void FillClusterInfoNew(size_t const& hit_iter,
			    hit_origin_t const& origin,
			    float const& charge,
			    std::vector< art::Ptr<recob::Cluster> > const& clusters_this_hit,
			    std::vector< std::vector< art::Ptr<anab::CosmicTag> > > const& tags_per_cluster,
			    std::vector<bool> & hitsAccounted_per_tag,
			    double & time);

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
  fMCModuleLabel            (pset.get< std::string >("MCModuleLabel")         ),
  fMCHitsModuleLabel        (pset.get< std::string >("MCHitsModuleLabel")         ),
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

  th0  = (TH1F*)tfs->make<TH1F>("h0_initTree","initTree;Time(s);Events",200,0.,2.);
  th1  = (TH1F*)tfs->make<TH1F>("h1_readHits","readHits;Time(s);Events",200,0.,2.);
  th2  = (TH1F*)tfs->make<TH1F>("h2_EveID","EveID;Time(s);Events",200,0.,2.);
  th3  = (TH1F*)tfs->make<TH1F>("h3_FillMCInfo","FillMCInfo;Time(s);Events",200,0.,2.);
  th4  = (TH1F*)tfs->make<TH1F>("h4_EveIDNew","EveIDNew;Time(s);Events",200,0.,2.);
  th5  = (TH1F*)tfs->make<TH1F>("h5_FillMCInfoNew","FillMCInfoNew;Time(s);Events",200,0.,2.);
  th6  = (TH1F*)tfs->make<TH1F>("h6_GrabAssociations","GrabAssociations;Time(s);Events",200,0.,2.);
  th7  = (TH1F*)tfs->make<TH1F>("h7_GrabAssociationsNew","GrabAssociationsNew;Time(s);Events",200,0.,2.);
  th12 = (TH1F*)tfs->make<TH1F>("h12_GetRecoInfo","GetRecoInfo;Time(s);Events",200,0.,2.);
  th17 = (TH1F*)tfs->make<TH1F>("h17_GetRecoInfoNew","GetRecoInfoNew;Time(s);Events",200,0.,2.);
  th18 = (TH1F*)tfs->make<TH1F>("h18_FillTree","FillTree;Time(s);Events",200,0.,2.);
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

  art::Handle< std::vector<sim::MCHitCollection> > mchitListHandle;
  evt.getByLabel(fMCHitsModuleLabel,mchitListHandle);
  std::vector<sim::MCHitCollection> const& mchitcolVector(*mchitListHandle);

  sw.Stop();
  times.push_back( std::make_pair("readHits",sw.RealTime()) );
  sw.Start();

  std::vector<hit_origin_t> hitOrigins(hitlist.size());
  double timeEveID=0;
  FillMCInfo(hitlist, hitOrigins, timeEveID);  
 
  sw.Stop();
  times.push_back( std::make_pair("EveID",timeEveID) );
  times.push_back( std::make_pair("FillMCInfo",sw.RealTime()) );
  float ei = sw.RealTime();
  sw.Start();
  
  std::vector<hit_origin_t> hitOriginsNew(hitlist.size());

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
  evt.getByLabel(fMCModuleLabel,mcParticleHandle);
  std::vector<simb::MCParticle> const& mcParticleVector(*mcParticleHandle);

  art::Handle< art::Assns<simb::MCParticle,simb::MCTruth> > assnMCParticleTruthHandle;
  evt.getByLabel(fMCModuleLabel,assnMCParticleTruthHandle);
  art::Assns<simb::MCParticle,simb::MCTruth> const& assnMCParticleTruth(*assnMCParticleTruthHandle);

  std::vector< const simb::MCTruth* > particle_to_truth(mcParticleVector.size());
  for(auto const& pair : assnMCParticleTruth)
    particle_to_truth.at(pair.first.key()) = &(*(pair.second));

  std::map<int,const simb::MCTruth* > trackIDToTruthMap;
  for(size_t p_iter=0; p_iter<mcParticleVector.size(); p_iter++)
    trackIDToTruthMap[ mcParticleVector[p_iter].TrackId() ] = particle_to_truth[p_iter];

  double timeEveIDNew=0;
  FillMCInfoNew(hitlist, hitOriginsNew, timeEveIDNew, mchitcolVector, trackIDToTruthMap);  
 
  sw.Stop();
  times.push_back( std::make_pair("EveIDNew",timeEveIDNew) );
  times.push_back( std::make_pair("FillMCInfoNew",sw.RealTime()) );
  float ei_new = sw.RealTime();
  
  size_t ndifferences=0;
  for(size_t hitlist_iter=0; hitlist_iter<hitOrigins.size(); hitlist_iter++){
    if( hitOrigins[hitlist_iter] != hitOriginsNew[hitlist_iter]){
      std::cout << "\tDIFFERNCE!!! " << hitOrigins[hitlist_iter] << " != " << hitOriginsNew[hitlist_iter] << std::endl;
      ndifferences++;
    }
  }
  std::cout << "Out of a total of " << hitOrigins.size() << " there were " << ndifferences << " differences." << std::endl;

  sw.Start();

  art::FindManyP<recob::Track> tracks_per_hit(hitlist, evt, fTrackModuleLabel);
  art::FindManyP<recob::Cluster> clusters_per_hit(hitlist, evt, fClusterModuleLabel);

  sw.Stop();
  times.push_back( std::make_pair("GrabAssociations",sw.RealTime()) );
  float ga = sw.RealTime();
  sw.Start();

  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);
  std::vector<art::Ptr<recob::Track> > tracklist;  
  art::fill_ptr_vector(tracklist, trackListHandle);

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);
  std::vector<art::Ptr<recob::Cluster> > clusterlist;  
  art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::Handle< art::Assns<recob::Hit,recob::Track> > assnHitTrackHandle;
  evt.getByLabel(fTrackModuleLabel,assnHitTrackHandle);
  art::Assns<recob::Hit,recob::Track> const& assnHitTrack(*assnHitTrackHandle);

  std::vector< std::vector<art::Ptr<recob::Track> > > tracks_per_hit_new(hitlist.size());
  for(auto const& pair : assnHitTrack)
    tracks_per_hit_new.at(pair.first.key()).push_back(pair.second);

  art::Handle< art::Assns<recob::Hit,recob::Cluster> > assnHitClusterHandle;
  evt.getByLabel(fClusterModuleLabel,assnHitClusterHandle);
  art::Assns<recob::Hit,recob::Cluster> const& assnHitCluster(*assnHitClusterHandle);

  std::vector< std::vector<art::Ptr<recob::Cluster> > > clusters_per_hit_new(hitlist.size());
  for(auto const& pair : assnHitCluster)
    clusters_per_hit_new.at(pair.first.key()).push_back(pair.second);

  sw.Stop();
  times.push_back( std::make_pair("GrabAssociationsNew",sw.RealTime()) );
  float ga_new = sw.RealTime();

  ndifferences=0;
  for(size_t hitlist_iter=0; hitlist_iter<tracks_per_hit_new.size(); hitlist_iter++){
    if( tracks_per_hit_new[hitlist_iter].size() != tracks_per_hit.at(hitlist_iter).size()){
      std::cout << "\tDIFFERNCE!!! " << tracks_per_hit.at(hitlist_iter).size() << " != " << tracks_per_hit_new[hitlist_iter].size() << std::endl;
      ndifferences++;
    }
  }
  std::cout << "Track assns: Out of a total of " << tracks_per_hit.size() << " there were " << ndifferences << " differences." << std::endl;

  ndifferences=0;
  for(size_t hitlist_iter=0; hitlist_iter<clusters_per_hit_new.size(); hitlist_iter++){
    if( clusters_per_hit_new[hitlist_iter].size() != clusters_per_hit.at(hitlist_iter).size()){
      std::cout << "\tDIFFERNCE!!! " << clusters_per_hit.at(hitlist_iter).size() << " != " << clusters_per_hit_new[hitlist_iter].size() << std::endl;
      ndifferences++;
    }
  }
  std::cout << "Cluster assns: Out of a total of " << clusters_per_hit.size() << " there were " << ndifferences << " differences." << std::endl;

  sw.Start();

  std::vector< std::vector<bool> > hitsAccounted(hitlist.size(),std::vector<bool>(fCosmicTagAssocLabel.size(),false));

  double timeTrack=0, timeCluster=0;
  double timeTrackTotal=0, timeClusterTotal=0;
  TStopwatch sw_inner;

  for(size_t hit_iter=0; hit_iter<hitlist.size(); hit_iter++){

    float charge = hitlist.at(hit_iter)->Charge();
    hit_origin_t origin = hitOrigins.at(hit_iter);

    sw_inner.Start();
    FillTrackInfo(hit_iter,
		  origin,
		  charge,
		  tracks_per_hit,
		  evt,
		  hitsAccounted.at(hit_iter),
		  timeTrack);
    sw_inner.Stop(); timeTrackTotal += sw_inner.RealTime();

    sw_inner.Start();
    FillClusterInfo(hit_iter,
		    origin,
		    charge,
		    clusters_per_hit,
		    evt,
		    hitsAccounted.at(hit_iter),
		    timeCluster);
    sw_inner.Stop(); timeClusterTotal += sw_inner.RealTime();

  }//end loop over all the hits

  sw.Stop();
  times.push_back( std::make_pair("TrackTime",timeTrack) );
  times.push_back( std::make_pair("TrackTimeTotal",timeTrackTotal) );
  times.push_back( std::make_pair("ClusterTime",timeCluster) );
  times.push_back( std::make_pair("ClusterTimeTotal",timeClusterTotal) );
  times.push_back( std::make_pair("GetRecoInfo",sw.RealTime()) );
  sw.Start();
  

  std::vector< art::Handle< std::vector<anab::CosmicTag> > > cosmicTagHandlesVector(fCosmicTagAssocLabel.size());
  std::vector< art::Handle< art::Assns<recob::Track,anab::CosmicTag> > > assnTrackTagHandlesVector(fCosmicTagAssocLabel.size());
  std::vector< std::vector< art::Ptr<anab::CosmicTag> > > tags_per_track(tracklist.size(), std::vector< art::Ptr<anab::CosmicTag> >(fCosmicTagAssocLabel.size()));
  std::vector< art::Handle< art::Assns<recob::Cluster,anab::CosmicTag> > > assnClusterTagHandlesVector(fCosmicTagAssocLabel.size());
  std::vector< std::vector< art::Ptr<anab::CosmicTag> > > tags_per_cluster(clusterlist.size(), std::vector< art::Ptr<anab::CosmicTag> >(fCosmicTagAssocLabel.size()));
  for(size_t label_i=0; label_i<fCosmicTagAssocLabel.size(); label_i++){
    try{
      evt.getByLabel(fCosmicTagAssocLabel[label_i],cosmicTagHandlesVector[label_i]);
      //std::cout << "Got the cosmics " << fCosmicTagAssocLabel[label_i] << std::endl;
    }
    catch(...){ continue; }
    try{ 
      evt.getByLabel(fCosmicTagAssocLabel[label_i],assnTrackTagHandlesVector[label_i]); 
      for(auto const& pair : *assnTrackTagHandlesVector[label_i])
	tags_per_track.at(pair.first.key())[label_i] = pair.second;
      std::cout << "Got the track associations " << fCosmicTagAssocLabel[label_i] << std::endl;
    }
    catch(...){}
    try{ 
      evt.getByLabel(fCosmicTagAssocLabel[label_i],assnClusterTagHandlesVector[label_i]); 
      for(auto const& pair : *assnClusterTagHandlesVector[label_i])
	tags_per_cluster.at(pair.first.key())[label_i] = pair.second;
      std::cout << "Got the cluster associations " << fCosmicTagAssocLabel[label_i] << std::endl;
    }
    catch(...){}
  }

  std::vector< std::vector<bool> > hitsAccounted_new(hitlist.size(),std::vector<bool>(fCosmicTagAssocLabel.size(),false));

  timeTrack=0, timeCluster=0;
  timeTrackTotal=0, timeClusterTotal=0;

  for(size_t hit_iter=0; hit_iter<hitlist.size(); hit_iter++){

    float charge = hitlist.at(hit_iter)->Charge();
    hit_origin_t origin = hitOrigins.at(hit_iter);

    sw_inner.Start();

    if(tracks_per_hit_new[hit_iter].size()!=0)
      FillTrackInfoNew(hit_iter,
		       origin,
		       charge,
		       tracks_per_hit_new[hit_iter],
		       tags_per_track,
		       hitsAccounted_new[hit_iter],
		       timeTrack);

    sw_inner.Stop(); timeTrackTotal += sw_inner.RealTime();

    sw_inner.Start();

    if(clusters_per_hit_new[hit_iter].size()!=0)
      FillClusterInfoNew(hit_iter,
			 origin,
			 charge,
			 clusters_per_hit_new[hit_iter],
			 tags_per_cluster,
			 hitsAccounted_new[hit_iter],
			 timeCluster);

    sw_inner.Stop(); timeClusterTotal += sw_inner.RealTime();
    
  }//end loop over all the hits

  sw.Stop();
  times.push_back( std::make_pair("TrackTimeNew",timeTrack) );
  times.push_back( std::make_pair("TrackTimeTotalNew",timeTrackTotal) );
  times.push_back( std::make_pair("ClusterTimeNew",timeCluster) );
  times.push_back( std::make_pair("ClusterTimeTotalNew",timeClusterTotal) );
  times.push_back( std::make_pair("GetRecoInfoNew",sw.RealTime()) );

  sw.Start();
  tEventTree->Fill();
  sw.Stop();
  times.push_back( std::make_pair("FillTree",sw.RealTime()) );
  sw.Start();


  for(size_t t_iter=0; t_iter<times.size(); t_iter++){
    auto const& time_check = times[t_iter];
    std::cout << time_check.first << ":" << time_check.second << std::endl;
  }
  th0->Fill(times[0].second);
  th1->Fill(times[1].second);
  th2->Fill(times[2].second);
  th3->Fill(times[3].second);
  th4->Fill(times[4].second);
  th5->Fill(times[5].second);
  th6->Fill(times[6].second);
  th7->Fill(times[7].second);
  th12->Fill(times[12].second);
  th17->Fill(times[17].second);
  th18->Fill(times[18].second);
  std::cout << "Comps: " << ei/ei_new << " and " << ga/ga_new << std::endl;
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
					       std::vector<hit_origin_t> & hitOrigins,
					       double & timeEveID){

  art::ServiceHandle<cheat::BackTracker> bt;
  TStopwatch sw; timeEveID=0;
  
  for(size_t itr=0; itr<hitlist.size(); itr++){
    
    art::Ptr<recob::Hit> const& hitptr = hitlist.at(itr);
    sw.Start();
    std::vector<sim::TrackIDE> eveIDs = bt->HitToEveID(hitptr);
    sw.Stop(); timeEveID += sw.RealTime();
    

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
// take in a list of hits, and determine the origin for those hits (and fill in the tree info)
void microboone::CosmicRemovalAna::FillMCInfoNew( std::vector< art::Ptr<recob::Hit> > const& hitlist,
						  std::vector<hit_origin_t> & hitOrigins,
						  double & timeEveID,
						  std::vector<sim::MCHitCollection> const& mchitCollectionVector,
						  std::map<int,const simb::MCTruth* > const& trackIdToTruthMap){

  art::ServiceHandle<util::TimeService> ts;
  TStopwatch sw; timeEveID=0;
  
  for(size_t itr=0; itr<hitlist.size(); itr++){
    
    art::Ptr<recob::Hit> const& hitptr = hitlist.at(itr);
    sw.Start();
    
    std::vector<int> trackIDs;
    std::vector<double> energy;

    //std::cout << "Hit we are checking is on ch " << hitptr->Channel() << " at time " << hitptr->PeakTime() << std::endl;

    for( auto const& mchit : mchitCollectionVector.at(hitptr->Channel()) ){
      //std::cout << "\tMCHit we are checking is at time " << ts->TPCTDC2Tick(mchit.PeakTime()) << std::endl;
      if( std::abs(ts->TPCTDC2Tick(mchit.PeakTime()) - hitptr->PeakTime()) < 4.){
	trackIDs.push_back(mchit.PartTrackId());
	energy.push_back(mchit.PartEnergy());
      }
    }

    sw.Stop(); timeEveID += sw.RealTime();
    

    if(trackIDs.size()==0){
      hitOrigins.at(itr) = hit_origin_Unknown;
      cEventVals.nHitsTotal_Unknown++;
      cEventVals.qTotal_Unknown += hitptr->Charge();
      continue;
    }
    
    float cosmic_energy=0;
    float non_cosmic_energy=0;

    for(size_t iter=0; iter<trackIDs.size(); iter++){ 
      int origin = trackIdToTruthMap.at(std::abs(trackIDs[iter]))->Origin();  
      if(origin == simb::kBeamNeutrino)
	non_cosmic_energy += energy[iter];
      else
	cosmic_energy += energy[iter];
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

}//end FillMCInfoNew

//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillTrackInfo(size_t const& hit_iter,
						 hit_origin_t const& origin,
						 float const& charge,
						 art::FindManyP<recob::Track> const& tracks_per_hit,
						 art::Event const& evt,
						 std::vector<bool> & hitsAccounted_per_tag,
						 double & time){

  TStopwatch sw;

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
  
  //sw.Stop(); time+= sw.RealTime();
  
  for(unsigned int nCT = 0; nCT < fCosmicTagAssocLabel.size(); nCT++){//<---This loops over the vector of cosmicTags in stored in the event
    
    if(hitsAccounted_per_tag.at(nCT)) continue;
    
    sw.Start();
    try{ 
      //sw.Start();
      art::FindManyP<anab::CosmicTag> cosmic_tags_per_track(tracks_this_hit,evt,fCosmicTagAssocLabel.at(nCT));
      //sw.Stop(); time+= sw.RealTime();
      //for(auto const& cosmic_tags_this_track : cosmic_tags_per_track){ 
      //sw.Start();
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
      //sw.Stop(); time+= sw.RealTime();
      
    }
    catch (...){       
      sw.Stop(); time+= sw.RealTime();
      return;
    }
    
    
  }
  
}//end FillTrackInfo


//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillTrackInfoNew(size_t const& hit_iter,
						    hit_origin_t const& origin,
						    float const& charge,
						    std::vector< art::Ptr<recob::Track> > const& tracks_this_hit,
						    std::vector< std::vector< art::Ptr<anab::CosmicTag> > > const& tags_per_track,
						    std::vector<bool> & hitsAccounted_per_tag,
						    double & time){

  TStopwatch sw;
  
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
  
  //sw.Stop(); time+= sw.RealTime();
  
  for(unsigned int nCT = 0; nCT < fCosmicTagAssocLabel.size(); nCT++){//<---This loops over the vector of cosmicTags in stored in the event
    if(hitsAccounted_per_tag.at(nCT)) continue;
    sw.Start();

    for(size_t track_iter=0; track_iter<tracks_this_hit.size(); track_iter++){
      if(!tags_per_track[track_iter][nCT]) continue;
      art::Ptr<anab::CosmicTag> const& currentTag(tags_per_track[track_iter][nCT]);
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
    

    sw.Stop(); time+= sw.RealTime(); 
  }
  
}//end FillTrackInfo

//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillClusterInfo(size_t const& hit_iter,
						   hit_origin_t const& origin,
						   float const& charge,
						   art::FindManyP<recob::Cluster> const& clusters_per_hit,
						   art::Event const& evt,
						   std::vector<bool> & hitsAccounted_per_tag,
						   double & time){

  TStopwatch sw;

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
      sw.Start();
      art::FindManyP<anab::CosmicTag> cosmic_tags_per_cluster(clusters_this_hit,evt,fCosmicTagAssocLabel.at(nCT));
      sw.Stop(); time += sw.RealTime();
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

//-------------------------------------------------------------------------------------------------------------------
void microboone::CosmicRemovalAna::FillClusterInfoNew(size_t const& hit_iter,
						    hit_origin_t const& origin,
						    float const& charge,
						    std::vector< art::Ptr<recob::Cluster> > const& clusters_this_hit,
						    std::vector< std::vector< art::Ptr<anab::CosmicTag> > > const& tags_per_cluster,
						    std::vector<bool> & hitsAccounted_per_tag,
						    double & time){

  TStopwatch sw;
  
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
  
  //sw.Stop(); time+= sw.RealTime();
  
  for(unsigned int nCT = 0; nCT < fCosmicTagAssocLabel.size(); nCT++){//<---This loops over the vector of cosmicTags in stored in the event
    if(hitsAccounted_per_tag.at(nCT)) continue;
    sw.Start();

    for(size_t cluster_iter=0; cluster_iter<clusters_this_hit.size(); cluster_iter++){
      if(!tags_per_cluster[cluster_iter][nCT]) continue;
      art::Ptr<anab::CosmicTag> const& currentTag(tags_per_cluster[cluster_iter][nCT]);
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
    

    sw.Stop(); time+= sw.RealTime(); 
  }
  
}//end FillClusterInfo


namespace microboone{
  
  DEFINE_ART_MODULE(CosmicRemovalAna)  
}








#endif
