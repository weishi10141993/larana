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

    //    double GetOverlapScore(std::vector<art::Ptr<recob::Hit> >& AllHits, std::vector<int>& HitsThisParticle, std::set<art::ProductID>&  TaggedHits);
    //    double GetOverlapScore(std::vector<art::Ptr<recob::Hit> >& AllHits, std::vector<int>& HitsThisParticle, std::set<art::Ptr<recob::Hit>::key_type>&  TaggedHits);
    double GetOverlapScore(std::vector<art::Ptr<recob::Hit> >& AllHits, std::map<int, std::vector<int> >& HitsThisParticle, std::set<art::Ptr<recob::Hit>::key_type>&  TaggedHits);


    int GetOrigin( std::vector<art::Ptr<recob::Hit> >& hitlist, const art::Event& evt, int &pdg, float &energy);

    void GetNewOverlapScore( std::vector< float > &scoreVector, 
			      std::vector< int > &hitsVector,
			     std::vector < std::vector<int> > &originVector, 
			     float &cOverA, float &dOverB );
			     
    void FillTrackTree( art::Ptr<anab::CosmicTag> &currentTag, float CosmicScore, const art::Event& evt);

    unsigned int nCosmicTags;
    
    std::vector<TH1D*> fCosmicScoresPerCT;
    std::vector<TH1D*> fFractionChargeTaggedPerCT_Cosmic;
    std::vector<TH1D*> fFractionChargeTaggedPerCT_NonCosmic;
    TH1D * fNAlgsRejected60_Cosmic;
    TH1D * fNAlgsRejected60_NonCosmic;
    TH1D * fNAlgsRejected80_Cosmic;
    TH1D * fNAlgsRejected80_NonCosmic;
    TH1D * fNAlgsRejected95_Cosmic;
    TH1D * fNAlgsRejected95_NonCosmic;
    TH1D * fTotalCharge_Cosmic;
    TH1D * fTotalCharge_NonCosmic;
    TH1D * fTrackedFraction_Cosmic;
    TH1D * fTrackedFraction_NonCosmic;
    TH1D * fClusteredFraction_Cosmic;
    TH1D * fClusteredFraction_NonCosmic;

    TH1F *fTPCCosmicType;

    TTree *tTreeTrack;
    TTree *tTreeEvent;

    std::string fGenieGenModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fHitsModuleLabel;
    std::string fClusterModuleLabel; 
    std::string fTrackModuleLabel;
    std::vector <std::string> fCosmicTagAssocLabel;
    std::vector <float> fCosmicScoreThresholds;
  };//<---End
}

typedef struct{
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
} cTTrack_t;
cTTrack_t cTrack;

typedef struct{
  int eventNumber;
  int nHitsTotal;
  int nHitsEvID;
  int nHitsTrack;
  int nHitsTrackEvID;
  float cOverA;  
  float dOverB;
  int A;
  int B;
  int C;
  int D;
} cTEvent_t;
cTEvent_t cEvent;

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

void microboone::CosmicRemovalAna::beginJob()
{
  
  nCosmicTags = fCosmicTagAssocLabel.size();
  
  art::ServiceHandle<art::TFileService> tfs;
  
  // ###################
  // ### Set up TH1s ###
  // ###################
  
  double TotalChargeLimit = 10000;
  
  fNAlgsRejected60_Cosmic = (TH1D*)tfs->make<TH1D>("NAlgsRejected60Cosmic","Number of algorithms rejecting 60% charge, per cosmic particle; N_algs; Particles", nCosmicTags, 0, nCosmicTags);
  fNAlgsRejected60_NonCosmic = (TH1D*)tfs->make<TH1D>("NAlgsRejected60NonCosmic","Number of algorithms rejecting 60% charge, per non-cosmic particle; N_algs; Particles", nCosmicTags, 0, nCosmicTags);

  fNAlgsRejected80_Cosmic = (TH1D*)tfs->make<TH1D>("NAlgsRejected80Cosmic","Number of algorithms rejecting 80% charge, per cosmic particle; N_algs; Particles", nCosmicTags,0, nCosmicTags);
  fNAlgsRejected80_NonCosmic = (TH1D*)tfs->make<TH1D>("NAlgsRejected80NonCosmic","Number of algorithms rejecting 80% charge, per non-cosmic particle; N_algs; Particles", nCosmicTags, 0, nCosmicTags);

  fNAlgsRejected95_Cosmic = (TH1D*)tfs->make<TH1D>("NAlgsRejected95Cosmic","Number of algorithms rejecting 95% charge, per cosmic particle; N_algs; Particles", nCosmicTags, 0, nCosmicTags);
  fNAlgsRejected95_NonCosmic = (TH1D*)tfs->make<TH1D>("NAlgsRejected95NonCosmic","Number of algorithms rejecting 95% charge, per non-cosmic particle; N_algs; Particles", nCosmicTags, 0, nCosmicTags);

  fTotalCharge_Cosmic = (TH1D*)tfs->make<TH1D>("TotalChargeCosmic", "Total Hit Charge for True Cosmic Particles; Charge; N", 100,0,TotalChargeLimit);
  fTotalCharge_NonCosmic = (TH1D*)tfs->make<TH1D>("TotalChargeNonCosmic", "Total Hit Charge for True NonCosmic Particles; Charge; N", 100,0,TotalChargeLimit);

  fTrackedFraction_Cosmic = (TH1D*)tfs->make<TH1D>("TrackedFractionCosmic", "Tracked Hit Charge Fraction for True Cosmic Particles; Charge; N", 100,0,1.5);
  fTrackedFraction_NonCosmic = (TH1D*)tfs->make<TH1D>("TrackedFractionNonCosmic", "Tracked Charge Fraction for True NonCosmic Particles; Charge; N", 100,0,1.5);

  fClusteredFraction_Cosmic = (TH1D*)tfs->make<TH1D>("ClusteredFractionCosmic", "Clustered Hit Charge Fraction for True Cosmic Particles; Charge; N", 100,0,1.5);
  fClusteredFraction_NonCosmic = (TH1D*)tfs->make<TH1D>("ClusteredFractionNonCosmic", "Clustered Charge Fraction for True NonCosmic Particles; Charge; N", 100,0, 1.5);

  fTPCCosmicType = (TH1F*)tfs->make<TH1F>("fTPCCosmicType","TPC Cosmic Type", 5, -0.5, 4.5);
  
  // #############################################################
  // ### Titling the TH1's uniquely for each CosmicRemoval Tag ###
  // #############################################################
  for(size_t i=0; i!=nCosmicTags; ++i)
    {
      std::stringstream sname, stitle;
      sname.str("");  sname.flush();
      stitle.str(""); stitle.flush();

      sname<<"CosmicScoresFor"<<fCosmicTagAssocLabel.at(i);
      stitle<<"Cosmic score per object for " << fCosmicTagAssocLabel.at(i)<<"; Score; N";
      fCosmicScoresPerCT.push_back( (TH1D*)tfs->make<TH1D>(sname.str().c_str(), stitle.str().c_str(), 101,0,1.01));

      sname.str("");  sname.flush();
      stitle.str(""); stitle.flush();
				    
      sname<<"FractionTaggedCosmicFor"<<fCosmicTagAssocLabel.at(i);
      stitle<<"Fraction of Cosmic Charge Tagged as Cosmic For " << fCosmicTagAssocLabel.at(i)<<"; Frac; N";
      fFractionChargeTaggedPerCT_Cosmic.push_back( (TH1D*)tfs->make<TH1D>(sname.str().c_str(), stitle.str().c_str(), 101,0,1.01));	  
      
      sname.str("");  sname.flush();
      stitle.str(""); stitle.flush();
      
      sname<<"FractionTaggedNonCosmicFor"<<fCosmicTagAssocLabel.at(i);
      stitle<<"Fraction of NonCosmic Charge Tagged as Cosmic For " << fCosmicTagAssocLabel.at(i)<<"; Frac; N";
      fFractionChargeTaggedPerCT_NonCosmic.push_back( (TH1D*)tfs->make<TH1D>(sname.str().c_str(), stitle.str().c_str(), 101,0,1.01));												   
    }
  
  // ##################################################
  // ### Setting up TTree for Track based CosmicTag ###
  // ##################################################
  tTreeTrack = tfs->make<TTree>("CosmicTree","CosmicTree");
  tTreeTrack->Branch("eventNumber", &cTrack.eventNumber, "eventNumber/I");
  tTreeTrack->Branch("tagType"	  , &cTrack.tagType    , "tagType/I");
  tTreeTrack->Branch("x0"	  , &cTrack.x0         , "x0/F");
  tTreeTrack->Branch("x1"	  , &cTrack.x1         , "x1/F");
  tTreeTrack->Branch("y0"	  , &cTrack.y0         , "y0/F");
  tTreeTrack->Branch("y1"	  , &cTrack.y1         , "y1/F");
  tTreeTrack->Branch("z0"	  , &cTrack.z0         , "z0/F");
  tTreeTrack->Branch("z1"	  , &cTrack.z1         , "z1/F");
  tTreeTrack->Branch("nHits"	  , &cTrack.nHits      , "nHits/I");
  tTreeTrack->Branch("nGoodHits"  , &cTrack.nGoodHits  , "nGoodHits/I");
  tTreeTrack->Branch("score"	  , &cTrack.score      , "score/F");
  tTreeTrack->Branch("origin" 	  , &cTrack.origin     , "origin/I");
  tTreeTrack->Branch("pdg" 	  , &cTrack.pdg        , "pdg/I");
  tTreeTrack->Branch("energy" 	  , &cTrack.energy     , "energy/F");

  tTreeEvent = tfs->make<TTree>("CosmicEvent","CosmicEvent");  
  tTreeEvent->Branch("eventNumber", &cEvent.eventNumber, "eventNumber/I");
  tTreeEvent->Branch("nHitsTotal" , &cEvent.nHitsTotal , "nHitsTotal/I");
  tTreeEvent->Branch("nHitsEvID"  , &cEvent.nHitsEvID  , "nHitsEvID/I");
  tTreeEvent->Branch("nHitsTrack" , &cEvent.nHitsTrack , "nHitsTrack/I");
  tTreeEvent->Branch("nHitsTrackEvID" , &cEvent.nHitsTrackEvID , "nHitsTrackEvID/I");
  tTreeEvent->Branch("cOverA" , &cEvent.cOverA , "cOverA/F");
  tTreeEvent->Branch("dOverB" , &cEvent.dOverB , "dOverB/F");
  tTreeEvent->Branch("A" , &cEvent.A , "A/I");
  tTreeEvent->Branch("B" , &cEvent.B , "B/I");
  tTreeEvent->Branch("C" , &cEvent.C , "C/I");
  tTreeEvent->Branch("D" , &cEvent.D , "D/I");





}


// =====================================================
// Event Loop
// =====================================================
void microboone::CosmicRemovalAna::analyze(const art::Event& evt)
{

  
  
  // ===========================================================================================================
  // ============================================ LOOKING AT MCTRUTH ===========================================
  // ===========================================================================================================
  
  // ######################################
  // ### Picking up BackTrackser Service ###
  // ######################################
  art::ServiceHandle<cheat::BackTracker> bt;
  
  // ###################################################################
  // ### Defining a std::map of trackIDEs and vector of hit indicies ###
  // ###################################################################
  std::map<int, std::vector<int> > trkIDEsHitIndex_Cosmic;
  std::map<int, std::vector<int> > trkIDEsHitIndex_NonCosmic;

  std::map<int, float> TotalChargePerParticle_Cosmic;
  std::map<int, float> TotalChargePerParticle_NonCosmic;

  
  // ##################################################################
  // ### Grabbing ALL HITS in the event to monitor the backtracking ###
  // ##################################################################
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hitlist, hitListHandle);}
  
  
  
 
  
  int counter = 0;
  int newCounter = 0;


  int ntimethrough = 0;
  int nBadHits=0;
  int nTaggedHits = 0;
  // ###########################################
  // ### Looping over hits to get TrackIDE's ###
  // ###########################################
  for ( auto const& itr : hitlist )
    {
      newCounter++;
      std::vector<cheat::TrackIDE> eveIDs    = bt->HitToEveID(itr);
      
      if(eveIDs.size() == 0)
	{
	  nBadHits++;
	  continue;
	}
      // ############################
      // ### Loop over eveIDs's ###
      // ############################

      for (size_t e = 0; e<eveIDs.size(); e++)
	{
	  // ###################################################
	  // ### Retrieve mcTruth information for this eveID ###
	  // ###################################################
	  art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(eveIDs[e].trackID);
	  
	  
	  // Origin == 1 Not Cosmic
	  // Origin == 2 Cosmic
	  int origin = mctruth->Origin();
	  
	  // ########################################################
	  // ### If the Origin of this track is a Cosmic (i.e. 2) ###
	  // ### then fill a map keyed by MCParticleID and the 
	  // ### the second entry is a vector of hit indicies 
	  // ### associated with that particle
	  // ### Origin == 1 Not Cosmic
	  // ### Origin == 2 Cosmic
	  // ########################################################
	  if(origin == 2)//<---Cosmic Hit
	    {
	      // Fill the map keyed by the pointer and push back on 
	      // the index of the current hit
	      trkIDEsHitIndex_Cosmic[eveIDs[e].trackID].push_back(counter);
	      // Fill the map keyed by the pointer and add the 
	      // hit "charge" (ADC's) to the map
	      TotalChargePerParticle_Cosmic[eveIDs[e].trackID]+=itr->Charge();
	    }
	  else if( origin == 1 )//<--Not cosmic hit
	    {
	      trkIDEsHitIndex_NonCosmic[eveIDs[e].trackID].push_back(counter);
	      TotalChargePerParticle_NonCosmic[eveIDs[e].trackID]+=itr->Charge();
	    }

	  //<---
	  ntimethrough++;
	}//<---End e loop
      
      counter++;
    }//<--End nhits loop


  std::cerr << "counters: " << counter << " " << newCounter << std::endl;

  // ##########################################################
  // ### Once we are done with all the MC Particles, fill a ###
  // ### histogram that has the total truth charge for all  ###
  // ###    Cosmic and Non-cosmic particles in the event    ###
  // ##########################################################
  for(auto it = TotalChargePerParticle_Cosmic.begin(); it!=TotalChargePerParticle_Cosmic.end(); ++it)
    fTotalCharge_Cosmic->Fill(it->second);

  for(auto it = TotalChargePerParticle_NonCosmic.begin(); it!=TotalChargePerParticle_NonCosmic.end(); ++it)
    fTotalCharge_NonCosmic->Fill(it->second);
  
   
  // #################################################
  // ### Picking up track information on the event ###
  // #################################################
  art::Handle< std::vector<recob::Track> > trackh; //<---Track Handle
  evt.getByLabel(fTrackModuleLabel,trackh);
  std::vector<art::Ptr<recob::Track> > tracklist;//<---Ptr Vector
  art::fill_ptr_vector(tracklist,trackh);//<---Fill the vector
  

  
  // ################################################### 
  // ### Picking up cluster information on the event ###
  // ###################################################	
  art::Handle< std::vector<recob::Cluster> > clusterh; //<---Cluster Handle
  evt.getByLabel(fClusterModuleLabel, clusterh); 
  std::vector<art::Ptr<recob::Cluster> > clusterlist;//<---Ptr Vector
  art::fill_ptr_vector(clusterlist,clusterh);//<---Fill the vector



  
  std::map<int, int> NRejected60_Cosmic;
  std::map<int, int> NRejected60_NonCosmic;
  std::map<int, int> NRejected80_Cosmic;
  std::map<int, int> NRejected80_NonCosmic;
  std::map<int, int> NRejected95_Cosmic;
  std::map<int, int> NRejected95_NonCosmic;
  
  // #################################################
  // ### Looping over the collection of CosmicTags ###
  // #################################################

  art::Ptr<anab::CosmicTag> currentTag;
  //std::set<art::ProductID> TaggedHitIDs; // <-- I think this is describing an ID for how something is made and not something like a hit id
  std::set<art::Ptr<recob::Hit>::key_type> TaggedHitIDs;
  std::set<art::Ptr<recob::Hit>::key_type> TrackedHitIDs;
  std::set<art::Ptr<recob::Hit>::key_type> ClusteredHitIDs;




  std::vector <float> cosmicTagValues;
	
  int nTotalHitsTrack =0;  	
  int nTotalHitsTrackEvID =0;  	



  // Get the list of tracked hits

  art::FindManyP<recob::Hit> TrkHit( tracklist, evt, fTrackModuleLabel);
 
  for(unsigned int trk = 0; trk < trackh->size(); trk++)
    for(size_t i=0; i!=TrkHit.at(trk).size();++i)
      TrackedHitIDs.insert(TrkHit.at(trk).at(i).key());


  // Get the list of clustered hits

  art::FindManyP<recob::Hit> CluHit(clusterlist, evt, fClusterModuleLabel);

  for(unsigned int clu = 0; clu < clusterh->size(); clu++)
    for(size_t i=0; i!=CluHit.at(clu).size();++i)
      ClusteredHitIDs.insert(CluHit.at(clu).at(i).key());


  // These get filled later, but they should be filled ~here in future - Ben J
  double TrackedFractionCosmic = 0;
  double TrackedFractionNonCosmic = 0;
  double ClusteredFractionCosmic = 0;
  double ClusteredFractionNonCosmic = 0;



  /// -------------------------------------------------


  for(unsigned int nCT = 0; nCT < nCosmicTags; nCT++)//<---This loops over the vector of cosmicTags in stored in the event
    {
      

      TaggedHitIDs.clear();
      
      try{ //<---Putting in a try/catch in case no tags are found

	// ### Getting current cosmic tag associations ###      
	art::FindManyP<anab::CosmicTag> cosmicTrackTag( tracklist, evt, fCosmicTagAssocLabel[nCT]);
	
	// ============================================================================================================
	// ============================================== LOOKING AT TRACKS ===========================================
	// ============================================================================================================

	
	// ###########################
	// ### Looping over tracks ###
	// ###########################
	for(unsigned int trk = 0; trk < trackh->size(); trk++)
	  {

	    std::vector< art::Ptr< recob::Hit> > HitVec = TrkHit.at(trk);
	    
	    // ############################################################
	    // ### Start assuming the track is not associated w/ cosmic ###
	    // ############################################################
	    bool TrkMatchToCosmicTag = false;
	    
	    if(cosmicTrackTag.at(nCT).size()>1) 
	      std::cerr << "Warning : more than one cosmic tag per track in module " << fCosmicTagAssocLabel[nCT] << ". Confused, but just taking the first one anyway.";
	    
	    if(cosmicTrackTag.at(nCT).size()==0) continue;
	    
	    currentTag = cosmicTrackTag.at(trk).at(0);
	    float Score = currentTag->CosmicScore();

	    if(Score > fCosmicScoreThresholds[nCT] ) TrkMatchToCosmicTag = true;
	    
	    //	    scoreVector.push_back( Score );

	    fCosmicScoresPerCT[nCT]->Fill(Score);
	    fTPCCosmicType->Fill( currentTag->CosmicType() );

	    ///////////////////////////////////////
	    /// sel, adding some track information
	    if( currentTag->CosmicType() == 2 ) {
	      std::cerr << "found one that crossed two y || z boundaries" << std::endl;
	    }




	  // ############################################
	  // ###   If Track is matched to a CosmicTag ###
	  // ### for the TPC then find associated hits###
	  // ############################################
	    if(TrkMatchToCosmicTag)
	      {
		for(size_t i=0; i!=TrkHit.at(trk).size(); ++ i)
		  {
		    //TaggedHitIDs.insert(TrkHit.at(trk).at(i).id());
		    TaggedHitIDs.insert(TrkHit.at(trk).at(i).key());
		  }
		
	      	
	      }//<---TrkMatchToCosmicTag

	    // IF YOU'RE WONDERING WHY I CHANGED HOW THE OVERLAP WAS DONE, UNCOMMENT THE FOLLOWING LINES:
	    //std::cerr << "debugging... " << TrkMatchToCosmicTag << ", hit ids are ";
	    //for(size_t i=0; i!=TrkHit.at(trk).size(); ++ i) std::cerr << " " << TrkHit.at(trk).at(i).id() << " " << TrkHit.at(trk).at(i).key();
	    //std::cerr << std::endl;
	    FillTrackTree(currentTag, Score, evt);
	    
	    

	  }//<---end trk loop
      }
      catch(...)
	{
	  std::cout<<"No tracks for tag " << nCT<<", moving on"<<std::endl;
	}

    
      // ============================================================================================================
      // ============================================ LOOKING AT CLUSTERS ===========================================
      // ============================================================================================================
      
      // ####################################################
      // ### Getting the hits associated with the cluster ###
      // ####################################################

      try{

	evt.getByLabel( fClusterModuleLabel, clusterh ); 
	std::vector<art::Ptr<recob::Cluster> > clusterVector2; //<---Ptr Vector
	art::fill_ptr_vector(clusterVector2,clusterh);        //<---Fill the vector


	art::FindOneP<anab::CosmicTag> cosmicClusterTag( clusterVector2, evt, fCosmicTagAssocLabel[nCT]); 


	for(unsigned int clu = 0; clu < clusterh->size(); clu++)
	  {
	    // ##############################################################
	    // ### Start assuming the cluster is not associated w/ cosmic ###
	    // ##############################################################
	    bool CluMatchToCosmicTag = false;

	    std::vector< art::Ptr<recob::Hit> > thisHitVec = CluHit.at(clu);
	    
	    currentTag = cosmicClusterTag.at(clu);

	    float CosmicScore = currentTag->CosmicScore();
	    if(CosmicScore > fCosmicScoreThresholds[nCT]) CluMatchToCosmicTag = true;	
	    
	    
	    // #############################################
	    // ### If Cluster is matched to a Flash then ###
	    // ###    find the hits associated with it   ###
	    // #############################################
	    if(CluMatchToCosmicTag)
	      {
		for(size_t i=0; i!=CluHit.at(clu).size(); ++ i)
		  {
		    //		    TaggedHitIDs.insert(CluHit.at(clu).at(i).id());
		    TaggedHitIDs.insert(CluHit.at(clu).at(i).key());
		    nTaggedHits++;
		  }
	      }
	    
	    cTrack.eventNumber = -1;
	    cTrack.tagType     = -1;
	    cTrack.x0          = -1;	   
	    cTrack.x1          = -1;
	    cTrack.y0          = -1;
	    cTrack.y1          = -1;
	    cTrack.z0          = -1;
	    cTrack.z1          = -1;
	    cTrack.nHits       = -1;
	    cTrack.nGoodHits   = -1;
	    cTrack.origin      = -1;
	    cTrack.score       = -1;
	    cTrack.pdg         = -999;
	    cTrack.energy      = -999;

	    cTrack.eventNumber = evt.event();
	    cTrack.tagType     = currentTag->CosmicType();
	    cTrack.x0          = currentTag->endPt1[0];
	    cTrack.x1          = currentTag->endPt2[0];
	    cTrack.y0          = currentTag->endPt1[1];
	    cTrack.y1          = currentTag->endPt2[1];
	    cTrack.z0          = currentTag->endPt1[2];
	    cTrack.z1          = currentTag->endPt2[2];
	    cTrack.nHits       = -1;//anotherCounter; 
	    cTrack.nGoodHits   = -1;//numberGoodEvID; 
//	    int count1s = count( tempOrigin.begin(), tempOrigin.end(), 1);
//	    int count2s = count( tempOrigin.begin(), tempOrigin.end(), 2);
//	    std::cerr << "----- COUNTS: nu: " << count1s << " cosmic: " << count2s << " position: " <<  currentTag->endPt1[0]<< " " << currentTag->endPt2[0] << " " ;
//	    int majorityOrigin = -1;
//	    if(tempOrigin.size()>0) majorityOrigin = count1s > count2s ? 1 : 2;
//	    std::cerr << "majority origin " << majorityOrigin << std::endl;
	    int pdg =-999;
	    float energy =-999;
	    float majorityOrigin = GetOrigin( thisHitVec, evt, pdg, energy);
	    cTrack.origin   = majorityOrigin;	    
	    cTrack.score = CosmicScore;
	    cTrack.pdg = pdg;
	    cTrack.energy = energy;
	    tTreeTrack->Fill();


	  }//<---End clu loop
      }
      catch(...)
	{
	  std::cout<<"No clusters for tag " << nCT << ", moving on"<<std::endl;
	}
    	
//       for(auto itCosmicParticle = trkIDEsHitIndex_Cosmic.begin(); itCosmicParticle!=trkIDEsHitIndex_Cosmic.end(); ++itCosmicParticle)
// 	{

// 	  //	  double OverlapScore = GetOverlapScore(hitlist, itCosmicParticle->second, TaggedHitIDs);
// 	  if(OverlapScore > 0.60) NRejected60_Cosmic[itCosmicParticle->first]++;
// 	  if(OverlapScore > 0.80) NRejected80_Cosmic[itCosmicParticle->first]++;
// 	  if(OverlapScore > 0.95) NRejected95_Cosmic[itCosmicParticle->first]++;
// 	  //fFractionChargeTaggedPerCT_Cosmic[nCT]->Fill(OverlapScore);
// 	}
      
//       for(auto itNonCosmicParticle = trkIDEsHitIndex_NonCosmic.begin(); itNonCosmicParticle!=trkIDEsHitIndex_NonCosmic.end(); ++itNonCosmicParticle)
// 	{

// 	  //double OverlapScore = GetOverlapScore(hitlist, itNonCosmicParticle->second, TaggedHitIDs);
// 	  if(OverlapScore > 0.60) NRejected60_NonCosmic[itNonCosmicParticle->first]++;
// 	  if(OverlapScore > 0.80) NRejected80_NonCosmic[itNonCosmicParticle->first]++;
// 	  if(OverlapScore > 0.95) NRejected95_NonCosmic[itNonCosmicParticle->first]++;
// 	  //fFractionChargeTaggedPerCT_NonCosmic[nCT]->Fill(OverlapScore);
// 	}



      
      // this needs moving outside the nCT loop, along with the truth stuff - Ben J

      if(nCT==0)
	{
	  TrackedFractionCosmic = GetOverlapScore(hitlist, trkIDEsHitIndex_Cosmic, TrackedHitIDs);
	  fTrackedFraction_Cosmic->Fill(TrackedFractionCosmic);

	  TrackedFractionNonCosmic = GetOverlapScore(hitlist, trkIDEsHitIndex_NonCosmic, TrackedHitIDs);

	  fTrackedFraction_NonCosmic->Fill(TrackedFractionNonCosmic);

	  ClusteredFractionCosmic = GetOverlapScore(hitlist, trkIDEsHitIndex_Cosmic, ClusteredHitIDs);

	  fClusteredFraction_Cosmic->Fill(ClusteredFractionCosmic);
	  
	  ClusteredFractionNonCosmic = GetOverlapScore(hitlist, trkIDEsHitIndex_NonCosmic, ClusteredHitIDs);

	  fClusteredFraction_NonCosmic->Fill(ClusteredFractionNonCosmic);
	}

      // Checking on an Efficiency Measure
      double OverlapScore = GetOverlapScore(hitlist, trkIDEsHitIndex_Cosmic, TaggedHitIDs);
      fFractionChargeTaggedPerCT_Cosmic[nCT]->Fill(OverlapScore);

      // Checking on a Purity Measure
      OverlapScore = GetOverlapScore(hitlist, trkIDEsHitIndex_NonCosmic, TaggedHitIDs);
      fFractionChargeTaggedPerCT_NonCosmic[nCT]->Fill(OverlapScore);

      std::cerr << "Report... Total hits in hit collection: " << hitlist.size() << std::endl;
      std::cerr << "Number of hits without an eveIDE: " << nBadHits << std::endl;
      int cosmicCounts=0;
      for( auto& x: trkIDEsHitIndex_Cosmic ) for(unsigned int i=0;i<x.second.size(); i++) cosmicCounts++;
      std::cerr << "Total truth cosmic hits: " << cosmicCounts << std::endl;
      int nuCounts=0;
      for( auto& x: trkIDEsHitIndex_NonCosmic ) for(unsigned int i=0;i<x.second.size(); i++) nuCounts++;
      std::cerr << "Total truth nu hits: " << nuCounts << std::endl;
      std::cerr << "Tagged Hits: " << TaggedHitIDs.size() << " " << nTaggedHits << std::endl;



      cEvent.eventNumber = -1;
      cEvent.nHitsTotal  = -1;
      cEvent.nHitsEvID   = -1;
      cEvent.nHitsTrack  = -1;  
      cEvent.nHitsTrackEvID = -1;
      cEvent.cOverA = -1;      
      cEvent.dOverB = -1;
      cEvent.A = -1;
      cEvent.B = -1;
      cEvent.C = -1;
      cEvent.D = -1;

      float cOverA = -1;
      float dOverB = -1;
//       GetNewOverlapScore(scoreVector, nHitsVector, originVector, cOverA, dOverB);
//       if(dOverB >0.5) std::cerr << "Large D/B: " << dOverB << " in Event: " << evt.event() << std::endl; 
//       fFractionChargeTaggedPerCT_Cosmic[nCT]->Fill( cOverA );
//       fFractionChargeTaggedPerCT_NonCosmic[nCT]->Fill( dOverB );


      cEvent.eventNumber = evt.event();
      cEvent.nHitsTotal = hitlist.size();
      cEvent.nHitsEvID  = counter;
      cEvent.nHitsTrack = nTotalHitsTrack;
      cEvent.nHitsTrackEvID = nTotalHitsTrackEvID;
      cEvent.cOverA = cOverA;      
      cEvent.dOverB = dOverB;

      tTreeEvent->Fill();

    }//<---End nCT (number of cosmic tag) loop


  for(auto it =NRejected60_Cosmic.begin(); it!=NRejected60_Cosmic.end(); ++it)
    fNAlgsRejected60_Cosmic->Fill(it->second);
  for(auto it =NRejected60_NonCosmic.begin(); it!=NRejected60_NonCosmic.end(); ++it)
    fNAlgsRejected60_NonCosmic->Fill(it->second);
  for(auto it =NRejected80_Cosmic.begin(); it!=NRejected80_Cosmic.end(); ++it)
    fNAlgsRejected80_Cosmic->Fill(it->second);
  for(auto it =NRejected80_NonCosmic.begin(); it!=NRejected80_NonCosmic.end(); ++it)
    fNAlgsRejected80_NonCosmic->Fill(it->second);
  for(auto it =NRejected95_Cosmic.begin(); it!=NRejected95_Cosmic.end(); ++it)
    fNAlgsRejected95_Cosmic->Fill(it->second);
  for(auto it =NRejected95_NonCosmic.begin(); it!=NRejected95_NonCosmic.end(); ++it)
    fNAlgsRejected95_NonCosmic->Fill(it->second);


  
}



//-------------------------------------------------------


//double microboone::CosmicRemovalAna::GetOverlapScore(std::vector<art::Ptr<recob::Hit> >& AllHits, std::vector<int>& HitsThisParticle, std::set<art::ProductID>& TaggedHits)
double microboone::CosmicRemovalAna::GetOverlapScore(std::vector<art::Ptr<recob::Hit> >& AllHits, 
						     std::map<int,std::vector<int> >& HitsThisParticle, 
						     std::set<art::Ptr<recob::Hit>::key_type>& TaggedHits)
{
  int CountAll    = 0;
  int CountTagged = 0;

  for( auto& x: HitsThisParticle ) {
    CountAll += (x.second).size();

    for( unsigned int j = 0; j < x.second.size(); j++ ) {
      int index = (x.second).at(j);
      if( TaggedHits.count(AllHits.at(index).key() ) > 0 ) CountTagged++;
    }

  }

  std::cerr << std::endl;
  std::cerr << "GetOverlapScore: count tagged: " << CountTagged << " CountAll: " << CountAll << " TaggedHits.size(): " << TaggedHits.size() << std::endl;

  return float(CountTagged) / float(CountAll);
}

//-------------------------------------------------------
int microboone::CosmicRemovalAna::GetOrigin( std::vector<art::Ptr<recob::Hit> >& hitlist, const art::Event& evt, int &pdg, float &en) {

  art::ServiceHandle<cheat::BackTracker> bt;

  //float cosmicOrg = 0;
  //float nuOrg = 0;
  int org = -1;

  for ( auto const& itr : hitlist ) {
    std::vector<cheat::TrackIDE> eveIDs    = bt->HitToEveID(itr);
      if(eveIDs.size() == 0) {
	//	nBadHits++;
	continue;
      }


      // MCPARTICLE IS BARFING RIGHT NOW
      // UNCOMMENT THE FOLLOWING LINES AFTER A NEW RELEASE IS CUT
      // THERE'S SOME SORT OF INCONSISTENCY RIGHT NOW WITH NUTOOLS
      /*
      double maxE = -999;
      size_t indx = -1;
      for (size_t e = 0; e<eveIDs.size(); e++) {
	const simb::MCParticle* mcparticle = bt->TrackIDToParticle( eveIDs[e].trackID );
	double energy = mcparticle->E();
	if( energy > maxE ) {
	  maxE = energy;
	  indx = e;
	}
      }

      const simb::MCParticle* mcparticle = bt->TrackIDToParticle( eveIDs[indx].trackID );
      pdg = mcparticle->PdgCode();
      en = (float)maxE;


      art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(eveIDs[indx].trackID);
      org = mctruth->Origin();
*/


//      for (size_t e = 0; e<eveIDs.size(); e++) {
//	art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(eveIDs[e].trackID);
//	const simb::MCParticle* mcparticle = bt->TrackIDToParticle( eveIDs[e].trackID );
//	// Origin == 1 Not Cosmic
//	// Origin == 2 Cosmic
//	int origin = mctruth->Origin();
//	if( origin==1 ) nuOrg += mcparticle->E();
//	if( origin==2 ) cosmicOrg += mcparticle->E();
//      }


  } // loop over hits
  
  //  return nuOrg > cosmicOrg ? 1 : 2;
  return org;

}

//-------------------------------------------------------
void microboone::CosmicRemovalAna::GetNewOverlapScore( std::vector< float > &scoreVector, 
							std::vector< int > &hitsVector,
						       std::vector< std::vector<int> > &originVector,
						       float &cOverA, float &dOverB ) {

  int countA = 0;
  int countB = 0;
  int countC = 0;
  int countD = 0;

  int prevHitVal=0;


  for( size_t s = 0; s < scoreVector.size(); s++ ) {
    int nHits = hitsVector.at(s);

    float score = scoreVector.at(s);

    // each track has nHits associated with it.  These are the hits we want to look at.
    for(int h = prevHitVal; h < prevHitVal+nHits; h++ ) {

      int cosmicOrg = std::count(originVector.at(h).begin(), originVector.at(h).end(),2); 
      int nuOrg = std::count(originVector.at(h).begin(), originVector.at(h).end(),1); 
      countA += cosmicOrg;
      countB += nuOrg;
      if(score>0.5) {
	countC += cosmicOrg;
	countD += nuOrg;
      }
    } // loop over origin vector
    prevHitVal = prevHitVal+nHits;
  }// loop over scoreVector
  
  std::cerr << "A: " << countA <<  " B: " << countB << " C: " << countC << " D: " << countD << std::endl;
  std::cerr << "-------------------------------     C/A: " << countC*1.0/countA << " D/B: " << countD*1.0/countB << std::endl;

  if( countA < 1 ) countA = 1;
  if( countB < 1 ) countB = 1;
  cOverA = countC*1.0/countA;
  dOverB = countD*1.0/countB;

  cEvent.A = countA;
  cEvent.B = countB;
  cEvent.C = countC;
  cEvent.D = countD;

}


void microboone::CosmicRemovalAna::FillTrackTree(art::Ptr<anab::CosmicTag> &currentTag, float CosmicScore, const art::Event& evt)
{
    cTrack.eventNumber = -1;
    cTrack.tagType     = -1;
    cTrack.x0          = -1;	   
    cTrack.x1          = -1;
    cTrack.y0          = -1;
    cTrack.y1          = -1;
    cTrack.z0          = -1;
    cTrack.z1          = -1;
    cTrack.nHits       = -1;
    cTrack.nGoodHits   = -1;
    cTrack.origin      = -1;
    cTrack.score       = -1;
    cTrack.pdg         = -999;
    cTrack.energy      = -999;

	    
    cTrack.eventNumber = evt.event();
    cTrack.tagType     = currentTag->CosmicType();
    cTrack.x0          = currentTag->endPt1[0];
    cTrack.x1          = currentTag->endPt2[0];
    cTrack.y0          = currentTag->endPt1[1];
    cTrack.y1          = currentTag->endPt2[1];
    cTrack.z0          = currentTag->endPt1[2];
    cTrack.z1          = currentTag->endPt2[2];
    cTrack.nHits       = -1;//anotherCounter; 
    cTrack.nGoodHits   = -1;//numberGoodEvID; 
//	    int count1s = count( tempOrigin.begin(), tempOrigin.end(), 1);
//	    int count2s = count( tempOrigin.begin(), tempOrigin.end(), 2);
//	    std::cerr << "----- COUNTS: nu: " << count1s << " cosmic: " << count2s << " position: " <<  currentTag->endPt1[0]<< " " << currentTag->endPt2[0] << " " ;
//	    int majorityOrigin = -1;
//	    if(tempOrigin.size()>0) majorityOrigin = count1s > count2s ? 1 : 2;
//	    std::cerr << "majority origin " << majorityOrigin << std::endl;
    cTrack.origin   = -1;//majorityOrigin;	    
    cTrack.score = CosmicScore;
    cTrack.pdg    = -999;
    cTrack.energy = -999;
    tTreeTrack->Fill();




}
//-------------------------------------------------------

namespace microboone{
  
  DEFINE_ART_MODULE(CosmicRemovalAna)  
}








#endif
