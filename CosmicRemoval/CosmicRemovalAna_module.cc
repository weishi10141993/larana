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

    double GetOverlapScore(std::vector<art::Ptr<recob::Hit> >& AllHits, std::vector<int>& HitsThisParticle, std::set<art::ProductID>&  TaggedHits);


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
      fCosmicScoresPerCT.push_back( (TH1D*)tfs->make<TH1D>(sname.str().c_str(), stitle.str().c_str(), 100,0,1));

      sname.str("");  sname.flush();
      stitle.str(""); stitle.flush();
				    
      sname<<"FractionTaggedCosmicFor"<<fCosmicTagAssocLabel.at(i);
      stitle<<"Fraction of Cosmic Charge Tagged as Cosmic For " << fCosmicTagAssocLabel.at(i)<<"; Frac; N";
      fFractionChargeTaggedPerCT_Cosmic.push_back( (TH1D*)tfs->make<TH1D>(sname.str().c_str(), stitle.str().c_str(), 100,0,1));	  
      
      sname.str("");  sname.flush();
      stitle.str(""); stitle.flush();
      
      sname<<"FractionTaggedNonCosmicFor"<<fCosmicTagAssocLabel.at(i);
      stitle<<"Fraction of NonCosmic Charge Tagged as Cosmic For " << fCosmicTagAssocLabel.at(i)<<"; Frac; N";
      fFractionChargeTaggedPerCT_NonCosmic.push_back( (TH1D*)tfs->make<TH1D>(sname.str().c_str(), stitle.str().c_str(), 100,0,1));												   
    }
  

  
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
  
  
  
 
  int nhits = hitlist.size();
  
  int counter = 0;
  
  //std::cout<<"nhits = "<<nhits<<std::endl;
  //std::cout<<std::endl;
  int ntimethrough = 0;
   // ###########################################
  // ### Looping over hits to get TrackIDE's ###
  // ###########################################
  for ( auto const& itr : hitlist )
    {
      std::vector<cheat::TrackIDE> eveIDs    = bt->HitToEveID(itr);
      
      if(eveIDs.size() == 0)
	{
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
	  else//<--Not cosmic hit
	    {
	      trkIDEsHitIndex_NonCosmic[eveIDs[e].trackID].push_back(counter);
	      TotalChargePerParticle_NonCosmic[eveIDs[e].trackID]+=itr->Charge();
	    }


	  //<---
	  ntimethrough++;
	}//<---End e loop
      
      counter++;
    }//<--End nhits loop

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
  std::set<art::ProductID> TaggedHitIDs;
	
  for(unsigned int nCT = 0; nCT < nCosmicTags; nCT++)//<---This loops over the vector of cosmicTags in stored in the event
    {
      
      try{ //<---Putting in a try/catch in case no tags are found
	
	// ### Getting current cosmic tag associations ###      
	art::FindManyP<anab::CosmicTag> cosmicTrackTag( tracklist, evt, fCosmicTagAssocLabel[nCT]);
	
	// ============================================================================================================
	// ============================================== LOOKING AT TRACKS ===========================================
	// ============================================================================================================
	// ### Getting hits associatied with tracks ###	
	art::FindManyP<recob::Hit> TrkHit( tracklist, evt, fHitsModuleLabel);
	
	// We will fill this set with the hit IDs of cosmic tagged hits
	
	// ###########################
	// ### Looping over tracks ###
	// ###########################
	for(unsigned int trk = 0; trk < trackh->size(); trk++)
	  {
	    // ############################################################
	    // ### Start assuming the track is not associated w/ cosmic ###
	    // ############################################################
	    bool TrkMatchToCosmicTag = false;
	    
	    if(cosmicTrackTag.at(nCT).size()>1) 
	      std::cerr << "Warning : more than one cosmic tag per track in module " << fCosmicTagAssocLabel[nCT] << ". Confused, but just taking the first one anyway.";
	    
	    if(cosmicTrackTag.at(nCT).size()==0) continue;
	    
	    currentTag = cosmicTrackTag.at(trk).at(0);
	    float Score = currentTag->CosmicScore();
	    std::cerr << "Score = " << Score << std::endl;
	    
	    fCosmicScoresPerCT[nCT]->Fill(Score);
	    
	  // ############################################
	  // ###   If Track is matched to a CosmicTag ###
	  // ### for the TPC then find associated hits###
	  // ############################################
	    if(TrkMatchToCosmicTag)
	      {
		for(size_t i=0; i!=TrkHit.at(trk).size(); ++ i)
		  {
		    TaggedHitIDs.insert(TrkHit.at(trk).at(i).id());
		  }
		
	      	
	      }//<---TrkMatchToCosmicTag
	    
	    
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

	art::FindManyP<recob::Hit> CluHit(clusterlist, evt, fHitsModuleLabel);
	
	art::FindManyP<anab::CosmicTag> cosmicClusterTag( clusterlist, evt, fCosmicTagAssocLabel[nCT]);
	
	for(unsigned int clu = 0; clu < clusterh->size(); clu++)
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
		for(size_t i=0; i!=CluHit.at(clu).size(); ++ i)
		  {
		    TaggedHitIDs.insert(CluHit.at(clu).at(i).id());
		  }
	      }
	    
	  }//<---End clu loop
      }
      catch(...)
	{
	  std::cout<<"No clusters for tag " << nCT << ", moving on"<<std::endl;
	}
	
      
      for(auto itCosmicParticle = trkIDEsHitIndex_Cosmic.begin(); itCosmicParticle!=trkIDEsHitIndex_Cosmic.end(); ++itCosmicParticle)
	{
	  double OverlapScore = GetOverlapScore(hitlist, itCosmicParticle->second, TaggedHitIDs);
	  if(OverlapScore > 0.60) NRejected60_Cosmic[itCosmicParticle->first]++;
	  if(OverlapScore > 0.80) NRejected80_Cosmic[itCosmicParticle->first]++;
	  if(OverlapScore > 0.95) NRejected95_Cosmic[itCosmicParticle->first]++;
	  fFractionChargeTaggedPerCT_Cosmic[nCT]->Fill(OverlapScore);
	}
      
      for(auto itNonCosmicParticle = trkIDEsHitIndex_NonCosmic.begin(); itNonCosmicParticle!=trkIDEsHitIndex_NonCosmic.end(); ++itNonCosmicParticle)
	{
	  double OverlapScore = GetOverlapScore(hitlist, itNonCosmicParticle->second, TaggedHitIDs);
	  if(OverlapScore > 0.60) NRejected60_NonCosmic[itNonCosmicParticle->first]++;
	  if(OverlapScore > 0.80) NRejected80_NonCosmic[itNonCosmicParticle->first]++;
	  if(OverlapScore > 0.95) NRejected95_NonCosmic[itNonCosmicParticle->first]++;
	  fFractionChargeTaggedPerCT_NonCosmic[nCT]->Fill(OverlapScore);
	}
  
      
      


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


double microboone::CosmicRemovalAna::GetOverlapScore(std::vector<art::Ptr<recob::Hit> >& AllHits, std::vector<int>& HitsThisParticle, std::set<art::ProductID>& TaggedHits)
{
  int CountAll = HitsThisParticle.size();
  int CountTagged = 0 ;
  for(size_t i=0; i!=HitsThisParticle.size(); ++i)
    {
      if( TaggedHits.count(AllHits.at(i).id() ) == 1 )
	{
	  CountTagged++;
	}
    }
  return float(CountTagged) / float(CountAll);
}


//-------------------------------------------------------


namespace microboone{
  
  DEFINE_ART_MODULE(CosmicRemovalAna)  
}








#endif
