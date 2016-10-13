// Some kinda description here, maybe
//
// It does optical stuff.



#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// ROOT includes.
#include <Rtypes.h>
#ifndef FlashClusterMatch_h
#define FlashClusterMatch_h 1


namespace recob{
  class OpFlash;
  class Cluster;
  class Hit;
}


namespace opdet {
  

  class FlashClusterMatch : public art::EDProducer{
  public:
    
    FlashClusterMatch(const fhicl::ParameterSet&);
    virtual ~FlashClusterMatch();
    
    void produce(art::Event&);
    void reconfigure(fhicl::ParameterSet const& p);
      
    
    void beginJob();
    
    
  private:
    
    std::vector<double>  GetLightHypothesis(std::vector<recob::SpacePoint> spts);
    bool                 CheckCompatibility(std::vector<double>& hypothesis, std::vector<double>& signal);

    
    trkf::SpacePointAlg       *  fSptalg;
    calo::CalorimetryAlg      *  fCaloAlg;

    std::string fClusterModuleLabel;
    std::string fFlashModuleLabel;
    int         fMinSptsForOverlap;
    double      fSingleChannelCut;
    double      fIntegralCut;
  };

}

#endif




////////////////////////////////////////////////////////////////////////
/// \file  FlashClusterMatch_module.cc
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(FlashClusterMatch)

}//end namespace opdet
////////////////////////////////////////////////////////////////////////


// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"


// FMWK includes
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>



namespace opdet {

  //-------------------------------------------------
  
  FlashClusterMatch::FlashClusterMatch(fhicl::ParameterSet const& pset)
  {
    produces< std::vector<anab::CosmicTag> >();
    produces< art::Assns<recob::Cluster, anab::CosmicTag> >();

    this->reconfigure(pset);
   }


  //-------------------------------------------------

  void FlashClusterMatch::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");   
    fFlashModuleLabel   = pset.get<std::string>("FlashModuleLabel");
    fMinSptsForOverlap  = pset.get<int>("MinSptsForOverlap");

    fSingleChannelCut   = pset.get<int>("SingleChannelCut");
    fIntegralCut        = pset.get<int>("IntegralCut");

    fSptalg             = new trkf::SpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
    fCaloAlg            = new calo::CalorimetryAlg(pset.get< fhicl::ParameterSet >("CaloAlg"));
    
  }


  //-------------------------------------------------

  void FlashClusterMatch::beginJob()
  {
  }



  //-------------------------------------------------

  FlashClusterMatch::~FlashClusterMatch()
  {
  }





  //-------------------------------------------------


  void FlashClusterMatch::produce(art::Event& evt)
  {

    int n=0;


  top:
    n++;

    // DO NOT REMOVE THIS LINE:
    mf::LogWarning("RecoBaseDefaultCtor") << "using default Hit ctor - should only ever"
      					  << " be done when getting hits out of an event"
      					  << " not when trying to produce new hits to store"
      					  << " in the event";

    std::cerr<< " Warning : you have disabled the RecoBaseDefaultCtor message."  ;
    std::cerr<< "  Should only ever be done when trying to avoid messages when getting hits out of an event, not when trying to produce new hits to store in the event."<< std::endl;
    
    ++n;
    if(n<10) goto top;


    
    std::unique_ptr< std::vector<anab::CosmicTag> > cosmic_tags ( new std::vector<anab::CosmicTag>);
    std::unique_ptr< art::Assns<recob::Cluster, anab::CosmicTag > > assn_tag( new art::Assns<recob::Cluster, anab::CosmicTag>);



    
    // Read in flashes from the event
    art::Handle< std::vector<recob::OpFlash> > flashh;
    evt.getByLabel(fFlashModuleLabel, flashh);
    std::vector<art::Ptr<recob::OpFlash> > Flashes;
    for(unsigned int i=0; i < flashh->size(); ++i)
      {
	art::Ptr<recob::OpFlash> flash(flashh,i);
        Flashes.push_back(flash);
      }

    // Read in clusters from the event
    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);
    std::vector<art::Ptr<recob::Cluster> >  Clusters;
    for(unsigned int i=0; i < clusterh->size(); ++i)
      {
	art::Ptr<recob::Cluster> cluster(clusterh,i);
	Clusters.push_back(cluster);
      }
    
    // Pull associated hits from event
    art::FindManyP<recob::Hit> hits(clusterh, evt, fClusterModuleLabel);


    // Extract flash shape info
    std::vector<std::vector<double> > FlashShapes;
    art::ServiceHandle<geo::Geometry> geom;
    size_t NOpDets = geom->NOpDets();

    for(size_t f=0; f!=Flashes.size(); ++f)
      {
	if(Flashes.at(f)->OnBeamTime())
          {
	    std::vector<double> ThisFlashShape(NOpDets,0);
            for (unsigned int c = 0; c < geom->NOpChannels(); c++){
              unsigned int o = geom->OpDetFromOpChannel(c);
	      ThisFlashShape[o]+=Flashes.at(f)->PE(c);
            }
            FlashShapes.push_back(ThisFlashShape);
          }
      }



    std::vector<std::vector<int> > SortedByViews(3);

    //    std::vector<int> MaxWire(Clusters.size(), 0);
    //    std::vector<int> MinWire(Clusters.size(), 10000);

    std::vector<int> MaxTime(Clusters.size(), 0);
    std::vector<int> MinTime(Clusters.size(), 10000);


    // Sort clusters by view
    for(size_t iClus=0; iClus!=Clusters.size(); ++iClus)
      {
	SortedByViews[Clusters.at(iClus)->View()].push_back(iClus);
	for(size_t iHit=0; iHit!=hits.at(iClus).size(); ++iHit)
	  {
	    double Time = hits.at(iClus).at(iHit)->PeakTime();
	    if(Time > MaxTime[iClus]) MaxTime[iClus] = Time;
	    if(Time < MinTime[iClus]) MinTime[iClus] = Time;
	    
	    //  Equivalent info for wires, maybe we want it later.
	    //	    int Wire = hits.at(iClus).at(iHit)->WireID().Wire;	 
	    //	    if(Wire < MinWire[iClus]) MinWire[iClus] = Wire;
	    //	    if(Wire > MaxWire[iClus]) MaxWire[iClus] = Wire;

	  }
      }

    std::vector<std::vector<double> > hypotheses;

    // Loop over sets of 3 clusters
    for(size_t nU=0; nU!=SortedByViews[0].size(); ++nU)
      for(size_t nV=0; nV!=SortedByViews[1].size(); ++nV)
	for(size_t nW=0; nW!=SortedByViews[2].size(); ++nW)
	  {
	    int indexU = SortedByViews[0][nU];
	    int indexV = SortedByViews[1][nV];
	    int indexW = SortedByViews[2][nW];

	    bool NoOverlap = false;
	    
	    // Skip over clusters with no time overlap
	    for(size_t v=0; v!=3; ++v)
	      {
		int v1 = (v+1)%3;
		int v2 = (v+2)%3;
		
		if(MinTime[v] > std::min(MaxTime[v1],MaxTime[v2]))
		  NoOverlap = true;

		if(MaxTime[v] < std::max(MinTime[v1],MinTime[v2]))
		  NoOverlap = true;
	      }
	    
	    if(NoOverlap) continue;

	    // Prepare flattened vector for space pointery
	    art::PtrVector<recob::Hit>  FlatHits;

	    FlatHits.insert(FlatHits.begin(), hits.at(indexU).begin(), hits.at(indexU).end());
	    FlatHits.insert(FlatHits.begin(), hits.at(indexV).begin(), hits.at(indexV).end());
	    FlatHits.insert(FlatHits.begin(), hits.at(indexW).begin(), hits.at(indexW).end());

	    // Make the spacepoints
	    std::vector<recob::SpacePoint> spts;
	    fSptalg->makeSpacePoints(FlatHits, spts);

	    if(int(spts.size()) < fMinSptsForOverlap) continue;

	    // Get light hypothesis for this collection
	    std::vector<double> hypothesis = GetLightHypothesis(spts);

	    bool IsCompatible = false;
	    
	    // Check for each flash, whether this subevent is compatible
	    for(size_t jFlash=0; jFlash!=FlashShapes.size(); ++jFlash)
	      {
		// It is compatible with the beam flash.
		if(CheckCompatibility(hypothesis,FlashShapes.at(jFlash)))
		  {
		    IsCompatible=true;
		  }	      
	      }

	    // If not compatible with any beam flash, throw out
	    if(!IsCompatible)
	      {
		cosmic_tags->push_back(anab::CosmicTag(1.));
		util::CreateAssn(*this, evt, *(cosmic_tags.get()), Clusters.at(indexU), *(assn_tag.get()), cosmic_tags->size()-1);
		util::CreateAssn(*this, evt, *(cosmic_tags.get()), Clusters.at(indexV), *(assn_tag.get()), cosmic_tags->size()-1);
		util::CreateAssn(*this, evt, *(cosmic_tags.get()), Clusters.at(indexW), *(assn_tag.get()), cosmic_tags->size()-1);
	      }
	    
	  }
    

    

    evt.put(std::move(cosmic_tags));
    evt.put(std::move(assn_tag));
  }



  //---------------------------------------------------------


  // Get a hypothesis for the light from a spacepoint collection
  std::vector<double> FlashClusterMatch::GetLightHypothesis(std::vector<recob::SpacePoint> spts)
  {
    art::ServiceHandle<geo::Geometry> geom;
    std::vector<double> ReturnVector(geom->NOpDets(),0);

    art::ServiceHandle<phot::PhotonVisibilityService> pvs;


    for (size_t s=0; s!=spts.size(); s++)
      {
	const art::PtrVector<recob::Hit>& assochits = fSptalg->getAssociatedHits(spts.at(s));
	
	double Charge     = 0;
	double WirePitch  = 0.3;
	
	for(size_t iHit=0; iHit!=assochits.size(); ++iHit)
	  if(assochits.at(iHit)->View()==2) Charge += WirePitch * fCaloAlg->dEdx_AMP(assochits.at(iHit), 1);
	
	
      	double xyz[3];
	
	for(size_t i=0; i!=3; ++i) xyz[i] = spts.at(s).XYZ()[i];
	
        const float* PointVisibility = pvs->GetAllVisibilities(xyz);
	if (!PointVisibility) continue; // point not covered by the service
        for(size_t OpDet =0; OpDet!=pvs->NOpChannels();  OpDet++)
          {
            ReturnVector.at(OpDet)+= PointVisibility[OpDet];
          }
      }
    double PhotonYield = 24000;
    double QE          = 0.01;
    
    for(size_t i=0; i!=ReturnVector.size(); ++i)
      {
	ReturnVector[i] *= QE * PhotonYield;
      }

    return ReturnVector;
  }


  //----------------------------------------------
  bool FlashClusterMatch::CheckCompatibility(std::vector<double>& hypothesis, std::vector<double>& signal)
  {
    double sigintegral=0, hypintegral=0;
    for(size_t i=0; i!=hypothesis.size(); ++i)
      {
	sigintegral+=signal.at(i);
	hypintegral+=hypothesis.at(i);
	double HypErr = pow(hypothesis.at(i),0.5);
	if(( (hypothesis.at(i) - signal.at(i)) / HypErr) > fSingleChannelCut) return false;
      }
    double HypIntErr= pow(hypintegral,0.5);

    if( ( (hypintegral - sigintegral)/HypIntErr) > fIntegralCut) return false;
    return true;
  }




}


