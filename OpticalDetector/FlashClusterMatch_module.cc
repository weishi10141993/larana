// Some kinda description here, maybe
//
// It does optical stuff.



#include "art/Framework/Core/EDProducer.h"
#include "AnalysisBase/FlashMatch.h"


// ROOT includes.
#include <Rtypes.h>
#ifndef FlashClusterMatch_h
#define FlashClusterMatch_h 1



namespace recob{
  class OpFlash;
  class Cluster;

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
    std::string fClusterModuleLabel;
    std::string fFlashModuleLabel;
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
#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/OpFlash.h"

// FMWK includes
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Utilities/LArProperties.h"

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>



namespace opdet {

  //-------------------------------------------------

  FlashClusterMatch::FlashClusterMatch(fhicl::ParameterSet const& pset)
  {
    produces< std::vector<anab::FlashMatch> >();
    produces< art::Assns<recob::Cluster, anab::FlashMatch> >();

    this->reconfigure(pset);
   }


  //-------------------------------------------------

  void FlashClusterMatch::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");   
    fFlashModuleLabel = pset.get<std::string>("FlashModuleLabel");

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

    // DO NOT REMOVE THIS LINE:
    mf::LogWarning("RecoBaseDefaultCtor") << "using default Hit ctor - should only ever"
      					  << " be done when getting hits out of an event"
      					  << " not when trying to produce new hits to store"
      					  << " in the event";

    
    
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

    std::unique_ptr< std::vector<anab::FlashMatch> > flash_matches ( new std::vector<anab::FlashMatch>);
    std::unique_ptr< art::Assns<recob::Cluster, anab::FlashMatch > > assn_cluster( new art::Assns<recob::Cluster, anab::FlashMatch>);


    evt.put(std::move(flash_matches));
    evt.put(std::move(assn_cluster));

    
  }


}


