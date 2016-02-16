// -*- mode: c++; c-basic-offset: 2; -*-
// Ben Jones, MIT, 2013
//
// This module finds periods of time-localized activity
// from the optical system, called Flashes, using OpHits as an input.
//
// Modified to make it more detector agnostic
// by Gleb Sinev, Duke, 2015
//


#ifndef OpFlashFinder_H
#define OpFlashFinder_H 1

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RecoBase/OpFlash.h"
#include "RecoBase/OpHit.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/TimeService.h"
#include "OpFlashAlg.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes

// C++ Includes
#include <map>
#include <vector>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

namespace opdet {
 
  class OpFlashFinder : public art::EDProducer{
  public:
 
    // Standard constructor and destructor for an ART module.
    explicit OpFlashFinder(const fhicl::ParameterSet&);
    virtual ~OpFlashFinder();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& pset);

    // The producer routine, called once per event. 
    void produce(art::Event&); 
    
  private:

    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpHit collection
    
    Int_t    fBinWidth;
    Float_t  fFlashThreshold;
    Float_t  fWidthTolerance;
    Double_t fTrigCoinc;

  };

} 

namespace opdet {
  DEFINE_ART_MODULE(OpFlashFinder)
}

#endif 

#include "TriggerAlgo/TriggerAlgoMicroBoone.h"

namespace opdet {

  //----------------------------------------------------------------------------
  // Constructor
  OpFlashFinder::OpFlashFinder(const fhicl::ParameterSet & pset)
  {

    reconfigure(pset);

    produces< std::vector< recob::OpFlash > >();
    produces< art::Assns< recob::OpFlash, recob::OpHit > >();

  }

  //----------------------------------------------------------------------------
  void OpFlashFinder::reconfigure(fhicl::ParameterSet const& pset)
  {

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get< std::string >("InputModule");
      
    fBinWidth       = pset.get< int >   ("BinWidth");
    fFlashThreshold = pset.get< float > ("FlashThreshold");
    fWidthTolerance = pset.get< float > ("WidthTolerance");
    fTrigCoinc      = pset.get< double >("TrigCoinc");

  }

  //----------------------------------------------------------------------------
  // Destructor
  OpFlashFinder::~OpFlashFinder() 
  {
  }
   
  //----------------------------------------------------------------------------
  void OpFlashFinder::beginJob()
  {
  }

  //----------------------------------------------------------------------------
  void OpFlashFinder::endJob()
  { 
  }

  //----------------------------------------------------------------------------
  void OpFlashFinder::produce(art::Event& evt) 
  {

    // These are the storage pointers we will put in the event
    std::unique_ptr< std::vector< recob::OpFlash > > 
                      FlashPtr(new std::vector< recob::OpFlash >);
    std::unique_ptr< art::Assns< recob::OpFlash, recob::OpHit > >  
                      AssnPtr(new art::Assns< recob::OpFlash, recob::OpHit >);

    // This will keep track of what flashes will assoc to what ophits
    // at the end of processing
    std::vector< std::vector< int > > AssocList;

    art::ServiceHandle< geo::Geometry > GeometryHandle;
    geo::Geometry const& Geometry(*GeometryHandle);

    art::ServiceHandle< util::TimeService > ts_ptr;
    util::TimeService const& ts(*ts_ptr);

    // Get OpHits from the event
    art::Handle< std::vector< recob::OpHit > > OpHitHandle;
    evt.getByLabel(fInputModule, OpHitHandle);

    RunFlashFinder(*OpHitHandle,
                   *FlashPtr,
                   AssocList,
                   fBinWidth,
                   Geometry,
                   fFlashThreshold,
                   fWidthTolerance,
                   ts,
                   fTrigCoinc);

    // Make the associations which we noted we need
    for (size_t i = 0; i != AssocList.size(); ++i)
      util::CreateAssn(*this, evt, *(AssnPtr.get()), i, 
                       AssocList[i].begin(), AssocList[i].end());
    
    // Store results into the event
    evt.put(std::move(FlashPtr));
    evt.put(std::move(AssnPtr));
    
  }

} // namespace opdet
