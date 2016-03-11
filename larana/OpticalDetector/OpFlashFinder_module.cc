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
#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoBase/OpFlash.h"
#include "lardata/RecoBase/OpHit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larana/OpticalDetector/OpFlashAlg.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// ROOT includes

// C++ Includes
#include <vector>
#include <string>
#include <memory>

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
                      flashPtr(new std::vector< recob::OpFlash >);
    std::unique_ptr< art::Assns< recob::OpFlash, recob::OpHit > >  
                      assnPtr(new art::Assns< recob::OpFlash, recob::OpHit >);

    // This will keep track of what flashes will assoc to what ophits
    // at the end of processing
    std::vector< std::vector< int > > assocList;

    auto const& geometry(*lar::providerFrom< geo::Geometry >());

    auto const& detectorClocks
       (*lar::providerFrom< detinfo::DetectorClocksService >());

    // Get OpHits from the event
    art::Handle< std::vector< recob::OpHit > > opHitHandle;
    evt.getByLabel(fInputModule, opHitHandle);

    RunFlashFinder(*opHitHandle,
                   *flashPtr,
                   assocList,
                   fBinWidth,
                   geometry,
                   fFlashThreshold,
                   fWidthTolerance,
                   detectorClocks,
                   fTrigCoinc);

    // Make the associations which we noted we need
    for (size_t i = 0; i != assocList.size(); ++i)
    {
      art::PtrVector< recob::OpHit > opHitPtrVector;
      for (int const& hitIndex : assocList.at(i))
      {
        art::Ptr< recob::OpHit > opHitPtr(opHitHandle, hitIndex);
        opHitPtrVector.push_back(opHitPtr);
      }

      util::CreateAssn(*this, evt, *flashPtr, opHitPtrVector, 
                                        *(assnPtr.get()), i);
    }
    
    // Store results into the event
    evt.put(std::move(flashPtr));
    evt.put(std::move(assnPtr));
    
  }

} // namespace opdet
