// -*- mode: c++; c-basic-offset: 2; -*-
// Ben Jones, MIT, 2013
//
// This module finds periods of time-localized activity
// from the optical system, called Flashes.
//
// Modified to make it more detector agnostic
// by Gleb Sinev, Duke, 2015
//


#ifndef OpFlashFinder_H
#define OpFlashFinder_H 1

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/OpDetGeo.h"
#include "RawData/OpDetWaveform.h"
#include "OpticalDetector/AlgoThreshold.h"
#include "OpticalDetector/AlgoSiPM.h"
#include "OpticalDetector/AlgoPedestal.h"
#include "OpticalDetector/PulseRecoManager.h"
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
    void produce (art::Event&); 
    
    std::map<int, int>  GetChannelMap();
    std::vector<double> GetSPEScales();


  private:

    // The parameters we'll read from the .fcl file.
    std::string fInputModule;              // Input tag for OpDet collection
    std::string fGenModule ;
    std::vector< std::string > fInputLabels;
    std::string fThreshAlgName;
    std::set<unsigned int> fChannelMasks;
    
    pmtana::PulseRecoManager  fPulseRecoMgr;
    pmtana::PMTPulseRecoBase* fThreshAlg;

    Int_t    fBinWidth;
    Float_t  fFlashThreshold;
    Float_t  fHitThreshold;
    Float_t  fWidthTolerance;
    Double_t fTrigCoinc;
    Bool_t   fAreaToPE;
    Float_t  fSPEArea;
    

    unsigned int fNplanes;
    unsigned int fNOpChannels;
    unsigned int fMaxOpChannel;
   
    std::vector<double> fSPESize;

  };


} 

namespace opdet {
  DEFINE_ART_MODULE(OpFlashFinder)
}

#endif 

#include "TriggerAlgo/TriggerAlgoMicroBoone.h"

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpFlashFinder::OpFlashFinder(const fhicl::ParameterSet & pset):
    fPulseRecoMgr()
  {

    reconfigure(pset);

    // Initialize the hit finder
    if      (fThreshAlgName == "Threshold") 
              fThreshAlg = new pmtana::AlgoThreshold();
    else if (fThreshAlgName == "SiPM") 
              fThreshAlg = new pmtana::AlgoSiPM(
                pset.get< fhicl::ParameterSet >("algo_threshold"));
    else throw art::Exception(art::errors::UnimplementedFeature)
                    << "Cannot find implementation for " 
                    << fThreshAlgName << " algorithm.\n";   

    produces<std::vector< recob::OpFlash> >();
    produces<std::vector< recob::OpHit> >();
    produces<art::Assns<recob::OpFlash, recob::OpHit> >();

    fPulseRecoMgr.AddRecoAlgo(fThreshAlg);
    fPulseRecoMgr.SetPedAlgo(pmtana::kHEAD);

  }

  //-----------------------------------------------------------------------
  void OpFlashFinder::reconfigure(fhicl::ParameterSet const& pset)
  {

    // Indicate that the Input Module comes from .fcl
    fInputModule    = pset.get<std::string>("InputModule");
    fGenModule      = pset.get<std::string>("GenModule");
    fInputLabels    = pset.get<std::vector<std::string> >("InputLabels");
    fThreshAlgName  = pset.get< fhicl::ParameterSet >("algo_threshold").get< std::string >("HitFinder");
      
    for(auto const& ch : pset.get<std::vector<unsigned int> >("ChannelMasks",std::vector<unsigned int>()))
      fChannelMasks.insert(ch);
    
    fBinWidth       = pset.get<int>          ("BinWidth");
    fFlashThreshold = pset.get<float>        ("FlashThreshold");
    fWidthTolerance = pset.get<float>        ("WidthTolerance");
    fTrigCoinc      = pset.get<double>       ("TrigCoinc");
    fHitThreshold   = pset.get<float>        ("HitThreshold");
    fAreaToPE       = pset.get<bool>         ("AreaToPE");
    fSPEArea        = pset.get<float>        ("SPEArea");

    art::ServiceHandle<geo::Geometry> geom;
    fNOpChannels  = geom->NOpChannels();
    fMaxOpChannel = geom->MaxOpChannel();
    fNplanes      = geom->Nplanes();
    
    fSPESize     = GetSPEScales();

  }

  //-----------------------------------------------------------------------
  // Destructor
  OpFlashFinder::~OpFlashFinder() 
  {
  
    delete fThreshAlg;

  }
   
  //-----------------------------------------------------------------------
  void OpFlashFinder::beginJob()
  {
  }

  //-----------------------------------------------------------------------
  void OpFlashFinder::endJob()
  { 
  }

  //-----------------------------------------------------------------------
  void OpFlashFinder::produce(art::Event& evt) 
  {

    // These are the storage pointers we will put in the event
    std::unique_ptr< std::vector< recob::OpHit > >   HitPtr (new std::vector<recob::OpHit >);
    std::unique_ptr< std::vector< recob::OpFlash > > FlashPtr (new std::vector<recob::OpFlash >);
    std::unique_ptr< art::Assns<recob::OpFlash, recob::OpHit > >  AssnPtr( new art::Assns<recob::OpFlash, recob::OpHit>);

    // This will keep track of what flashes will assoc to what ophits
    //  at the end of processing
    std::vector<std::vector<int> > AssocList;

    std::vector<const sim::BeamGateInfo*> beamGateArray;
    try { evt.getView(fGenModule, beamGateArray); }
    catch ( art::Exception const& err ){ 
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;
    }

    art::ServiceHandle<geo::Geometry> GeometryHandle;
    geo::Geometry const& Geometry(*GeometryHandle);

    art::ServiceHandle<util::TimeService> ts_ptr;
    util::TimeService const& ts(*ts_ptr);

    //
    // Get the pulses from the event
    //

    // Reserve a large enough array
    int totalsize = 0;
    for ( auto label : fInputLabels) {
      art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
      evt.getByLabel(fInputModule, label, wfHandle);
      if (!wfHandle.isValid()) continue; // Skip non-existent collections
      totalsize += wfHandle->size();
    }

    // Load pulses into WaveformVector
    std::vector< raw::OpDetWaveform > WaveformVector;
    WaveformVector.reserve(totalsize);
    for ( auto label : fInputLabels) {
      art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
      evt.getByLabel(fInputModule, label, wfHandle);
      if (!wfHandle.isValid()) continue; // Skip non-existent collections
      //WaveformVector.insert(WaveformVector.end(), wfHandle->begin(), wfHandle->end());
      for(auto const& wf : *wfHandle) {
	if(fChannelMasks.find(wf.ChannelNumber()) != fChannelMasks.end()) continue;
	WaveformVector.push_back(wf);
      }
    }

    RunFlashFinder(WaveformVector,
                   *HitPtr,
                   *FlashPtr,
                   AssocList,
                   fBinWidth,
                   fPulseRecoMgr,
                   *fThreshAlg,
                   Geometry,
                   fHitThreshold,
                   fFlashThreshold,
                   fWidthTolerance,
                   ts,
                   fSPESize,
                   fAreaToPE,
                   fTrigCoinc);

    // Make the associations which we noted we need
    for(size_t i=0; i!=AssocList.size(); ++i)
      util::CreateAssn(*this, evt, *(AssnPtr.get()), i, AssocList[i].begin(), AssocList[i].end());
    
    // Store results into the event
    evt.put(std::move(FlashPtr));
    evt.put(std::move(HitPtr));
    evt.put(std::move(AssnPtr));
    
  }


  //-----------------------------------------------------------------------
  std::vector<double> OpFlashFinder::GetSPEScales()
  {
    // This will eventually interface to some kind of gain service
    // or database. For now all SPE scales are set to 20 ADC.
    // Alternatively all SPE scales are set to fSPEArea 
    // if hit area is used to calculate number of PEs.
    
    if (fAreaToPE) return std::vector<double>(fNOpChannels,fSPEArea);
    else           return std::vector<double>(fMaxOpChannel+1,20); // temp fix while we work out the expeiment-agnostic service that provides this info.
  }



} // namespace opdet

