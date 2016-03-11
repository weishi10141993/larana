// -*- mode: c++; c-basic-offset: 2; -*-
// Gleb Sinev, Duke, 2016
//
// This module finds periods of time-localized activity
// on each channel of the optical system and creates OpHits.
//
// Split from OpFlashFinder_module.cc
// by Ben Jones, MIT, 2013
//


#ifndef OpHitFinder_H
#define OpHitFinder_H 1

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSlidingWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoFixedWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoCFD.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoRollingMean.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "lardata/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
#include "larsim/Simulation/BeamGateInfo.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/Exception.h"

// ROOT includes

// C++ Includes
#include <map>
#include <vector>
#include <string>
#include <memory>

namespace opdet {
 
  class OpHitFinder : public art::EDProducer{
  public:
 
    // Standard constructor and destructor for an ART module.
    explicit OpHitFinder(const fhicl::ParameterSet&);
    virtual ~OpHitFinder();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& pset);

    // The producer routine, called once per event. 
    void produce(art::Event&); 
    
    std::map< int, int >  GetChannelMap();
    std::vector< double > GetSPEScales();

  private:

    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpDetWaveform collection
    std::string fGenModule;
    std::vector< std::string > fInputLabels;
    std::set< unsigned int > fChannelMasks;
    
    pmtana::PulseRecoManager  fPulseRecoMgr;
    pmtana::PMTPulseRecoBase* fThreshAlg;
    pmtana::PMTPedestalBase*  fPedAlg;

    Float_t  fHitThreshold;
    Bool_t   fAreaToPE;
    Float_t  fSPEArea;

    unsigned int fMaxOpChannel;
   
    std::vector< double > fSPESize;

  };

} 

namespace opdet {
  DEFINE_ART_MODULE(OpHitFinder)
}

#endif 

namespace opdet {

  //----------------------------------------------------------------------------
  // Constructor
  OpHitFinder::OpHitFinder(const fhicl::ParameterSet & pset):
    fPulseRecoMgr()
  {

    reconfigure(pset);

    // Initialize the hit finder algorithm
    auto const hit_alg_pset = pset.get< fhicl::ParameterSet >("HitAlgoPset");
    std::string threshAlgName = hit_alg_pset.get< std::string >("Name");
    if      (threshAlgName == "Threshold") 
      fThreshAlg = new pmtana::AlgoThreshold(hit_alg_pset);
    else if (threshAlgName == "SiPM") 
      fThreshAlg = new pmtana::AlgoSiPM(hit_alg_pset);
    else if (threshAlgName == "SlidingWindow")
      fThreshAlg = new pmtana::AlgoSlidingWindow(hit_alg_pset);
    else if (threshAlgName == "FixedWindow")
      fThreshAlg = new pmtana::AlgoFixedWindow(hit_alg_pset);
    else if (threshAlgName == "CFD" )
      fThreshAlg = new pmtana::AlgoCFD(hit_alg_pset);
    else throw art::Exception(art::errors::UnimplementedFeature)
                    << "Cannot find implementation for " 
                    << threshAlgName << " algorithm.\n";   

    auto const ped_alg_pset = pset.get< fhicl::ParameterSet >("PedAlgoPset");
    std::string pedAlgName = ped_alg_pset.get< std::string >("Name");
    if      (pedAlgName == "Edges")
      fPedAlg = new pmtana::PedAlgoEdges(ped_alg_pset);
    else if (pedAlgName == "RollingMean")
      fPedAlg = new pmtana::PedAlgoRollingMean(ped_alg_pset);
    else if (pedAlgName == "UB"   )
      fPedAlg = new pmtana::PedAlgoUB(ped_alg_pset);
    else throw art::Exception(art::errors::UnimplementedFeature)
                    << "Cannot find implementation for " 
                    << pedAlgName << " algorithm.\n";   

    produces< std::vector< recob::OpHit > >();

    fPulseRecoMgr.AddRecoAlgo(fThreshAlg);
    fPulseRecoMgr.SetDefaultPedAlgo(fPedAlg);

  }

  //----------------------------------------------------------------------------
  void OpHitFinder::reconfigure(fhicl::ParameterSet const& pset)
  {

    // Indicate that the Input Module comes from .fcl
    fInputModule   = pset.get< std::string >("InputModule");
    fGenModule     = pset.get< std::string >("GenModule");
    fInputLabels   = pset.get< std::vector< std::string > >("InputLabels");
      
    for (auto const& ch : pset.get< std::vector< unsigned int > >
                            ("ChannelMasks", std::vector< unsigned int >()))
      fChannelMasks.insert(ch);
    
    fHitThreshold = pset.get< float >("HitThreshold");
    fAreaToPE     = pset.get< bool > ("AreaToPE");
    fSPEArea      = pset.get< float >("SPEArea");

    auto const& geometry(*lar::providerFrom< geo::Geometry >());
    fMaxOpChannel = geometry.MaxOpChannel();
    
    fSPESize = GetSPEScales();

  }

  //----------------------------------------------------------------------------
  // Destructor
  OpHitFinder::~OpHitFinder() 
  {
  
    delete fThreshAlg;
    delete fPedAlg;

  }
   
  //----------------------------------------------------------------------------
  void OpHitFinder::beginJob()
  {
  }

  //----------------------------------------------------------------------------
  void OpHitFinder::endJob()
  { 
  }

  //----------------------------------------------------------------------------
  void OpHitFinder::produce(art::Event& evt) 
  {

    // These is the storage pointer we will put in the event
    std::unique_ptr< std::vector< recob::OpHit > >
                      HitPtr(new std::vector< recob::OpHit >);

    std::vector< const sim::BeamGateInfo* > beamGateArray;
    try 
    { 
      evt.getView(fGenModule, beamGateArray); 
    }
    catch (art::Exception const& err) 
    { 
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;
    }

    auto const& geometry(*lar::providerFrom< geo::Geometry >());
    auto const& detectorClocks
                  (*lar::providerFrom< detinfo::DetectorClocksService >());

    //
    // Get the pulses from the event
    //

    // Reserve a large enough array
    int totalsize = 0;
    for (auto label : fInputLabels) 
    {
      art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
      evt.getByLabel(fInputModule, label, wfHandle);
      if (!wfHandle.isValid()) continue; // Skip non-existent collections
      totalsize += wfHandle->size();
    }

    // Load pulses into WaveformVector
    std::vector< raw::OpDetWaveform > WaveformVector;
    WaveformVector.reserve(totalsize);
    for (auto label : fInputLabels) 
    {
      art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
      evt.getByLabel(fInputModule, label, wfHandle);
      if (!wfHandle.isValid()) continue; // Skip non-existent collections

      //WaveformVector.insert(WaveformVector.end(), 
      //                      wfHandle->begin(), wfHandle->end());
      for(auto const& wf : *wfHandle) 
      {
        if (fChannelMasks.find(wf.ChannelNumber()) 
            != fChannelMasks.end()) continue;
        WaveformVector.push_back(wf);
      }
    }

    RunHitFinder(WaveformVector,
                 *HitPtr,
                 fPulseRecoMgr,
                 *fThreshAlg,
                 geometry,
                 fHitThreshold,
                 detectorClocks,
                 fSPESize,
                 fAreaToPE);

    // Store results into the event
    evt.put(std::move(HitPtr));
    
  }

  //----------------------------------------------------------------------------
  std::vector< double > OpHitFinder::GetSPEScales()
  {
    // This will eventually interface to some kind of gain service
    // or database. For now all SPE scales are set to 20 ADC.
    // Alternatively all SPE scales are set to fSPEArea 
    // if hit area is used to calculate number of PEs
    
    if (fAreaToPE) return std::vector< double >(fMaxOpChannel + 1, fSPEArea);
    // temp fix while we work out the experiment-agnostic service 
    // that provides this info
    else           return std::vector< double >(fMaxOpChannel + 1, 20); 
  }

} // namespace opdet
