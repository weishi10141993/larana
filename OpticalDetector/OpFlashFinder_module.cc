// Ben Jones, MIT, 2013
//
// This module finds periods of time-localized activity
// from the optical system, called Flashes.


#ifndef OpFlashFinder_H
#define OpFlashFinder_H 1

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/OpDetGeo.h"
#include "OpticalDetectorData/FIFOChannel.h"
#include "OpticalDetector/AlgoThreshold.h"
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
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes

// C++ Includes
#include <map>
#include <vector>
#include <iostream>
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

    
    pmtana::PulseRecoManager  fPulseRecoMgr;
    pmtana::AlgoThreshold     fThreshAlg;

    Int_t   fChannelMapMode;
    Int_t   fBinWidth;
    Float_t fFlashThreshold;
    Float_t fHitThreshold;
    Float_t fWidthTolerance;
    Double_t fTrigCoinc;
    

    unsigned int fNplanes;
    unsigned int fNOpChannels;
   
    std::vector<double> fSPESize;
    std::map<int, int>  fChannelMap;
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
    fPulseRecoMgr(),
    fThreshAlg()
  {

    reconfigure(pset);

    produces<std::vector< recob::OpFlash> >();
    produces<std::vector< recob::OpHit> >();
    produces<art::Assns<recob::OpFlash, recob::OpHit> >();

    fPulseRecoMgr.AddRecoAlgo(&fThreshAlg);
    fPulseRecoMgr.SetPedAlgo(pmtana::kHEAD);


  }

  //---------------------------------------------

  void OpFlashFinder::reconfigure(fhicl::ParameterSet const& pset)
  {
    // Indicate that the Input Module comes from .fcl
    fInputModule    = pset.get<std::string>("InputModule");
    fGenModule      = pset.get<std::string>("GenModule");

    fChannelMapMode = pset.get<int>          ("ChannelMapMode");

    fBinWidth       = pset.get<int>          ("BinWidth");
    fFlashThreshold = pset.get<float>        ("FlashThreshold");
    fWidthTolerance = pset.get<float>        ("WidthTolerance");
    fTrigCoinc      = pset.get<double>       ("TrigCoinc");
    fHitThreshold   = pset.get<float>        ("HitThreshold");
    

    art::ServiceHandle<geo::Geometry> geom;
    fNOpChannels = geom->NOpChannels();
    fNplanes     = geom->Nplanes();
    
    fSPESize     = GetSPEScales();
    fChannelMap  = GetChannelMap();

  }

  //-----------------------------------------------------------------------
  // Destructor
  OpFlashFinder::~OpFlashFinder() 
  {}
   
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
    std::unique_ptr<std::vector< recob::OpHit > >   HitPtr (new std::vector<recob::OpHit >);
    std::unique_ptr<std::vector< recob::OpFlash > > FlashPtr (new std::vector<recob::OpFlash >);
    std::unique_ptr< art::Assns<recob::OpFlash, recob::OpHit > >  AssnPtr( new art::Assns<recob::OpFlash, recob::OpHit>);

    // This will keep track of what flashes will assoc to what ophits
    //  at the end of processing
    std::vector<std::vector<int> > AssocList;

    // Temporary - needs to be gotten from somewhere
    //  and not hard coded
    art::ServiceHandle<trigger::TriggerAlgoMicroBoone> trig_mod;

    art::ServiceHandle<util::TimeService> ts_sptr;
    util::TimeService const& ts = *ts_sptr;

    std::vector<const sim::BeamGateInfo*> beamGateArray;
    try { evt.getView(fGenModule, beamGateArray); }
    catch ( art::Exception const& err ){ 
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;
    }

    // Get the pulses from the event
    art::Handle< std::vector< optdata::FIFOChannel > > FIFOChannelHandle;
    evt.getByLabel(fInputModule, FIFOChannelHandle);
    std::vector<optdata::FIFOChannel> const& FIFOChannelVector(*FIFOChannelHandle);

    art::ServiceHandle<geo::Geometry> GeometryHandle;
    geo::Geometry const& Geometry(*GeometryHandle);

    RunFlashFinder(FIFOChannelVector,
		   *HitPtr,
		   *FlashPtr,
		   AssocList,
		   fBinWidth,
		   fPulseRecoMgr,
		   fThreshAlg,
		   fChannelMap,
		   Geometry,
		   fHitThreshold,
		   fFlashThreshold,
		   fWidthTolerance,
		   ts,
		   fSPESize,
		   fTrigCoinc);


    // Make the associations which we noted we need
    for(size_t i=0; i!=AssocList.size(); ++i)
      for(size_t j=0; j!=AssocList.at(i).size(); ++j)
	{
	  util::CreateAssn(*this, evt, *(FlashPtr), *(HitPtr), *(AssnPtr.get()), AssocList[i][j], AssocList[i][j], i);
	}
    
    // Store results into the event
    evt.put(std::move(FlashPtr));
    evt.put(std::move(HitPtr));
    evt.put(std::move(AssnPtr));
    
  }


  //--------------------------------------


  std::vector<double> OpFlashFinder::GetSPEScales()
  {
    // This will eventually interface to some kind of gain service
    //  or database.  For now all SPE scales are set to 20 ADC.
  
    return std::vector<double>(fNOpChannels,20);
  }

  //-------------------------------------
  std::map<int, int>  OpFlashFinder::GetChannelMap()
  {
    // This will eventually interface to a shaper map database
    // For now just hacked in by hand.
    // 
    // Two modes at present:
    //  Mode0 : direct map for MC, 1<->1, 2<->2, etc
    //  Mode1 : shaper map as used in PMT pre-com
    //
    std::map<int,int> ReturnMap;

    if(fChannelMapMode==0)
      {
	for(size_t i=0; i!=32; ++i)
	  ReturnMap[i] = i;
      }
    else if(fChannelMapMode==1)
      {
	ReturnMap[0] = 5; 
	ReturnMap[1] = 4;
	ReturnMap[2] = 7;
	ReturnMap[3] = 6;  
	ReturnMap[4] = 2;
	ReturnMap[5] = 1;
	ReturnMap[6] = 0;
	ReturnMap[7] = 3;
	ReturnMap[8] = 13;
	ReturnMap[9] = 12;
	ReturnMap[10]= 15;
	ReturnMap[11]= 14;
	ReturnMap[12]= 10;
	ReturnMap[13]= 9;
	ReturnMap[14]= 8;
	ReturnMap[15]= 11;
	ReturnMap[16]= 21;
	ReturnMap[17]= 20;
	ReturnMap[18]= 23;
	ReturnMap[19]= 22;
	ReturnMap[20]= 18;
	ReturnMap[21]= 17;
	ReturnMap[22]= 16;
	ReturnMap[23]= 19;
	ReturnMap[24]= 29;
	ReturnMap[25]= 28;
	ReturnMap[26]= -1;  // eventually 31
	ReturnMap[27]= -1;  // eventually 30
	ReturnMap[28]= 26;
	ReturnMap[29]= 25;
	ReturnMap[30]= 24;
	ReturnMap[31]= 27;
	ReturnMap[32]= -1;
	ReturnMap[33]= -1;
	ReturnMap[34]= -1;
	ReturnMap[35]= -1;
	ReturnMap[36]= -1;
	ReturnMap[37]= -1; // eventually paddle1
	ReturnMap[38]= -1; // eventually paddle2
	ReturnMap[39]= -1; // eventually paddle3
	ReturnMap[40]= -1; // eventually paddle4
      }
    else
      {
	mf::LogError("OpFlashFinder") << "Error : Unknown channel map mode!";
      }
    return ReturnMap;
  }


} // namespace opdet

