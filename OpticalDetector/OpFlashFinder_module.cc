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
    
    void GetTriggerTime(art::Event& evt, unsigned int& Frame, unsigned short& Sample);

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
    Float_t fTrigCoinc;
    

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
      fPulseRecoMgr(pset.get<fhicl::ParameterSet>("reco_man")),
      fThreshAlg(pset.get<fhicl::ParameterSet>("algo_threshold"))    
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
    fTrigCoinc      = pset.get<float>        ("TrigCoinc");
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
    uint32_t fTimeSlicesPerFrame=110000;
    
    // Find out when the trigger happened 
    unsigned int   TrigFrame = 0;
    unsigned short TrigTime = 0;
    GetTriggerTime(evt, TrigFrame, TrigTime);
    

    // Get the pulses from the event
    art::Handle< std::vector< optdata::FIFOChannel > > FIFOChannelHandle;
    evt.getByLabel(fInputModule, FIFOChannelHandle);
    
    // First we organize all the pulses we saw by frame
    std::map<unsigned int, std::vector<art::Ptr<optdata::FIFOChannel> > > FIFOChanByFrame;
    for(size_t i=0; i!=FIFOChannelHandle->size(); ++i)
      {
	art::Ptr<optdata::FIFOChannel> ptr(FIFOChannelHandle,i);
	FIFOChanByFrame[FIFOChannelHandle->at(i).Frame()].push_back(ptr);
      }

    // Now proceed to flash-find frame by frame
    for(auto itFrame :  FIFOChanByFrame)
      {

	// OpHits and Flashes we find in this frame
	//  (will get transferred into storage containers later)
	std::vector<recob::OpHit>   OpHitsThisFrame;
	std::vector<recob::OpFlash> OpFlashesThisFrame;

	// These are the accumulators which will hold broad-binned light yields
	std::vector<double>  Binned1((fTimeSlicesPerFrame + fBinWidth)/fBinWidth);
	std::vector<double>  Binned2((fTimeSlicesPerFrame + fBinWidth)/fBinWidth);
	
	// These will keep track of which pulses put activity in each bin
	std::vector<std::vector<int> > Contributors1((fTimeSlicesPerFrame + fBinWidth)/fBinWidth);
	std::vector<std::vector<int> > Contributors2((fTimeSlicesPerFrame + fBinWidth)/fBinWidth);
	  
	// These will keep track of where we have met the flash condition
	//  (in order to prevent second pointless loop)
	std::vector<int> FlashesInAccumulator1;
	std::vector<int> FlashesInAccumulator2;
	
	uint32_t Frame = itFrame.first;

	for(size_t i=0; i!=itFrame.second.size(); ++i)
	  {
	    const optdata::FIFOChannel* fifo_ptr = ( (itFrame.second).at(i) ).get() ;
	
	    fPulseRecoMgr.RecoPulse(fifo_ptr);

	    size_t NPulses = fThreshAlg.GetNPulse();
	
	    int Channel   = fChannelMap[FIFOChannelHandle->at(i).ChannelNumber()];
	    uint32_t TimeSlice = FIFOChannelHandle->at(i).TimeSlice();
	    
	    if((Channel<0)||(Channel > int(fNOpChannels)))
	      {
		mf::LogError("OpFlashFinder")<<"Error! unrecognized channel number " << Channel<<". Ignoring pulse";
		continue;
	      }

	    if((TimeSlice<0)||(TimeSlice>fTimeSlicesPerFrame))
	      {
		mf::LogError("OpFlashFinder")<<"This slice " << TimeSlice<< "is outside the countable region - skipping";
		continue;
	      }

	    for(size_t k=0; k!=NPulses; ++k)
	      {
		{
		  double Peak  = fThreshAlg.GetPulse(k)->peak;

		  if(Peak>fHitThreshold)
		    {
		      double TMax  = fThreshAlg.GetPulse(k)->t_max;
		      double Width = 
			fThreshAlg.GetPulse(k)->t_end - 
			fThreshAlg.GetPulse(k)->t_start;
		      double Area  = fThreshAlg.GetPulse(k)->area;
		      
		      double AbsTime = TMax + TimeSlice;
		      double RelTime = TMax + TimeSlice - TrigTime + 
			(Frame - TrigFrame)*fTimeSlicesPerFrame;
		      
		      double PE = Peak/fSPESize.at(Channel);
		      
		      OpHitsThisFrame.push_back( recob::OpHit( Channel,
							       AbsTime,
							       RelTime,
							       Frame,
							       Width,
							       Area,
							       Peak,
							       PE,
							       0) ) ;
		  
		      size_t HitIndex = OpHitsThisFrame.size()-1;
		      
		      size_t Accum1Index = int(TimeSlice  / fBinWidth);
		      size_t Accum2Index = int((TimeSlice + fBinWidth/2)/fBinWidth);
		      
		      Contributors1[Accum1Index].push_back(HitIndex);
		      Contributors2[Accum2Index].push_back(HitIndex);
		      
		      double NewVal1 = Binned1[Accum1Index] += PE; 
		      double NewVal2 = Binned2[Accum2Index] += PE;  
		      
		      
		      if(NewVal1>=fFlashThreshold)
			{
			  
			  // If this wasn't a flash already, add it to the list
			  if( (NewVal1-PE)<fFlashThreshold)
			    FlashesInAccumulator1.push_back(Accum1Index);
			}
		      
		      if(NewVal2>=fFlashThreshold)
			{
			  if( (NewVal2-PE)<fFlashThreshold)
			    FlashesInAccumulator2.push_back(Accum2Index);
			}
		    }
		}
	      }
	    
	  } // end loop over hits this frame
    
	// Sort all the flashes we found by size.  The structure is
	//      FlashesBySize[flash size][accumulator_num] = [flash_index1, flash_index2...]     
	std::map<double, std::map<int, std::vector<int > > > FlashesBySize;
      
	// This will track which hits go with each flash
	//      HitsPerFlash[flash_index] = [hit_index1, hit_index2...]
	std::vector<std::vector<int> > HitsPerFlash;

	// Sort the flashes by size using map, and count them

	size_t NHits    = itFrame.second.size();
	size_t NFlashes = 0;
      
 	for(size_t i1=0; i1!=FlashesInAccumulator1.size(); ++i1)
	  {
	    int    BinIndex  = FlashesInAccumulator1[i1];
	    double FlashSize = Binned1[BinIndex];
	    FlashesBySize[FlashSize][1].push_back(BinIndex);
	    NFlashes++;
	  }
	for(size_t i2=0; i2!=FlashesInAccumulator2.size(); ++i2)
	  {
	    int    BinIndex  = FlashesInAccumulator2[i2];
	    double FlashSize = Binned2[BinIndex];
	    FlashesBySize[FlashSize][2].push_back(BinIndex);
	    NFlashes++;
	  }
    
	// This keeps track of which hits are claimed by which flash
	std::vector<int >      HitClaimedByFlash (NHits,-1) ;
	
	
	// Walk from largest to smallest, claiming hits. biggest flash always gets dibbs,
	//  but we keep track of overlaps for re-merging later

	// Walk through flashes in size order, largest to smallest
	for(auto itFlash = FlashesBySize.rbegin(); itFlash!=FlashesBySize.rend(); ++itFlash)
	  {
	    // If several with same size, walk walk through accumulators
	    for(auto itAcc = itFlash->second.begin(); itAcc!=itFlash->second.end(); ++itAcc)
	      {
		int Accumulator = itAcc->first;

		// Walk through flash-tagged bins in this accumulator
		for(size_t iBin=0; iBin!=itAcc->second.size(); ++iBin)
		  {

		    std::vector<int>   HitsThisFlash;
		    std::vector<int>   Overlaps;


		    int Bin = itAcc->second.at(iBin);
		    if(Accumulator==1)
		      // for each hit in the flash
		      for(size_t iHit = 0; iHit!=Contributors1[Bin].size(); ++iHit)
			{
			  int HitIndex = Contributors1[Bin][iHit];
			  
			  // if unclaimed, claim it
			  if(HitClaimedByFlash[HitIndex]==-1)
			    {
			      HitsThisFlash.push_back(HitIndex);
			    }
			  else
			    // if claimed, note the overlap to deal with later
			    Overlaps.push_back(HitClaimedByFlash[HitIndex]);
			}
		    else if (Accumulator==2)
		      for(size_t iHit =0; iHit!=Contributors2[Bin].size(); ++iHit)
			{
			  int HitIndex = Contributors2[Bin][iHit];
			  if(HitClaimedByFlash[HitIndex]==-1)
			    {
			      // if unclaimed, claim it
			      HitsThisFlash.push_back(HitIndex);
			    }
			  else
			    // if already claimed, note the overlap to deal with later
			    Overlaps.push_back(HitClaimedByFlash[HitIndex]);
			  
			}
		  
		    // Were there any unclaimed hits this flash?
		    if(HitsThisFlash.size()>0)
		      {
			// if there were some, check if this guy still gets above threshold
			double PE = 0;
			for(size_t iHit=0; iHit!=HitsThisFlash.size(); ++iHit)
			  PE += OpHitsThisFrame.at(HitsThisFlash.at(iHit)).PE();
			
			// if it still gets over threshold
			if(PE >= fFlashThreshold)
			  {
			    // add the flash to the list
			    HitsPerFlash.push_back(HitsThisFlash);
			    
			    // and claim all the hits
			    for(size_t ih=0; ih!=HitsThisFlash.size(); ++ih)
			      {
				if(HitClaimedByFlash[HitsThisFlash[ih]]==-1)
				  HitClaimedByFlash[HitsThisFlash[ih]]=HitsPerFlash.size()-1;
			      }
			  }
		      }
		    HitsThisFlash.clear();
		  }
	      }
	  } // end of loops over sorted flashes

	// Now we do the fine grained part.  Subdivide each flash into sub-flashes with overlaps within hit widths (assumed wider than photon travel time)
	std::vector<std::vector<int> > RefinedHitsPerFlash;
       
	for(size_t iFlash=0; iFlash!=HitsPerFlash.size(); ++iFlash)
	  {
	    // Sort the hits in each flash by their size using map
	    //   HitsBySize[HitSize] = [hit1, hit2 ...]
	    std::map<double, std::vector<int> > HitsBySize;
	    for(size_t iHit=0; iHit!=HitsPerFlash.at(iFlash).size(); ++iHit)
	      {
		size_t HitID = HitsPerFlash.at(iFlash).at(iHit); 
		HitsBySize[OpHitsThisFrame.at(HitID).PE()].push_back(HitID);
	      }

	    // Heres what we do:
	    //  1.Start with the biggest remaining hit
	    //  2.Look for any within one width of this hit
	    //  3.Find the new upper and lower bounds of the flash
	    //  4.Collect again
	    //  5.Repeat until no new hits collected
	    //  6.Remove these hits from consideration and repeat

	    std::vector<bool> HitsUsed(OpHitsThisFrame.size(),false);
	   	    
	    bool HitsLeft = true;
	    while(HitsLeft)
	      {
		HitsLeft = false;
	
		double PEAccumulated=0;
		
		double FlashMaxTime=0, FlashMinTime=0;
		std::vector<int> HitsThisRefinedFlash;
		
		// find the largest remaining hit
		for(auto itHit = HitsBySize.rbegin(); itHit!=HitsBySize.rend(); ++itHit)
		  {
		    for(size_t iHit=0; iHit!=itHit->second.size(); ++iHit)
		      {
			int HitID = itHit->second.at(iHit);
			if(!HitsUsed[HitID])
			  {
			    HitsLeft = true;
			    double FlashTime = OpHitsThisFrame.at(HitID).PeakTime();
			    double FlashWidth = OpHitsThisFrame.at(HitID).Width();
			    PEAccumulated +=OpHitsThisFrame.at(HitID).PE();
			    FlashMaxTime = FlashTime + fWidthTolerance * FlashWidth;
			    FlashMinTime = FlashTime - fWidthTolerance * FlashWidth;
			    HitsThisRefinedFlash.push_back(HitID);
			    HitsUsed[HitID]=true;
			    
			    break;
			  }
		      }
		    if(HitsLeft==true) break;
		  }

		// Now iteratively collect hits within the flash width
		bool StillCollecting = true;
		while(StillCollecting)		
		  {
		    StillCollecting = false;
		    for(size_t iHit=0; iHit!=HitsPerFlash.at(iFlash).size(); ++iHit)
		      {
			size_t HitID = HitsPerFlash.at(iFlash).at(iHit);
			
			if(!HitsUsed[HitID])
			  {
			    double HitTime  =   OpHitsThisFrame.at(HitID).PeakTime();
			    double HitWidth =   OpHitsThisFrame.at(HitID).Width();
			    double FlashTime =  0.5*(FlashMaxTime + FlashMinTime);
			    double FlashWidth = 0.5*(FlashMaxTime - FlashMinTime);
			    if( ( fabs(HitTime-FlashTime) < fWidthTolerance*(HitWidth + FlashWidth) ) )
			      {
				HitsThisRefinedFlash.push_back(HitID);
				HitsUsed[HitID] = true;
				FlashMaxTime = std::max(FlashMaxTime, HitTime + HitWidth);
				FlashMinTime = std::min(FlashMinTime, HitTime - HitWidth);
				PEAccumulated += OpHitsThisFrame.at(HitID).PE();
				StillCollecting = true;
				
			      }
			  }
		      }
		  }
		// We did our collecting, now check if the flash is
		//  still good and push back
		if(PEAccumulated >= fFlashThreshold)
		  {
		    RefinedHitsPerFlash.push_back(HitsThisRefinedFlash);
		  }
	      }
	  }
	

	
	// We may have lost some flashes along the way, 
	//  so update NFlashes
	NFlashes = RefinedHitsPerFlash.size();

	// Loop over those flashes we want to save
	for(size_t i=0; i!=NFlashes; ++i)
	  {
	    double MaxTime = -100, MinTime = fTimeSlicesPerFrame;
	    std::vector<double> PEs(fNOpChannels);
	    std::vector<double> sumw(fNplanes,0), sumw2(fNplanes,0);
	    double TotalPE=0, AveTime=0, FastToTotal=0, sumy=0, sumz=0, sumy2=0, sumz2=0;
	    
	    // For each hit in the flash
	    for(size_t ih=0; ih!=RefinedHitsPerFlash.at(i).size(); ++ih)
	      {
		int HitID = RefinedHitsPerFlash[i][ih];
		double FastToTotalThisHit  = OpHitsThisFrame.at(HitID).FastToTotal();
		double PEThisHit           = OpHitsThisFrame.at(HitID).PE();
		double TimeThisHit         = OpHitsThisFrame.at(HitID).PeakTime();
		double ChannelThisHit      = OpHitsThisFrame.at(HitID).OpChannel();
		if(TimeThisHit > MaxTime) MaxTime = TimeThisHit;
		if(TimeThisHit < MinTime) MinTime = TimeThisHit;

		
		// These quantities for the flash are defined as the weighted averages
		//   over the hits
		FastToTotal += FastToTotalThisHit *PEThisHit;
		AveTime     += TimeThisHit        * PEThisHit;
		
		// These are totals
		TotalPE     += PEThisHit;
		PEs[ChannelThisHit]+=PEThisHit;
	    

		// now dig into geometrical stuff...
	        art::ServiceHandle<geo::Geometry> geom;
		
		unsigned int o=0, c=0;
		geom->OpChannelToCryoOpDet(ChannelThisHit,o,c);
		
		double xyz[3];
		geom->Cryostat(c).OpDet(o).GetCenter(xyz);
		
		for(size_t p=0; p!=fNplanes; p++)
		  {
		    unsigned int w = geom->NearestWire(xyz,p);
		    sumw[p]  += w *        PEThisHit;
		    sumw2[p] += pow(w,2) * PEThisHit;
		  }             
		
		sumy+=xyz[1] * PEThisHit; sumy2+=pow(xyz[1],2) * PEThisHit;
		sumz+=xyz[2] * PEThisHit; sumz2+=pow(xyz[2],2) * PEThisHit;
	      } 
	    
	    AveTime     /= TotalPE;
	    FastToTotal /= TotalPE;
	    
	    double meany = sumy / TotalPE;
	    double meanz = sumz / TotalPE;
	
	    double widthy = pow(sumy2 -  pow(sumy,2)/TotalPE, 0.5) / pow(TotalPE,0.5);
	    double widthz = pow(sumz2 -  pow(sumz,2)/TotalPE, 0.5) / pow(TotalPE,0.5);
	
	    std::vector<double> WireCenters(fNplanes,0);
	    std::vector<double> WireWidths(fNplanes,0);
	    
	    for(size_t p=0; p!=fNplanes; ++p)
	      {
		WireCenters[p] = sumw[p]/TotalPE;
		WireWidths[p]  = pow(sumw2[p] - pow(sumw[p],2)/TotalPE, 0.5) / pow(TotalPE,0.5);
	      }

	    bool InBeamFrame = (Frame==TrigFrame);
	    double RelTime   = AveTime - TrigTime +
	      (Frame - TrigFrame) * fTimeSlicesPerFrame;
	    double TimeWidth = (MaxTime-MinTime)/2.;
	    
	    int OnBeamTime =0; 
	    if((InBeamFrame) && (fabs(RelTime) < fTrigCoinc)) OnBeamTime=1;


	    // Make the flash
	    OpFlashesThisFrame.push_back( recob::OpFlash( RelTime,
							  TimeWidth,
							  AveTime,
							  Frame,
							  PEs, 
							  InBeamFrame,
							  OnBeamTime,
							  FastToTotal,
							  meany, 
							  widthy, 
							  meanz, 
							  widthz, 
							  WireCenters, 
							  WireWidths ));
	    
	    
	  }


	// We now toss out flashes which are consistent with being late light of another flash
	//  The expected light in a trailing flash is A0 exp(-dt / tau) * w1/w0
	// We allow 2sigma fluctuations around this also.
	
	std::vector<bool> MarkedForRemoval(OpFlashesThisFrame.size(),false);
	
	for(size_t iFlash=0; iFlash!=OpFlashesThisFrame.size(); ++iFlash)
	  {
	    double iTime = OpFlashesThisFrame.at(iFlash).Time();
	    double iPE   = OpFlashesThisFrame.at(iFlash).TotalPE();
	    double iWidth= OpFlashesThisFrame.at(iFlash).TimeWidth();
	    
	    for(size_t jFlash=0; jFlash!=OpFlashesThisFrame.size(); ++jFlash)
	      {
		double jTime = OpFlashesThisFrame.at(jFlash).Time();
		
		if(jTime <= iTime) continue;
		
		double jPE   = OpFlashesThisFrame.at(jFlash).TotalPE();
		double jWidth= OpFlashesThisFrame.at(jFlash).TimeWidth();
		
		// Calculate hypothetical PE if this were actually a late flash from i
		//  Argon time const is 1600, so 100 samples.
		double HypPE = iPE * jWidth / iWidth * exp(-(jTime-iTime)/100.);
		double nsigma = (jPE-HypPE)/pow(HypPE,0.5);
		
		// If smaller than, or within 2sigma of expectation,
		//  attribute to late light and toss out
		if( nsigma < 3. ) 
		  {
		    MarkedForRemoval[jFlash]=true;
		  }
	      }
	    
	  }
	
	for(int iFlash=OpFlashesThisFrame.size()-1.; iFlash!=-1; --iFlash)
	  {
	    if(MarkedForRemoval[iFlash])
	      {
		MarkedForRemoval.erase(MarkedForRemoval.begin()       + iFlash);
		RefinedHitsPerFlash.erase(RefinedHitsPerFlash.begin() + iFlash);
		OpFlashesThisFrame.erase(OpFlashesThisFrame.begin()   + iFlash);
	      }
	  }


	// Now add into the main storage containers, and note where we're gona want associations
	//  (making the assocs prematurely leads to a segfault, it turns out)

	size_t FirstHit = HitPtr->size();
	
	for(size_t i=0; i!=OpFlashesThisFrame.size(); ++i)
	  FlashPtr->push_back(OpFlashesThisFrame.at(i));	  
	
	for(size_t i=0; i!=OpHitsThisFrame.size(); ++i)
	  HitPtr->push_back(OpHitsThisFrame.at(i));

	for(size_t i=0; i!=RefinedHitsPerFlash.size(); ++i)
	  {
	    for(size_t j=0; j!=RefinedHitsPerFlash.at(i).size(); ++j)
	      RefinedHitsPerFlash[i][j]+=FirstHit;
	    AssocList.push_back(RefinedHitsPerFlash.at(i));
	  }


	// Clean up
	FlashesBySize.clear();
	HitsPerFlash.clear();
	RefinedHitsPerFlash.clear();
	OpHitsThisFrame.clear();
	OpFlashesThisFrame.clear();
	HitClaimedByFlash.clear();
	Binned1.clear();
	Binned2.clear();
	Contributors1.clear();
	Contributors2.clear();
	FlashesInAccumulator1.clear();
	FlashesInAccumulator2.clear();

      }
  


	


    

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


    // Clean up
    AssocList.clear();    
    for(size_t i=0; i!=FIFOChanByFrame.size(); ++i)
      FIFOChanByFrame.at(i).clear();
    FIFOChanByFrame.clear();

  }


  //--------------------------------------

  void OpFlashFinder::GetTriggerTime(art::Event& evt, unsigned int& Frame, unsigned short& Sample)
  {
    // This code gets the trigger time from the BeamGateInfo
    //  Eventually it will be replaced with code to look in a 
    //  corresponding reco object


    //
    // Read in BeamGateInfo array ... handle the case if there's no BeamGateInfo 
    //
    std::vector<const sim::BeamGateInfo*> beamGateArray;
    try{
      evt.getView(fGenModule, beamGateArray);
    }
    catch ( art::Exception const& err ) {
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;
    }

    art::ServiceHandle<trigger::TriggerAlgoMicroBoone> _trig_mod;
    // 
    for(size_t index=0; index < beamGateArray.size(); ++index){

      const sim::BeamGateInfo* trig(beamGateArray.at(index));

      Sample = _trig_mod->ConvertTime(trig->Start());

      Frame  = Sample/(_trig_mod->FrameSizeTrigger());
      
      Sample = Sample - (Frame * (_trig_mod->FrameSizeTrigger()));

      if(index>0)

	mf::LogError("OpFlashFinder")<<"Found more than 1 beam gate info!";
    }
   
  }

  //-------------------------------------

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

