/*!
 * Title:   OpFlash Algorithims
 * Author:  Ben Jones, MIT (Edited by wketchum@lanl.gov)
 *
 * Description:
 * These are the algorithms used by OpFlashFinder to produce flashes.
 */

#include <algorithm>
#include <functional>

#include "OpFlashAlg.h"
#include "RecoBase/OpHit.h"
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/OpDetGeo.h"
#include "TH1D.h"
#include "TFile.h"

namespace opdet{

//-------------------------------------------------------------------------------------------------
  void GetTriggerTime(std::vector<const sim::BeamGateInfo*> const& beamGateArray,
		      double const& opdigi_TimeBegin,
		      double const& opdigi_SampleFreq,
		      optdata::TimeSlice_t const& trigger_frame_size,
		      unsigned int& Frame, double& TrigTime)
  {
    // This code gets the trigger time from the BeamGateInfo
    //  Eventually it will be replaced with code to look in a 
    //  corresponding reco object
    
    if(beamGateArray.size()>1)
      mf::LogError("OpFlashFinder")<<"Found more than 1 beam gate info! " << beamGateArray.size();
    
    for(auto trig : beamGateArray){
      
      double start_time = 1.e-3 * trig->Start() - opdigi_TimeBegin;
      
      if(start_time < 0) 
	throw cet::exception(__FUNCTION__) << "Found beam time (=" 
					   << trig->Start()
					   << ") before discrete clock count start (=" 
					   << opdigi_TimeBegin*1.e3 << ")\n";
      
      unsigned int ticks = (unsigned int)(start_time * opdigi_SampleFreq);
      Frame  = ticks / trigger_frame_size;
      //Sample = (unsigned short)(ticks - (Frame * trigger_frame_size));

      TrigTime = start_time; //trigger time

    }

    

  }
  
  //-------------------------------------------------------------------------------------------------
  void RunFlashFinder(std::vector<optdata::FIFOChannel> const& FIFOChannelVector,
		      std::vector<recob::OpHit>& HitVector,
		      std::vector<recob::OpFlash>& FlashVector,
		      std::vector< std::vector<int> >& AssocList,
		      optdata::TimeSlice_t const& trigger_frame_size,
		      int const& BinWidth,
		      pmtana::PulseRecoManager const& PulseRecoMgr,
		      pmtana::AlgoThreshold const& ThreshAlg,
		      std::map<int,int> const& ChannelMap,
		      geo::Geometry const& geom,
		      float const& HitThreshold,
		      float const& FlashThreshold,
		      float const& WidthTolerance,
		      unsigned int const& TrigFrame,
		      double const& TrigTimeAbs,
		      double const& opdigi_SampleFreq,
		      std::vector<double> const& SPESize,
		      float const& TrigCoinc)
  {
    
    std::map<unsigned short, std::vector<const optdata::FIFOChannel*> > FIFOChanByFrame;
    for(auto const& fifochannel : FIFOChannelVector)
      FIFOChanByFrame[fifochannel.Frame()].push_back(&fifochannel);
    
    for(auto fifoframe : FIFOChanByFrame)
      ProcessFrame(fifoframe.first,
		   fifoframe.second,
		   HitVector,
		   FlashVector,
		   AssocList,
		   trigger_frame_size,
		   BinWidth,
		   PulseRecoMgr,
		   ThreshAlg,
		   ChannelMap,
		   geom,
		   HitThreshold,
		   FlashThreshold,
		   WidthTolerance,
		   TrigFrame,
		   TrigTimeAbs,
		   opdigi_SampleFreq,
		   SPESize,
		   TrigCoinc);
    
  }
  
  //-------------------------------------------------------------------------------------------------
  void writeHistogram(std::vector<double> const& binned){

    TH1D *binned_histogram = new TH1D("binned_histogram","Collection of All OpHits;Time (ms);PEs",binned.size(),0,binned.size());
    for(size_t i=0; i<binned.size(); i++)
      binned_histogram->SetBinContent(i,binned.at(i));

    TFile f_out("output_hist.root","RECREATE");
    binned_histogram->Write();
    f_out.Close();

    delete binned_histogram;
  }

  //-------------------------------------------------------------------------------------------------
  void checkOnBeamFlash(std::vector<recob::OpFlash> const& FlashVector){
    for(auto const& flash : FlashVector){
      if(flash.OnBeamTime()==1) std::cout << "OnBeamFlash with time " <<  flash.Time() << std::endl;;
    }
  }

  //-------------------------------------------------------------------------------------------------
  void ProcessFrame(unsigned short Frame,
		    std::vector<const optdata::FIFOChannel*> const& FIFOChannelFramePtrVector,
		    std::vector<recob::OpHit>& HitVector,
		    std::vector<recob::OpFlash>& FlashVector,
		    std::vector< std::vector<int> >& AssocList,
		    optdata::TimeSlice_t const& TimeSlicesPerFrame,
		    int const& BinWidth,
		    pmtana::PulseRecoManager const& PulseRecoMgr,
		    pmtana::AlgoThreshold const& ThreshAlg,
		    std::map<int,int> const& ChannelMap,
		    geo::Geometry const& geom,
		    float const& HitThreshold,
		    float const& FlashThreshold,
		    float const& WidthTolerance,
		    unsigned int const& TrigFrame,
		    double const& TrigTimeAbs,
		    double const& opdigi_SampleFreq,
		    std::vector<double> const& SPESize,
		    float const& TrigCoinc)

  {

    //The +3000 here is microboone specific, to account for the beam-gate window size.

    // These are the accumulators which will hold broad-binned light yields
    std::vector<double>  Binned1((TimeSlicesPerFrame + 3000 + BinWidth)/BinWidth);
    std::vector<double>  Binned2((TimeSlicesPerFrame + 3000 + BinWidth)/BinWidth);
    
    // These will keep track of which pulses put activity in each bin
    std::vector<std::vector<int> > Contributors1((TimeSlicesPerFrame + 3000 + BinWidth)/BinWidth);
    std::vector<std::vector<int> > Contributors2((TimeSlicesPerFrame + 3000 + BinWidth)/BinWidth);
    
    // These will keep track of where we have met the flash condition
    //  (in order to prevent second pointless loop)
    std::vector<int> FlashesInAccumulator1;
    std::vector<int> FlashesInAccumulator2;
    
    const size_t NHits_prev = HitVector.size();
    unsigned int NOpChannels = geom.NOpChannels();

    for(auto const& fifo_ptr : FIFOChannelFramePtrVector){

      const int Channel = ChannelMap.at((int)fifo_ptr->ChannelNumber());
      const uint32_t TimeSlice = fifo_ptr->TimeSlice();

      if( Channel<0 || Channel > int(NOpChannels) ) {
	mf::LogError("OpFlashFinder")<<"Error! unrecognized channel number " << Channel<<". Ignoring pulse";
	continue;
      }
      
      if( TimeSlice>TimeSlicesPerFrame ){
	mf::LogError("OpFlashFinder")<<"This slice " << TimeSlice<< "is outside the countable region - skipping";
	continue;
      }
      
      PulseRecoMgr.RecoPulse(&(*fifo_ptr));
      
      const size_t NPulses = ThreshAlg.GetNPulse();
      for(size_t k=0; k<NPulses; ++k){
	
	ConstructHit( HitThreshold,
		      Channel,
		      TimeSlice,
		      Frame,
		      ThreshAlg.GetPulse(k),
		      TimeSlicesPerFrame,
		      opdigi_SampleFreq,
		      TrigTimeAbs,
		      SPESize.at(Channel),
		      HitVector );

	unsigned int AccumIndex1 = GetAccumIndex(ThreshAlg.GetPulse(k)->t_max, 
						TimeSlice, 
						BinWidth, 
						0);
	FillAccumulator(AccumIndex1,
			HitVector.size()-1,
			HitVector[HitVector.size()-1].PE(),
			FlashThreshold,
			Binned1,
			Contributors1,
			FlashesInAccumulator1);

	unsigned int AccumIndex2 = GetAccumIndex(ThreshAlg.GetPulse(k)->t_max, 
						TimeSlice, 
						BinWidth, 
						BinWidth/2);
	FillAccumulator(AccumIndex2,
			HitVector.size()-1,
			HitVector[HitVector.size()-1].PE(),
			FlashThreshold,
			Binned2,
			Contributors2,
			FlashesInAccumulator2);
  
      }
      
    }//end loop over FIFO channels in frame

    //Now start to create flashes
    //First, need vector to keep track of which hits belong to which flashes
    std::vector< std::vector<int> > HitsPerFlash;
    size_t NHitsThisFrame = HitVector.size() - NHits_prev;
    
    //if(Frame==1) writeHistogram(Binned1);

    AssignHitsToFlash(FlashesInAccumulator1,
		      FlashesInAccumulator2,
		      Binned1,
		      Binned2,
		      Contributors1,
		      Contributors2,
		      NHitsThisFrame,
		      HitVector,
		      HitsPerFlash,
		      FlashThreshold);

    // Now we do the fine grained part.  
    // Subdivide each flash into sub-flashes with overlaps within hit widths (assumed wider than photon travel time)
    std::vector<std::vector<int> > RefinedHitsPerFlash;
    RefineHitsToFlash(HitsPerFlash,
		      HitVector,
		      RefinedHitsPerFlash,
		      WidthTolerance,
		      FlashThreshold);

    //Now we have all our hits assigned to a flash. Make the recob::OpFlash objects.
    ConstructFlashes(RefinedHitsPerFlash,
		     HitVector,
		     FlashVector,
		     TimeSlicesPerFrame,
		     geom,
		     TrigFrame,
		     Frame,
		     TrigCoinc);


    RemoveLateLight(FlashVector,
    		    RefinedHitsPerFlash);

    //checkOnBeamFlash(FlashVector);

    //Finally, write the association list
    //The transform adds a constant offset to the elements of each vector in RefinedHitsPerFlash
    //back_inserter tacks the result onto the end of AssocList
    for(auto & HitIndicesThisFlash : RefinedHitsPerFlash){
      for(auto & HitIndex : HitIndicesThisFlash)
	HitIndex += NHits_prev;
      AssocList.push_back(HitIndicesThisFlash);
    }
    
  }//end ProcessFrame

  //-------------------------------------------------------------------------------------------------
  void ConstructHit( float const& HitThreshold,
		     int const& Channel,
		     uint32_t const& TimeSlice,
		     unsigned short const& Frame,
		     const pmtana::pulse_param* pulse,
		     optdata::TimeSlice_t const& TimeSlicesPerFrame,
		     double const& opdigi_SampleFreq,
		     double const& TrigTimeAbs,
		     double const& SPESize,
		     std::vector<recob::OpHit>& HitVector)
  {
    if( pulse->peak<HitThreshold ) return;
    
    double AbsTime = (pulse->t_max + TimeSlice + Frame*TimeSlicesPerFrame)/opdigi_SampleFreq;
    double RelTime = AbsTime - TrigTimeAbs;
    double PE = pulse->peak/SPESize;
    
    HitVector.emplace_back( Channel,
			    RelTime,
			    AbsTime,
			    Frame,
			    pulse->t_end - pulse->t_start,
			    pulse->area,
			    pulse->peak,
			    PE,
			    0.);
  }


  //-------------------------------------------------------------------------------------------------
  unsigned int GetAccumIndex(double const& TMax, 
			     uint32_t const& TimeSlice, 
			     int const& BinWidth, 
			     double const& BinOffset){
    return ( (TMax + TimeSlice) + BinOffset) / BinWidth;
  }

  //-------------------------------------------------------------------------------------------------
  void FillAccumulator(unsigned int const& AccumIndex,
		       unsigned int const& HitIndex,
		       double const& PE,
		       float const& FlashThreshold,
		       std::vector<double> & Binned,
		       std::vector< std::vector<int> > & Contributors,
		       std::vector<int> & FlashesInAccumulator)
  {
  
      
    (Contributors.at(AccumIndex)).push_back(HitIndex);
    
    Binned.at(AccumIndex) += PE; 
    
    // If this wasn't a flash already, add it to the list
    if( Binned.at(AccumIndex)>=FlashThreshold &&
	(Binned.at(AccumIndex)-PE)<FlashThreshold )
      FlashesInAccumulator.push_back(AccumIndex);
    
  }

  //-------------------------------------------------------------------------------------------------
  void FillFlashesBySizeMap(std::vector<int> const& FlashesInAccumulator,
			    std::vector<double> const& BinnedPE,
			    int const& Accumulator,
			    std::map<double, std::map<int,std::vector<int> >, std::greater<double> > & FlashesBySize){
    for( auto const& flash : FlashesInAccumulator)
      FlashesBySize[BinnedPE.at(flash)][Accumulator].push_back(flash);
  }

  //-------------------------------------------------------------------------------------------------
  void FillHitsThisFlash(std::vector< std::vector<int> > const& Contributors,
			 int const& Bin,
			 size_t const& NHits_prev,
			 std::vector<int> const& HitClaimedByFlash,
			 std::vector<int> & HitsThisFlash){
    
    // for each hit in the flash
    for(auto const& HitIndex : Contributors.at(Bin)){
      
      // if unclaimed, claim it
      if(HitClaimedByFlash.at(HitIndex-NHits_prev)==-1)
	HitsThisFlash.push_back(HitIndex);
    }
  }

  //-------------------------------------------------------------------------------------------------
  void ClaimHits(std::vector<recob::OpHit> const& HitVector,
		 std::vector<int> const& HitsThisFlash,
		 float const& FlashThreshold,
		 std::vector< std::vector<int> > & HitsPerFlash,
		 size_t const& NHits_prev,
		 std::vector<int> & HitClaimedByFlash){

    //Check for newly claimed hits
    double PE = 0;
    for(auto const& Hit : HitsThisFlash)
      PE += HitVector.at(Hit).PE();
    
    if(PE < FlashThreshold) return;
    
    // add the flash to the list
    HitsPerFlash.push_back(HitsThisFlash);
    
    // and claim all the hits
    for(auto const& Hit : HitsThisFlash){
      if(HitClaimedByFlash.at(Hit-NHits_prev)==-1)
	HitClaimedByFlash.at(Hit-NHits_prev)=HitsPerFlash.size()-1;
    }//end loop over hits in this flash
    
  }

  //-------------------------------------------------------------------------------------------------
  void AssignHitsToFlash( std::vector<int> const& FlashesInAccumulator1,
			  std::vector<int> const& FlashesInAccumulator2,
			  std::vector<double> const& Binned1,
			  std::vector<double> const& Binned2,
			  std::vector< std::vector<int> > const& Contributors1,
			  std::vector< std::vector<int> > const& Contributors2,
			  size_t const& NHits,
			  std::vector<recob::OpHit> const& HitVector,
			  std::vector< std::vector<int> >& HitsPerFlash,
			  float const& FlashThreshold)
  {

    size_t NHits_prev = HitVector.size() - NHits;

    // Sort all the flashes found by size. The structure is:
    // FlashesBySize[flash size][accumulator_num] = [flash_index1, flash_index2...]     
    std::map<double, std::map<int,std::vector<int> >, std::greater<double> > FlashesBySize;
      
    // Sort the flashes by size using map
    FillFlashesBySizeMap(FlashesInAccumulator1,
			 Binned1,
			 1,
			 FlashesBySize);
    FillFlashesBySizeMap(FlashesInAccumulator2,
			 Binned2,
			 2,
			 FlashesBySize);
    
  
    // This keeps track of which hits are claimed by which flash
    std::vector<int > HitClaimedByFlash(NHits,-1);

    // Walk from largest to smallest, claiming hits. The biggest flash always gets dibbs,
    // but we keep track of overlaps for re-merging later (do we? ---WK)
    for(auto const& itFlash : FlashesBySize){

      // If several with same size, walk walk through accumulators
      for(auto const& itAcc : itFlash.second){

	  int Accumulator = itAcc.first;
	  
	  // Walk through flash-tagged bins in this accumulator
	  for(auto const& Bin : itAcc.second){

	    std::vector<int>   HitsThisFlash;

	    if(Accumulator==1)
	      FillHitsThisFlash(Contributors1,
				Bin,
				NHits_prev,
				HitClaimedByFlash,
				HitsThisFlash);	    
	    else if(Accumulator==2)
	      FillHitsThisFlash(Contributors2,
				Bin,
				NHits_prev,
				HitClaimedByFlash,
				HitsThisFlash);

	    ClaimHits(HitVector,
		      HitsThisFlash,
		      FlashThreshold,
		      HitsPerFlash,
		      NHits_prev,
		      HitClaimedByFlash);
	    
	  }//end loop over this accumulator

      }//end loops over accumulators

    } // end of loops over sorted flashes
    
  }//end AssignHitsToFlash


  //-------------------------------------------------------------------------------------------------
  void RefineHitsToFlash(std::vector< std::vector<int> > const& HitsPerFlash,
			 std::vector<recob::OpHit> const& HitVector,
			 std::vector< std::vector<int> >& RefinedHitsPerFlash,
			 float const& WidthTolerance,
			 float const& FlashThreshold){

    for(size_t iFlash=0; iFlash!=HitsPerFlash.size(); ++iFlash){}
    for(auto const& flash : HitsPerFlash){

      // Sort the hits in each flash by their size using map
      // HitsBySize[HitSize] = [hit1, hit2 ...]
      std::map<double, std::vector<int> > HitsBySize;
      for(auto const& HitID : flash)
	HitsBySize[HitVector.at(HitID).PE()].push_back(HitID);

      // Heres what we do:
      //  1.Start with the biggest remaining hit
      //  2.Look for any within one width of this hit
      //  3.Find the new upper and lower bounds of the flash
      //  4.Collect again
      //  5.Repeat until no new hits collected
      //  6.Remove these hits from consideration and repeat
      
      std::vector<bool> HitsUsed(HitVector.size(),false);
      bool HitsLeft = true, StillCollecting = true;
      double PEAccumulated, FlashMaxTime, FlashMinTime;
      std::vector<int> HitsThisRefinedFlash;
      while(HitsLeft){

	HitsLeft = false;
	PEAccumulated = 0; FlashMaxTime = 0; FlashMinTime = 0;

	for(auto const& itHit : HitsBySize){
	  for(auto const& HitID : itHit.second){

	    if(!HitsUsed.at(HitID)){
	      HitsLeft = true;
	      PEAccumulated += HitVector.at(HitID).PE();
	      FlashMaxTime = HitVector.at(HitID).PeakTime() + 
		WidthTolerance * HitVector.at(HitID).Width();
	      FlashMinTime = HitVector.at(HitID).PeakTime() - 
		WidthTolerance * HitVector.at(HitID).Width();
	      HitsThisRefinedFlash.push_back(HitID);
	      HitsUsed.at(HitID)=true; 
	      break;
	    }

	  }//end loop over inner vector
	  
	  if(HitsLeft==true) break;
	
	}// end loop over HitsBySize map


	// Now iteratively collect hits within the flash width
	StillCollecting = true;
	while(StillCollecting){
	  
	  StillCollecting = false;
	  for(auto const& HitID : flash){
	    
	    if(!HitsUsed.at(HitID)){
	      double HitTime  =   HitVector.at(HitID).PeakTime();
	      double HitWidth =   HitVector.at(HitID).Width();
	      double FlashTime =  0.5*(FlashMaxTime + FlashMinTime);
	      double FlashWidth = 0.5*(FlashMaxTime - FlashMinTime);
	      if( fabs(HitTime-FlashTime) < WidthTolerance*(HitWidth + FlashWidth) ){
		HitsThisRefinedFlash.push_back(HitID);
		HitsUsed.at(HitID) = true;
		FlashMaxTime = std::max(FlashMaxTime, HitTime + HitWidth);
		FlashMinTime = std::min(FlashMinTime, HitTime - HitWidth);
		PEAccumulated += HitVector.at(HitID).PE();
		StillCollecting = true;
	      }//end if hit matches flash
	    }//end if hit unused
	  }// end loop over hits
	}//end while collecting hits into flashes
	

	// We did our collecting, now check if the flash is
	// still good and push back
	if(PEAccumulated >= FlashThreshold)
	    RefinedHitsPerFlash.push_back(HitsThisRefinedFlash);
		
      }//end while there are hits left
    }//end loop over hits per flash

  }//end RefineHitsToFlash


  //-------------------------------------------------------------------------------------------------
  void ConstructFlashes(std::vector< std::vector<int> > const& RefinedHitsPerFlash,
			std::vector<recob::OpHit> const& HitVector,
			std::vector<recob::OpFlash>& FlashVector,
			uint32_t const& TimeSlicesPerFrame,
			geo::Geometry const& geom,
			unsigned int const& TrigFrame,
			unsigned short const& Frame,
			float const& TrigCoinc){

    for(auto const& HitsPerFlashVec : RefinedHitsPerFlash){

      double MaxTime = -1e9, MinTime = 1e9;

      std::vector<double> PEs(geom.NOpChannels());
      unsigned int Nplanes = geom.Nplanes();
      std::vector<double> sumw(Nplanes,0), sumw2(Nplanes,0);

      double TotalPE=0, AveTime=0, AveAbsTime=0, FastToTotal=0, sumy=0, sumz=0, sumy2=0, sumz2=0;

      for(auto const& HitID : HitsPerFlashVec){

	double FastToTotalThisHit   = HitVector.at(HitID).FastToTotal();
	double PEThisHit            = HitVector.at(HitID).PE();
	double TimeThisHit          = HitVector.at(HitID).PeakTime();
	double AbsTimeThisHit       = HitVector.at(HitID).PeakTimeAbs();
	unsigned int ChannelThisHit = HitVector.at(HitID).OpChannel();
	if(TimeThisHit > MaxTime) MaxTime = TimeThisHit;
	if(TimeThisHit < MinTime) MinTime = TimeThisHit;

	// These quantities for the flash are defined as the weighted averages
	//   over the hits
	FastToTotal += FastToTotalThisHit *PEThisHit;
	AveTime     += TimeThisHit        *PEThisHit;
	AveAbsTime  += AbsTimeThisHit     *PEThisHit;
	
	// These are totals
	TotalPE     += PEThisHit;
	PEs.at(ChannelThisHit)+=PEThisHit;
	
	unsigned int o=0, c=0;
	geom.OpChannelToCryoOpDet(ChannelThisHit,o,c);
	
	double xyz[3];
	geom.Cryostat(c).OpDet(o).GetCenter(xyz);
	
	for(size_t p=0; p!=Nplanes; p++){
	  unsigned int w = geom.NearestWire(xyz,p);
	  sumw.at(p)  += w*PEThisHit;
	  sumw2.at(p) += w*w*PEThisHit;
	}             
		
	sumy+=xyz[1]*PEThisHit; sumy2+=xyz[1]*xyz[1]*PEThisHit;
	sumz+=xyz[2]*PEThisHit; sumz2+=xyz[2]*xyz[2]*PEThisHit;
      
      }//end loop over hits in (refined) flash

	    
      AveTime     /= TotalPE;
      AveAbsTime  /= TotalPE;
      FastToTotal /= TotalPE;
      
      double meany = sumy / TotalPE;
      double meanz = sumz / TotalPE;
      
      double widthy = std::sqrt(sumy2*TotalPE - sumy*sumy)/TotalPE;
      double widthz = std::sqrt(sumz2*TotalPE - sumz*sumz)/TotalPE;
 
      std::vector<double> WireCenters(Nplanes,0);
      std::vector<double> WireWidths(Nplanes,0);
      
      for(size_t p=0; p!=Nplanes; ++p){
	WireCenters.at(p) = sumw.at(p)/TotalPE;
	WireWidths.at(p)  = std::sqrt(sumw2.at(p)*TotalPE - sumw.at(p)*sumw.at(p))/TotalPE;
      }
      
      bool InBeamFrame = (Frame==TrigFrame);
      double TimeWidth = (MaxTime-MinTime)/2.;
      
      int OnBeamTime =0; 
      if( std::abs(AveTime) < TrigCoinc ) OnBeamTime=1;
      /*
      std::cout << "Time summary for this flash:"
		<< "\n\tAveTime: " << AveTime
		<< "\n\tAveAbsTime: " << AveAbsTime
		<< "\n\tTotalPE: " << TotalPE
		<< "\n\tmeany,widthy: " << meany << "," << widthy
		<< "\n\tmeanz,widthz: " << meanz << "," << widthz
		<< "\n\tFrame: " << Frame
		<< "\n\tInBeamFrame: " << InBeamFrame
		<< "\n\tOnBeamTime: " << OnBeamTime << std::endl;
      */
      FlashVector.emplace_back( AveTime,
				TimeWidth,
				AveAbsTime,
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
				WireWidths);
      
      
    }//end loop over all flashes

  }//end ConstructFlashes


  //-------------------------------------------------------------------------------------------------
  void RemoveLateLight(std::vector<recob::OpFlash>& FlashVector,
		       std::vector< std::vector<int> >& RefinedHitsPerFlash){

    std::vector<bool> MarkedForRemoval(RefinedHitsPerFlash.size(),false);

    size_t BeginFlash = FlashVector.size() - RefinedHitsPerFlash.size();

    recob::OpFlashSortByTime sort_flash_by_time;
    std::sort(FlashVector.begin()+BeginFlash,
	      FlashVector.end(),
	      sort_flash_by_time);

    for(size_t iFlash=BeginFlash; iFlash!=FlashVector.size(); ++iFlash){

      double iTime = FlashVector.at(iFlash).Time();
      double iPE   = FlashVector.at(iFlash).TotalPE();
      double iWidth= FlashVector.at(iFlash).TimeWidth();
      
      for(size_t jFlash=iFlash+1; jFlash!=FlashVector.size(); ++jFlash){

	if(MarkedForRemoval.at(jFlash-BeginFlash)) continue;

	double jTime = FlashVector.at(jFlash).Time();
	double jPE   = FlashVector.at(jFlash).TotalPE();
	double jWidth= FlashVector.at(jFlash).TimeWidth();
	
	// Calculate hypothetical PE if this were actually a late flash from i
	//  Argon time const is 1600 ns, so 1.6.
	double HypPE = iPE * jWidth / iWidth * std::exp(-(jTime-iTime)/1.6);
	double nsigma = (jPE-HypPE)/std::sqrt(HypPE);
	
	// If smaller than, or within 2sigma of expectation,
	//  attribute to late light and toss out
	if( nsigma < 3. ) 
	  MarkedForRemoval.at(jFlash-BeginFlash)=true;
      }	    
    }

    for(int iFlash=MarkedForRemoval.size()-1; iFlash!=-1; --iFlash){
      if(MarkedForRemoval.at(iFlash)){
	RefinedHitsPerFlash.erase(RefinedHitsPerFlash.begin()+iFlash);
	FlashVector.erase(FlashVector.begin()+BeginFlash+iFlash);
      }
    }


  }

}//end namespace opdet

