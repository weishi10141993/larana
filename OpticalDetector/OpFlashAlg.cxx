/*!
 * Title:   OpFlash Algorithims
 * Author:  Ben Jones, MIT (Edited by wketchum@lanl.gov)
 *
 * Description:
 * These are the algorithms used by OpFlashFinder to produce flashes.
 */

#include "OpFlashAlg.h"
#include "RecoBase/OpHit.h"
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace opdet{

//-------------------------------------------------------------------------------------------------
  void GetTriggerTime(std::vector<const sim::BeamGateInfo*> const& beamGateArray,
		      double const& opdigi_TimeBegin,
		      double const& opdigi_SampleFreq,
		      optdata::TimeSlice_t const& trigger_frame_size,
		      unsigned int& Frame, unsigned short& Sample)
  {
    // This code gets the trigger time from the BeamGateInfo
    //  Eventually it will be replaced with code to look in a 
    //  corresponding reco object
    
    if(beamGateArray.size()>0)
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
      Sample = (unsigned short)(ticks - (Frame * trigger_frame_size));

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
		      unsigned int const& NOpChannels,
		      unsigned int const& Nplanes,
		      float const& HitThreshold,
		      float const& FlashThreshold,
		      float const& WidthTolerance,
		      unsigned int const& TrigFrame,
		      unsigned int const& TrigSample,
		      std::vector<double> const& SPESize)
  {
    
    std::map<unsigned short, std::vector<const optdata::FIFOChannel*> > FIFOChanByFrame;
    for(auto fifochannel : FIFOChannelVector)
      FIFOChanByFrame[fifochannel.Frame()].push_back(&fifochannel);
    
    for(auto fifoframe : FIFOChanByFrame)
      ProcessFrame(fifoframe.first,
		   fifoframe.second,
		   HitVector,
		   FlashVector,
		   trigger_frame_size,
		   BinWidth,
		   PulseRecoMgr,
		   ThreshAlg,
		   ChannelMap,
		   NOpChannels,
		   HitThreshold,
		   FlashThreshold,
		   WidthTolerance,
		   TrigFrame,
		   TrigSample,
		   SPESize);
    
  }
  
  //-------------------------------------------------------------------------------------------------
  void ProcessFrame(unsigned short Frame,
		    std::vector<const optdata::FIFOChannel*> const& FIFOChannelFramePtrVector,
		    std::vector<recob::OpHit>& HitVector,
		    std::vector<recob::OpFlash>& FlashVector,
		    optdata::TimeSlice_t const& TimeSlicesPerFrame,
		    int const& BinWidth,
		    pmtana::PulseRecoManager const& PulseRecoMgr,
		    pmtana::AlgoThreshold const& ThreshAlg,
		    std::map<int,int> const& ChannelMap,
		    unsigned int const& NOpChannels,
		    float const& HitThreshold,
		    float const& FlashThreshold,
		    float const& WidthTolerance,
		    unsigned int const& TrigFrame,
		    unsigned int const& TrigTime,
		    std::vector<double> const& SPESize)
  {

    // These are the accumulators which will hold broad-binned light yields
    std::vector<double>  Binned1((TimeSlicesPerFrame + BinWidth)/BinWidth);
    std::vector<double>  Binned2((TimeSlicesPerFrame + BinWidth)/BinWidth);
    
    // These will keep track of which pulses put activity in each bin
    std::vector<std::vector<int> > Contributors1((TimeSlicesPerFrame + BinWidth)/BinWidth);
    std::vector<std::vector<int> > Contributors2((TimeSlicesPerFrame + BinWidth)/BinWidth);
    
    // These will keep track of where we have met the flash condition
    //  (in order to prevent second pointless loop)
    std::vector<int> FlashesInAccumulator1;
    std::vector<int> FlashesInAccumulator2;
    
    size_t NHits_prev = HitVector.size();

    for(auto & fifo_ptr : FIFOChannelFramePtrVector){

	const int Channel = ChannelMap.at(fifo_ptr->ChannelNumber());
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

	ConstructHits(Channel,
		      TimeSlice,
		      Frame,
		      ThreshAlg,
		      HitVector,
		      TimeSlicesPerFrame,
		      BinWidth,
		      HitThreshold,
		      FlashThreshold,
		      TrigFrame,
		      TrigTime,
		      SPESize.at(Channel),
		      Binned1, Binned2,
		      Contributors1, Contributors2,
		      FlashesInAccumulator1, FlashesInAccumulator2);

    }//end loop over FIFO channels in frame
    
    //Now start to create flashes
    //First, need vector to keep track of which hits belong to which flashes
    std::vector< std::vector<int> > HitsPerFlash;
    size_t NHits = HitVector.size() - NHits_prev;
    
    AssignHitsToFlash(FlashesInAccumulator1,
		      FlashesInAccumulator2,
		      Binned1,
		      Binned2,
		      Contributors1,
		      Contributors2,
		      NHits,
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
    
  }//end ProcessFrame

  //-------------------------------------------------------------------------------------------------
  void ConstructHits(int const& Channel,
		     uint32_t const& TimeSlice,
		     unsigned short const& Frame,
		     pmtana::AlgoThreshold const& ThreshAlg,
		     std::vector<recob::OpHit>& HitVector,
		     optdata::TimeSlice_t const& TimeSlicesPerFrame,
		     int const& BinWidth,
		     float const& HitThreshold,
		     float const& FlashThreshold,
		     unsigned int const& TrigFrame,
		     unsigned int const& TrigTime,
		     double const& SPESize,
		     std::vector<double> & Binned1,
		     std::vector<double> & Binned2,
		     std::vector< std::vector<int> > & Contributors1,
		     std::vector< std::vector<int> > & Contributors2,
		     std::vector<int> & FlashesInAccumulator1,
		     std::vector<int> & FlashesInAccumulator2)
  {
    const size_t NPulses = ThreshAlg.GetNPulse();
    for(size_t k=0; k<NPulses; ++k){
      
      double Peak  = ThreshAlg.GetPulse(k)->peak;
      if( Peak<HitThreshold ) continue;
      
      //note, these times are in units of pmt ticks. 
      //These need to be translated to real times once the time service gets up and running.
      double TMax  = ThreshAlg.GetPulse(k)->t_max;
      double Width = ThreshAlg.GetPulse(k)->t_end - ThreshAlg.GetPulse(k)->t_start;
      double Area  = ThreshAlg.GetPulse(k)->area;
      double AbsTime = TMax + TimeSlice;
      double RelTime = TMax + (double)TimeSlice - (double)TrigTime;
      RelTime += ((double)Frame - (double)TrigFrame)*TimeSlicesPerFrame;
      double PE = Peak/SPESize;
      
      HitVector.emplace_back( Channel,
			      AbsTime,
			      RelTime,
			      Frame,
			      Width,
			      Area,
			      Peak,
			      PE,
			      0.);
      
      int HitIndex = HitVector.size()-1;
      
      size_t Accum1Index = int(TimeSlice  / BinWidth);
      size_t Accum2Index = int((TimeSlice + BinWidth/2)/BinWidth);
      
      (Contributors1.at(Accum1Index)).push_back(HitIndex);
      (Contributors2.at(Accum2Index)).push_back(HitIndex);
      
      Binned1.at(Accum1Index) += PE; 
      Binned2.at(Accum2Index) += PE;  
      
      // If this wasn't a flash already, add it to the list
      if( Binned1.at(Accum1Index)>=FlashThreshold &&
	  (Binned1.at(Accum1Index)-PE)<FlashThreshold )
	FlashesInAccumulator1.push_back(Accum1Index);
      
      if( Binned2.at(Accum1Index)>=FlashThreshold &&
	  (Binned2.at(Accum1Index)-PE)<FlashThreshold )
	FlashesInAccumulator2.push_back(Accum1Index);
      
    }//end loop over pulses
    
  } // end ConstructHits
  

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

    // Sort all the flashes found by size. The structure is:
    // FlashesBySize[flash size][accumulator_num] = [flash_index1, flash_index2...]     
    std::map<double, std::map<int,std::vector<int> > > FlashesBySize;
      
    // Sort the flashes by size using map
    for( auto const& flash : FlashesInAccumulator1)
      FlashesBySize[Binned1.at(flash)][1].push_back(flash);
    for( auto const& flash : FlashesInAccumulator2)
      FlashesBySize[Binned2.at(flash)][2].push_back(flash);
  
    // This keeps track of which hits are claimed by which flash
    std::vector<int > HitClaimedByFlash(NHits,-1);

    // Walk from largest to smallest, claiming hits. The biggest flash always gets dibbs,
    // but we keep track of overlaps for re-merging later

    // Walk through flashes in size order, largest to smallest
    for(auto const& itFlash : FlashesBySize){

      // If several with same size, walk walk through accumulators
      for(auto const& itAcc : itFlash.second){

	  int Accumulator = itAcc.first;
	  
	  // Walk through flash-tagged bins in this accumulator
	  for(auto const& Bin : itAcc.second){

	    std::vector<int>   HitsThisFlash;

	    if(Accumulator==1)

	      // for each hit in the flash
	      for(auto const& HitIndex : Contributors1.at(Bin)){

		// if unclaimed, claim it
		if(HitClaimedByFlash.at(HitIndex)==-1)
		  HitsThisFlash.push_back(HitIndex);
	      }
	    else if(Accumulator==2)

	      // for each hit in the flash
	      for(auto const& HitIndex : Contributors2.at(Bin)){
		
		// if unclaimed, claim it
		if(HitClaimedByFlash.at(HitIndex)==-1)
		  HitsThisFlash.push_back(HitIndex);
	      }
	    
	    //Check for newly claimed hits
	    double PE = 0;
	    for(auto const& Hit : HitsThisFlash)
	      PE += HitVector.at(Hit).PE();
	    
	    // if it still gets over threshold
	    if(PE >= FlashThreshold){
	      
	      // add the flash to the list
	      HitsPerFlash.push_back(HitsThisFlash);
	      
	      // and claim all the hits
	      for(auto const& Hit : HitsThisFlash){
		if(HitClaimedByFlash.at(Hit)==-1)
		  HitClaimedByFlash.at(Hit)=HitsPerFlash.size()-1;
	      }//end loop over hits in this flash

	    }//end if PE above threshold

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


}//end namespace opdet

