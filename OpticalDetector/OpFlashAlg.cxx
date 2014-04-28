/*!
 * Title:   OpFlash Algorithims
 * Author:  Ben Jones, MIT (Edited by wketchum@lanl.gov)
 *
 * Description:
 * These are the algorithms used by OpFlashFinder to produce flashes.
 */

#include "OpFlashAlg.h"
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
		      unsigned int const& TrigFrame,
		      unsigned int const& TrigSample,
		      std::vector<double> const& SPESize)
  {
    
    std::map<unsigned int, std::vector<const optdata::FIFOChannel*> > FIFOChanByFrame;
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
		   TrigFrame,
		   TrigSample,
		   SPESize);
    
  }
  
  //-------------------------------------------------------------------------------------------------
  void ProcessFrame(unsigned int Frame,
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
    

    for(std::vector<const optdata::FIFOChannel*>::const_iterator i_chan=FIFOChannelFramePtrVector.begin();
	i_chan!=FIFOChannelFramePtrVector.end();
	i_chan++)
      {

	const optdata::FIFOChannel* fifo_ptr = i_chan;
	//const std::vector<uint16_t> channel_data_vector = fifo_ptr;

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
	
	bool result = PulseRecoMgr.RecoPulse(fifo_ptr);

	ConstructHits(Channel,
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
		      Contributors1, Contributors2,
		      Binned1, Binned2,
		      FlashesInAccumulator1, FlashesInAccumulator2);

      }//end loop over FIFO channels in frame
  }//end ProcessFrame

  //-------------------------------------------------------------------------------------------------
  void ConstructHits(int const& Channel,
		     unsigned int const& Frame,
		     pmtana::AlgoThreshold const& ThreshAlg,
		     std::vector<recob::OpHit>& HitVector,
		     optdata::TimeSlice_t const& TimeSlicesPerFrame,
		     int const& BinWidth,
		     float const& HitThreshold,
		     float const& FlashThreshold,
		     unsigned int const& TrigFrame,
		     unsigned int const& TrigTime,
		     double const& SPESize,
		     std::vector<double> const& Contributors1,
		     std::vector<double> const& Contributors1,
		     std::vector< std::vector<int> > const& Binned1,
		     std::vector< std::vector<int> > const& Binned2,
		     std::vector<int> const& FlashesInAccumulator1,
		     std::vector<int> const& FlashesInAccumulator2)
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
			      0);
      
      size_t HitIndex = OpHitsThisFrame.size()-1;
      
      size_t Accum1Index = int(TimeSlice  / BinWidth);
      size_t Accum2Index = int((TimeSlice + BinWidth/2)/BinWidth);
      
      Contributors1.at(Accum1Index).push_back(HitIndex);
      Contributors2.at(Accum2Index).push_back(HitIndex);
      
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
  
}
