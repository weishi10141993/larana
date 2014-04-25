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

//-------------------------------------------------------------------------------------------------
void opdet::GetTriggerTime(std::vector<const sim::BeamGateInfo*> const& beamGateArray,
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
		    std::vector< std::vector<int> >& AssocList)
{
}
