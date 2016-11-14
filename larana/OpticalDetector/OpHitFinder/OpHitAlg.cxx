// -*- mode: c++; c-basic-offset: 2; -*-
/*!
 * Title:   OpHit Algorithims
 * Authors, editors:  Ben Jones, MIT
 *                    Wes Ketchum wketchum@lanl.gov
 *                    Gleb Sinev  gleb.sinev@duke.edu
 *                    Alex Himmel ahimmel@fnal.gov
 *                    Kevin Wood  kevin.wood@stonybrook.edu
 *
 * Description:
 * These are the algorithms used by OpHitFinder to produce optical hits.
 */

#include "OpHitAlg.h"

namespace opdet{
  //----------------------------------------------------------------------------
  void RunHitFinder(std::vector< raw::OpDetWaveform > const& 
                                                    opDetWaveformVector,
                    std::vector< recob::OpHit >&    hitVector,
                    pmtana::PulseRecoManager const& pulseRecoMgr,
                    pmtana::PMTPulseRecoBase const& threshAlg,
                    geo::GeometryCore const&        geometry,
                    float                           hitThreshold,
                    detinfo::DetectorClocks const&  detectorClocks,
                    std::vector< double > const&    SPESize,
                    bool                            areaToPE,
		    std::vector< double > const&    SPEShiftPerChan ) {

    for (auto const& waveform : opDetWaveformVector) {

      const int channel = static_cast< int >(waveform.ChannelNumber());

      if (!geometry.IsValidOpChannel(channel)) {
        mf::LogError("OpHitFinder") << "Error! unrecognized channel number " 
                         << channel << ". Ignoring pulse";
        continue;
      }

      pulseRecoMgr.Reconstruct(waveform);
      
      // Get the result
      auto const& pulses = threshAlg.GetPulses();

      const double timeStamp = waveform.TimeStamp();


      for (auto const& pulse : pulses)
        ConstructHit(hitThreshold,
                     channel,
                     timeStamp,
                     pulse,
                     detectorClocks,
                     SPESize.at(channel),
                     areaToPE,
                     hitVector,
		     SPEShiftPerChan.at(channel) );
      
    }

  }
  //----------------------------------------------------------------------------
  // For backward compatibility
  void RunHitFinder(std::vector< raw::OpDetWaveform > const&
                                                    opDetWaveformVector,
                    std::vector< recob::OpHit >&    hitVector,
                    pmtana::PulseRecoManager const& pulseRecoMgr,
                    pmtana::PMTPulseRecoBase const& threshAlg,
                    geo::GeometryCore const&        geometry,
                    float                           hitThreshold,
                    detinfo::DetectorClocks const&  detectorClocks,
                    std::vector< double > const&    SPESize,
                    bool                            areaToPE) {

    // if no SPEShiftPerChan vec is given, use one with no shift for all channels
    std::vector< double > noSPEShift(geometry.NOpChannels() , 0. );

    RunHitFinder(opDetWaveformVector,
                 hitVector,
                 pulseRecoMgr,
                 threshAlg,
                 geometry,
                 hitThreshold,
                 detectorClocks,
                 SPESize,
                 areaToPE,
                 noSPEShift);

  }
  //----------------------------------------------------------------------------
  void ConstructHit(float                          hitThreshold,
                    int                            channel,
                    double                         timeStamp,
                    pmtana::pulse_param const&     pulse,
                    detinfo::DetectorClocks const& detectorClocks,
                    double                         SPESize,
                    bool                           areaToPE,
                    std::vector< recob::OpHit >&   hitVector,
		    double                         SPEShift=0.) {

    if (pulse.peak < hitThreshold) return;

    double absTime = timeStamp 
                   + pulse.t_max*detectorClocks.OpticalClock().TickPeriod();

    double relTime = absTime - detectorClocks.TriggerTime();
    
    int frame = detectorClocks.OpticalClock().Frame(timeStamp);

    double PE = 0.0;
    if (areaToPE) PE = pulse.area/SPESize + SPEShift;
    else          PE = pulse.peak/SPESize + SPEShift;
    
    double width = (pulse.t_end - pulse.t_start)
                     *detectorClocks.OpticalClock().TickPeriod();

    hitVector.emplace_back(channel,
                           relTime,
                           absTime,
                           frame,
                           width,
                           pulse.area,
                           pulse.peak,
                           PE,
                           0.0);

  }

} // End namespace opdet
