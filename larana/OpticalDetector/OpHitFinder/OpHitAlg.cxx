// -*- mode: c++; c-basic-offset: 2; -*-
/*!
 * Title:   OpHit Algorithims
 * Authors, editors:  Ben Jones, MIT
 *                    Wes Ketchum wketchum@lanl.gov
 *                    Gleb Sinev  gleb.sinev@duke.edu
 *                    Alex Himmel ahimmel@fnal.gov
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
                    float const&                    hitThreshold,
                    detinfo::DetectorClocks const&  detectorClocks,
                    std::vector< double > const&    SPESize,
                    bool const&                     areaToPE) {

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
                     hitVector);
      
    }

  }

  //----------------------------------------------------------------------------
  void ConstructHit(float const&                   hitThreshold,
                    int const&                     channel,
                    double const&                  timeStamp,
                    pmtana::pulse_param const&     pulse,
                    detinfo::DetectorClocks const& detectorClocks,
                    double const&                  SPESize,
                    bool const&                    areaToPE,
                    std::vector< recob::OpHit >&   hitVector) {

    if (pulse.peak < hitThreshold) return;

    double absTime = timeStamp 
                   + pulse.t_max*detectorClocks.OpticalClock().TickPeriod();

    double relTime = absTime - detectorClocks.BeamGateTime();
    if (detectorClocks.BeamGateTime() < 0.0) 
      relTime = absTime - detectorClocks.TriggerTime();
    
    int frame = detectorClocks.OpticalClock().Frame(timeStamp);

    double PE = 0.0;
    if (areaToPE) PE = pulse.area/SPESize;
    else          PE = pulse.peak/SPESize;
    
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
