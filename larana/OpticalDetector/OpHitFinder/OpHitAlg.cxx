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
                                                    OpDetWaveformVector,
                    std::vector< recob::OpHit >&    HitVector,
                    pmtana::PulseRecoManager const& PulseRecoMgr,
                    pmtana::PMTPulseRecoBase const& ThreshAlg,
                    geo::Geometry const&            geom,
                    float const&                    HitThreshold,
                    detinfo::DetectorClocks const&  ts,
                    std::vector< double > const&    SPESize,
                    bool const&                     AreaToPE) {

    for (auto const& wf_ptr : OpDetWaveformVector) {

      const int    Channel   = static_cast< int >(wf_ptr.ChannelNumber());
      const double TimeStamp = wf_ptr.TimeStamp();

      if (!geom.IsValidOpChannel(Channel)) {
        mf::LogError("OpHitFinder") << "Error! unrecognized channel number " 
                         << Channel << ". Ignoring pulse";
        continue;
      }
      
      PulseRecoMgr.Reconstruct(wf_ptr);
      
      // Get the result
      auto const& pulses = ThreshAlg.GetPulses();

      for (auto const& pulse : pulses)
        ConstructHit(HitThreshold,
                     Channel,
                     TimeStamp,
                     pulse,
                     ts,
                     SPESize.at(Channel),
                     AreaToPE,
                     HitVector);
      
    }

  }

  //----------------------------------------------------------------------------
  void ConstructHit(float const&                   HitThreshold,
                    int const&                     Channel,
                    double const&                  TimeStamp,
                    pmtana::pulse_param const&     pulse,
                    detinfo::DetectorClocks const& ts,
                    double const&                  SPESize,
                    bool const&                    AreaToPE,
                    std::vector< recob::OpHit >&   HitVector) {

    if (pulse.peak < HitThreshold) return;

    double AbsTime = TimeStamp + pulse.t_max*ts.OpticalClock().TickPeriod();

    double RelTime = AbsTime - ts.BeamGateTime();
    if (ts.BeamGateTime() < 0.0) RelTime = AbsTime - ts.TriggerTime();
    
    int Frame = ts.OpticalClock().Frame(TimeStamp);

    double PE = 0.0;
    if (AreaToPE) PE = pulse.area/SPESize;
    else          PE = pulse.peak/SPESize;
    
    double width = (pulse.t_end - pulse.t_start)*ts.OpticalClock().TickPeriod();

    HitVector.emplace_back(Channel,
                           RelTime,
                           AbsTime,
                           Frame,
                           width,
                           pulse.area,
                           pulse.peak,
                           PE,
                           0.0);

  }
    
} // End namespace opdet
