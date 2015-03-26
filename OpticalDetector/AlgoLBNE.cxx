//===================================
// AlgoLBNE.cxx
//===================================

#include "AlgoLBNE.h"

//#include <iostream>

#ifndef ALGOLBNE_CXX
#define ALGOLBNE_CXX

namespace pmtana {

  //---------------------------------------------------------------------
  AlgoLBNE::AlgoLBNE()
  {

    _adc_thres = 15;

    _nsigma = 5;

    Reset();

    // Creating single photoelectron waveform
    _voltage_to_adc = 151.5;
    _sample_freq    = 0.15;
    _pulse_length   = 400.0;
    _offset         = _pulse_length/4;
    _max_amplitude  = 0.12*_voltage_to_adc;
    _front_time     = 9.0*_sample_freq;
    _back_time      = 476.0*_sample_freq;
    CreateStandardPulse();

  }

  //---------------------------------------------------------------------
  AlgoLBNE::~AlgoLBNE()
  {}

  //---------------------------------------------------------------------
  void AlgoLBNE::Reset()
  {

    PMTPulseRecoBase::Reset();

  }

  //---------------------------------------------------------------------
  bool AlgoLBNE::RecoPulse(const std::vector< uint16_t > &wf)
  {
    
    double threshold    = (_adc_thres > (_nsigma*_ped_rms) ? _adc_thres : (_nsigma*_ped_rms));
    int    counter      = 0;
    double abs_max      = 0.0;
    int    t_max        = 0;
    bool   subtractable = true;

    Reset();

    std::vector< double > waveform(wf.begin(), wf.end());

    while(subtractable)
    {
      counter = 0;
      t_max   = 0;
      abs_max = 0.0;

      for (auto const &value : waveform)
      {
      
        if (value > abs_max) 
        {
          abs_max = value;
          t_max   = counter;
        }

        counter++;
      }

      if (abs_max > threshold) 
      {
        //std::cout << "Waveform" << "\n";
        SubtractPulse(waveform, t_max);
        //std::cout << abs_max << "\n";

        _pulse.t_start = t_max - _offset;
        _pulse.t_end   = t_max + _pulse_length - _offset;
        _pulse.peak    = abs_max;
        _pulse.area    = 0;
        _pulse.t_max   = t_max;

        _pulse_v.push_back(_pulse);

      }
      else subtractable = false;
    }

    return true;

  }

  //---------------------------------------------------------------------
  double AlgoLBNE::Pulse1PE(double time)
  {

    if (time < 0.0) return (_max_amplitude*std::exp(time/_front_time));
    else return (_max_amplitude*std::exp(-(time)/_back_time));

  }

  //---------------------------------------------------------------------
  void AlgoLBNE::CreateStandardPulse()
  {

    int length = _pulse_length;
    _standard_pulse.resize(length);
    for (size_t tick = 0; tick != _standard_pulse.size(); tick++)
    {
      _standard_pulse.at(tick) = Pulse1PE(double(tick) - _offset);
    }

  }

  //---------------------------------------------------------------------
  void AlgoLBNE::SubtractPulse(std::vector< double > &waveform, int tick)
  {

    size_t pulseStart = tick - _offset;
    if (pulseStart < 0) pulseStart = 0;
    size_t pulseEnd   = tick + _pulse_length - _offset;
    if (pulseEnd > waveform.size()) pulseEnd = waveform.size();
    
    for (size_t i = pulseStart; i < pulseEnd; i++)
    {
      if (pulseStart) waveform.at(i) -= _standard_pulse.at(i - pulseStart);
      else            waveform.at(i) -= _standard_pulse.at(i - tick + _offset);
//      std::cout << waveform.at(i) << "\n";
    }

  }

}

#endif
