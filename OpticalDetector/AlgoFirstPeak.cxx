//=======================
// AlgoFirstPeak.cxx
//=======================

#include "AlgoFirstPeak.h"

#ifndef ALGOFIRSTPEAK_CXX
#define ALGOFIRSTPEAK_CXX

namespace pmtana {

  //---------------------------------------------------------------------------
  AlgoFirstPeak::AlgoFirstPeak()
  {

    //_adc_thres = 3;
    //_adc_thres = 13;
    _adc_thres = 5;
    _min_width = 10;

    _nsigma = 5;

    Reset();

  }

  //---------------------------------------------------------------------------
  AlgoFirstPeak::~AlgoFirstPeak()
  {}

  //---------------------------------------------------------------------------
  void AlgoFirstPeak::Reset()
  {

    PMTPulseRecoBase::Reset();

  }

  //---------------------------------------------------------------------------
  bool AlgoFirstPeak::RecoPulse(const std::vector< short > &wf)
  {
    
    bool   fire        = false;
    bool   first_found = false;
    int    counter     = 0;
//    double threshold   = (_adc_thres > (_nsigma*_ped_rms) ? _adc_thres 
//                                                    : (_nsigma*_ped_rms));
    double threshold = _adc_thres;
    threshold += _ped_mean;

    Reset();

    for (short const &value : wf) {

      if (!fire && (double(value) >= threshold)) {
        
        // Found a new pulse
        fire           = true;
        first_found    = false;
        _pulse.t_start = counter;

      }

      if (fire && (double(value) < threshold)) {

        // Found the end of a pulse
        fire = false;
        _pulse.t_end = counter - 1;
        if ((_pulse.t_end - _pulse.t_start) >= _min_width) 
          _pulse_v.push_back(_pulse);
        _pulse.reset_param();

      }

      if (fire) {

        // Add this ADC count to the integral
        _pulse.area += (double(value) - double(_ped_mean));

        if (!first_found && 
            (_pulse.peak < (double(value) - double(_ped_mean)))) {

          // Found a new maximum
          _pulse.peak  = (double(value) - double(_ped_mean));
          _pulse.t_max = counter;

        }
        else if (!first_found)
          // Found the first peak
          first_found = true;

      }

      counter++;
    
    }

    if (fire) {

      // Take care of a pulse that did not finish within the readout window
      fire = false;
      _pulse.t_end = counter - 1;
      if ((_pulse.t_end - _pulse.t_start) >= _min_width) 
        _pulse_v.push_back(_pulse);
      _pulse.reset_param();

    }

    return true;

  }

}

#endif
