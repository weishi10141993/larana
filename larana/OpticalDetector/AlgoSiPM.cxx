//=======================
// AlgoSiPM.cxx
//=======================

#include "AlgoSiPM.h"

#ifndef ALGOSIPM_CXX
#define ALGOSIPM_CXX

namespace pmtana {

  //---------------------------------------------------------------------------
  AlgoSiPM::AlgoSiPM(const fhicl::ParameterSet &pset)
  {

    _adc_thres = pset.get< float >("ADCThreshold"   );
    _min_width = pset.get< float >("MinWidth"       );
    _2nd_thres = pset.get< float >("SecondThreshold");
    _pedestal  = pset.get< float >("Pedestal"       );

//    _nsigma = 5;

    Reset();

  }

  //---------------------------------------------------------------------------
  AlgoSiPM::~AlgoSiPM()
  {}

  //---------------------------------------------------------------------------
  void AlgoSiPM::Reset()
  {

    PMTPulseRecoBase::Reset();

  }

  //---------------------------------------------------------------------------
  bool AlgoSiPM::RecoPulse(const std::vector< short > &wf)
  {
    
    bool   fire          = false;
    bool   first_found   = false;
    bool   record_hit    = false;
    int    counter       = 0;
//    double threshold   = (_2nd_thres > (_nsigma*_ped_rms) ? _2nd_thres 
//                                                    : (_nsigma*_ped_rms));
    double pedestal      = _pedestal;
    double threshold     = _adc_thres;
    threshold           += pedestal;
    double pre_threshold = _2nd_thres;
    pre_threshold       += pedestal;

    Reset();

    for (short const &value : wf) {

      if (!fire && (double(value) >= pre_threshold)) {
        
        // Found a new pulse
        fire           = true;
        first_found    = false;
        record_hit     = false;
        _pulse.t_start = counter;

      }

      if (fire && (double(value) < pre_threshold)) {

        // Found the end of a pulse
        fire = false;
        _pulse.t_end = counter - 1;
        if (record_hit && ((_pulse.t_end - _pulse.t_start) >= _min_width)) 
        {
          _pulse_v.push_back(_pulse);
          record_hit = false;
        }
        _pulse.reset_param();

      }

      if (fire) {

        // We want to record the hit only if _adc_thres is reached
        if (!record_hit && (double(value) >= threshold)) record_hit = true;

        // Add this ADC count to the integral
        _pulse.area += (double(value) - double(pedestal));

        if (!first_found && 
            (_pulse.peak < (double(value) - double(pedestal)))) {

          // Found a new maximum
          _pulse.peak  = (double(value) - double(pedestal));
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
      if (record_hit && ((_pulse.t_end - _pulse.t_start) >= _min_width)) 
      {
        _pulse_v.push_back(_pulse);
        record_hit = false;
      }
      _pulse.reset_param();

    }

    return true;

  }

}

#endif
