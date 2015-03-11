//============================================
// AlgoLBNE.h
// This is a hit finding algorithm 
// for LBNE photon detectors.
//
// Gleb Sinev, Duke, 2015
// Based on AlgoTreshold.h
//============================================

#ifndef ALGOLBNE_H
#define ALGOLBNE_H

#include "PMTPulseRecoBase.h"

#include <vector>

namespace pmtana {

  class AlgoLBNE : public PMTPulseRecoBase {

    public:

      // Default constructor
      AlgoLBNE();

      // Default destructor
      virtual ~AlgoLBNE();

      // Implementation of PMTPulseRecoBase::RecoPulse() method
      virtual bool RecoPulse(const std::vector< uint16_t > &wf);

      // Implementation of PMTPulseRecoBase::Reset() method
      virtual void Reset();

      // A method to set user-defined ADC threshold value
      void SetADCThreshold(double v) {_adc_thres = v;};

      // A method to set a multiplication factor to the pedestal standard deviation
      // which is used as one of two input values to define a threshold
      void SetNSigma(double v) {_nsigma = v;};

    protected:

      // A variable holder for a user-defined absolute ADC threshold value
      double _adc_thres;

      // A variable holder for a multiplicative factor for the pedestal standard deviation to difine the threshold
      double _nsigma;

      // Functional response to one photoelectron (time in ticks)
      double Pulse1PE(double time);

      // A vector containing the single photoelectron shape
      std::vector< double > _standard_pulse;

      // A function to initialize the standard pulse vector
      void CreateStandardPulse();

      // One photoelectron pulse parameters
      double _pulse_length;   // 1PE pulse length in ticks
      double _offset;         // How many ticks to the left of the peak we want cover
      double _max_amplitude;  // Maximum amplitude of the pulse in ADC counts
      double _front_time;     // Constant in the exponential function in ticks
      double _back_time;      // Constant in the exponential function in ticks
      double _sample_freq;    // Frequency in GHz (number of ticks in one ns)
      double _voltage_to_adc; // Conversion factor mV to ADC counts

      // A function to subtract one pulse from the waveform
      void SubtractPulse(std::vector< double > &waveform, int tick);

  };

}

#endif
