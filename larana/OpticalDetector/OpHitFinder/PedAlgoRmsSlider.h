/**
 * \file PedAlgoRmsSlider.h
 *
 * \ingroup PulseReco
 * 
 * \brief Class definition file of PedAlgoRmsSlider
 *
 * @author Marco - Oxford 2017
 */

/** \addtogroup PulseReco
    
@{*/

#ifndef larana_OPTICALDETECTOR_PEDALGORMSSLIDER_H
#define larana_OPTICALDETECTOR_PEDALGORMSSLIDER_H

// STL
#include "PMTPedestalBase.h"
#include "fhiclcpp/ParameterSet.h"

#include <fstream>

namespace pmtana
{

  /**
   \class PedAlgoRmsSlider
   A class that calculates pedestal mean & standard deviation (here and elsewhere called as "RMS").   
  */
  class PedAlgoRmsSlider : public PMTPedestalBase{

  public:

    /// Default constructor
    PedAlgoRmsSlider(const std::string name="PedRmsSlider");

    /// Alternative ctor
    PedAlgoRmsSlider(const fhicl::ParameterSet &pset,const std::string name="PedRmsSlider");

    /// Default destructor
    virtual ~PedAlgoRmsSlider();

    /// Print settings
    void PrintInfo();


  protected:

    /// Method to compute a pedestal of the input waveform using "nsample" ADC samples from "start" index.
    bool ComputePedestal( const pmtana::Waveform_t& wf,
			  pmtana::PedestalMean_t&   mean_v,
			  pmtana::PedestalSigma_t&  sigma_v);
    
  private:

    size_t _sample_size;  ///< How many samples are used to calculate local rms and mean
    double _threshold;    ///< Threshold applied to local rms to claim a pulse

    float _max_sigma;     ///< Max sigma to consider adc as 'sane'
    float _ped_range_max; ///< Max value of adc to consider adc as 'sane'
    float _ped_range_min; ///< Min value of adc to consider adc as 'sane'

    bool _verbose;        ///< For debugging
    int _n_wf_to_csvfile; ///< If greater than zero saves firsts waveforms with pedestal to csv file
    int _wf_saved = 0;
    std::ofstream _csvfile;    

    /// Returns the mean of the elements of the vector from start to start+nsample
    double CalcMean(const std::vector<double>& wf, size_t start, size_t nsample);

    /// Returns the std of the elements of the vector from start to start+nsample
    double CalcStd(const std::vector<double>& wf, const double ped_mean, size_t start, size_t nsample);

    /// Checks the sanity of the estimated pedestal, returns false if not sane
    bool CheckSanity(pmtana::PedestalMean_t& mean_v, pmtana::PedestalSigma_t& sigma_v);
  };
}
#endif

/** @} */ // end of doxygen group
