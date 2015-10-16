////////////////////////////////////////////////////////////////////////
//
//  AlgoSlidingWindow source
//
////////////////////////////////////////////////////////////////////////

#ifndef ALGOSLIDINGWINDOW_CXX
#define ALGOSLIDINGWINDOW_CXX

#include "AlgoSlidingWindow.h"

namespace pmtana{

  //############################
  AlgoSlidingWindow::AlgoSlidingWindow(const fhicl::ParameterSet &pset)
  //############################
  {

    _adc_thres = pset.get<float>("ADCThreshold");
  
    _nsigma = pset.get<float>("NSigmaThreshold");

    _min_wf_size = pset.get<size_t>("MinWindowSize");
    
    _max_sigma = pset.get<float>("MaxSigma");

    Reset();

  }

  //***************************************************************
  AlgoSlidingWindow::~AlgoSlidingWindow()
  //***************************************************************
  {}

  //***************************************************************
  void AlgoSlidingWindow::Reset()
  //***************************************************************
  {
    PMTPulseRecoBase::Reset();
  }

  //**********************************************************************
  bool AlgoSlidingWindow::ConstructPedestal(const std::vector<short> &wf)
  //**********************************************************************
  {
    // parameters
    if(wf.size()<_min_wf_size) return false;
    
    //
    // Compute pedestal by itself
    //
    _local_mean.resize(wf.size(),0);
    _local_sigma.resize(wf.size(),100);

    for(size_t i=0; i<wf.size(); ++i) {

      if( (i+_min_wf_size) == wf.size() ){
	float last_mean = _local_mean.at(i-1);
	float last_sigma = _local_sigma.at(i-1);
	for(size_t j=i; j<wf.size(); ++j) {
	  _local_mean.at(j) = last_mean;
	  _local_sigma.at(j) = last_sigma;
	}
	break;
      }

      double mean = 0;
      double sigma = 0;
      for(size_t j=0; j<_min_wf_size; ++j) mean += wf.at(i+j);

      mean /= ((float)_min_wf_size);

      for(size_t j=0; j<_min_wf_size; ++j) sigma += pow( (wf.at(i+j) - mean), 2 );

      sigma = sqrt(sigma/((float)(_min_wf_size)));

      _local_mean.at(i)=mean;
      _local_sigma.at(i)=sigma;
    }

    float min_sigma = 1e9;
    size_t min_sigma_index = 0;
    size_t num_good_sigma = 0;
    for(size_t i=0; i<_local_sigma.size(); ++i) {
      auto const& sigma = _local_sigma.at(i);
      if(sigma<min_sigma) {
	min_sigma = sigma;
	min_sigma_index = i;
      }
      if(sigma < _max_sigma) num_good_sigma += 1;
    }
    
    // If no good mean, use the best guess within this waveform
    if(min_sigma > _max_sigma || num_good_sigma<3) {
      for(size_t i=0; i<_local_mean.size(); ++i) {
	_local_mean.at(i) = _local_mean.at(min_sigma_index);
	_local_sigma.at(i) = _local_sigma.at(min_sigma_index);
      }
      //std::cout<<"Only "<<num_good_sigma<<" with min = "<<min_sigma<<std::endl;
      return true;
    }

    // Else do extrapolation
    int last_good_index=-1;
    //size_t good_ctr = 0;
    for(size_t i=0; i<wf.size(); ++i) {
      if(_local_sigma.at(i) <= _max_sigma) {

	if(last_good_index<0) {
	  last_good_index = (int)i;
	  continue;
	}

	if( (last_good_index+1) < (int)i ) {

	  float slope = (_local_mean.at(i) - _local_mean.at(last_good_index)) / (float(i - last_good_index));

	  for(size_t j=last_good_index+1; j<i; ++j) {
	    _local_mean.at(j) = slope * (float(j-last_good_index)) + _local_mean.at(last_good_index);
	    _local_sigma.at(j) = _max_sigma;
	  }
	}
	last_good_index = (int)i;
      }
    }
    //std::cout<<"Used "<<good_ctr<<std::endl;
    // Next do extrapolation to the first and end
    if(_local_sigma.front() > _max_sigma) {
      int first_index=-1;
      int second_index=-1;
      for(size_t i=0; i<wf.size(); ++i) {
	if(_local_sigma.at(i)<_max_sigma) {
	  if(first_index<0) first_index = (int)i;
	  else if(second_index<0) {
	    second_index = (int)i;
	    break;
	  }
	}
      }
      if(first_index<0 || second_index<0) throw std::exception();

      float slope = (_local_mean.at(second_index) - _local_mean.at(first_index)) / (float(second_index - first_index));
      for(int i=0; i<first_index; ++i) {
	_local_mean.at(i) = _local_mean.at(first_index) - slope * (first_index - i);
	_local_sigma.at(i) = _max_sigma;
      }
    }
    
    if(_local_sigma.back() > _max_sigma) {
      int first_index=-1;
      int second_index=-1;
      for(int i=wf.size()-1; i>=0; --i) {
	if(_local_sigma.at(i)<_max_sigma) {
	  if(second_index<0) second_index = (int)i;
	  else if(first_index<0) {
	    first_index = (int)i;
	    break;
	  }
	}
      }
      float slope = (_local_mean.at(second_index) - _local_mean.at(first_index)) / (float(second_index - first_index));
      for(int i=second_index+1; i<int(wf.size()); ++i) {
	_local_mean.at(i) = _local_mean.at(second_index) + slope * (i-second_index);
	_local_sigma.at(i) = _max_sigma;
      }
    }
    return true;
  }


  //***************************************************************
  bool AlgoSlidingWindow::RecoPulse(const std::vector<short> &wf)
  //***************************************************************
  {

    if(!ConstructPedestal(wf)) return false;

    bool fire = false;
    
    double counter=0;

    //double threshold = ( _adc_thres > (_nsigma * _ped_rms) ? _adc_thres : (_nsigma * _ped_rms) );

    //threshold += _ped_mean;

    Reset();

    for(size_t i=0; i<wf.size(); ++i) {

      auto const& value = wf[i];

      float threshold = _local_mean.at(i);
      if(_local_sigma.at(i) * _nsigma < _adc_thres) threshold += _adc_thres;
      else threshold += _local_sigma.at(i) * _nsigma;
      
      //std::cout << "SlidingWindow=" << threshold << ", value=" << value << ", counter=" << counter << std::endl;

      if( !fire && ((double)value) >= threshold ){

	// Found a new pulse

	fire = true;

	_pulse.t_start = counter;

      }
    
      if( fire && ((double)value) < threshold ){
      
	// Found the end of a pulse

	fire = false;

	_pulse.t_end = counter - 1;
      
	_pulse_v.push_back(_pulse);

	_pulse.reset_param();

      }


      //std::cout << "\tFire=" << fire << std::endl;

      if(fire){

	// Add this adc count to the integral

	_pulse.area += ((double)value - (double)_local_mean.at(i));

	if(_pulse.peak < ((double)value - (double)_local_mean.at(i))) {

	  // Found a new maximum
	  
	  _pulse.peak = ((double)value - (double)_local_mean.at(i));

	  _pulse.t_max = counter;

	}

      }
    
      counter++;
    }

    if(fire){

      // Take care of a pulse that did not finish within the readout window.
    
      fire = false;
    
      _pulse.t_end = counter - 1;
    
      _pulse_v.push_back(_pulse);
    
      _pulse.reset_param();

    }

    return true;

  }

}

#endif
