////////////////////////////////////////////////////////////////////////
//
//  AlgoPedestal source
//
////////////////////////////////////////////////////////////////////////

#ifndef ALGOPEDESTAL_CC
#define ALGOPEDESTAL_CC

#include "AlgoPedestal.h"

namespace pmtana{

  //#########################################################
  AlgoPedestal::AlgoPedestal(fhicl::ParameterSet const& /* pset */)
  //#########################################################
  {

    _mean  = -1;

    _sigma = -1;

  }

  //#########################################################
  AlgoPedestal::~AlgoPedestal()
  //#########################################################
  {}

  //*************************************************************************************************
  void AlgoPedestal::ComputePedestal(const std::vector<uint16_t>* wf, size_t start, size_t nsample)
  //*************************************************************************************************
  {  
    _mean  = -1;
    _sigma = -1;
  
    if( (start + nsample) > wf->size() ){
      mf::LogWarning(__PRETTY_FUNCTION__)
	<<Form("Wavelength too short (%zu ADC samples) to compute pedestal! (minimum %zu)",
	       wf->size(),(start + nsample));
      return;
    }

    for(size_t index=start; index < (start + nsample); ++index)

      _mean += wf->at(index);

    _mean = _mean / ((double)nsample);

    for(size_t index=0; index < (start+nsample); ++index)
    
      _sigma += pow( (wf->at(index) - _mean), 2 );

    _sigma = sqrt(_sigma/((double)(nsample)));

    return;

  }

}

#endif
