#ifndef OPFLASHALG_H
#define OPFLASHALG_H
/*!
 * Title:   OpFlash Algorithims
 * Author:  Ben Jones, MIT (Edited by wketchum@lanl.gov)
 *
 * Description:
 * These are the algorithms used by OpFlashFinder to produce flashes.
 */

#include <functional>
#include "Simulation/BeamGateInfo.h"
#include "OpticalDetectorData/OpticalTypes.h"
#include "OpticalDetectorData/FIFOChannel.h"
#include "OpticalDetector/AlgoThreshold.h"
#include "OpticalDetector/PulseRecoManager.h"
#include "OpticalDetector/PMTPulseRecoBase.h"
#include "RecoBase/OpHit.h"
#include "RecoBase/OpFlash.h"
#include "Geometry/Geometry.h"

namespace opdet{

  void GetTriggerTime(std::vector<const sim::BeamGateInfo*> const&,
		      double const&,
		      double const&,
		      optdata::TimeSlice_t const&,
		      unsigned int&, double&);

  void RunFlashFinder(std::vector<optdata::FIFOChannel> const&,
		      std::vector<recob::OpHit>&,
		      std::vector<recob::OpFlash>&,
		      std::vector< std::vector<int> >&,
		      optdata::TimeSlice_t const&,
		      int const&,
		      pmtana::PulseRecoManager const&,
		      pmtana::AlgoThreshold const&,
		      std::map<int,int> const&,
		      geo::Geometry const&,
		      float const&,
		      float const&,
		      float const&,
		      unsigned int const&,
		      double const&,
		      double const&,
		      std::vector<double> const&,
		      float const&);
  
  void ProcessFrame(unsigned short,
		    std::vector<const optdata::FIFOChannel*> const&,
		    std::vector<recob::OpHit>&,
		    std::vector<recob::OpFlash>&,
		    std::vector< std::vector<int> >&,
		    optdata::TimeSlice_t const&,
		    int const&,
		    pmtana::PulseRecoManager const&,
		    pmtana::AlgoThreshold const&,
		    std::map<int,int> const&,
		    geo::Geometry const&,
		    float const&,
		    float const&,
		    float const&,
		    unsigned int const&,
		    double const&,
		    double const&,
		    std::vector<double> const&,
		    float const&);

  void ConstructHits(int const&,
		     uint32_t const&,
		     unsigned short const&,
		     pmtana::AlgoThreshold const&,
		     std::vector<recob::OpHit>&,
		     optdata::TimeSlice_t const&,
		     int const&,
		     float const&,
		     float const&,
		     unsigned int const&,
		     double const&,
		     double const&,
		     double const&,
		     std::vector<double> &,
		     std::vector<double> &,
		     std::vector< std::vector<int> > &,
		     std::vector< std::vector<int> > &,
		     std::vector<int> &,
		     std::vector<int> &);


  void ConstructHit( float const&, 
		     int const&,
		     uint32_t const&,
		     unsigned short const&,
		     const pmtana::pulse_param*,
		     optdata::TimeSlice_t const&,
		     double const&,
		     double const&,
		     double const&,
		     std::vector<recob::OpHit>&);

  unsigned int GetAccumIndex(double const& TMax, 
			     uint32_t const& TimeSlice, 
			     int const& BinWidth, 
			     double const& BinOffset);

  void FillAccumulator(unsigned int const& AccumIndex,
		       unsigned int const& HitIndex,
		       double const& PE,
		       float const& FlashThreshold,
		       std::vector<double> & Binned,
		       std::vector< std::vector<int> > & Contributors,
		       std::vector<int> & FlashesInAccumulator);

  void AssignHitsToFlash( std::vector<int> const&,
			  std::vector<int> const&,
			  std::vector<double> const&,
			  std::vector<double> const&,
			  std::vector< std::vector<int> > const&,
			  std::vector< std::vector<int> > const&,
			  size_t const&,
			  std::vector<recob::OpHit> const&,
			  std::vector< std::vector<int> >&,
			  float const&);

  void FillFlashesBySizeMap(std::vector<int> const& FlashesInAccumulator,
			    std::vector<double> const& BinnedPE,
			    int const& Accumulator,
			    std::map<double, std::map<int,std::vector<int> > > & FlashesBySize);

  void RefineHitsToFlash(std::vector< std::vector<int> > const&,
			 std::vector<recob::OpHit> const&,
			 std::vector< std::vector<int> >&,
			 float const&,
			 float const&);

  void ConstructFlashes(std::vector< std::vector<int> > const&,
			std::vector<recob::OpHit> const&,
			std::vector<recob::OpFlash>&,
			uint32_t const&,
			geo::Geometry const&,
			unsigned int const&,
			unsigned short const&,
			float const&);

  void RemoveLateLight(std::vector<recob::OpFlash>&,
		       std::vector< std::vector<int> >&);

}//end opdet namespace

#endif
