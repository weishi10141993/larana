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
#include "RecoBase/OpHit.h"
#include "RecoBase/OpFlash.h"

namespace opdet{

  void GetTriggerTime(std::vector<const sim::BeamGateInfo*> const&,
		      double const&,
		      double const&,
		      optdata::TimeSlice_t const&,
		      unsigned int&, unsigned short&);

  void RunFlashFinder(std::vector<optdata::FIFOChannel> const&,
		      std::vector<recob::OpHit>&,
		      std::vector<recob::OpFlash>&,
		      std::vector< std::vector<int> >&);
  
}//end opdet namespace

#endif
