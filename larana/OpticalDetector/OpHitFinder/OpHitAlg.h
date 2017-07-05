// -*- mode: c++; c-basic-offset: 2; -*-
#ifndef OPHITALG_H
#define OPHITALG_H
/*!
 * Title:   OpHit Algorithims
 * Author:  Ben Jones, MIT ( Edited by wketchum@lanl.gov, gleb.sinev@duke.edu
 *                           and kevin.wood@stonybrook.edu )
 *
 * Description:
 * These are the algorithms used by OpHit to produce optical hits.
 */

#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfo/DetectorClocks.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>

namespace opdet{

  void RunHitFinder(std::vector< raw::OpDetWaveform > const&,
                    std::vector< recob::OpHit >&,
                    pmtana::PulseRecoManager const&,
                    pmtana::PMTPulseRecoBase const&,
                    geo::GeometryCore const&,
                    float,
                    detinfo::DetectorClocks const&,
                    std::vector< double > const&,
                    bool,
		    std::vector< double > const&);

  // For backwards compatibility
  void RunHitFinder(std::vector< raw::OpDetWaveform > const&,
                    std::vector< recob::OpHit >&,
                    pmtana::PulseRecoManager const&,
                    pmtana::PMTPulseRecoBase const&,
                    geo::GeometryCore const&,
                    float,
                    detinfo::DetectorClocks const&,
                    std::vector< double > const&,
                    bool);

  void ConstructHit(float, 
                    int,
                    double,
                    pmtana::pulse_param const&,
                    detinfo::DetectorClocks const&,
                    double,
                    bool,
                    std::vector< recob::OpHit >&,
		    double);

} // End opdet namespace

#endif
