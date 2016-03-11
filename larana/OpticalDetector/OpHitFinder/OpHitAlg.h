// -*- mode: c++; c-basic-offset: 2; -*-
#ifndef OPHITALG_H
#define OPHITALG_H
/*!
 * Title:   OpHit Algorithims
 * Author:  Ben Jones, MIT (Edited by wketchum@lanl.gov and gleb.sinev@duke.edu)
 *
 * Description:
 * These are the algorithms used by OpHit to produce optical hits.
 */

#include "lardata/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "lardata/RecoBase/OpHit.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfo/DetectorClocks.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>

namespace opdet{

  void RunHitFinder(std::vector< raw::OpDetWaveform > const&,
                    std::vector< recob::OpHit >&,
                    pmtana::PulseRecoManager const&,
                    pmtana::PMTPulseRecoBase const&,
                    geo::GeometryCore const&,
                    float const&,
                    detinfo::DetectorClocks const&,
                    std::vector< double > const&,
                    bool const&);

  void ConstructHit(float const&, 
                    int const&,
                    double const&,
                    pmtana::pulse_param const&,
                    detinfo::DetectorClocks const&,
                    double const&,
                    bool const&,
                    std::vector< recob::OpHit >&);

} // End opdet namespace

#endif
