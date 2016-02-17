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

#include "RawData/OpDetWaveform.h"
#include "OpticalDetector/PulseRecoManager.h"
#include "OpticalDetector/PMTPulseRecoBase.h"
#include "RecoBase/OpHit.h"
#include "Geometry/Geometry.h"
#include "Utilities/TimeService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>

namespace opdet{

  void RunHitFinder(std::vector< raw::OpDetWaveform > const&,
                    std::vector< recob::OpHit >&,
                    pmtana::PulseRecoManager const&,
                    pmtana::PMTPulseRecoBase const&,
                    geo::Geometry const&,
                    float const&,
                    util::TimeService const&,
                    std::vector< double > const&,
                    bool const&);

  void ConstructHit(float const&, 
                    int const&,
                    double const&,
                    pmtana::pulse_param const&,
                    util::TimeService const&,
                    double const&,
                    bool const&,
                    std::vector< recob::OpHit >&);

} // End opdet namespace

#endif
