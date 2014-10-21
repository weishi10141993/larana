#ifndef TRACKCALORIMETRYALG_H
#define TRACKCALORIMETRYALG_H
/*!
 * Title:   Track Calorimetry Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on code the Calorimetry_module
 *
 * Description: Algorithm that produces a calorimetry object given a track
 * Input:       recob::Track, Assn<recob::Spacepoint,recob::Track>, Assn<recob::Hit,recob::Track>
 * Output:      anab::Calorimetry, (and Assn<anab::Calorimetry,recob::Track>) 
*/
#include <iostream>

#include "fhiclcpp/ParameterSet.h"

#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Track.h"
#include "Filters/ChannelFilter.h"
#include "AnalysisBase/Calorimetry.h"

#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

#include "TTree.h"

namespace calo{
  class TrackCalorimetryAlg;
}

class calo::TrackCalorimetryAlg{
 public:
  TrackCalorimetryAlg(fhicl::ParameterSet const& p);
  void reconfigure(fhicl::ParameterSet const& p);

  void ExtractCalorimetry(std::vector<recob::Track> const&,
			  std::vector<recob::Hit> const&,
			  std::vector< std::vector<size_t> > const&,
			  std::vector<recob::SpacePoint> const&,
			  std::vector< std::vector<size_t> > const&,
			  filter::ChannelFilter const&,
			  std::vector<anab::Calorimetry>&,
			  std::vector<size_t>&,
			  geo::Geometry const&,
			  util::LArProperties const&,
			  util::DetectorProperties const&);

 private:


};

#endif
