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

#include "AnalysisAlg/CalorimetryAlg.h"

#include "TTree.h"
#include "TVector3.h"

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
			  util::DetectorProperties &);

 private:

  CalorimetryAlg caloAlg;

  std::vector<float>    fQVector;
  std::vector<float>    fdQdxVector;
  std::vector<float>    fdEdxVector;
  std::vector<float>    fPitchVector;
  std::vector<TVector3> fXYZVector;
  void ClearInternalVectors()
  { fQVector.clear(); fdQdxVector.clear(); fdEdxVector.clear(); fPitchVector.clear(); fXYZVector.clear(); }
  void ReserveInternalVectors(size_t s)
  { fQVector.reserve(s); fdQdxVector.reserve(s); fdEdxVector.reserve(s); 
    fPitchVector.reserve(s); fXYZVector.reserve(s); }


  void AnalyzeHit(recob::Hit const&,
		  recob::Track const&,
		  std::vector< std::pair<geo::WireID,float> > const&, 
		  geo::Geometry const&);

};

#endif
