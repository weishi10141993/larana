/*!
 * Title:   Track Calorimetry Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on code the Calorimetry_module
 *
 * Description: Algorithm that produces a calorimetry object given a track
 * Input:       recob::Track, Assn<recob::Spacepoint,recob::Track>, Assn<recob::Hit,recob::Track>
 * Output:      anab::Calorimetry, (and Assn<anab::Calorimetry,recob::Track>) 
*/

#include "TrackCalorimetryAlg.h"
#include <limits>

calo::TrackCalorimetryAlg::TrackCalorimetryAlg(fhicl::ParameterSet const& p)
{
  this->reconfigure(p);
}

void calo::TrackCalorimetryAlg::reconfigure(fhicl::ParameterSet const& p){
}

void calo::TrackCalorimetryAlg::ExtractCalorimetry(std::vector<recob::Track> const& trackVector,
						   std::vector<recob::Hit> const& hitVector,
						   std::vector< std::vector<size_t> > const& hit_indices_per_track,
						   std::vector<recob::SpacePoint> const& spptVector,
						   std::vector< std::vector<size_t> > const& sppt_indices_per_track,
						   filter::ChannelFilter const& chanFilt,
						   std::vector<anab::Calorimetry>& caloVector,
						   std::vector<size_t>& assnTrackCaloVector,
						   geo::Geometry const& geom,
						   util::LArProperties const& larp,
						   util::DetectorProperties const& detprop)
{

  for(size_t i_track=0; i_track<trackVector.size(); i_track++){
  }

}
