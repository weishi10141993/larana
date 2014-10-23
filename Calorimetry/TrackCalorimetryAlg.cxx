/*!
 * Title:   Track Calorimetry Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on code in the Calorimetry_module
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

  //loop over the track list
  for(size_t i_track=0; i_track<trackVector.size(); i_track++){

    recob::Track const& track = trackVector[i_track];

    //sort hits into each plane
    std::vector< std::vector<size_t> > hit_indices_per_plane(geom.Nplanes());
    for(size_t i_plane=0; i_plane<geom.Nplanes(); i_plane++)
      for(auto const& i_hit : hit_indices_per_track[i_track])
	hit_indices_per_plane[hitVector[i_hit].WireID().Plane].push_back(i_hit);
    
    //loop over the planes
    for(size_t i_plane=0; i_plane<geom.Nplanes(); i_plane++){

      ClearInternalVectors();

      //project down the track into wire/tick space for this plane
      std::vector< std::pair<geo::WireID,float> > traj_points_in_plane(track.NumberTrajectoryPoints());
      for(size_t i_trjpt=0; i_trjpt<track.NumberTrajectoryPoints(); i_trjpt++)
	traj_points_in_plane.push_back(std::make_pair(geom.NearestWireID(track.LocationAtPoint(i_trjpt),i_plane),
						      detprop.ConvertXToTicks((double)track.LocationAtPoint(i_trjpt).X(),i_plane,0,0) ));

      

      //now loop through hits
      for(auto const& i_hit : hit_indices_per_plane[i_plane])
	AnalyzeHit(hitVector[i_hit],
		   traj_points_in_plane,
		   geom);

    }//end loop over planes

  }//end loop over tracks

}//end ExtractCalorimetry


void calo::TrackCalorimetryAlg::AnalyzeHit(recob::Hit const& hit,
					   std::vector< std::pair<geo::WireID,float> > const& traj_points_in_plane,
					   geo::Geometry const& geom){
}
