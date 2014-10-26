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

calo::TrackCalorimetryAlg::TrackCalorimetryAlg(fhicl::ParameterSet const& p):
  caloAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  this->reconfigure(p);
}

void calo::TrackCalorimetryAlg::reconfigure(fhicl::ParameterSet const& p){
  caloAlg.reconfigure(p.get<fhicl::ParameterSet>("CalorimetryAlg"));
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
						   util::DetectorProperties & detprop)
{

  //loop over the track list
  for(size_t i_track=0; i_track<trackVector.size(); i_track++){

    recob::Track const& track = trackVector[i_track];
    std::vector<float> path_length_fraction_vec(CreatePathLengthFractionVector(track));

    //sort hits into each plane
    std::vector< std::vector<size_t> > hit_indices_per_plane(geom.Nplanes());
    for(size_t i_plane=0; i_plane<geom.Nplanes(); i_plane++)
      for(auto const& i_hit : hit_indices_per_track[i_track])
	hit_indices_per_plane[hitVector[i_hit].WireID().Plane].push_back(i_hit);
    
    //loop over the planes
    for(size_t i_plane=0; i_plane<geom.Nplanes(); i_plane++){

      ClearInternalVectors();
      ReserveInternalVectors(hit_indices_per_plane[i_plane].size());

      //project down the track into wire/tick space for this plane
      std::vector< std::pair<geo::WireID,float> > traj_points_in_plane(track.NumberTrajectoryPoints());
      for(size_t i_trjpt=0; i_trjpt<track.NumberTrajectoryPoints(); i_trjpt++){
	double x_pos = track.LocationAtPoint(i_trjpt).X();
	float tick = detprop.ConvertXToTicks(x_pos,(int)i_plane,0,0);
	traj_points_in_plane[i_trjpt] = std::make_pair(geom.NearestWireID(track.LocationAtPoint(i_trjpt),i_plane),
						       tick);
      }
      

      //now loop through hits
      for(auto const& i_hit : hit_indices_per_plane[i_plane])
	AnalyzeHit(hitVector[i_hit],
		   track,
		   traj_points_in_plane,
		   path_length_fraction_vec,
		   geom);

    }//end loop over planes

  }//end loop over tracks

}//end ExtractCalorimetry


class dist_projected{
public:
  dist_projected(recob::Hit const& h, geo::Geometry const& g):
    hit(h), geom(g){}
  bool operator() (std::pair<geo::WireID,float> i, std::pair<geo::WireID,float> j)
  {
    float dw_i = ((int)(i.first.Wire) - (int)(hit.WireID().Wire))*geom.WirePitch(0,1,i.first.Plane);
    float dw_j = ((int)(j.first.Wire) - (int)(hit.WireID().Wire))*geom.WirePitch(0,1,j.first.Plane);
    float dt_i = i.second - hit.PeakTime();
    float dt_j = j.second - hit.PeakTime();

    return (std::sqrt(dw_i*dw_i + dt_i*dt_i) < std::sqrt(dw_j*dw_j + dt_j*dt_j));
  }
private:
  recob::Hit const& hit;
  geo::Geometry const& geom;

};

std::vector<float> calo::TrackCalorimetryAlg::CreatePathLengthFractionVector(recob::Track const& track){

  std::vector<float> trk_path_length_frac_vec(track.NumberTrajectoryPoints());

  float cumulative_path_length=0;
  const float total_path_length = track.Length();
  for(size_t i_trj=1; i_trj<track.NumberTrajectoryPoints()-1; i_trj++){    
    cumulative_path_length+=(track.LocationAtPoint(i_trj)-track.LocationAtPoint(i_trj-1)).Mag();
    trk_path_length_frac_vec[i_trj]=cumulative_path_length/total_path_length;
  }

  return trk_path_length_frac_vec;
}

void calo::TrackCalorimetryAlg::AnalyzeHit(recob::Hit const& hit,
					   recob::Track const& track,
					   std::vector< std::pair<geo::WireID,float> > const& traj_points_in_plane,
					   std::vector<float> const& path_length_fraction_vec,
					   geo::Geometry const& geom){

  size_t traj_iter = std::distance(traj_points_in_plane.begin(),
				   std::min_element(traj_points_in_plane.begin(),
						    traj_points_in_plane.end(),
						    dist_projected(hit,geom)));
  float pitch = track.PitchInView(geom.View(hit.WireID().Plane),traj_iter);
  fHitPropertiesMultiset.emplace(hit.Charge(false),
				 hit.Charge(false)/pitch,
				 caloAlg.dEdx_AREA(hit,pitch),
				 pitch,
				 track.LocationAtPoint(traj_iter),
				 path_length_fraction_vec[traj_iter]);
  

}
