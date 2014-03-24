//Analyzer for flash<--->track matching (testing against MC)
#include "OpticalRecoAna.h"

// LArSoft includes
#include "OpticalDetector/OpDigiProperties.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opreco {
  DEFINE_ART_MODULE(OpticalRecoAna)
}


// ART includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// ROOT includes
#include "TTree.h"

// C++ Includes
#include <iostream>
#include <sstream>

  // Constructor
opreco::OpticalRecoAna::OpticalRecoAna(fhicl::ParameterSet const& pset): EDAnalyzer(pset)
{
  
  // Indicate that the Input Module comes from .fcl
  fFlashModuleLabel  = pset.get<std::string>("FlashModule");
  fTrackModuleLabel  = pset.get<std::string>("TrackModule");
  fKineticEnergyMin  = pset.get<float>("KineticEnergyMin");
  fPEMin             = pset.get<float>("PEMin");
  fTimeMatchMax      = pset.get<float>("TimeMatchMax");
  
  fFlash_match_vector.clear();
  fTrack_match_vector.clear();
  fParticle_match_vector.clear();
}



// Destructor
opreco::OpticalRecoAna::~OpticalRecoAna() 
{}



// Do something here to setup the file (like make a TTree)
void opreco::OpticalRecoAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTimeDiff = tfs->make<TH1F>("htdiff","Time difference between particles and flashes; t_diff (ns); flash/particle pairs",1e3,-5e6,5e6);
  fTimeDiff_fine = 
    tfs->make<TH1F>("htdiff_fine","Time difference between particles and flashes; t_diff (ns); flash/particle pairs",100,-1000,1000);
}



// The analyzer itself
void opreco::OpticalRecoAna::analyze(const art::Event& evt)
{
  
  //get flag for MC
  const bool is_MC= !(evt.isRealData());
  
  //first get out track, flash, and particle list handles
  art::Handle< std::vector<recob::OpFlash> > flash_handle;
  evt.getByLabel(fFlashModuleLabel, flash_handle);
  std::vector<recob::OpFlash> const& flash_vector(*flash_handle);
  fFlash_match_vector.resize(flash_vector.size());
  
  LOG_INFO ("OpticalRecoAna")  
    << "Number of flashes is " << flash_vector.size() << std::flush;
  
  art::Handle< std::vector<recob::Track> > track_handle;
  evt.getByLabel(fTrackModuleLabel, track_handle);
  std::vector<recob::Track> const& track_vector(*track_handle);
  fTrack_match_vector.resize(track_vector.size());
  
  LOG_INFO ("OpticalRecoAna")  
    << "Number of tracks is " << track_vector.size() << std::flush;
  
  match_flashes_to_tracks(flash_vector, track_vector);
  
  //all this for the MC matching
  if(is_MC){
    
    art::ServiceHandle<cheat::BackTracker> bt;
    std::vector<simb::MCParticle> particle_list;
    get_MC_particle_list(bt->ParticleList(),particle_list);
    fParticle_match_vector.resize(particle_list.size());
    
    mf::LogInfo("OpticalRecoAna")  
      << "Number of MC particles is " << particle_list.size();
    

    art::ServiceHandle<opdet::OpDigiProperties> odp;
    const float ns_per_PMT_tick = ( 1e3 / odp->SampleFreq()) ; //SampleFreq is in MHz
    art::ServiceHandle<geo::Geometry> geometry_handle;

    match_flashes_to_particles(flash_vector,
			       particle_list,
			       ns_per_PMT_tick,
			       *geometry_handle);
    //check_flash_matches();

    //match_tracks_to_particles(track_handle,particle_list);
    
  }
  
}



void opreco::OpticalRecoAna::get_MC_particle_list(sim::ParticleList const& plist,std::vector<simb::MCParticle> & particle_vector) {
  
  for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
    
    const simb::MCParticle* particle = (*ipart).second;
    
    //do not store if it's below our energy cut
    if( particle->E() < 0.001*particle->Mass() +  fKineticEnergyMin ) continue;
    
    //check to see if it's a charged particle we expect to leave ionization
    const int pdg = std::abs(particle->PdgCode());
    if( pdg==11 // electron
	|| pdg==13 //muon
	|| pdg==15 //tau
	|| pdg==211 //charged pion
	|| pdg==321 //charged kaon
	|| pdg==2212 //proton
	|| pdg==2214 //delta
	|| pdg==2224 //delta
	|| pdg==3222 //sigma
	|| pdg==3112 //sigma
	|| pdg==3312 //xi
	|| pdg==3334 //omega
	) {
      particle_vector.emplace_back(*particle);
      //std::cout << "Particle " << particle_vector.size() << " is " << pdg << " with K_energy " << part.E()-0.001*part.Mass() << std::endl;
    }
  }
  
}

float opreco::OpticalRecoAna::update_MC_particle_time(simb::MCParticle const& particle, bool & pass_check, geo::Geometry const& geometry){

  pass_check = false;

  unsigned int tpc   = 0;
  unsigned int cstat = 0;

  const size_t numtraj = particle.NumberTrajectoryPoints();
  size_t t_iter = 0;
  while(t_iter < numtraj){
    try{
      // check if the particle is inside a TPC                                                                                         
      double pos[3] = {particle.Vx(t_iter), particle.Vy(t_iter), particle.Vz(t_iter)};
      geometry.PositionToTPC(pos, tpc, cstat);
    }
    catch(cet::exception &e){
      t_iter++;
      continue;
    }
    break;
  }

  if(t_iter == numtraj)
    return 0;

  pass_check = true;
  return particle.T(t_iter);
}

void opreco::OpticalRecoAna::match_flashes_to_tracks(std::vector<recob::OpFlash> const& flash_vector, 
						     std::vector<recob::Track> const& track_vector){
  bool matching=false;
  
  for(size_t i_flash=0; i_flash < flash_vector.size(); i_flash++){

    recob::OpFlash const& my_flash( flash_vector.at(i_flash) );
    if(my_flash.TotalPE() < fPEMin) continue;

    for(size_t i_track=0; i_track < track_vector.size(); i_track++){
      
      recob::Track const& my_track( track_vector.at(i_track) );
      
      matching=false;
      compare_track_and_flash(my_track,my_flash,matching);
      if(matching){
	fFlash_match_vector.at(i_flash).track_indices.push_back(i_track);
	fTrack_match_vector.at(i_track).flash_indices.push_back(i_flash);
      }

    }//end inner loop over tracks
    
  }//end loop over flashes
  
}//end match_flashes_to_tracks

void opreco::OpticalRecoAna::compare_track_and_flash(recob::Track const& track,
						     recob::OpFlash const& flash,
						     bool & matching){
}

void opreco::OpticalRecoAna::match_flashes_to_particles(std::vector<recob::OpFlash> const& flash_vector, 
							std::vector<simb::MCParticle> const& particle_vector,
							float const& ns_per_PMT_tick,
							geo::Geometry const& geometry){

  for(size_t i_particle=0; i_particle < particle_vector.size(); i_particle++){
    
    simb::MCParticle const& my_particle( particle_vector.at(i_particle) );
    bool pass_check = false;
    const float particle_time = update_MC_particle_time(my_particle,pass_check,geometry);
    if(!pass_check) continue;
    
    std::cout << "Particle " << i_particle << " (" << my_particle.PdgCode() << ")" << particle_time << std::endl;
    
    for(size_t i_flash=0; i_flash < flash_vector.size(); i_flash++){
      
      recob::OpFlash const& my_flash( flash_vector.at(i_flash) );
      if(my_flash.TotalPE() < fPEMin) continue;
      
      const float flash_time = my_flash.Time()*ns_per_PMT_tick;

      std::cout << "\tFlash " << i_flash << " time is " << flash_time << std::endl;

      fTimeDiff->Fill(particle_time-flash_time);
      fTimeDiff_fine->Fill(particle_time-flash_time);
      
      if( std::abs(particle_time - flash_time ) < fTimeMatchMax){
	fFlash_match_vector.at(i_flash).particle_indices.push_back(i_particle);
	fParticle_match_vector.at(i_particle).flash_indices.push_back(i_flash);
      }

    }//end inner loop over particles
    
  }//end loop over flashes
  
}//end match_flashes_to_particles

//currently unused...
void opreco::OpticalRecoAna::compare_particle_and_flash(simb::MCParticle const& particle,
							recob::OpFlash const& flash,
							bool & matching,
							float const& ns_per_PMT_tick,
							geo::Geometry const& geometry){
} 

void opreco::OpticalRecoAna::match_flashes_to_particles(std::vector<recob::Track> const& track_vector, 
							std::vector<simb::MCParticle> const& particle_vector,
							geo::Geometry const& geometry){
  bool matching=false;
  
  for(size_t i_track=0; i_track < track_vector.size(); i_track++){

    recob::Track const& my_track( track_vector.at(i_track) );

    for(size_t i_particle=0; i_particle < particle_vector.size(); i_particle++){
      
      simb::MCParticle const& my_particle( particle_vector.at(i_particle) );
      
      matching=false;
      compare_particle_and_track(my_particle,my_track,matching,geometry);
      if(matching){
	fTrack_match_vector.at(i_track).particle_indices.push_back(i_track);
	fParticle_match_vector.at(i_particle).track_indices.push_back(i_track);
      }

    }//end inner loop over particles
    
  }//end loop over tracks
  
}//end match_tracks_to_particles

void opreco::OpticalRecoAna::compare_particle_and_track(simb::MCParticle const& particle,
							recob::Track const& track,
							bool & matching,
							geo::Geometry const& geometry){
}
