//Analyzer for flash<--->track matching (testing against MC)
#include "OpticalRecoAna.h"

// LArSoft includes
#include "Geometry/Geometry.h"

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
}



// The analyzer itself
void opreco::OpticalRecoAna::analyze(const art::Event& evt)
{
  
  fFlash_match_vector.clear();
  fTrack_match_vector.clear();
  fParticle_match_vector.clear();
  
  //get flag for MC
  bool is_MC= !(evt.isRealData());
  
  //first get out track, flash, and particle list handles
  art::Handle< std::vector<recob::OpFlash> > flash_handle;
  evt.getByLabel(fFlashModuleLabel, flash_handle);
  fFlash_match_vector.resize(flash_handle->size());
  
  LOG_INFO ("OpticalRecoAna")  
    << "Number of flashes is " << fFlash_match_vector.size() << std::flush;
  
  art::Handle< std::vector<recob::Track> > track_handle;
  evt.getByLabel(fTrackModuleLabel, track_handle);
  fTrack_match_vector.resize(track_handle->size());
  
  LOG_INFO ("OpticalRecoAna")  
    << "Number of tracks is " << fTrack_match_vector.size() << std::flush;
  
  match_flashes_to_tracks(flash_handle, track_handle);
  
  //all this for the MC matching
  if(is_MC){
    
    art::ServiceHandle<cheat::BackTracker> bt;
    std::vector<simb::MCParticle> particle_list;
    get_MC_particle_list(bt->ParticleList(),particle_list);
    fParticle_match_vector.resize(particle_list.size());
    
    mf::LogInfo("OpticalRecoAna")  
      << "Number of MC particles is " << fParticle_match_vector.size();
    
    match_flashes_to_particles(flash_handle,particle_list);
    check_flash_matches();

    //match_tracks_to_particles(track_handle,particle_list);
    
  }
  
}



void opreco::OpticalRecoAna::get_MC_particle_list(sim::ParticleList plist,std::vector<simb::MCParticle> & particle_vector) {
  
  for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
    
    const simb::MCParticle* part_ptr = (*ipart).second;
    simb::MCParticle part(*part_ptr);
    
    //do not store if it's below our energy cut
    if( part.E() < 0.001*part.Mass() +  fKineticEnergyMin ) continue;
    
    //check to see if it's a charged particle we expect to leave ionization
    int pdg = std::abs(part.PdgCode());
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
      particle_vector.push_back(part);
      std::cout << "Particle " << particle_vector.size() << " is " << pdg << std::endl;
    }
  }
  
}

void opreco::OpticalRecoAna::match_flashes_to_tracks(art::Handle< std::vector<recob::OpFlash> > flash_handle, 
						     art::Handle< std::vector<recob::Track> >   track_handle){
  bool matching=false;
  
  for(size_t i_flash=0; i_flash < flash_handle->size(); i_flash++){
    art::Ptr<recob::OpFlash> my_flash(flash_handle, i_flash);
    if(!fFlash_match_vector.at(i_flash).flash) 
      fFlash_match_vector.at(i_flash).flash = my_flash;
    
    for(size_t i_track=0; i_track < track_handle->size(); i_track++){
      
      art::Ptr<recob::Track> my_track(track_handle, i_track);
      
      matching=false;
      //put matching code here?
      if(matching)
	fFlash_match_vector.at(i_flash).tracks.push_back(my_track);
      
    }//end inner loop over tracks
    
  }//end loop over flashes
  
}//end match_falshes_to_tracks

void opreco::OpticalRecoAna::match_flashes_to_particles(art::Handle< std::vector<recob::OpFlash> > flash_handle, 
							std::vector<simb::MCParticle>   particle_list){
  bool matching=false;
  
  for(size_t i_flash=0; i_flash < flash_handle->size(); i_flash++){
    
    std::cout << "Processing flash " << i_flash << std::endl;
    
    art::Ptr<recob::OpFlash> my_flash(flash_handle, i_flash);
    if(!fFlash_match_vector.at(i_flash).flash) 
      fFlash_match_vector.at(i_flash).flash = my_flash;
    
    std::cout << "Got ptr to flash " << i_flash << std::endl;
    std::cout << "Total PE in this falsh is " << my_flash->TotalPE() << std::endl;

    for(size_t i_particle=0; i_particle < particle_list.size(); i_particle++){
      const simb::MCParticle my_particle = particle_list.at(i_particle);
      
      std::cout << "\tChecking particle " << i_particle << std::endl;
      
      matching=false;
      
      fTimeDiff->Fill(my_particle.T() - my_flash->Time()*16);

      if( std::abs(my_particle.T() - my_flash->Time()*16 ) < fTimeMatchMax)
	matching=true;
      
      std::cout << "\t\t(Flash,Particle) = (" << my_flash->Time()*16 << "," << my_particle.T() << ") match? " << matching << std::endl;
      
      if(matching)
	fFlash_match_vector.at(i_flash).particles.push_back(my_particle);
      
    }//end inner loop over particles
    
  }//end loop over flashes
  
  
  
  for(size_t i_particle=0; i_particle < particle_list.size(); i_particle++){
    const simb::MCParticle my_particle = particle_list.at(i_particle);
    
    for(size_t i_flash=0; i_flash < flash_handle->size(); i_flash++){
      art::Ptr<recob::OpFlash> my_flash(flash_handle, i_flash);
      if(!fFlash_match_vector.at(i_flash).flash) 
	fFlash_match_vector.at(i_flash).flash = my_flash;
      
      matching=false;
      
      if( std::abs(my_particle.T() - my_flash->Time() ) < fTimeMatchMax)
	matching=true;
      
      if(matching)
	fParticle_match_vector.at(i_particle).flashes.push_back(my_flash);
      
    }//end inner loop over flashes
    
  }//end loop over particles
  
}//end match_flashes_to_particles

void opreco::OpticalRecoAna::check_flash_matches(){
  
  for( std::vector<opreco::flash_match>::iterator f_iter=fFlash_match_vector.begin(); 
       f_iter<fFlash_match_vector.end();
       f_iter++ ){

    size_t particles_per_flash = (*f_iter).particles.size();
    size_t tracks_per_flash    = (*f_iter).tracks.size();

    std::cout << "Flash " << std::distance(fFlash_match_vector.begin(),f_iter) 
	      << " (particles,tracks)=(" << particles_per_flash << "," << tracks_per_flash << ")" << std::endl;
   
  }
  
}
  
void opreco::OpticalRecoAna::match_tracks_to_particles(art::Handle< std::vector<recob::Track> > track_handle, 
						 std::vector<simb::MCParticle>   particle_list){
    bool matching=false;
    
    for(size_t i_track=0; i_track < track_handle->size(); i_track++){
      art::Ptr<recob::Track> my_track(track_handle, i_track);
      if(!fTrack_match_vector.at(i_track).track) 
	fTrack_match_vector.at(i_track).track = my_track;
      
      for(size_t i_particle=0; i_particle < particle_list.size(); i_particle++){
	const simb::MCParticle my_particle = particle_list.at(i_particle);
		
	matching=false;
	//put matching code here?
	if(matching)
	  fTrack_match_vector.at(i_track).particles.push_back(my_particle);
	
      }//end inner loop over particles
      
    }//end loop over tracks
    
  }//end match_tracks_to_particles

bool compare_flashes(const recob::OpFlash& f1, const recob::OpFlash& f2){
  return (f1.Time() < f2.Time());
}

void opreco::OpticalRecoAna::sort_and_print_flashes(std::vector<recob::OpFlash> flash_vector){
  
  std::sort(flash_vector.begin(), flash_vector.end(), compare_flashes);
  for ( auto const& flash : flash_vector)
    std::cout << "Flash time is " << flash.Time();

}
