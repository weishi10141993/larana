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

namespace opreco {

  // Constructor
  OpticalRecoAna::OpticalRecoAna(fhicl::ParameterSet const& pset)
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
  OpticalRecoAna::~OpticalRecoAna() 
  {}
   


  // Do something here to setup the file (like make a TTree)
  void OpticalRecoAna::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
  }
   


  // The analyzer itself
  void OpticalRecoAna::analyze(const art::Event& evt) 
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
      std::vector<const simb::MCParticle*> particle_list = get_MC_particle_list(bt->ParticleList());
      fParticle_match_vector.resize(particle_list.size());
 
      mf::LogInfo("OpticalRecoAna")  
	<< "Number of MC particles is " << fParticle_match_vector.size();

      match_flashes_to_particles(flash_handle,particle_list);

      match_tracks_to_particles(track_handle,particle_list);

    }

  }
  


  std::vector<const simb::MCParticle*> OpticalRecoAna::get_MC_particle_list(sim::ParticleList plist) {

    std::vector<const simb::MCParticle*> my_list;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      
      const simb::MCParticle* part = (*ipart).second;

      //do not store if it's below our energy cut
      if( part->E() < 0.001*part->Mass() +  fKineticEnergyMin ) continue;
      
      //check to see if it's a charged particle we expect to leave ionization
      int pdg = std::abs(part->PdgCode());
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
	  ) my_list.push_back(part);
    }

    return my_list;

  }
  
  void OpticalRecoAna::match_flashes_to_tracks(art::Handle< std::vector<recob::OpFlash> > flash_handle, 
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

 void OpticalRecoAna::match_flashes_to_particles(art::Handle< std::vector<recob::OpFlash> > flash_handle, 
						 std::vector<const simb::MCParticle*>   particle_list){
    bool matching=false;
    
    for(size_t i_flash=0; i_flash < flash_handle->size(); i_flash++){
      art::Ptr<recob::OpFlash> my_flash(flash_handle, i_flash);
      if(!fFlash_match_vector.at(i_flash).flash) 
	fFlash_match_vector.at(i_flash).flash = my_flash;
      
      for(size_t i_particle=0; i_particle < particle_list.size(); i_particle++){
	const simb::MCParticle* my_particle = particle_list.at(i_particle);
		
	matching=false;

	if( std::abs(my_particle->T() - my_flash->Time() ) < fTimeMatchMax)
	  matching=true;

	if(matching)
	  fFlash_match_vector.at(i_flash).particles.push_back(my_particle);
	
      }//end inner loop over particles
      
    }//end loop over flashes
    


    for(size_t i_particle=0; i_particle < particle_list.size(); i_particle++){
      const simb::MCParticle* my_particle = particle_list.at(i_particle);
      
      for(size_t i_flash=0; i_flash < flash_handle->size(); i_flash++){
	art::Ptr<recob::OpFlash> my_flash(flash_handle, i_flash);
	if(!fFlash_match_vector.at(i_flash).flash) 
	  fFlash_match_vector.at(i_flash).flash = my_flash;
		
	matching=false;

	if( std::abs(my_particle->T() - my_flash->Time() ) < fTimeMatchMax)
	  matching=true;

	if(matching)
	  fParticle_match_vector.at(i_particle).flashes.push_back(my_flash);
	
      }//end inner loop over flashes
      
    }//end loop over particles

  }//end match_flashes_to_particles

 void OpticalRecoAna::match_tracks_to_particles(art::Handle< std::vector<recob::Track> > track_handle, 
						 std::vector<const simb::MCParticle*>   particle_list){
    bool matching=false;
    
    for(size_t i_track=0; i_track < track_handle->size(); i_track++){
      art::Ptr<recob::Track> my_track(track_handle, i_track);
      if(!fTrack_match_vector.at(i_track).track) 
	fTrack_match_vector.at(i_track).track = my_track;
      
      for(size_t i_particle=0; i_particle < particle_list.size(); i_particle++){
	const simb::MCParticle* my_particle = particle_list.at(i_particle);
		
	matching=false;
	//put matching code here?
	if(matching)
	  fTrack_match_vector.at(i_track).particles.push_back(my_particle);
	
      }//end inner loop over particles
      
    }//end loop over tracks
    
  }//end match_tracks_to_particles

} // namespace opdet


