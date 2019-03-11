//Analyzer for flash<--->track matching (testing against MC)
#include "OpticalRecoAna.h"

// LArSoft includes
#include "larana/OpticalDetector/OpDigiProperties.h"

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
  fTimeDiff_fine = tfs->make<TH1F>("htdiff_fine","Time difference between particles and flashes; t_diff (ns); flash/particle pairs",100,-1000,1000);
  fNBeamFlashes = tfs->make<TH1I>("hNBeamFlashes","Number of flashes OnBeamTime per event; N_{Flashes}; Events",5,0,5);

  fMatchTree_PF = tfs->make<TTree>("MatchTree_PF","MatchTree_PF");
  fMatchTree_PF->Branch("Run",&fRun,"Run/I");
  fMatchTree_PF->Branch("Event",&fEvent,"Event/I");
  fMatchTree_PF->Branch("p_id",&fParticleID,"p_id/I");
  fMatchTree_PF->Branch("p_time",&fParticleTime,"p_time/F");
  fMatchTree_PF->Branch("p_vx",&fParticleVx,"p_vx/F");
  fMatchTree_PF->Branch("p_vy",&fParticleVy,"p_vy/F");
  fMatchTree_PF->Branch("p_vz",&fParticleVz,"p_vz/F");
  fMatchTree_PF->Branch("p_mid",&fParticleMother,"p_mid/I");
  fMatchTree_PF->Branch("p_trackid",&fParticleTrackID,"p_trackid/I");
  fMatchTree_PF->Branch("FlashID",&fFlashID,"flash_id/I");
  fMatchTree_PF->Branch("f_time",&fFlashTime,"f_time/F");
  fMatchTree_PF->Branch("f_timewidth",&fFlashTimeWidth,"f_timewidth/F");
  fMatchTree_PF->Branch("f_y",&fFlashY,"f_y/F");
  fMatchTree_PF->Branch("f_z",&fFlashZ,"f_z/F");
  fMatchTree_PF->Branch("f_ywidth",&fFlashYWidth,"f_ywidth/F");
  fMatchTree_PF->Branch("f_zwidth",&fFlashYWidth,"f_zwidth/F");
  fMatchTree_PF->Branch("f_onbeamtime",&fFlashOnBeamTime,"f_onbeamtime/I");
  fMatchTree_PF->Branch("f_pe",&fFlashPE,"f_pe/F");
  fMatchTree_PF->Branch("matchIndex",&fMatchIndex,"matchIndex/I");
  
  fMatchTree_PF_NotNu = tfs->make<TTree>("MatchTree_PF_NotNu","MatchTree_PF_NotNu");
  fMatchTree_PF_NotNu->Branch("Run",&fRun,"Run/I");
  fMatchTree_PF_NotNu->Branch("Event",&fEvent,"Event/I");
  fMatchTree_PF_NotNu->Branch("p_id",&fParticleID,"p_id/I");
  fMatchTree_PF_NotNu->Branch("p_time",&fParticleTime,"p_time/F");
  fMatchTree_PF_NotNu->Branch("p_vx",&fParticleVx,"p_vx/F");
  fMatchTree_PF_NotNu->Branch("p_vy",&fParticleVy,"p_vy/F");
  fMatchTree_PF_NotNu->Branch("p_vz",&fParticleVz,"p_vz/F");
  fMatchTree_PF_NotNu->Branch("p_mid",&fParticleMother,"p_mid/I");
  fMatchTree_PF_NotNu->Branch("p_trackid",&fParticleTrackID,"p_trackid/I");
  fMatchTree_PF_NotNu->Branch("FlashID",&fFlashID,"flash_id/I");
  fMatchTree_PF_NotNu->Branch("f_time",&fFlashTime,"f_time/F");
  fMatchTree_PF_NotNu->Branch("f_timewidth",&fFlashTimeWidth,"f_timewidth/F");
  fMatchTree_PF_NotNu->Branch("f_y",&fFlashY,"f_y/F");
  fMatchTree_PF_NotNu->Branch("f_z",&fFlashZ,"f_z/F");
  fMatchTree_PF_NotNu->Branch("f_ywidth",&fFlashYWidth,"f_ywidth/F");
  fMatchTree_PF_NotNu->Branch("f_zwidth",&fFlashYWidth,"f_zwidth/F");
  fMatchTree_PF_NotNu->Branch("f_onbeamtime",&fFlashOnBeamTime,"f_onbeamtime/I");
  fMatchTree_PF_NotNu->Branch("f_pe",&fFlashPE,"f_pe/F");
  fMatchTree_PF_NotNu->Branch("matchIndex",&fMatchIndex_NotNu,"matchIndex/I");

}



// The analyzer itself
void opreco::OpticalRecoAna::analyze(const art::Event& evt)
{
  
  //get flag for MC
  const bool is_MC= !(evt.isRealData());

  fRun = evt.id().run();
  fEvent = evt.id().event();
  
  //first get out track, flash, and particle list handles
  art::Handle< std::vector<recob::OpFlash> > flash_handle;
  evt.getByLabel(fFlashModuleLabel, flash_handle);
  std::vector<recob::OpFlash> const& flash_vector(*flash_handle);
  fFlash_match_vector.resize(flash_vector.size());
  
  MF_LOG_INFO ("OpticalRecoAna")  
    << "Number of flashes is " << flash_vector.size() << std::flush;
  
  //art::Handle< std::vector<recob::Track> > track_handle;
  //evt.getByLabel(fTrackModuleLabel, track_handle);
  //std::vector<recob::Track> const& track_vector(*track_handle);
  //fTrack_match_vector.resize(track_vector.size());
  
  //MF_LOG_INFO ("OpticalRecoAna")  
  //<< "Number of tracks is " << track_vector.size() << std::flush;
  
  //match_flashes_to_tracks(flash_vector, track_vector);

  //all this for the MC matching
  if(is_MC){
    
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::vector<simb::MCParticle> particle_list;
    get_MC_particle_list(pi_serv->ParticleList(),particle_list);
    fParticle_match_vector.resize(particle_list.size());
    
    mf::LogInfo("OpticalRecoAna")  
      << "Number of MC particles is " << particle_list.size();
    

    art::ServiceHandle<opdet::OpDigiProperties> odp;
    const float ns_per_PMT_tick = 1e3;// already converted to microseconds//( 1e3 / odp->SampleFreq()) ; //SampleFreq is in MHz
    art::ServiceHandle<geo::Geometry> geometry_handle;

    match_flashes_to_particles(flash_vector,
			       particle_list,
			       ns_per_PMT_tick,
			       *geometry_handle);
    /*
    FillMatchTree_PF(flash_vector,
		     particle_list,
		     ns_per_PMT_tick,
		     *geometry_handle);
		     */
    int n_flashes = 0;
    for(auto const& my_flash : flash_vector){
      if(my_flash.OnBeamTime()==1) n_flashes++;
    }
    fNBeamFlashes->Fill(n_flashes);

    //check_flash_matches();

    //match_tracks_to_particles(track_handle,particle_list);
    
  }
  
}

simb::Origin_t opreco::OpticalRecoAna::get_MC_particle_origin(simb::MCParticle const& particle){
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  return (pi_serv->TrackIdToMCTruth_P(particle.TrackId()))->Origin();
}

void opreco::OpticalRecoAna::get_MC_particle_list(sim::ParticleList const& plist,std::vector<simb::MCParticle> & particle_vector) {
  
  for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
    
    const simb::MCParticle* particle = (*ipart).second;
    /*
    std::cout << "Particle " << particle_vector.size() << " (" << (*ipart).first << ") is " << particle->PdgCode() << " (mother=" << particle->Mother() 
		<< ") with K_energy " << particle->E()-0.001*particle->Mass() << std::endl;
      std::cout << "\t(X,Y,Z,T)=(" << particle->Vx() << "," << particle->Vy() << "," << particle->Vz() << "," << particle->T() << ")" << std::endl;
    */
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
      //std::cout << "Particle " << particle_vector.size() << " (" << (*ipart).first << ") is " << pdg << " (status=" << particle->StatusCode() 
      //	<< ") with K_energy " << particle->E()-0.001*particle->Mass() << std::endl;
      //std::cout << "\t(X,Y,Z,T)=(" << particle->Vx() << "," << particle->Vy() << "," << particle->Vz() << "," << particle->T() << ")" << std::endl;
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

  fMatchIndex=0;
  fMatchIndex_NotNu=0;
  int lastFlashID = -1;

  for(size_t i_particle=0; i_particle < particle_vector.size(); i_particle++){
    
    simb::MCParticle const& my_particle( particle_vector.at(i_particle) );
    bool pass_check = false;
    const float particle_time = update_MC_particle_time(my_particle,pass_check,geometry);
    if(!pass_check) continue;
    
    //std::cout << "Particle " << i_particle << " (" << my_particle.PdgCode() << ")" << particle_time << std::endl;
    
    for(size_t i_flash=0; i_flash < flash_vector.size(); i_flash++){
      
      recob::OpFlash const& my_flash( flash_vector.at(i_flash) );
      if(my_flash.TotalPE() < fPEMin) continue;
      
      const float flash_time = my_flash.Time()*ns_per_PMT_tick;

      //if(my_flash.OnBeamTime()) 
	//std::cout << "\tFlash " << i_flash << " time is " << flash_time << std::endl;

      fTimeDiff->Fill(particle_time-flash_time);
      fTimeDiff_fine->Fill(particle_time-flash_time);
      
      fParticleID = i_particle;
      fParticleTime = particle_time;
      fParticleVx = my_particle.Vx(0);
      fParticleVy = my_particle.Vy(0);
      fParticleVz = my_particle.Vz(0);
      fParticleMother = my_particle.Mother();
      fParticleTrackID = my_particle.TrackId();
      fFlashID = i_flash;
      fFlashTime = flash_time;
      fFlashTimeWidth = my_flash.TimeWidth()*ns_per_PMT_tick;
      fFlashY = my_flash.YCenter();
      fFlashZ = my_flash.ZCenter();
      fFlashYWidth = my_flash.YWidth();
      fFlashZWidth = my_flash.ZWidth();
      fFlashPE = my_flash.TotalPE();
      fFlashOnBeamTime = my_flash.OnBeamTime();

      bool beam_match = ((std::abs(particle_time - flash_time )/fFlashTimeWidth)<fTimeMatchMax &&
			 get_MC_particle_origin(my_particle)==simb::kBeamNeutrino);
      bool nonbeam_match = ((std::abs(particle_time - flash_time )/fFlashTimeWidth)<fTimeMatchMax && 
			    get_MC_particle_origin(my_particle)!=simb::kBeamNeutrino);

      if( beam_match ){	
	fMatchTree_PF->Fill();
	fMatchIndex++;

	//std::cout << "\t\tParticle (x,y,z)=(" << my_particle.Vx() << "," << my_particle.Vy() << "," << my_particle.Vz() << std::endl;
	//std::cout << "\t\tParticle PDG = " << my_particle.PdgCode() << " mother=" << my_particle.Mother() << std::endl;
	//std::cout << "\t\tFlash (y,z) = (" 
	//	  << my_flash.YCenter() << " +/- " << my_flash.YWidth() << ","
	//	  << my_flash.ZCenter() << " +/- " << my_flash.ZWidth() << std::endl;

	fFlash_match_vector.at(i_flash).particle_indices.push_back(i_particle);
	fParticle_match_vector.at(i_particle).flash_indices.push_back(i_flash);
      }
      else if( nonbeam_match ){	
	if(lastFlashID!=fFlashID){
	  fMatchTree_PF_NotNu->Fill();
	  fMatchIndex_NotNu++;
	  lastFlashID = fFlashID;
	}
      }


    }//end inner loop over particles
    
  }//end loop over flashes
  
}//end match_flashes_to_particles

void opreco::OpticalRecoAna::FillMatchTree_PF(std::vector<recob::OpFlash> const& flash_vector, 
					      std::vector<simb::MCParticle> const& particle_vector,
					      float const& ns_per_PMT_tick,
					      geo::Geometry const& geometry){

    for(size_t i_flash=0; i_flash < flash_vector.size(); i_flash++){

      recob::OpFlash const& my_flash( flash_vector.at(i_flash) );
      const float flash_time = my_flash.Time()*ns_per_PMT_tick;

      for(auto i_particle : fFlash_match_vector.at(i_flash).particle_indices){

	simb::MCParticle const& my_particle( particle_vector.at(i_particle) );
	bool pass_check = false;
	const float particle_time = update_MC_particle_time(my_particle,pass_check,geometry);
	
	if(my_particle.Mother()==0 && my_particle.TrackId()>1e4){

	  fParticleID = i_particle;
	  fParticleTime = particle_time;
	  fParticleVx = my_particle.Vx(0);
	  fParticleVy = my_particle.Vy(0);
	  fParticleVz = my_particle.Vz(0);
	  fFlashID = i_flash;
	  fFlashTime = flash_time;
	  fFlashY = my_flash.YCenter();
	  fFlashZ = my_flash.ZCenter();
	  fFlashYWidth = my_flash.YWidth();
	  fFlashZWidth = my_flash.ZWidth();

	  fMatchTree_PF->Fill();

	  return;
	}

      }

  }

}

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
