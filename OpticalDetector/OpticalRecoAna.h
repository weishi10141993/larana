//Analyzer for flash <----> track matching
#ifndef OpticalRecoAna_H
#define OpticalRecoAna_H

#include <string>
#include <vector>

// ART includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"

//LArSoft includes
#include "MCCheater/BackTracker.h"
#include "SimulationBase/MCParticle.h"
#include "RecoBase/Track.h"
#include "RecoBase/OpFlash.h"

#include "TH1.h"

namespace opreco {

  //bookkeeping on all the matches
  struct flash_match{
    art::Ptr<recob::OpFlash> flash;
    std::vector<simb::MCParticle> particles;
    std::vector< art::Ptr<recob::Track> > tracks;
  };
  struct track_match{
    art::Ptr<recob::Track> track;
    std::vector< art::Ptr<recob::OpFlash> > flashes;
    std::vector<simb::MCParticle> particles;
  };
  struct particle_match{
    const simb::MCParticle particle;
    std::vector< art::Ptr<recob::OpFlash> > flashes;
    std::vector< art::Ptr<recob::Track> > tracks;
  };

  class OpticalRecoAna : public art::EDAnalyzer{
  public:
 
    OpticalRecoAna(const fhicl::ParameterSet&);
    virtual ~OpticalRecoAna();

    void beginJob();
    void analyze (const art::Event&); 

  private:

    //fcl parameters here
    std::string fFlashModuleLabel;
    std::string fTrackModuleLabel;
    float       fKineticEnergyMin;
    float       fTimeMatchMax;

    std::vector<flash_match>    fFlash_match_vector;
    std::vector<track_match>    fTrack_match_vector;
    std::vector<particle_match> fParticle_match_vector;

    // function to get sublist of mc_particles that could leave a flash
    void get_MC_particle_list(sim::ParticleList parent_list,std::vector<simb::MCParticle>&);

    //matching functions
    void match_flashes_to_tracks(art::Handle<std::vector<recob::OpFlash>>, art::Handle<std::vector<recob::Track>>);

    void match_flashes_to_particles(art::Handle<std::vector<recob::OpFlash>>, std::vector<simb::MCParticle>);
    void match_tracks_to_particles(art::Handle<std::vector<recob::Track>>, std::vector<simb::MCParticle>);
    
    void check_flash_matches();
    void sort_and_print_flashes(std::vector<recob::OpFlash>);
    //bool compare_flashes(const recob::OpFlash&, const recob::OpFlash&);

    TH1F *fTimeDiff;
  };

} 




#endif
