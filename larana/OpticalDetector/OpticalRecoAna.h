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
#include "canvas/Persistency/Common/Ptr.h"

//LArSoft includes
//#include "larsim/MCCheater/BackTrackerService.h" //Not used
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcore/Geometry/Geometry.h"

#include "TH1.h"
#include "TTree.h"

namespace opreco {

  //bookkeeping on all the matches
  struct flash_match{
    std::vector<size_t> particle_indices;
    std::vector<size_t> track_indices;
  };
  struct track_match{
    std::vector<size_t> particle_indices;
    std::vector<size_t> flash_indices;
  };
  struct particle_match{
    std::vector<size_t> flash_indices;
    std::vector<size_t> track_indices;
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
    float       fPEMin;
    float       fTimeMatchMax;

    std::vector<flash_match>    fFlash_match_vector;
    std::vector<track_match>    fTrack_match_vector;
    std::vector<particle_match> fParticle_match_vector;

    void get_MC_particle_list(sim::ParticleList const& ,std::vector<simb::MCParticle> & );
    simb::Origin_t get_MC_particle_origin(simb::MCParticle const& );
    float update_MC_particle_time(simb::MCParticle const&, bool& ,geo::Geometry const&);

    //matching functions
    void match_flashes_to_tracks(std::vector<recob::OpFlash> const& , 
				 std::vector<recob::Track> const&);
    void compare_track_and_flash(recob::Track const&,
				 recob::OpFlash const&,
				 bool &);

    void match_flashes_to_particles(std::vector<recob::OpFlash> const&, 
				    std::vector<simb::MCParticle> const&,
				    float const&,
				    geo::Geometry const&);
    void compare_particle_and_flash(simb::MCParticle const&,
				    recob::OpFlash const&,
				    bool &,
				    float const&,
				    geo::Geometry const&);

    void match_flashes_to_particles(std::vector<recob::Track> const&, 
				    std::vector<simb::MCParticle> const&,
				    geo::Geometry const&);
    void compare_particle_and_track(simb::MCParticle const&,
				    recob::Track const&,
				    bool &,
				    geo::Geometry const&);    

    void FillMatchTree_PF(std::vector<recob::OpFlash> const&, 
			  std::vector<simb::MCParticle> const&,
			  float const&,
			  geo::Geometry const&);


    TH1F *fTimeDiff;
    TH1F *fTimeDiff_fine;
    TH1I *fNBeamFlashes;
    
    TTree *fMatchTree_PF;
    TTree *fMatchTree_PF_NotNu;
    int fRun;
    int fEvent;
    int fParticleID;
    int fParticleMother;
    int fParticleTrackID;
    float fParticleTime;
    float fParticleVx;
    float fParticleVy;
    float fParticleVz;
    int fFlashID;
    float fFlashTime;
    float fFlashTimeWidth;
    float fFlashY;
    float fFlashZ;
    float fFlashYWidth;
    float fFlashZWidth;
    float fFlashPE;
    int fFlashOnBeamTime;
    int fMatchIndex;
    int fMatchIndex_NotNu;

  };

} 




#endif
