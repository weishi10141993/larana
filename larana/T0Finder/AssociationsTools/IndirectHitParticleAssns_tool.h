
#include "larana/T0Finder/AssociationsTools/IHitParticleAssociations.h"

#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"

#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <algorithm>

namespace t0
{
////////////////////////////////////////////////////////////////////////
//
// Class:       IndirectHitParticleAssns
// Module Type: art tool
// File:        IndirectHitParticleAssns.h
//
//              This provides MC truth information by using output
//              reco Hit <--> MCParticle associations
//
// Configuration parameters:
//
// TruncMeanFraction     - the fraction of waveform bins to discard when
//
// Created by Tracy Usher (usher@slac.stanford.edu) on November 21, 2017
//
////////////////////////////////////////////////////////////////////////

class IndirectHitParticleAssns : virtual public IHitParticleAssociations
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit IndirectHitParticleAssns(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~IndirectHitParticleAssns();
    
    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset) override;
    
    /**
     *  @brief This rebuilds the internal maps
     */
    void CreateHitParticleAssociations(art::Event&, HitParticleAssociations*) override;

private:

    art::InputTag fHitModuleLabel;
    art::InputTag fHitPartAssnsModuleLabel;
    
};
}
