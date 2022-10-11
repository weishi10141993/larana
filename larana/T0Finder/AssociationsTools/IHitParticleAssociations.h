////////////////////////////////////////////////////////////////////////
///
/// \file  IMCTruthMatching.h
/// \brief This provides an interface which defines truth matching functions
///        made available to downstream analysis code
///
/// \author  usher@slac.stanford.edu
///
////////////////////////////////////////////////////////////////////////
#ifndef IMCTRUTHMATCHING_H
#define IMCTRUTHMATCHING_H

namespace art {
  class Event;
}
namespace fhicl {
  class ParameterSet;
}
namespace recob {
  class Hit;
}
namespace simb {
  class MCParticle;
}

#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

///code to link reconstructed objects back to the MC truth information
namespace t0 {

  using HitParticleAssociations =
    art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>;

  class IHitParticleAssociations {
  public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IHitParticleAssociations() noexcept = default;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;

    /**
     *  @brief This rebuilds the internal maps
     */
    virtual void CreateHitParticleAssociations(art::Event&, HitParticleAssociations*) = 0;
  };

} // namespace
#endif // IMCTRUTHMATCHING_H
