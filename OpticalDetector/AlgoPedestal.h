/**
 * \file AlgoPedestal.h
 *
 * \ingroup PulseReco
 * 
 * \brief Class definition file of AlgoPedestal
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PulseReco
    
@{*/

#ifndef ALGOPEDESTAL_H
#define ALGOPEDESTAL_H


#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArSoft
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoAlg/SpacePointAlg.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/Geometry.h"
// STL
#include <set>
#include <vector>
#include <cmath>
#include <functional>
#include <numeric>


// ROOT
#include <TString.h>
#include <TTree.h>

namespace pmtana
{


  /**
   \class AlgoPedestal
   A class that calculates pedestal mean & standard deviation (here and elsewhere called as "RMS").   
  */
  class AlgoPedestal {

  public:

    /// Default constructor
    AlgoPedestal(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~AlgoPedestal();

    /// Method to compute a pedestal of the input waveform using "nsample" ADC samples from "start" index.
    void ComputePedestal(const std::vector<uint16_t>* wf, size_t start, size_t nsample) const;

    /// Getter of the pedestal mean value
    double Mean() const {return _mean;};

    /// Getter of the pedestal standard deviation
    double Sigma() const {return _sigma;};

  protected:

    /// A variable holder for pedestal mean value
    mutable double _mean;

    /// A variable holder for pedestal standard deviation
    mutable double _sigma;


  };
}
#endif

/** @} */ // end of doxygen group
