////////////////////////////////////////////////////////////////////////
// \file OpDetResponse.h
//
// \brief service containing information about the response of optical detectors
//
// \author ahimmel@phy.duke.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef OPDET_RESPONSE_INTERFACE_H
#define OPDET_RESPONSE_INTERFACE_H

// LArSoft includes
//#include "OpticalDetectorData/OpticalTypes.h"
#include "Geometry/Geometry.h"
#include "Simulation/SimPhotons.h"

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "CLHEP/Random/RandFlat.h"


namespace opdet
{
    class OpDetResponseInterface {
    public:

        virtual ~OpDetResponseInterface() = default;

        virtual void reconfigure(fhicl::ParameterSet const& p);
        virtual bool detected(int OpChannel, const sim::OnePhoton& Phot) const;
        virtual bool detectedLite(int OpChannel) const;

        virtual float wavelength(double energy) const;

    private:
        virtual void doReconfigure(fhicl::ParameterSet const& p) = 0;
        virtual bool doDetected(int OpChannel, const sim::OnePhoton& Phot) const = 0;
        virtual bool doDetectedLite(int OpChannel) const = 0;


    }; // class OpDetResponse

    inline void OpDetResponseInterface::reconfigure(fhicl::ParameterSet const& p) { doReconfigure(p); }
    inline bool OpDetResponseInterface::detected(int OpChannel, const sim::OnePhoton& Phot) const { return doDetected(OpChannel, Phot); }
    inline bool OpDetResponseInterface::detectedLite(int OpChannel) const { return doDetectedLite(OpChannel); }

    inline float OpDetResponseInterface::wavelength(double energy) const { return (2.0*3.142)*0.000197/energy; }

} //namespace opdet


DECLARE_ART_SERVICE_INTERFACE(opdet::OpDetResponseInterface, LEGACY)

#endif //OPDET_RESPONSE_H
