////////////////////////////////////////////////////////////////////////
// \file DefaultOpDetResponse.h
//
// \brief service containing information about the response of optical detectors, which by default does nothing
//
// \author ahimmel@phy.duke.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef DEFAULT_OPDET_RESPONSE_H
#define DEFAULT_OPDET_RESPONSE_H

// LArSoft includes
#include "Simulation/SimPhotons.h"
#include "OpticalDetector/OpDetResponseInterface.h"



namespace opdet
{
    class DefaultOpDetResponse : public opdet::OpDetResponseInterface {
    public:

        DefaultOpDetResponse(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
        ~DefaultOpDetResponse() throw();



    private:

        virtual void doReconfigure(fhicl::ParameterSet const& p);
        virtual bool doDetected(int OpChannel, const sim::OnePhoton& Phot, int &newOpChannel) const;
        virtual bool doDetectedLite(int OpChannel, int &newOpChannel) const;

    }; // class DefaultOpDetResponse

    
} //namespace opdet


DECLARE_ART_SERVICE_INTERFACE_IMPL(opdet::DefaultOpDetResponse, opdet::OpDetResponseInterface, LEGACY)

#endif //DEFAULT_OPDET_RESPONSE_H
