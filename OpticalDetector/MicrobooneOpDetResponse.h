////////////////////////////////////////////////////////////////////////
// \file MicrobooneOpDetResponse.h
//
// \brief service containing information about the response of optical detectors in Microboone
//
// \author ahimmel@phy.duke.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MICROBOONE_OPDET_RESPONSE_H
#define MICROBOONE_OPDET_RESPONSE_H

// LArSoft includes
#include "Simulation/SimPhotons.h"
#include "OpticalDetector/OpDetResponseInterface.h"



namespace opdet
{
    class MicrobooneOpDetResponse : public opdet::OpDetResponseInterface {
    public:

        MicrobooneOpDetResponse(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
        ~MicrobooneOpDetResponse() throw();



    private:

        virtual void doReconfigure(fhicl::ParameterSet const& p);
        virtual bool doDetected(int OpChannel, const sim::OnePhoton& Phot, int &newOpChannel) const;
        virtual bool doDetectedLite(int OpChannel, int &newOpChannel) const;

        float fQE;                     // Quantum efficiency of tube
        
        float fWavelengthCutLow;       // Sensitive wavelength range 
        float fWavelengthCutHigh;      // 
        


    }; // class MicrobooneOpDetResponse

    
} //namespace opdet


DECLARE_ART_SERVICE_INTERFACE_IMPL(opdet::MicrobooneOpDetResponse, opdet::OpDetResponseInterface, LEGACY)

#endif //MICROBOONE_OPDET_RESPONSE_H
