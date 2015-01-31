////////////////////////////////////////////////////////////////////////
//
//  \file MicrobooneOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////


#include "OpticalDetector/MicrobooneOpDetResponse.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Random/RandFlat.h"


namespace opdet{


    //--------------------------------------------------------------------
    MicrobooneOpDetResponse::MicrobooneOpDetResponse(fhicl::ParameterSet const& pset, 
                                         art::ActivityRegistry &/*reg*/)
    {
        this->doReconfigure(pset);
    }
    
    //--------------------------------------------------------------------
    MicrobooneOpDetResponse::~MicrobooneOpDetResponse() throw()
    { }


    //--------------------------------------------------------------------
    void MicrobooneOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
    {
        fQE=                       pset.get<double>("QuantumEfficiency");
        fWavelengthCutLow=         pset.get<double>("WavelengthCutLow");
        fWavelengthCutHigh=        pset.get<double>("WavelengthCutHigh");
    }


    //--------------------------------------------------------------------
    bool MicrobooneOpDetResponse::doDetected(int OpChannel, const sim::OnePhoton& Phot, int &newOpChannel) const
    {
        newOpChannel = OpChannel;
        
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        double wavel = wavelength(Phot.Energy);
        // Check wavelength acceptance
        if (wavel < fWavelengthCutLow) return false;
        if (wavel > fWavelengthCutHigh) return false;

        return true;
    }
    
    //--------------------------------------------------------------------
    bool MicrobooneOpDetResponse::doDetectedLite(int OpChannel, int &newOpChannel) const
    {
        newOpChannel = OpChannel;
        
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        return true;
    }



} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::MicrobooneOpDetResponse, opdet::OpDetResponseInterface)

