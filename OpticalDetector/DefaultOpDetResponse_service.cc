////////////////////////////////////////////////////////////////////////
//
//  \file DefaultOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////


#include "OpticalDetector/DefaultOpDetResponse.h"


namespace opdet{


    //--------------------------------------------------------------------
    DefaultOpDetResponse::DefaultOpDetResponse(fhicl::ParameterSet const& pset, 
                                         art::ActivityRegistry &/*reg*/)
    {
        this->doReconfigure(pset);
    }
    
    //--------------------------------------------------------------------
    DefaultOpDetResponse::~DefaultOpDetResponse() throw()
    { }


    //--------------------------------------------------------------------
    void DefaultOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
    { }


    //--------------------------------------------------------------------
    bool DefaultOpDetResponse::doDetected(int /*OpChannel*/, const sim::OnePhoton& /*Phot*/) const
    {
        return true;
    }
    
    //--------------------------------------------------------------------
    bool DefaultOpDetResponse::doDetectedLite(int /*OpChannel*/) const
    {
        return true;
    }



} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::DefaultOpDetResponse, opdet::OpDetResponseInterface)

