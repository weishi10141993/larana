////////////////////////////////////////////////////////////////////////
//
//  \file DefaultOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////


#include "OpticalDetector/DefaultOpDetResponse.h"
#include "Utilities/LArProperties.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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
    {
        art::ServiceHandle<util::LArProperties>   LarProp;

        if ( LarProp->ScintPreScale() < 1 ) {
            mf::LogWarning("DefaultOpDetResponse_service") << "A prescale of " << LarProp->ScintPreScale() << " has been applied during optical MC production, "
                                                           << "but DefaultOpDetResponse does not include any QE so this effect is not being corrected out.";
            assert(false);
        }

    }


    //--------------------------------------------------------------------
    bool DefaultOpDetResponse::doDetected(int OpChannel, const sim::OnePhoton& /*Phot*/, int &newOpChannel) const
    {
        newOpChannel = OpChannel;
        return true;
    }
    
    //--------------------------------------------------------------------
    bool DefaultOpDetResponse::doDetectedLite(int OpChannel, int &newOpChannel) const
    {
        newOpChannel = OpChannel;
        return true;
    }



} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::DefaultOpDetResponse, opdet::OpDetResponseInterface)

