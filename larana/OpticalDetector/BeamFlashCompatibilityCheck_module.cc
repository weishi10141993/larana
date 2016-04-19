//  This module makes a set of hypotheses for track light profiles
//   based on beizer tracks in the event.
//
//  This module will likely be subsumed into the OpFlashFinder module
//   eventually. But at present they are being developped in parallel.
//
//
// \file BeamFlashCompatabilityCheck.cxx
// \author Ben Jones and Christie Chiu, MIT 2010
//         Heavily revised by Ben Jones and Wes Ketchum, 2013
//



#include "art/Framework/Core/EDProducer.h"
//#include "AnalysisBase/FlashMatch.h"
#include "lardata/AnalysisBase/CosmicTag.h"


// ROOT includes.
#include <Rtypes.h>
#ifndef BeamFlashCompatabilityCheck_h
#define BeamFlashCompatabilityCheck_h 1



namespace trkf{
  class BezierTrack;
}

namespace recob{
  class OpFlash;
  class Track;

}


namespace opdet {
  
  bool BeamFlashCompatabilityCheck_tracksort(art::Ptr<recob::Track> t1, art::Ptr<recob::Track> t2);

  class BeamFlashCompatabilityCheck : public art::EDProducer{
  public:
    
    BeamFlashCompatabilityCheck(const fhicl::ParameterSet&);
    virtual ~BeamFlashCompatabilityCheck();
    
    void produce(art::Event&);
    void reconfigure(fhicl::ParameterSet const& p);
      
    std::vector<double>               GetMIPHypotheses(trkf::BezierTrack* BTrack, double XOffset=0);
    bool                              CheckCompatibility(std::vector<double>& hypothesis, std::vector<double>& signal);
    
    void beginJob();
    
    
  private:
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    int         fBezierResolution;
    double      fSingleChannelCut;
    double      fIntegralCut;
  };

  

}

#endif




////////////////////////////////////////////////////////////////////////
/// \file  BeamFlashCompatabilityCheck_module.cc
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(BeamFlashCompatabilityCheck)

}//end namespace opdet
////////////////////////////////////////////////////////////////////////



//  This module makes a set of hypotheses for track light profiles
//   based on beizer tracks in the event.
//
//  This module will likely be subsumed into the OpFlashFinder module
//   eventually. But at present they are being developped in parallel.
//
//
// \file BeamFlashCompatabilityCheck.cxx
// \author Ben Jones and Christie Chiu, MIT 2010
//
//

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "lardata/RecoObjects/BezierTrack.h"
#include "lardata/RecoBase/OpFlash.h"

// FMWK includes
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "TH2D.h"

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// Debug flag; only used during code development.
const bool debug = true;

namespace opdet {

  //-------------------------------------------------

  BeamFlashCompatabilityCheck::BeamFlashCompatabilityCheck(fhicl::ParameterSet const& pset)
  {

    produces< std::vector<anab::CosmicTag> >();
    produces< art::Assns<recob::Track, anab::CosmicTag> >();

    this->reconfigure(pset);
   }


  //-------------------------------------------------

  void BeamFlashCompatabilityCheck::reconfigure(fhicl::ParameterSet const& pset)
  {
    fTrackModuleLabel = pset.get<std::string>("TrackModuleLabel");   
    fFlashModuleLabel = pset.get<std::string>("FlashModuleLabel");
    fBezierResolution = pset.get<int>("BezierResolution");
    fSingleChannelCut = pset.get<int>("SingleChannelCut");
    fIntegralCut      = pset.get<int>("IntegralCut");
  }


  //-------------------------------------------------

  void BeamFlashCompatabilityCheck::beginJob()
  {
  }



  //-------------------------------------------------

  BeamFlashCompatabilityCheck::~BeamFlashCompatabilityCheck()
  {
  }







  // Get a hypothesis for the light collected for a bezier track
  std::vector<double> BeamFlashCompatabilityCheck::GetMIPHypotheses(trkf::BezierTrack* Btrack, double XOffset)
  {
    art::ServiceHandle<geo::Geometry> geom;
    std::vector<double> ReturnVector(geom->NOpDets(),0);
    
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;

    float TrackLength = Btrack->GetLength();

    double xyz[3];
    for (int b=0; b!=fBezierResolution; b++)
      {
	float s               = float(b) / float(fBezierResolution);
		
	double MIPYield   = 24000;
	double QE         = 0.01;
	double MIPdQdx    = 2.1;
	double PromptFrac = 0.25;
	double PromptMIPScintYield = MIPYield * QE * MIPdQdx * PromptFrac;

	
	Btrack->GetTrackPoint(s,xyz);
	xyz[0]+=XOffset;
	const float* PointVisibility = pvs->GetAllVisibilities(xyz);
	if (!PointVisibility) continue; // point not covered by the service
	
	float LightAmount = PromptMIPScintYield*TrackLength/float(fBezierResolution);
	
	for(size_t OpDet =0; OpDet!=pvs->NOpChannels();  OpDet++)
	  {
	    ReturnVector.at(OpDet)+= PointVisibility[OpDet] * LightAmount;
	  }
      }
    return ReturnVector;
  }

  //-------------------------------------------------


  void BeamFlashCompatabilityCheck::produce(art::Event& evt)
  {
    
    //    int EventID = evt.id().event();
    
    
    // Read in flashes from the event
    art::Handle< std::vector<recob::OpFlash> > flashh;
    evt.getByLabel(fFlashModuleLabel, flashh);
    std::vector<art::Ptr<recob::OpFlash> > Flashes;
    for(unsigned int i=0; i < flashh->size(); ++i)
      {
	art::Ptr<recob::OpFlash> flash(flashh,i);
        Flashes.push_back(flash);
      }

    // Read in tracks from the event
    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);
    std::vector<art::Ptr<recob::Track> >  Tracks;
    for(unsigned int i=0; i < trackh->size(); ++i)
      {
	art::Ptr<recob::Track> track(trackh,i);
	Tracks.push_back(track);
      }


    // Use these to produce Bezier tracks
    std::vector<trkf::BezierTrack*> BTracks;
    BTracks.clear();
    for(size_t i=0; i!=Tracks.size(); i++)
      BTracks.push_back(new trkf::BezierTrack(*Tracks.at(i)));
        
      
    std::vector<std::vector<double> > TrackHypotheses;
    std::vector<std::vector<double> > FlashShapes;
       

    // For each track
    for (size_t i=0; i!=BTracks.size(); ++i)
      {
	TrackHypotheses.push_back(GetMIPHypotheses(BTracks.at(i)));
      }
    
    art::ServiceHandle<geo::Geometry> geom;
    size_t NOpDets = geom->NOpDets();
    size_t NOpChannels = geom->NOpChannels();
    
    std::vector<bool> Compatible(TrackHypotheses.size(),false);

    for(size_t f=0; f!=Flashes.size(); ++f)
      {
	if(Flashes.at(f)->OnBeamTime())
	  {
	    std::vector<double> ThisFlashShape(NOpDets,0);
	    for(size_t i=0; i!=NOpDets; ++i)
	      ThisFlashShape[i] =  0;
            for(size_t i=0; i < NOpChannels; ++i) {
              int opdet = geom->OpDetFromOpChannel(i);
              ThisFlashShape[opdet] += Flashes.at(f)->PE(i);
            }
	    FlashShapes.push_back(ThisFlashShape);

	    for(size_t i=0; i!=TrackHypotheses.size(); ++i)
	      {
		if(CheckCompatibility(TrackHypotheses.at(i),ThisFlashShape)) 
		  {
		    Compatible[i]=true;
		  }
	      }
	  }
      }

    std::unique_ptr< std::vector<anab::CosmicTag> > CosmicTagPtr ( new std::vector<anab::CosmicTag>);
    std::vector<anab::CosmicTag> & CosmicTagVector(*CosmicTagPtr);

    std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag > > assn_track( new art::Assns<recob::Track, anab::CosmicTag>);

    float cosmic_score;
    double xyz_begin[3];
    double xyz_end[3];
    for(size_t itrack=0; itrack<BTracks.size(); itrack++){
      

      if(Compatible.at(itrack)) cosmic_score = 0; //not a cosmic
      else if(!Compatible.at(itrack)) cosmic_score = 1; //is a cosmic
      
      BTracks.at(itrack)->GetTrackPoint(0,xyz_begin); //load in beginning point
      std::vector<float> endPt1 = { (float)xyz_begin[0], (float)xyz_begin[1], (float)xyz_begin[2] };
      BTracks.at(itrack)->GetTrackPoint(1,xyz_end); //load in ending point
      std::vector<float> endPt2 = { (float)xyz_end[0], (float)xyz_end[1], (float)xyz_end[2] };

      CosmicTagVector.emplace_back(endPt1,endPt2,cosmic_score,anab::CosmicTagID_t::kFlash_BeamIncompatible);
      util::CreateAssn(*this, evt, *(CosmicTagPtr.get()), Tracks.at(itrack), *(assn_track.get()), itrack); 

    }


    evt.put(std::move(CosmicTagPtr));
    evt.put(std::move(assn_track));

    
  }


  //---------------------------------------
  //  Check whether a hypothesis can be accomodated in a flash
  //   Flashes fail if 1 bin is far in excess of the observed signal 
  //   or if the whole flash intensity is much too large for the hypothesis.
  //  MIP dEdx is assumed for now.  Accounting for real dQdx will 
  //   improve performance of this algorithm.
  //
  bool BeamFlashCompatabilityCheck::CheckCompatibility(std::vector<double>& hypothesis, std::vector<double>& signal)
  {
    double sigintegral=0, hypintegral=0;
    for(size_t i=0; i!=hypothesis.size(); ++i)
      {
	sigintegral+=signal.at(i);
	hypintegral+=hypothesis.at(i);
	double HypErr = pow(hypothesis.at(i),0.5);
	if(( (hypothesis.at(i) - signal.at(i)) / HypErr) > fSingleChannelCut) return false;
      }
    double HypIntErr= pow(hypintegral,0.5);
    
    if( ( (hypintegral - sigintegral)/HypIntErr) > fIntegralCut) return false;
    return true;
  }


}


