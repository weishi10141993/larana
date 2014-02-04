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
#include "AnalysisBase/FlashMatch.h"


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
#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "RecoObjects/BezierTrack.h"
#include "RecoBase/OpFlash.h"

// FMWK includes
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Utilities/LArProperties.h"
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

    produces< std::vector<anab::FlashMatch> >();
    produces< art::Assns<recob::Track, anab::FlashMatch> >();

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
    std::vector<double> ReturnVector(geom->NOpDet(),0);
    
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
	std::vector<float>* PointVisibility = pvs->GetAllVisibilities(xyz);
	float LightAmount = PromptMIPScintYield*TrackLength/float(fBezierResolution);
	
	for(size_t OpDet =0; OpDet!=PointVisibility->size();  OpDet++)
	  {
	    ReturnVector.at(OpDet)+= PointVisibility->at(OpDet) * LightAmount;
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
    size_t NOpDets = geom->NOpDet();
    
    for(size_t f=0; f!=Flashes.size(); ++f)
      {
	if(Flashes.at(f)->OnBeamTime())
	  {
	    std::vector<double> ThisFlashShape(NOpDets,0);
	    for(size_t i=0; i!=NOpDets; ++i)
	      ThisFlashShape[i]=Flashes.at(f)->PE(i);
	    FlashShapes.push_back(ThisFlashShape);
	  }
      }


    // This will hold the list of matches
    std::vector<anab::FlashMatch> Matches;

    
    std::vector<bool> Compatible(TrackHypotheses.size(),false);
    for(size_t i=0; i!=TrackHypotheses.size(); ++i)
      {
	for(size_t j=0; j!=FlashShapes.size(); ++j)	    
	  {
	    if(CheckCompatibility(TrackHypotheses.at(i),FlashShapes.at(j))) 
	      {
		Compatible[i]=true;
	      }
	  }
	if(!Compatible[i])  Matches.push_back(anab::FlashMatch(-1,-1,i,false)); 
      }
    

    std::unique_ptr< std::vector<anab::FlashMatch> > flash_matches ( new std::vector<anab::FlashMatch>);
    std::unique_ptr< art::Assns<recob::Track, anab::FlashMatch > > assn_track( new art::Assns<recob::Track, anab::FlashMatch>);

    for(size_t i=0; i!=Matches.size(); ++i)
      {
	std::cout<<"Storing negative match for " << Matches.at(i).SubjectID()<<std::endl;
	
	flash_matches->push_back(Matches.at(i));
	
	util::CreateAssn(*this, evt, *(flash_matches.get()), Tracks.at(Matches.at(i).SubjectID()), *(assn_track.get()), i); 

      }
    

    evt.put(std::move(flash_matches));
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


