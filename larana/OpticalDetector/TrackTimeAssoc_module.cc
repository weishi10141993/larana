//  This module makes a set of hypotheses for track light profiles
//   based on beizer tracks in the event.
//
//  This module will likely be subsumed into the OpFlashFinder module
//   eventually. But at present they are being developped in parallel.
//
//
// \file TrackTimeAssoc.cxx
// \author Ben Jones and Christie Chiu, MIT 2010
//         Heavily revised by Ben Jones and Wes Ketchum, 2013
//



#include "art/Framework/Core/EDProducer.h"
#include "AnalysisBase/FlashMatch.h"


// ROOT includes.
#include <Rtypes.h>
#ifndef TrackTimeAssoc_h
#define TrackTimeAssoc_h 1



namespace trkf{
  class BezierTrack;
}

namespace recob{
  class OpFlash;
  class Track;

}


namespace opdet {
  
  bool TrackTimeAssoc_tracksort(art::Ptr<recob::Track> t1, art::Ptr<recob::Track> t2);

  class TrackTimeAssoc : public art::EDProducer{
  public:
    
    TrackTimeAssoc(const fhicl::ParameterSet&);
    virtual ~TrackTimeAssoc();
    
    void produce(art::Event&);
    void reconfigure(fhicl::ParameterSet const& p);
      
    std::vector<double>               GetMIPHypotheses(trkf::BezierTrack* BTrack, double XOffset=0);
    std::vector<std::vector<double> > ScanMIPHypotheses(trkf::BezierTrack * Btrack);
    void                              PrintHypotheses(std::vector<std::vector<double> > TrackHypotheses);
    double                            GetChi2(std::vector<double> signal, std::vector<double> hypothesis, double UpperLim=0);
    double                            GetMinChi2(std::vector<std::vector<double> > ScannedHypotheses, std::vector<double> FlashShape);

    
    void StoreFlashMatches(std::vector<art::Ptr<recob::Track> >& Tracks, std::vector<art::Ptr<recob::OpFlash> >& Flashes, std::vector<anab::FlashMatch>& Matches, art::Event& evt);

    
    void beginJob();
    
    
  private:
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    int         fBezierResolution;
    int         fPairingMode;
    double      fLengthCut;
    double      fPECut;
  };

  

}

#endif




////////////////////////////////////////////////////////////////////////
/// \file  TrackTimeAssoc_module.cc
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(TrackTimeAssoc)

}//end namespace opdet
////////////////////////////////////////////////////////////////////////



//  This module makes a set of hypotheses for track light profiles
//   based on beizer tracks in the event.
//
//  This module will likely be subsumed into the OpFlashFinder module
//   eventually. But at present they are being developped in parallel.
//
//
// \file TrackTimeAssoc.cxx
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

  TrackTimeAssoc::TrackTimeAssoc(fhicl::ParameterSet const& pset)
  {

    produces< std::vector<anab::FlashMatch> >();
    produces< art::Assns<recob::Track, anab::FlashMatch> >();
    produces< art::Assns<recob::OpFlash, anab::FlashMatch> >();

    this->reconfigure(pset);
   }


  //-------------------------------------------------

  void TrackTimeAssoc::reconfigure(fhicl::ParameterSet const& pset)
  {
    fTrackModuleLabel = pset.get<std::string>("TrackModuleLabel");   
    fFlashModuleLabel = pset.get<std::string>("FlashModuleLabel");
    fBezierResolution = pset.get<int>("BezierResolution");
    fLengthCut        = pset.get<double>("LengthCut");
    fPECut            = pset.get<double>("PECut");
    fPairingMode      = pset.get<int>("PairingMode");
  

  }


  //-------------------------------------------------

  void TrackTimeAssoc::beginJob()
  {
  }



  //-------------------------------------------------

  TrackTimeAssoc::~TrackTimeAssoc()
  {
  }


  //-------------------------------------------------
  std::vector<std::vector<double> > TrackTimeAssoc::ScanMIPHypotheses(trkf::BezierTrack * Btrack)
  {

    double MinX   = 0;    //temporary
    double MaxX   = 250;  //temporary
    size_t XSteps = 50;  //temporary

    art::ServiceHandle<geo::Geometry> geom;

    std::vector<std::vector<double> > ReturnVector(XSteps);
    for(size_t i=0; i!=XSteps; ++i)
      ReturnVector[i].resize(geom->NOpDets());
    
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;

    float TrackLength = Btrack->Length();
    float OldVertex   = Btrack->Vertex()[0];
    
    std::vector<bool> ValidTrajectory(XSteps, true);
    
    double xyz[3];
    for (int b=0; b!=fBezierResolution; b++)
      {
	art::ServiceHandle<util::LArProperties> larp;

	double MIPYield   = larp->ScintYield();
	double QE         = 0.01;
	double MIPdQdx    = 2.1;
	double PromptFrac = 0.25;
	double PromptMIPScintYield = MIPYield * QE * MIPdQdx * PromptFrac;
	
	//	std::cout<<"check : " << PromptMIPScintYield<<std::endl;
	
	float s               = float(b) / float(fBezierResolution);
	float LightAmount     = PromptMIPScintYield * TrackLength/float(fBezierResolution);
	
	
	for(size_t i=0; i!=XSteps; ++i)
	  {
	    if(ValidTrajectory[i])
	      {
		float NewVertex = MinX + float(i)/float(XSteps)*(MaxX-MinX);	
		Btrack->GetTrackPoint(s,xyz);
		xyz[0] += NewVertex-OldVertex;
		
		if((xyz[0] > MaxX) || (xyz[0] < MinX) ) ValidTrajectory[i]=false; 
		
		const std::vector<float>* PointVisibility = pvs->GetAllVisibilities(xyz);
		
		for(size_t OpDet =0; OpDet!=PointVisibility->size();  OpDet++)
		  {
		    ReturnVector.at(i).at(OpDet) += PointVisibility->at(OpDet) * LightAmount;
		  }
	      }
	  }
      }
    for(size_t i=0; i!=ReturnVector.size(); ++i)
      if(!ValidTrajectory[i]) ReturnVector.at(i).clear();
    return ReturnVector;
  }

   



  //-------------------------------------------------
  // This method gets the minimum chi2 achievable by marginalizing over the track x position
  // The track dQdx is allowed to vary from the MIP value and upwards
  //
  double TrackTimeAssoc::GetMinChi2(std::vector<std::vector<double> > ScannedHypotheses, std::vector<double> FlashShape)
  {
    double MinChi2  = 10000;
    if(FlashShape.size()==0) return MinChi2;
    
    for(size_t i=0; i!=ScannedHypotheses.size(); ++i)
      {
	if(ScannedHypotheses.at(i).size()>0)
	  {
	    double Chi2 = GetChi2(FlashShape, ScannedHypotheses.at(i), MinChi2);
	    if(Chi2 < MinChi2)
	      {
		MinChi2 = Chi2;
	      }
	  }
      }
    return MinChi2;
  }

  //-------------------------------------------------


  // Get a hypothesis for the light collected for a bezier track
  std::vector<double> TrackTimeAssoc::GetMIPHypotheses(trkf::BezierTrack* Btrack, double XOffset)
  {
    art::ServiceHandle<geo::Geometry> geom;
    std::vector<double> ReturnVector(geom->NOpDets(),0);
    
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;

    float TrackLength = Btrack->GetLength();

    double xyz[3];
    for (int b=0; b!=fBezierResolution; b++)
      {
	float s               = float(b) / float(fBezierResolution);
	float dQdx            = 2.1;    // Assume MIP value
	
	Btrack->GetTrackPoint(s,xyz);
	xyz[0]+=XOffset;
	const std::vector<float>* PointVisibility = pvs->GetAllVisibilities(xyz);
	float LightAmount = dQdx*TrackLength/float(fBezierResolution);
	
	for(size_t OpDet =0; OpDet!=PointVisibility->size();  OpDet++)
	  {
	    ReturnVector.at(OpDet)+= PointVisibility->at(OpDet) * LightAmount;
	  }
      }
    return ReturnVector;
  }

  //-------------------------------------------------


  void TrackTimeAssoc::produce(art::Event& evt)
  {
    
    //    int EventID = evt.id().event();
    
    
    // Read in flashes from the event
    art::Handle< std::vector<recob::OpFlash> > flashh;
    evt.getByLabel(fFlashModuleLabel, flashh);
    std::vector<art::Ptr<recob::OpFlash> > Flashes;
    for(unsigned int i=0; i < flashh->size(); ++i)
      {
	art::Ptr<recob::OpFlash> flash(flashh,i);
        if(flash->TotalPE()>fPECut) Flashes.push_back(flash);
      }

    // Read in tracks from the event
    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);
    std::vector<art::Ptr<recob::Track> >  Tracks;
    for(unsigned int i=0; i < trackh->size(); ++i)
      {
	art::Ptr<recob::Track> track(trackh,i);
      	if(track->Length()>fLengthCut) Tracks.push_back(track);
      }

    std::sort(Tracks.begin(), Tracks.end(), TrackTimeAssoc_tracksort);

    // Use these to produce Bezier tracks
    std::vector<bool> TracksToCut(Tracks.size(),false);
    std::vector<trkf::BezierTrack*> BTracks;
    BTracks.clear();
    for(size_t i=0; i!=Tracks.size(); i++)
      BTracks.push_back(new trkf::BezierTrack(*Tracks.at(i)));
        
    art::ServiceHandle<geo::Geometry> geom;
    size_t NOpDets = geom->NOpDets();
    
    std::map<int, bool> OnBeamFlashes;
    
    std::vector<std::vector<std::vector<double> > > TrackHypotheses;
    std::vector<std::vector<double> > FlashShapes;
       

    // For each track
    for (size_t i=0; i!=BTracks.size(); ++i)
      {
	TrackHypotheses.push_back(ScanMIPHypotheses(BTracks.at(i)));
      }
    

    for(size_t f=0; f!=Flashes.size(); ++f)
      {

	std::vector<double> ThisFlashShape(NOpDets,0);
	    
	//	if(Flashes.at(f)->InBeamFrame())
	//	  {
        for (unsigned int c = 0; c < geom->NOpChannels(); c++){
          unsigned int o = geom->OpDetFromOpChannel(c);
	  ThisFlashShape[o]+=Flashes.at(f)->PE(c);
        }
	if(Flashes.at(f)->OnBeamTime()) OnBeamFlashes[f]=true;
	//	  }
	FlashShapes.push_back(ThisFlashShape);

      }

    // This map stores the Chi2 for every match:
    //    Chi2Map[Track][Flash] = chi2
    std::map<int, std::map<int, double> > Chi2Map;

    // This map sorts the preferences of each track for each flash
    //    SortedPrefs[TrackID][Chi2] = {flashid1, flashid2, ...}
    std::vector<std::map<double, std::vector<int> > > SortedPrefs(Tracks.size());
    

    for(size_t i=0; i!=TrackHypotheses.size(); ++i)
      {
	for(size_t j=0; j!=FlashShapes.size(); ++j)	    
	  {
	    double Chi2 = GetMinChi2(TrackHypotheses.at(i), FlashShapes.at(j));
	    
	    Chi2Map[i][j]=Chi2;
	    SortedPrefs[i][Chi2].push_back(j);

	  }
      }


    // This will hold the list of matches
    std::vector<anab::FlashMatch> Matches;
	


    if(fPairingMode==0)
      {
	// In pairing mode 0, store all combinatoric matches and let the user
	//   deal with Chi2s later
	for(size_t i=0; i!=BTracks.size(); ++i)
	  for(size_t j=0; j!=Flashes.size(); ++j)
 	    Matches.push_back( anab::FlashMatch(Chi2Map[i][j], j, i, (Flashes.at(j)->OnBeamTime()>0)));
      }
    
    else if(fPairingMode==1)
      {
 	// In pairing mode 1, use the stable marriage algorithm to make a guess
 	//   at good 1<->1 pairings
	
	bool StillPairing =true;
 	std::vector<int> FlashesPaired(Flashes.size(),-1);
 	std::vector<int> TracksPaired(BTracks.size(),-1);
	
        // If we made a new match in the last round, don't stop
	while(StillPairing)
	  {
	    StillPairing=false;
	    for(size_t i=0; i!=BTracks.size(); ++i)
	      {
		// If this track still to be paired
		if(TracksPaired[i]<0)
		  {
		    // Find the flash with best remaining chi2 
		    bool MadeMatch  = false;
		    for(auto itPref = SortedPrefs[i].begin(); itPref!=SortedPrefs[i].end(); ++itPref)
		      {
			for(size_t iflash =0; iflash!=itPref->second.size(); ++iflash)
			  {
			    int FlashID           = itPref->second.at(iflash);
			    int FlashExistingPref = FlashesPaired[FlashID];
			    if(FlashExistingPref < 0)
			      {
				// This flash is available - make pairing
				TracksPaired[i]         = FlashID;
				
				// if the flash is on beam, claim to be 
				//  satisfied, but don't occupy flash
				// if flash is cosmic, claim it.
				if(!OnBeamFlashes[FlashID])
				  FlashesPaired[FlashID]  = i;
				
				StillPairing = true;
				MadeMatch    = true;
			      }
			    else
			      {
				// This flash taken - flash gets to vote
				if(Chi2Map[i][FlashID] < Chi2Map[FlashExistingPref][FlashID])
				  {
				    // If the flash prefers the new guy, switch
				    FlashesPaired[FlashID]          = i;
				    TracksPaired[i]                 = FlashID;
				    TracksPaired[FlashExistingPref] = -1;
				    MadeMatch    = true;
				    StillPairing = true;
				    break;
				  }
				// or else just roll on...			    
			      }
			  }
			if(MadeMatch)  break;
		      } // end loop over chi2s
		  } // end if unpaired track
	      } // end loop over tracks
	  } // end loop until no more pairing
	
		
	for(size_t i=0; i!=BTracks.size(); ++i)
	  {
	    if(TracksPaired[i]>0)
	      {
		int TrackID = i;
		int FlashID = TracksPaired[i];
		
		Matches.push_back( anab::FlashMatch(Chi2Map[TrackID][FlashID], FlashID, TrackID, Flashes.at(FlashID)->OnBeamTime() ));
	      }
	  }
      }


    StoreFlashMatches(Tracks, Flashes, Matches, evt);
    
  }


  //--------------------------------------------------
  // Get the chi2 between a flash hypothesis and a track 
  //  Optional : stop counting if Chi2 bigger than UpperLim
  //
  double TrackTimeAssoc::GetChi2(std::vector<double> signal, std::vector<double> hypothesis, double UpperLim)
  {
       
    double SignalIntegral = 0;
    double HypoIntegral = 0;

    for(size_t i=0; i!=signal.size(); ++i)
      {
	SignalIntegral+=signal.at(i);
	HypoIntegral+=hypothesis.at(i);
      }

    // If light yield indicates >1 MIP light yield, 
    //   Normalize the hypothesis to the signal size.
    //   We do not allow <1 MIP hypotheses.

    double NormFactor = SignalIntegral/HypoIntegral;
    if(NormFactor > 1) 
      {
	for(size_t i=0; i!=hypothesis.size(); ++i)
	  {
	    hypothesis.at(i) *= NormFactor;
	  }
      }
    
    double Chi2=0;
    
    for(size_t i=0; i!=signal.size(); ++i)
      {    
	// We assume width = sqrt(hypothesis)

	if(hypothesis.at(i)!=0)
	  Chi2 += pow(hypothesis.at(i) - signal.at(i),2)/(hypothesis.at(i));	
	
	if(Chi2>UpperLim)
	  break;
      }
    return Chi2;
  }
  

  //--------------------------------------------------

  void TrackTimeAssoc::PrintHypotheses(std::vector<std::vector<double> > TrackHypotheses)
  {
    // List the light hypotheses per track, per PMT
    for (size_t i=0; i!=TrackHypotheses.size(); ++i)
      {
	mf::LogVerbatim("TrackTimeAssoc")<< "Visbility for track " << i <<std::endl;

	for(size_t j=0; j!=TrackHypotheses.at(i).size(); ++j)
	  {
	    mf::LogVerbatim("TrackTimeAssoc") << "Signal at PMT " << j << ", "  << TrackHypotheses.at(i).at(j)<<std::endl;
	  }
      }
    
  }



  //--------------------------------------------------

  void TrackTimeAssoc::StoreFlashMatches(std::vector<art::Ptr<recob::Track> > & Tracks, std::vector<art::Ptr<recob::OpFlash> > & Flashes, std::vector<anab::FlashMatch>& Matches, art::Event& evt)
  {
    std::unique_ptr< std::vector<anab::FlashMatch> > flash_matches ( new std::vector<anab::FlashMatch>);
    std::unique_ptr< art::Assns<recob::Track, anab::FlashMatch > > assn_track( new art::Assns<recob::Track, anab::FlashMatch>);
    std::unique_ptr< art::Assns<recob::OpFlash, anab::FlashMatch > > assn_flash( new art::Assns<recob::OpFlash, anab::FlashMatch>);

    for(size_t i=0; i!=Matches.size(); ++i)
      {
	flash_matches->push_back(Matches.at(i));
	
	util::CreateAssn(*this, evt, *(flash_matches.get()), Tracks.at(Matches.at(i).SubjectID()), *(assn_track.get()), i); 

	util::CreateAssn(*this, evt, *(flash_matches.get()), Flashes.at(Matches.at(i).FlashID()), *(assn_flash.get()), i); 
      }
    

    evt.put(std::move(flash_matches));
    evt.put(std::move(assn_track));
    evt.put(std::move(assn_flash));
  }


  //-----------------------------------
  // Helper function for performing track length sort
  bool TrackTimeAssoc_tracksort(art::Ptr<recob::Track> t1, art::Ptr<recob::Track> t2)
  {
    return (t1->Length()>t2->Length());
    
  }


}


