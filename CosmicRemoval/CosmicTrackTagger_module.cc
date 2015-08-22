////////////////////////////////////////////////////////////////////////
// Class:       CosmicTagger
// Module Type: producer
// File:        CosmicTrackTagger_module.cc
//
// Generated at Mon Sep 24 18:21:00 2012 by Sarah Lockwitz using artmod
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CosmicTrackTagger
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "Geometry/Geometry.h"
#include "Geometry/geo.h"

#include "RecoBase/SpacePoint.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Track.h"

#include "AnalysisBase/CosmicTag.h"
#include "RecoAlg/SpacePointAlg.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"

#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include "TTree.h"
#include "TH1.h"
#include "TStopwatch.h"

class TTree;
class TH1;

namespace cosmic {
  class CosmicTrackTagger;
  class SpacePoint;
  class Track;
}


class cosmic::CosmicTrackTagger : public art::EDProducer {
public:
  explicit CosmicTrackTagger(fhicl::ParameterSet const & p);
  virtual ~CosmicTrackTagger();

  void produce(art::Event & e) override;

  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void endJob() override;


private:

  // Length of reconstructed track, trajectory by trajectory.
  double length(art::Ptr<recob::Track> track);

  //  float fTotalBoundaryLimit; // 15
  //  float f3DSpillDistance;    // 12
  //  int   fSpillVetoCtr;       // 2
  ////  int   fdTLimit;            // 8
  //int   fdWLimit;            // 8
  std::string fTrackModuleLabel;
  //  std::string fTrackAssocToClusterModuleLabel;
  //  std::string fClusterModuleLabel;
  //  int fDoTrackCheck;
  //  int fDoClusterCheck;
  //  int fClusterAssociatedToTracks;
  int fDetectorWidthTicks;
  float fTPCXBoundary, fTPCYBoundary, fTPCZBoundary;
  float fDetHalfHeight, fDetWidth, fDetLength;
};


cosmic::CosmicTrackTagger::CosmicTrackTagger(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);

  // Call appropriate Produces<>() functions here.
  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();
}

cosmic::CosmicTrackTagger::~CosmicTrackTagger() {
  // Clean up dynamic memory and other resources here.
}

void cosmic::CosmicTrackTagger::produce(art::Event & e) {
  // Implementation of required member function here.

  std::unique_ptr< std::vector< anab::CosmicTag > >              cosmicTagTrackVector( new std::vector<anab::CosmicTag> );
  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag > >  assnOutCosmicTagTrack( new art::Assns<recob::Track, anab::CosmicTag>);

  TStopwatch ts;

  art::Handle<std::vector<recob::Track> > Trk_h;
  e.getByLabel( fTrackModuleLabel, Trk_h );
  std::vector<art::Ptr<recob::Track> > TrkVec;
  art::fill_ptr_vector(TrkVec, Trk_h);

  /////////////////////////////////
  // LOOPING OVER INSPILL TRACKS
  /////////////////////////////////
  
    art::FindManyP<recob::Hit>        hitsSpill   (Trk_h, e, fTrackModuleLabel);
    //    art::FindManyP<recob::Cluster>    ClusterSpill(Trk_h, e, fTrackModuleLabel);

    for( unsigned int iTrack=0; iTrack<Trk_h->size(); iTrack++ ) {

        int isCosmic    =  0;
        anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kNotTagged;

        art::Ptr<recob::Track>              tTrack  = TrkVec.at(iTrack);
        std::vector<art::Ptr<recob::Hit> >  HitVec  = hitsSpill.at(iTrack);
        
        if (iTrack != tTrack.key())
        {
            std::cout << "Mismatch in track index/key" << std::endl;
        }

        // A BETTER WAY OF FINDING END POINTS:
        TVector3 tVector1 = tTrack->Vertex();
        TVector3 tVector2 = tTrack->End();

        float trackEndPt1_X = tVector1[0];
        float trackEndPt1_Y = tVector1[1];
        float trackEndPt1_Z = tVector1[2];
        float trackEndPt2_X = tVector2[0];
        float trackEndPt2_Y = tVector2[1];
        float trackEndPt2_Z = tVector2[2];


        if( trackEndPt1_X != trackEndPt1_X ||
            trackEndPt1_Y != trackEndPt1_Y ||
            trackEndPt1_Z != trackEndPt1_Z ||
            trackEndPt2_X != trackEndPt2_X ||
            trackEndPt2_Y != trackEndPt2_Y ||
            trackEndPt2_Z != trackEndPt2_Z ) {
                std::cerr << "!!! FOUND A PROBLEM... the length is: " << tTrack->Length() <<
                    " np: " << tTrack->NumberTrajectoryPoints() << " id: " << tTrack->ID() << " " << tTrack << std::endl;
                //for( size_t hh=0; hh<tTrack->NumberTrajectoryPoints(); hh++) {
                //  std::cerr << hh << " " << tTrack->LocationAtPoint(hh)[0] << ", " <<
                //    tTrack->LocationAtPoint(hh)[1] << ", " <<
                //    tTrack->LocationAtPoint(hh)[2] << std::endl;
                //}
                std::vector<float> tempPt1, tempPt2;
                tempPt1.push_back(-999);
                tempPt1.push_back(-999);
                tempPt1.push_back(-999);
                tempPt2.push_back(-999);
                tempPt2.push_back(-999);
                tempPt2.push_back(-999);
                cosmicTagTrackVector->emplace_back( tempPt1, tempPt2, -999, tag_id);
                util::CreateAssn(*this, e, *cosmicTagTrackVector, tTrack, *assnOutCosmicTagTrack );
                continue; // I don't want to deal with these "tracks"
        }



        /////////////////////////////////////
        // Getting first and last ticks
        /////////////////////////////////////
        float tick1 =  9999;
        float tick2 = -9999;
        bool dumpMe(false);
      
        for ( unsigned int p = 0; p < HitVec.size(); p++) {
            if (dumpMe)
            {
                std::cout << "###>> Hit key: " << HitVec[p].key() << ", peak - RMS: " << HitVec[p]->PeakTimeMinusRMS() << ", peak + RMS: " << HitVec[p]->PeakTimePlusRMS() << std::endl;
            }
            if( HitVec[p]->PeakTimeMinusRMS() < tick1 ) tick1 =  HitVec[p]->PeakTimeMinusRMS();
            if( HitVec[p]->PeakTimePlusRMS()  > tick2 ) tick2 =  HitVec[p]->PeakTimePlusRMS();
        }
      

        /////////////////////////////////////////////////////////
        // Are any of the ticks outside of the ReadOutWindow ?
        /////////////////////////////////////////////////////////
        if(tick1 < fDetectorWidthTicks || tick2 > 2*fDetectorWidthTicks ) {
            isCosmic = 1;
            tag_id = anab::CosmicTagID_t::kOutsideDrift_Partial;
        }



        /////////////////////////////////
        // Now check Y & Z boundaries:
        /////////////////////////////////
        int nBdY = 0, nBdZ=0;
        if( isCosmic == 0 ) {

            // Checking lower side of TPC
            if( fabs( fDetHalfHeight + trackEndPt1_Y ) < fTPCYBoundary ||
               fabs( fDetHalfHeight + trackEndPt2_Y ) < fTPCYBoundary ||
               trackEndPt1_Y < -fDetHalfHeight ||
               trackEndPt2_Y < -fDetHalfHeight ) nBdY++;
	
            // Checking upper side of TPC
            if( fabs( fDetHalfHeight - trackEndPt1_Y ) < fTPCYBoundary ||
               fabs( fDetHalfHeight - trackEndPt2_Y ) < fTPCYBoundary ||
               trackEndPt1_Y > fDetHalfHeight ||
               trackEndPt2_Y > fDetHalfHeight ) nBdY++;
	
            if(fabs(trackEndPt1_Z - fDetLength)<fTPCZBoundary || fabs(trackEndPt2_Z - fDetLength) < fTPCZBoundary ) nBdZ++;
            if(fabs(trackEndPt1_Z )< fTPCZBoundary || fabs(trackEndPt2_Z )< fTPCZBoundary ) nBdZ++;
            if( (nBdY+nBdZ)>1 ) {
                isCosmic = 2;
                if(nBdY>1) tag_id = anab::CosmicTagID_t::kGeometry_YY;
                else if(nBdZ>1) tag_id = anab::CosmicTagID_t::kGeometry_ZZ;
                else tag_id = anab::CosmicTagID_t::kGeometry_YZ;
            }
            else if( (nBdY+nBdZ)==1) {
                isCosmic = 3 ;
                if(nBdY==1) tag_id = anab::CosmicTagID_t::kGeometry_Y;
                else if(nBdZ==1) tag_id = anab::CosmicTagID_t::kGeometry_Z;
            }
        }

        std::vector<float> endPt1;
        std::vector<float> endPt2;
        endPt1.push_back( trackEndPt1_X );
        endPt1.push_back( trackEndPt1_Y );
        endPt1.push_back( trackEndPt1_Z );
        endPt2.push_back( trackEndPt2_X );
        endPt2.push_back( trackEndPt2_Y );
        endPt2.push_back( trackEndPt2_Z );

        float cosmicScore = isCosmic > 0 ? 1 : 0;
        if( isCosmic==3 ) cosmicScore = 0.5;


        ///////////////////////////////////////////////////////
        // Doing a very basic check on X boundaries
        // this gets the types of tracks that go through both X boundaries of the detector
        if( fabs( trackEndPt1_X - trackEndPt2_X ) > fDetWidth-fTPCXBoundary ) {
            cosmicScore = 1 ;
            isCosmic = 4;
            tag_id = anab::CosmicTagID_t::kGeometry_XX;
        }

        /*
         //////////////////////////////
         // Now check for X boundary
         //////////////////////////////
         if( isCosmic==0 ) {
         //anab::CosmicTag cctt = anab::CosmicTag(endPt1, endPt2, tag_id, isCosmic );
         int nBdX =0;
         float xBnd1 = -9999; //cctt.getXInteraction(endPt1[0], 2.0*fDetWidth, fReadOutWindowSize, trackTime, std::floor(tick1) );
         float xBnd2 = -9999; //cctt.getXInteraction(endPt1[0], 2.0*fDetWidth, fReadOutWindowSize, trackTime, std::floor(tick2) );
         if(xBnd1 < fTPCXBoundary || xBnd2 < fTPCXBoundary) nBdX++;
         if( ( fDetWidth - xBnd1 < fTPCXBoundary ) || ( fDetWidth - xBnd1 < fTPCXBoundary ) ) nBdX++;
         if(  nBdX+nBdY+nBdZ>1 && 0 ) isCosmic = 3; // THIS ISN'T SETUP YET -- NEED A HANDLE TO TIME INFO
         if( nBd >0 ) {isCosmic=3; cosmicScore = 0.5;}
         }
         */

        cosmicTagTrackVector->emplace_back( endPt1,
					  endPt2,
					  cosmicScore,
					  tag_id);

//      std::cerr << "The IsCosmic value is "<< isCosmic << " end pts " 
//		<< trackEndPt1_X<<","<< trackEndPt1_Y << "," << trackEndPt1_Z<< " | | " 
//		<< trackEndPt2_X<< ","<< trackEndPt2_Y <<"," << trackEndPt2_Z << std::endl;
      
      //outTracksForTags->push_back( *tTrack );

        util::CreateAssn(*this, e, *cosmicTagTrackVector, tTrack, *assnOutCosmicTagTrack );
        //util::CreateAssn(*this, e, *cosmicTagTrackVector, HitVec, *assnOutCosmicTagHit);
    }
    // END OF LOOPING OVER INSPILL TRACKS
    

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //////TAGGING DELTA RAYS (and other stub) ASSOCIATED TO A ALREADY TAGGED COSMIC TRACK//////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    float         dE=0, dS=0, temp=0, IScore=0;
    unsigned int  IndexE = 0, iTrk1=0, iTrk=0;
    anab::CosmicTagID_t IType=anab::CosmicTagID_t::kNotTagged;
        
    for(iTrk=0; iTrk<Trk_h->size(); iTrk++ ){
        art::Ptr<recob::Track> tTrk  = TrkVec.at(iTrk);
        if ((*cosmicTagTrackVector)[iTrk].CosmicScore()==0){
            TVector3 tStart = tTrk->Vertex();
            TVector3 tEnd   = tTrk->End();
            unsigned int l=0;
            for(iTrk1=0; iTrk1<Trk_h->size(); iTrk1++ ){
                art::Ptr<recob::Track> tTrk1  = TrkVec.at(iTrk1);
                float getScore = (*cosmicTagTrackVector)[iTrk1].CosmicScore();
                if (getScore == 1 || getScore == 0.5){
                    anab::CosmicTagID_t getType = (*cosmicTagTrackVector)[iTrk1].CosmicType();
                    TVector3 tStart1 = tTrk1->Vertex();
                    TVector3 tEnd1   = tTrk1->End();
                    TVector3 NumE    = (tEnd-tStart1).Cross(tEnd-tEnd1);
                    TVector3 DenE    = tEnd1-tStart1;
                    dE = NumE.Mag()/DenE.Mag();
                    if (l==0){
                        temp = dE;
                        IndexE = iTrk1;
                        IScore = getScore;
                        IType  = getType;
                    }
                    if (dE<temp){
                        temp = dE;
                        IndexE = iTrk1;
                        IScore = getScore;
                        IType  = getType;
                    }
                    l++;
                }
            }//End Trk1 loop
            art::Ptr<recob::Track> tTrkI = TrkVec.at(IndexE);
            TVector3 tStartI = tTrkI->Vertex();
            TVector3 tEndI   = tTrkI->End();
            TVector3 NumS    = (tStart-tStartI).Cross(tStart-tEndI);
            TVector3 DenS    = tEndI-tStartI;
            dS = NumS.Mag()/DenS.Mag();
            if (((dS<5 && temp<5) || (dS<temp && dS<5)) && (length(tTrk)<60)){
                (*cosmicTagTrackVector)[iTrk].CosmicScore() = IScore-0.05;
                (*cosmicTagTrackVector)[iTrk].CosmicType()  = IType;
          //util::CreateAssn(*this, e, *cosmicTagTrackVector, tTrk, *assnOutCosmicTagTrack, iTrk);
            }
        }//end cosmicScore==0 loop
    }//end iTrk loop	  
 
 /*std::cout<<"\n"<<Trk_h->size()<<"\t"<<(*cosmicTagTrackVector).size();
   for(unsigned int f=0;f<Trk_h->size();f++){
   	std::cout<<"\n\t"<<f<<"\t"<<(*cosmicTagTrackVector)[f].CosmicScore()<<"\t"<<(*cosmicTagTrackVector)[f].CosmicType();
   }*/
 
  // e.put( std::move(outTracksForTags) );
  e.put( std::move(cosmicTagTrackVector) );
  e.put( std::move(assnOutCosmicTagTrack) );

  TrkVec.clear();

} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////


// Length of reconstructed track, trajectory by trajectory.
double cosmic::CosmicTrackTagger::length(art::Ptr<recob::Track> track){
  double result = 0.;
  TVector3 disp = track->LocationAtPoint(0);
  int n = track->NumberTrajectoryPoints();
  for(int i = 1; i < n; ++i) {
    const TVector3& pos = track->LocationAtPoint(i);
    disp -= pos;
    result += disp.Mag();
    disp = pos;
  }
  return result;
}


void cosmic::CosmicTrackTagger::beginJob(){
}


void cosmic::CosmicTrackTagger::reconfigure(fhicl::ParameterSet const & p) {
  // Implementation of optional member function here.
  
  ////////  fSptalg  = new cosmic::SpacePointAlg(p.get<fhicl::ParameterSet>("SpacePointAlg"));
  art::ServiceHandle<util::DetectorProperties> detp;
  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<geo::Geometry> geo;

  fDetHalfHeight = geo->DetHalfHeight();
  fDetWidth      = 2.*geo->DetHalfWidth();
  fDetLength     = geo->DetLength();

  float fSamplingRate = detp->SamplingRate();

  fTrackModuleLabel = p.get< std::string >("TrackModuleLabel", "track");

  fTPCXBoundary = p.get< float >("TPCXBoundary", 5);
  fTPCYBoundary = p.get< float >("TPCYBoundary", 5);
  fTPCZBoundary = p.get< float >("TPCZBoundary", 5);

  const double driftVelocity = larp->DriftVelocity( larp->Efield(), larp->Temperature() ); // cm/us

  //std::cerr << "Drift velocity is " << driftVelocity << " cm/us.  Sampling rate is: "<< fSamplingRate << " detector width: " <<  2*geo->DetHalfWidth() << std::endl;
  fDetectorWidthTicks = 2*geo->DetHalfWidth()/(driftVelocity*fSamplingRate/1000); // ~3200 for uB
}

void cosmic::CosmicTrackTagger::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(cosmic::CosmicTrackTagger)
