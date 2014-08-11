////////////////////////////////////////////////////////////////////////
// Class:       CosmicTagger
// Module Type: producer
// File:        CosmicTagger_module.cc
//
// Generated at Mon Sep 24 18:21:00 2012 by Sarah Lockwitz using artmod
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CosmicTagger
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
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
  class CosmicClusterTagger;
  class SpacePoint;
  class Track;
}



class cosmic::CosmicClusterTagger : public art::EDProducer {
public:
  explicit CosmicClusterTagger(fhicl::ParameterSet const & p);
  virtual ~CosmicClusterTagger();

  void produce(art::Event & e) override;

  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void endJob() override;
  void fillOutOfSpillClusters();
  void doClusterCheck(  art::FindManyP<recob::Hit> tHits, art::FindManyP<recob::SpacePoint> tSpacePoints   );// override;
  void doTrackClusterCheck( std::vector< art::Ptr< recob::Cluster> > ClusterVect, std::vector<double> &t1Times, 
			    std::vector<double> &t2Times, std::vector<int> &fail );
  void doSomeSpacePointStuff(art::FindManyP<recob::SpacePoint> sptsSpill , int iTrack, int &fa, int &fb);


private:

  // Declare member data here.

  int   fReadOutWindowSize;
  //  int fNumberTimeSamples;
  float fSamplingRate;

  // stuff to set in the fcl file
  float fTotalBoundaryLimit; // 15
  float f3DSpillDistance;    // 12
  int   fSpillVetoCtr;       // 2
  int   fdTLimit;            // 8
  int   fdWLimit;            // 8
  std::string fTrackModuleLabel;
  std::string fTrackAssocToClusterModuleLabel;
  std::string fClusterModuleLabel;
  int fDoTrackCheck;
  int fDoClusterCheck;
  int fClusterAssociatedToTracks;
  int fDetectorWidthTicks;




  std::vector<double> tdDist_window_post, tdDist_window_pre;

  int failClusterCheck;








};


int nEvent;

std::vector<float> all_x, all_y, all_z; 




cosmic::CosmicClusterTagger::CosmicClusterTagger(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  this->reconfigure(p);

  // Call appropriate Produces<>() functions here.

  //  produces< std::vector<recob::Track> >();
  //  produces< std::vector<recob::Cluster> >();
  produces< std::vector<anab::CosmicTag> >();
  //  produces< std::vector<anab::CosmicTag> >("tracks");
  produces< art::Assns<anab::CosmicTag, recob::Track> >();


  produces< art::Assns<anab::CosmicTag, recob::Cluster> >();
  //  produces< art::Assns<anab::CosmicTag, recob::Hit> >();


}

cosmic::CosmicClusterTagger::~CosmicClusterTagger() {
  // Clean up dynamic memory and other resources here.
}

void cosmic::CosmicClusterTagger::produce(art::Event & e) {
  // Implementation of required member function here.




  //  std::unique_ptr< std::vector<recob::Track > >      outTracksForTags( new std::vector<recob::Track>);
  //  std::unique_ptr< std::vector<recob::Cluster> >     outClusters( new std::vector<recob::Cluster> );
  std::unique_ptr< std::vector< anab::CosmicTag > > cosmicTagTrackVector( new std::vector<anab::CosmicTag> );
  std::unique_ptr< std::vector< anab::CosmicTag > > cosmicTagClusterVector( new std::vector<anab::CosmicTag> );
  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag > >    assnOutCosmicTagTrack( new art::Assns<recob::Track, anab::CosmicTag>);
  std::unique_ptr< art::Assns<recob::Cluster, anab::CosmicTag > >    assnOutCosmicTagCluster( new art::Assns<recob::Cluster, anab::CosmicTag>);
  //  std::unique_ptr< art::Assns<anab::CosmicTag, recob::Hit > >    assnOutCosmicTagHit( new art::Assns<anab::CosmicTag, recob::Hit>);




  TStopwatch ts;





  art::Handle<std::vector<recob::Track> > Trk_h;

  if(fDoTrackCheck) e.getByLabel( fTrackModuleLabel, Trk_h );
  std::vector<art::Ptr<recob::Track> > TrkVec;
  art::fill_ptr_vector(TrkVec, Trk_h);



  

  failClusterCheck = 0;
  tdDist_window_post.clear();
  tdDist_window_pre .clear();

  for( int hh=0; hh<3; hh++ ) {
    tdDist_window_post.push_back(-999);
    tdDist_window_pre .push_back(-999);
  }





  art::ServiceHandle<geo::Geometry> geo;

  art::Handle<std::vector<recob::Cluster> > Cluster_h;
  e.getByLabel( fClusterModuleLabel, Cluster_h );
  std::vector<art::Ptr<recob::Cluster> > ClusterVec;
  art::fill_ptr_vector(ClusterVec, Cluster_h);

  /////////////////////////////////
  // LOOP OVER CLUSTERS
  /////////////////////////////////



  //  art::FindManyP<recob::Track> tracksFromCluster (Cluster_h, e, fTrackAssocToClusterModuleLabel);



  for( unsigned int iCluster = 0; iCluster < Cluster_h->size(); iCluster++ ) {

    int isCosmic = 0;

    art::Ptr<recob::Cluster> tCluster = ClusterVec.at(iCluster);
    art::Ptr<recob::Track> tTrack;

     std::vector<float> endPt1;
     std::vector<float> endPt2;

//    /////////////////////////////////////////////////////////////
//    // Check to see if the cluster is associated to a track
//    /////////////////////////////////////////////////////////////
//     if( tracksFromCluster.isValid() ) {
//       std::cerr << "How many tracks are associated to the cluster? " << tracksFromCluster.at(iCluster).size() << std::endl;
//       if( tracksFromCluster.at(iCluster).size() > 0 ) {
//	 
//	 if( tracksFromCluster.at(iCluster).size() > 1 ) { 
//	   std::cerr << "Well this is unexpected... The cluster is associated to " << tracksFromCluster.at(iCluster).size() << " tracks." << std::endl; 
//	 }
//	 
//	 std::vector< TVector3 > tempEP;
//	 
//	 for( unsigned int iTrk=0; iTrk < tracksFromCluster.at(iCluster).size(); iTrk++ ) { 
//	   
//	   tTrack = tracksFromCluster.at(iCluster).at(iTrk);
//	   //     if( tracksFromCluster.at(iCluster) ) {
//	   //       tTrack = tracksFromCluster.at(iCluster);
//	   
//	   
//	   TVector3 tVector1 = tTrack->Vertex();
//	   TVector3 tVector2 = tTrack->End();
//	   
//	   tempEP.push_back( tTrack->Vertex() );
//	   tempEP.push_back( tTrack->End() );
//	   
//	   
//	   //	 float trackEndPt1_X = tVector1[0]; 
//	   //	 float trackEndPt1_Y = tVector1[1]; 
//	   //	 float trackEndPt1_Z = tVector1[2]; 
//	   //	 float trackEndPt2_X = tVector2[0]; 
//	   //	 float trackEndPt2_Y = tVector2[1]; 
//	   //	 float trackEndPt2_Z = tVector2[2]; 
//	   
//	   
//	   if( tVector1[0] != tVector1[0] ){
//	     std::cerr << "!!! FOUND A PROBLEM... the length is: " << tTrack->Length() << 
//	       " np: " << tTrack->NumberTrajectoryPoints() << " id: " << tTrack->ID() << " " << tTrack << std::endl;
//	     for( size_t hh=0; hh<tTrack->NumberTrajectoryPoints(); hh++) {
//	       std::cerr << hh << " " << tTrack->LocationAtPoint(hh)[0] << ", " <<
//		 tTrack->LocationAtPoint(hh)[1] << ", " <<
//		 tTrack->LocationAtPoint(hh)[2] << std::endl;
//	     }
//	     continue; // I don't want to deal with these "tracks"
//	   }
//	   
//	 }// end of loop over tracks
//	 
//	 float maxDistance = -999;
//	 int indxEP1 = -1, indxEP2 =-2;
//	 for( unsigned int pp=0; pp<tempEP.size()-1; pp++ ) {
//	   for( unsigned int ww=pp+1; ww<tempEP.size(); ww++ ) {
//	     float distance = sqrt( pow( tempEP.at(pp)[0] - tempEP.at(ww)[0], 2 ) +
//				    pow( tempEP.at(pp)[1] - tempEP.at(ww)[1], 2 ) +
//				    pow( tempEP.at(pp)[2] - tempEP.at(ww)[2], 2 ) );
//	     if( distance > maxDistance ) {
//	       indxEP1 = tempEP.at(pp)[1] > tempEP.at(ww)[1] ? pp : ww;
//	       indxEP2 = tempEP.at(pp)[1] > tempEP.at(ww)[1] ? ww : pp;
//	       maxDistance = distance;
//	     }
//	   }
//	 }
//	 
//	 
//	 
//	 endPt1.push_back( tempEP.at(indxEP1)[0] );
//	 endPt1.push_back( tempEP.at(indxEP1)[1] );
//	 endPt1.push_back( tempEP.at(indxEP1)[2] );
//	 endPt2.push_back( tempEP.at(indxEP2)[0] );
//	 endPt2.push_back( tempEP.at(indxEP2)[1] );
//	 endPt2.push_back( tempEP.at(indxEP2)[2] );
//	 
//	 /////////////////////////////////
//	 // Now check Y & Z boundaries:
//	 /////////////////////////////////
//	 int nBd = 0;
//	 float bndDist = 5;
//	 if(isCosmic==0 ) {
//	   if(fabs( tempEP.at(indxEP1)[1] - geo->DetHalfHeight())<bndDist ) nBd++;
//	   if(fabs( tempEP.at(indxEP2)[1] + geo->DetHalfHeight())<bndDist ) nBd++;
//	   if(fabs( tempEP.at(indxEP1)[2] - geo->DetLength())<bndDist || fabs( tempEP.at(indxEP2)[2]- geo->DetLength())<bndDist ) nBd++;
//	   if(fabs( tempEP.at(indxEP1)[2])<bndDist || fabs( tempEP.at(indxEP2)[2] ) < bndDist ) nBd++;
//	   
//	   if( nBd>1 ) isCosmic = 2;
//	   if( nBd==1 ) isCosmic = 3; // only 1 boundary w/o time information
//	 }
//	 
//	 
//       }
//     } //     if( tracksFromCluster.isValid() )
//       /////////////////////////////////////////////////////////////
//       // End of track check
//       /////////////////////////////////////////////////////////////


    // Doing some checks on the cluster to determine if it's a cosmic
     bool failClusterTickCheck = false;
     
     int timeLimit = 0;//5;
     
     double t0 = tCluster->StartPos()[1] < tCluster->EndPos()[1] ? tCluster->StartPos()[1] : tCluster->EndPos()[1];
     double t1 = tCluster->StartPos()[1] > tCluster->EndPos()[1] ? tCluster->StartPos()[1] : tCluster->EndPos()[1]; 
     if( t0+timeLimit < fDetectorWidthTicks ) { // This is into the pre-spill window
       failClusterTickCheck = true;
     }
     if( t1-timeLimit > 2*fDetectorWidthTicks ) { // This is into the post-spill window
       failClusterTickCheck = true;
     }
     
     if(failClusterTickCheck) isCosmic = 4;
     


     float cosmicScore = isCosmic > 0 ? 1 : 0;
     if(isCosmic == 3 ) cosmicScore = 0.5;

     if( endPt1.size()<1 ) {
       for(int s=0; s<3; s++ ) {
	 endPt1.push_back( -999 );
	 endPt2.push_back( -999 );
       }
     }
     
     
     // Making stuff to save!
     std::cerr << "Cosmic Score, isCosmic, t0, t1: " << cosmicScore << " " << isCosmic << " t's: " << t0 << " " << t1 << " | " << fReadOutWindowSize<< " | " << fDetectorWidthTicks << std::endl;
     cosmicTagClusterVector->push_back( anab::CosmicTag(endPt1,
							endPt2,
							cosmicScore,
							isCosmic
							) );
     
     util::CreateAssn(*this, e, *cosmicTagClusterVector, tCluster, *assnOutCosmicTagCluster );
     
  }

  /////////////////////////////////
  // END OF CLUSTER LOOP
  /////////////////////////////////








  //  e.put( std::move(outTracksForTags) );
  //  e.put( std::move(outClusters) );
  //  e.put( std::move(cosmicTagTrackVector) ,"tracks");
  e.put( std::move(cosmicTagClusterVector) );
  e.put( std::move(assnOutCosmicTagTrack) );
  e.put( std::move(assnOutCosmicTagCluster) );
  //  e.put( std::move(assnOutCosmicTagHit) );



  TrkVec.clear();



} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////







void cosmic::CosmicClusterTagger::doTrackClusterCheck( std::vector< art::Ptr< recob::Cluster> > ClusterVect, 
					      std::vector<double> &t1Times, std::vector<double> &t2Times, 
					      std::vector<int> &fail  ) {
  // Let's have this return the earliest and latest t0 times
  // Also some sort of boolean for if these are outside of the
  // spill window

  t1Times.clear();
  t2Times.clear();
  fail.clear();
  int timeLimit = 0;//5;

  for( unsigned int k = 0; k < ClusterVect.size(); k++ ) {
    fail.push_back(0);

    double t0 = ClusterVect[k]->StartPos()[1] < ClusterVect[k]->EndPos()[1] ? ClusterVect[k]->StartPos()[1] : ClusterVect[k]->EndPos()[1];
    double t1 = ClusterVect[k]->StartPos()[1] > ClusterVect[k]->EndPos()[1] ? ClusterVect[k]->StartPos()[1] : ClusterVect[k]->EndPos()[1]; 

    t1Times.push_back(t0);
    t2Times.push_back(t1);

    if( t0+timeLimit < fReadOutWindowSize ) { // This is into the pre-spill window
      fail.at(k)=1;
    }
    if( t1-timeLimit > 2*fReadOutWindowSize ) { // This is into the post-spill window
      fail.at(k)=1;
    }

    std::cerr << "----------- Times are: " << t0 << ", " << t1 << " fail? " << fail.at(k) << " " << fReadOutWindowSize << std::endl;

  }

  return;
}






void cosmic::CosmicClusterTagger::doSomeSpacePointStuff(art::FindManyP<recob::SpacePoint> sptsSpill , int iTrack, int &face1, int &face2) {


  std::vector<art::Ptr<recob::SpacePoint> > spts = sptsSpill.at(iTrack);
  
  std::cerr << " SPACE POINTS SIZE IS : " << spts.size() << std::endl;
  //      if(spts.size()<3 ) continue;
  
  art::ServiceHandle<geo::Geometry> geo;
  
  
  if(spts.size()>3) {
    
    all_x.clear();
    all_y.clear();
    all_z.clear();
    
    
    double testXYZ[3]={0.};
    unsigned int testTPC=0, testCRYO=0;
    
    for(unsigned int j=0; j<spts.size(); j++) {
      try{
	testXYZ[0]=spts[j]->XYZ()[0];
	testXYZ[1]=spts[j]->XYZ()[1];
	testXYZ[2]=spts[j]->XYZ()[2];
	geo->PositionToTPC(testXYZ, testTPC, testCRYO);
      }
      catch(cet::exception &e) {
	continue;
      }
      all_x.push_back( spts[j]->XYZ()[0] );
      all_y.push_back( spts[j]->XYZ()[1] );
      all_z.push_back( spts[j]->XYZ()[2] );
    }
    
    // This is a way to see how "wide" a track is
    TVectorD xV;
    float centroidX=0, centroidY=0, centroidZ=0;
    TMatrixD mCluster(all_x.size(), 3);
    
    centroidX = accumulate(all_x.begin(),all_x.end(),0);
    centroidY = accumulate(all_y.begin(),all_y.end(),0);
    centroidZ = accumulate(all_z.begin(),all_z.end(),0);
    centroidX /= 1.0*all_x.size();
    centroidY /= 1.0*all_y.size();
    centroidZ /= 1.0*all_z.size();
    
    
    for(unsigned int ii=0; ii<all_x.size(); ii++) {
      double vals[3];
      vals[0] = all_x[ii] - centroidX;
      vals[1] = all_y[ii] - centroidY;
      vals[2] = all_z[ii] - centroidZ;
      xV.Use(3, vals);
      TMatrixDRow(mCluster, ii) = xV;
    }
    
    
    TMatrixD Aw = mCluster;
    std::cerr << "nrows: " << Aw.GetNrows() << " ncols: " << Aw.GetNcols() << std::endl;
    if(Aw.GetNrows()<3) return;
    TDecompSVD svd(Aw);
    TMatrixD m_svdV = svd.GetV();
    TVectorD c_svd = svd.GetSig();
    
    
    
    double xyz[3], dxyz[3];//, xyzout[6];
    xyz[0] = centroidX;
    xyz[1] = centroidY;
    xyz[2] = centroidZ;
    dxyz[0] = m_svdV(0,0);
    dxyz[1] = m_svdV(1,0);
    dxyz[2] = m_svdV(2,0);
    
    

    /*
    
    // see if the track should pass two y planes
    bool shouldPassTwoYPlanes=true;
    float t = (geo->DetHalfHeight()-xyz[1])/dxyz[1];
    float testX = xyz[0]+dxyz[0]*t;
    float testZ = xyz[2]+dxyz[2]*t;
    if( testX<0 || testX>2.*geo->DetHalfWidth() ) shouldPassTwoYPlanes=false;
    if( testZ<0 || testZ>geo->DetLength() )       shouldPassTwoYPlanes=false;
    t = (geo->DetHalfHeight()-xyz[1])/dxyz[1];
    testX = xyz[0]+dxyz[0]*t;
    testZ = xyz[2]+dxyz[2]*t;
    if( testX<0 || testX>2.*geo->DetHalfWidth() ) shouldPassTwoYPlanes=false;
    if( testZ<0 || testZ>geo->DetLength() )       shouldPassTwoYPlanes=false;
    */
    
    
    double xyzout1[3], xyzout2[3];
    geo::ProjectToBoxEdge( xyz, dxyz, 0, 2.*geo->DetHalfWidth(), -1.*geo->DetHalfHeight(), geo->DetHalfHeight(), 0, geo->DetLength(), xyzout1 );
    double ndxyz[3]; 
    ndxyz[0] = -dxyz[0];
    ndxyz[1] = -dxyz[1];
    ndxyz[2] = -dxyz[2];
    geo::ProjectToBoxEdge( xyz, ndxyz, 0, 2.*geo->DetHalfWidth(), -1.*geo->DetHalfHeight(), geo->DetHalfHeight(), 0, geo->DetLength(), xyzout2 );
    mf::LogInfo("CosmicClusterTaggerAna Info") << "1 project to box edges: " << xyzout1[0] << " " << xyzout1[1] << " " << xyzout1[2] ;
    mf::LogInfo("CosmicClusterTaggerAna Info") << "2 project to box edges: " << xyzout2[0] << " " << xyzout2[1] << " " << xyzout2[2] ;
    if( xyzout1[0]<2 ) face1=2;
    else if( fabs( 2.*geo->DetHalfWidth()-xyzout1[0] ) < 2 ) face1=5;
    else if( fabs( geo->DetHalfHeight() - xyzout1[1] ) < 2 ) face1=3;
    else if( fabs( xyzout1[1] +geo->DetHalfHeight() ) < 2 ) face1=4;
    else if( fabs( xyzout1[2] ) < 2 ) face1=1;
    else if( fabs( geo->DetLength() - xyzout1[2] ) < 2 ) face1=6;
    
    if( xyzout2[0]<2 ) face2=2;
    else if( fabs( 2.*geo->DetHalfWidth()-xyzout2[0] ) < 2 ) face2=5;
    else if( fabs( geo->DetHalfHeight() - xyzout2[1] ) < 2 ) face2=3;
    else if( fabs( xyzout2[1] + geo->DetHalfHeight() ) < 2 ) face2=4;
    else if( fabs( xyzout2[2] ) < 2 ) face2=1;
    else if( fabs( geo->DetLength() - xyzout2[2] ) < 2 ) face2=6;
    
    
    
    if( face1<0 || face2<0 ) {
      mf::LogInfo("CosmicClusterTagger Info") << "!!!!!!!!!!!! Problem projecting to box edge!!!! " ;
      mf::LogInfo("CosmicClusterTagger Info") << xyzout1[0]                                ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( 2.*geo->DetHalfWidth()-xyzout1[0] ) ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( geo->DetHalfHeight() - xyzout1[1] ) ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( xyzout1[1] + geo->DetHalfHeight() ) ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( xyzout1[2] )			    ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( geo->DetLength() - xyzout1[2] )     ;
      
      mf::LogInfo("CosmicClusterTagger Info") << "---- " ;
      mf::LogInfo("CosmicClusterTagger Info") << xyzout2[0]                                ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( 2.*geo->DetHalfWidth()-xyzout2[0] ) ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( geo->DetHalfHeight() - xyzout2[1] ) ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( xyzout2[1] + geo->DetHalfHeight() ) ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( xyzout2[2] )			    ;
      mf::LogInfo("CosmicClusterTagger Info") << fabs( geo->DetLength() - xyzout2[2] )     ;
      
      
    }
  }// spts > 3 so we can to the matrix stuff
  
  
  
} // doSomeSpacePointStuff







void cosmic::CosmicClusterTagger::beginJob() {



}

void cosmic::CosmicClusterTagger::reconfigure(fhicl::ParameterSet const & p) {
  // Implementation of optional member function here.
  
  ////////  fSptalg                = new cosmic::SpacePointAlg(p.get<fhicl::ParameterSet>("SpacePointAlg"));


  art::ServiceHandle<util::DetectorProperties> detp;

  //fReadOutWindowSize = detp->ReadOutWindowSize();
  //fNumberTimeSamples = detp->NumberTimeSamples();
  fSamplingRate = detp->SamplingRate();
  fTotalBoundaryLimit = p.get<float>("TPCBoundaryLimit", 15);
  f3DSpillDistance    = p.get<float>("SpillDistance",12); 

  fSpillVetoCtr = p.get<int>("SpillVetoCounter", 2);
  fdTLimit      = p.get<int>("dTLimit",8);
  fdWLimit      = p.get<int>("dWLimit",8);

  fClusterModuleLabel = p.get< std::string >("ClusterModuleLabel", "cluster");
  fTrackAssocToClusterModuleLabel =p.get< std::string >("TrackAssocToClusterModuleLabel", "track");
  fTrackModuleLabel   = p.get< std::string >("TrackModuleLabel", "track");
  fDoTrackCheck = p.get< int >("DoTrackCheck", 0);
  fDoClusterCheck = p.get< int >("DoClusterCheck", 1);
  fClusterAssociatedToTracks = p.get< int >("ClustersAssociatedToTracks",1);





  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<geo::Geometry> geo;
  const double driftVelocity = larp->DriftVelocity( larp->Efield(), larp->Temperature() ); // cm/us

  std::cerr << "Drift velocity is " << driftVelocity << " cm/us.  Sampling rate is: "<< fSamplingRate << " detector width: " <<  2*geo->DetHalfWidth() << std::endl;
  fDetectorWidthTicks = 2*geo->DetHalfWidth()/(driftVelocity*fSamplingRate/1000); // ~3200 for uB
  std::cerr << fDetectorWidthTicks<< std::endl;



}

void cosmic::CosmicClusterTagger::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(cosmic::CosmicClusterTagger)
