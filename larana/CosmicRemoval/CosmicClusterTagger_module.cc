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

#include <string>
#include <vector>
#include <memory> // std::unique_ptr<>
#include <utility> // std::pair<>, std::move()
#include <algorithm> // std::minmax() ...

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/geo.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "larreco/RecoAlg/SpacePointAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

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
  void reconfigure(fhicl::ParameterSet const & p) ;
  void endJob() override;




private:

  // Declare member data here.

//  int         fReadOutWindowSize;
  float       fSamplingRate;
  std::string fClusterModuleLabel;
  int         fDetectorWidthTicks;
  int         fTickLimit;
  int         fMinTickDrift;
  int         fMaxTickDrift;
};






cosmic::CosmicClusterTagger::CosmicClusterTagger(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  this->reconfigure(p);

  // Call appropriate Produces<>() functions here.
  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Cluster, anab::CosmicTag> >();



}

cosmic::CosmicClusterTagger::~CosmicClusterTagger() {
  // Clean up dynamic memory and other resources here.
}

void cosmic::CosmicClusterTagger::produce(art::Event & e) {
  // Implementation of required member function here.


  std::unique_ptr< std::vector< anab::CosmicTag > > cosmicTagClusterVector( new std::vector<anab::CosmicTag> );
  std::unique_ptr< art::Assns<recob::Cluster, anab::CosmicTag> >    assnOutCosmicTagCluster( new art::Assns<recob::Cluster, anab::CosmicTag>);



  art::Handle<std::vector<recob::Cluster> > Cluster_h;
  e.getByLabel( fClusterModuleLabel, Cluster_h );
  std::vector<art::Ptr<recob::Cluster> > ClusterVec;
  art::fill_ptr_vector(ClusterVec, Cluster_h);

  /////////////////////////////////
  // LOOP OVER CLUSTERS
  /////////////////////////////////

  for( unsigned int iCluster = 0; iCluster < Cluster_h->size(); iCluster++ ) {

    float cosmicScore = 0;
    anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kUnknown;

    art::Ptr<recob::Cluster> tCluster = ClusterVec.at(iCluster);
    art::Ptr<recob::Track> tTrack;

     std::vector<float> endPt1;
     std::vector<float> endPt2;


    // Doing some checks on the cluster to determine if it's a cosmic
     bool failClusterTickCheck = false;
     
     //int timeLimit = 0;//5;
     const std::pair<double, double> t_minmax
       = std::minmax(tCluster->StartTick(), tCluster->EndTick());
     double t0 = t_minmax.first; // minimum
     double t1 = t_minmax.second; // maximum
     //     if( t0+fTickLimit < fDetectorWidthTicks ) { // This is into the pre-spill window
     if( t0+fTickLimit < fMinTickDrift ) { // This is into the pre-spill window
       failClusterTickCheck = true;
     }
     //if( t1-fTickLimit > 2*fDetectorWidthTicks ) { // This is into the post-spill window
     if( t1-fTickLimit > fMaxTickDrift ) { // This is into the post-spill window
       failClusterTickCheck = true;
     }
     
     if(failClusterTickCheck) {
       cosmicScore=1.;
       tag_id = anab::CosmicTagID_t::kOutsideDrift_Partial;
     }

     if( endPt1.size()<1 ) {
       for(int s=0; s<3; s++ ) {
	 endPt1.push_back( -999 );
	 endPt2.push_back( -999 );
       }
     }
     
     
     // Making stuff to save!
     //std::cerr << "Cosmic Score, isCosmic, t0, t1: " << cosmicScore << " " << isCosmic << " t's: " << t0 << " " << t1 << " | " << fReadOutWindowSize<< " | " << fDetectorWidthTicks << std::endl;
     cosmicTagClusterVector->emplace_back( endPt1,
					   endPt2,
					   cosmicScore,
					   tag_id);
     
     util::CreateAssn(*this, e, *cosmicTagClusterVector, tCluster, *assnOutCosmicTagCluster );
     
  }

  /////////////////////////////////
  // END OF CLUSTER LOOP
  /////////////////////////////////


  e.put( std::move(cosmicTagClusterVector) );
  e.put( std::move(assnOutCosmicTagCluster) );




} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////










void cosmic::CosmicClusterTagger::beginJob() {

}



void cosmic::CosmicClusterTagger::reconfigure(fhicl::ParameterSet const & p) {
  // Implementation of optional member function here.
 
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const* geo = lar::providerFrom<geo::Geometry>();

  fSamplingRate = detp->SamplingRate();
  fClusterModuleLabel = p.get< std::string >("ClusterModuleLabel", "cluster");
  fTickLimit = p.get< int >("TickLimit", 0);
  const double driftVelocity = detp->DriftVelocity( detp->Efield(), detp->Temperature() ); // cm/us

  //  std::cerr << "Drift velocity is " << driftVelocity << " cm/us.  Sampling rate is: "<< fSamplingRate << " detector width: " <<  2*geo->DetHalfWidth() << std::endl;
  fDetectorWidthTicks = 2*geo->DetHalfWidth()/(driftVelocity*fSamplingRate/1000); // ~3200 for uB
  //  std::cerr << fDetectorWidthTicks<< std::endl;
  fMinTickDrift = p.get("MinTickDrift", 3200);
  fMaxTickDrift = fMinTickDrift + fDetectorWidthTicks;


}

void cosmic::CosmicClusterTagger::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(cosmic::CosmicClusterTagger)
