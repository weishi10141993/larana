/////////////////////////////////////////////////////////////////
//  \fileMVAAlg.h
//  m.haigh@warwick.ac.uk
////////////////////////////////////////////////////////////////////
#ifndef MVAAlg_H
#define MVAAlg_H
#include <vector>
#include <map>
//#include <cmath>
#include <iostream>
//#include <stdint.h>

#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h" 

//#include "SimpleTypesAndConstants/geo_types.h"
//#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
//#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
//#include "RecoAlg/APAGeometryAlg.h"
#include "AnalysisAlg/CalorimetryAlg.h"
//#include "MCCheater/BackTracker.h"

//#include "TMatrixD.h"
//#include "TVectorD.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include <Math/Vector3D.h>
#include "TMVA/Reader.h"

#include "AnalysisBase/MVAPIDResult.h"


namespace mvapid{



  //--------------------------------------------------------------- 
  class MVAAlg {
  public:

struct SortedTrack{
  const art::Ptr<recob::Track> track;
  TVector3 trackStart, trackEnd, trackDir;
  double trackLength;
  std::map<double,const art::Ptr<recob::Hit> > hitMap;
};

struct SumDistance2 {
  // the TGraph is a data member of the object
 TGraph2D* fGraph;
  
 SumDistance2(TGraph2D* g) : fGraph(g) {}
  
  // implementation of the function to be minimized
  double operator() (const double * p) {
 
    ROOT::Math::XYZVector x0(p[0], p[2], 0. ); 
    ROOT::Math::XYZVector u(p[1],p[3], 1. ); 
    u=u.Unit();
    double * x = fGraph->GetX();
    double * y = fGraph->GetY();
    double * z = fGraph->GetZ();
    int npoints = fGraph->GetN();
    double sum = 0;
    for (int i  = 0; i < npoints; ++i) {
      ROOT::Math::XYZVector xp(x[i],y[i],z[i]); 
      sum += ((xp-x0).Cross(u)) .Mag2();
    }
    return sum;
  }
};
        
    MVAAlg(fhicl::ParameterSet const& pset, const art::EDProducer* parentModule);
    
    virtual ~MVAAlg();
    
    void reconfigure(fhicl::ParameterSet const& p);

    void RunPID(art::Event& evt,std::vector<anab::MVAPIDResult>& result,
		art::Assns<recob::Track, anab::MVAPIDResult, void>& assns);

  private:
 
    void PrepareEvent(const art::Event& event);
    
    void FitAndSortTrack(art::Ptr<recob::Track> track,
			 SortedTrack& sortedTrack);
    
    void _Var_EValRatio(const art::Ptr<recob::Track>& track,double& eValRatio);
    
    void _Var_Shape(const SortedTrack& track,
		    double& coreHaloRatio,double& concentration,
		    double& conicalness);
    
    //void _Var_Calo(const SortedTrack& track);			       
    
    double CalcSegmentdEdxFrac(const SortedTrack& track,double start,double end);
    
    double CalcSegmentdEdxDist(const SortedTrack& track,double start,double end);
    
    int LinFit(const art::Ptr<recob::Track> track,TVector3& trackPoint,TVector3& trackDir);
    
    const calo::CalorimetryAlg fCaloAlg;
    
    double fEventT0;

    const art::EDProducer* fParentModule;

    std::string fTrackLabel;
    std::string fHitLabel;
    std::string fSpacePointLabel;

    std::vector<art::Ptr<recob::Track> > fTracks;
    std::vector<art::Ptr<recob::SpacePoint> > fSpacePoints;
    std::vector<art::Ptr<recob::Hit> > fHits;
  
    std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<recob::Hit> > > fTracksToHits;
    std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<recob::SpacePoint> > > fTracksToSpacePoints;
    std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint> > fHitsToSpacePoints;
    std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit > > fSpacePointsToHits;

    anab::MVAPIDResult fResHolder;

    TMVA::Reader fReader;
    
    std::vector<std::string> fMVAMethods;
    std::vector<std::string> fWeightFiles;

  }; // class MVAAlg

} // namespace mvapid

#endif // ifndef MVAAlg_H
