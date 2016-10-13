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
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h" 

//#include "SimpleTypesAndConstants/geo_types.h"
//#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
//#include "RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
//#include "RecoAlg/APAGeometryAlg.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
//#include "MCCheater/BackTracker.h"

//#include "TMatrixD.h"
//#include "TVectorD.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include <Math/Vector3D.h>
#include "TMVA/Reader.h"

#include "lardataobj/AnalysisBase/MVAPIDResult.h"


namespace mvapid{



  //--------------------------------------------------------------- 
  class MVAAlg {
  public:

struct SortedObj{
  TVector3 start, end, dir;
  double length;
  std::map<double,const art::Ptr<recob::Hit> > hitMap;
};

struct SumDistance2 {
  // the TGraph is a data member of the object
 TGraph2D* fGraph;
  
 SumDistance2(TGraph2D* g) : fGraph(g) {}
  
  // implementation of the function to be minimized
  double operator() (const double * p) {
    
    ROOT::Math::XYZVector x0(p[0], p[2], p[4] );
    ROOT::Math::XYZVector u(p[1],p[3], p[5] );

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
		art::Assns<recob::Track, anab::MVAPIDResult, void>& trackAssns,
		art::Assns<recob::Shower, anab::MVAPIDResult, void>& showerAssns);

  private:

    int IsInActiveVol(const TVector3& pos);
 
    void PrepareEvent(const art::Event& event);
    
    void FitAndSortTrack(art::Ptr<recob::Track> track,int& isStoppingReco,
			 SortedObj& sortedObj);

    //void SortShower(art::Ptr<recob::Shower> shower,TVector3 dir,int& isStoppingReco,
    //		    mvapid::MVAAlg::SortedObj& sortedShower);
    void SortShower(art::Ptr<recob::Shower> shower,int& isStoppingReco,mvapid::MVAAlg::SortedObj& sortedShower);

    void RunPCA(std::vector< art::Ptr<recob::Hit> >& hits,std::vector<double>& eVals,std::vector<double>& eVecs);

    void _Var_Shape(const SortedObj& track,
		    double& coreHaloRatio,double& concentration,
		    double& conicalness);
    
    double CalcSegmentdEdxFrac(const SortedObj& track,double start,double end);
    
    double CalcSegmentdEdxDist(const SortedObj& track,double start,double end);

    double CalcSegmentdEdxDistAtEnd(const mvapid::MVAAlg::SortedObj& track,double distAtEnd);
    
    int LinFit(const art::Ptr<recob::Track> track,TVector3& trackPoint,TVector3& trackDir);

    int LinFitShower(const art::Ptr<recob::Shower> shower,TVector3& showerPoint,TVector3& showerDir);
    
    const calo::CalorimetryAlg fCaloAlg;
    
    double fEventT0;

    const art::EDProducer* fParentModule;

    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fHitLabel;
    std::string fSpacePointLabel;
    std::string fTrackingLabel;

    std::vector<art::Ptr<recob::Track> > fTracks;
    std::vector<art::Ptr<recob::Shower> > fShowers;
    std::vector<art::Ptr<recob::SpacePoint> > fSpacePoints;
    std::vector<art::Ptr<recob::Hit> > fHits;
  
    std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<recob::Hit> > > fTracksToHits;
    std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<recob::SpacePoint> > > fTracksToSpacePoints;
    std::map<art::Ptr<recob::Shower>,std::vector<art::Ptr<recob::Hit> > > fShowersToHits;
    std::map<art::Ptr<recob::Shower>,std::vector<art::Ptr<recob::SpacePoint> > > fShowersToSpacePoints;
    std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint> > fHitsToSpacePoints;
    std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit > > fSpacePointsToHits;

    anab::MVAPIDResult fResHolder;

    TMVA::Reader fReader;
    
    std::vector<std::string> fMVAMethods;
    std::vector<std::string> fWeightFiles;

    bool fCheatVertex;

    TLorentzVector fVertex4Vect;

  }; // class MVAAlg

} // namespace mvapid

#endif // ifndef MVAAlg_H
