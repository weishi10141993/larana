/////////////////////////////////////////////////////////////////
//  \fileMVAAlg.cxx
//  m.haigh@warwick.ac.uk
////////////////////////////////////////////////////////////////////

#include "MVAAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "TPrincipal.h"
#include "TFile.h"
//#include "TF1.h"
//#include "TF2.h"
#include "TMatrix.h"
#include "TVectorD.h"
#include <Math/Functor.h>
#include <Fit/Fitter.h>
//#include "TVirtualFitter.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/FindOneP.h"
#include "lardata/Utilities/AssociationUtil.h"

mvapid::MVAAlg::MVAAlg(fhicl::ParameterSet const& pset, const art::EDProducer* parentModule):
  fCaloAlg(pset), fParentModule(parentModule),fReader(""){

  fHitLabel="hit35t";
  fTrackLabel="particlestitcher";
  fSpacePointLabel="pandora";

  fReader.AddVariable("evalRatio",&fResHolder.evalRatio);
  fReader.AddVariable("concentration",&fResHolder.concentration);
  fReader.AddVariable("coreHaloRatio",&fResHolder.coreHaloRatio);
  fReader.AddVariable("conicalness",&fResHolder.conicalness);
  fReader.AddVariable("dEdxStart",&fResHolder.dEdxStart);
  fReader.AddVariable("dEdxEnd",&fResHolder.dEdxEnd);
  fReader.AddVariable("dEdxPenultimate",&fResHolder.dEdxPenultimate);

  fMVAMethods=pset.get<std::vector<std::string> >("MVAMethods");
  fWeightFiles=pset.get<std::vector<std::string> >("WeightFiles");

  if(fMVAMethods.size()!=fWeightFiles.size()){
    std::cerr<<"Mismatch in number of MVA methods and weight files!"<<std::endl;
    exit(1);
  }

  for(unsigned int iMethod=0;iMethod!=fMVAMethods.size();++iMethod){
    fReader.BookMVA(fMVAMethods[iMethod], fWeightFiles[iMethod]);
  }
}

mvapid::MVAAlg::~MVAAlg(){}
    
void mvapid::MVAAlg::reconfigure(fhicl::ParameterSet const& p){}
  
void mvapid::MVAAlg::RunPID(art::Event& evt,std::vector<anab::MVAPIDResult>& result,
			    art::Assns<recob::Track, anab::MVAPIDResult, void>& assns){

  //Need to get these from geometry really
  /*const double activeVolMinX = -35.18;
  const double activeVolMaxX = 222.46;
  const double activeVolMinY = -84.22;
  const double activeVolMaxY = 115.09;
  const double activeVolMinZ = -2.04;
  const double activeVolMaxZ = 156.78;
  */

  this->PrepareEvent(evt);

  for(auto trackIter=fTracks.begin();trackIter!=fTracks.end();++trackIter){

    mvapid::MVAAlg::SortedTrack sortedTrack;
    
    this->FitAndSortTrack(*trackIter,sortedTrack);
    double evalRatio;
    this->_Var_EValRatio(*trackIter,evalRatio);
    double coreHaloRatio,concentration,conicalness;
    this->_Var_Shape(sortedTrack,coreHaloRatio,concentration,conicalness);
    double dEdxStart = CalcSegmentdEdxFrac(sortedTrack,0.,0.2);
    double dEdxEnd = CalcSegmentdEdxFrac(sortedTrack,0.8,1.0);
    double dEdxPenultimate = CalcSegmentdEdxFrac(sortedTrack,0.6,0.8);

    std::cout<<dEdxStart<<" "<<dEdxEnd<<" "<<dEdxPenultimate<<std::endl;

    fResHolder.nSpacePoints=sortedTrack.hitMap.size();
    fResHolder.trackID=(*trackIter)->ID();

    fResHolder.evalRatio=evalRatio;
    fResHolder.concentration=concentration;
    fResHolder.coreHaloRatio=coreHaloRatio;
    fResHolder.conicalness=conicalness;
    fResHolder.dEdxStart=dEdxStart;
    fResHolder.dEdxEnd=dEdxEnd;
    fResHolder.dEdxPenultimate=dEdxPenultimate;

    for(auto methodIter=fMVAMethods.begin();methodIter!=fMVAMethods.end();++methodIter){
      fResHolder.mvaOutput[*methodIter]=fReader.EvaluateMVA(*methodIter);
    }

    result.push_back(fResHolder);
    util::CreateAssn(*fParentModule, evt, result, *trackIter, assns);
  }

  /*  if(recoTrackEnd.X() > (activeVolMinX + fiducialDist) && recoTrackEnd.X() < (activeVolMaxX - fiducialDist)
     && recoTrackEnd.Y() > (activeVolMinY + fiducialDist) && recoTrackEnd.Y() < (activeVolMaxY - fiducialDist)
     && recoTrackEnd.Z() > (activeVolMinZ + fiducialDist) && recoTrackEnd.Z() < (activeVolMaxZ - fiducialDist))
    outputPtr->IsStoppingReco = true;
  else
  outputPtr->IsStoppingReco = false;*/
}

void mvapid::MVAAlg::PrepareEvent(const art::Event& evt){

  fHits.clear();
  fSpacePoints.clear();
  fTracks.clear();
  fSpacePointsToHits.clear();
  fHitsToSpacePoints.clear();
  fTracksToHits.clear();
  fTracksToSpacePoints.clear();

  auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fEventT0=detProp->TriggerOffset();

  art::Handle< std::vector<recob::Hit> > hitsHandle;
  evt.getByLabel(fHitLabel, hitsHandle);

  for (unsigned int iHit = 0; iHit < hitsHandle->size(); ++iHit){
    const art::Ptr<recob::Hit> hit(hitsHandle, iHit);
    fHits.push_back(hit);
  }

  art::Handle< std::vector<recob::Track> > tracksHandle;
  evt.getByLabel(fTrackLabel, tracksHandle);

  for (unsigned int iTrack = 0; iTrack < tracksHandle->size(); ++iTrack){
    const art::Ptr<recob::Track> track(tracksHandle, iTrack);
    fTracks.push_back(track);
  }

  art::Handle< std::vector<recob::SpacePoint> > spHandle;
  evt.getByLabel(fSpacePointLabel, spHandle);

  for (unsigned int iSP = 0; iSP < spHandle->size(); ++iSP){
    const art::Ptr<recob::SpacePoint> spacePoint(spHandle, iSP);
    fSpacePoints.push_back(spacePoint);
  }

  art::FindManyP<recob::Hit> findTracksToHits(fTracks,evt,fTrackLabel);
  art::FindOneP<recob::Hit> findSPToHits(fSpacePoints,evt,fSpacePointLabel);
  
  for(unsigned int iSP = 0; iSP < fSpacePoints.size(); ++iSP)
    {
      const art::Ptr<recob::SpacePoint> spacePoint=fSpacePoints.at(iSP);

      const art::Ptr<recob::Hit> hit = findSPToHits.at(iSP);
      fSpacePointsToHits[spacePoint]=hit;
      fHitsToSpacePoints[hit]=spacePoint;
  }	

  
  for(unsigned int iTrack = 0; iTrack < fTracks.size(); ++iTrack){
    const art::Ptr<recob::Track> track=fTracks.at(iTrack);
    
    const std::vector< art::Ptr<recob::Hit> > trackHits = findTracksToHits.at(iTrack);
    for (unsigned int iHit=0; iHit<trackHits.size(); ++iHit)
    {
      const art::Ptr<recob::Hit> hit = trackHits.at(iHit);          
      fTracksToHits[track].push_back(hit);
      if(fHitsToSpacePoints.count(hit)){
	fTracksToSpacePoints[track].push_back(fHitsToSpacePoints.at(hit));
      }
    }
  }
}

void mvapid::MVAAlg::FitAndSortTrack(art::Ptr<recob::Track> track,
				     mvapid::MVAAlg::SortedTrack& sortedTrack){
  
  sortedTrack.hitMap.clear();
  TVector3 trackPoint,trackDir;
  this->LinFit(track,trackPoint,trackDir);

  TVector3 nearestPointStart = trackPoint+trackDir*(trackDir.Dot(track->Vertex()-trackPoint)/trackDir.Mag2());
  TVector3 nearestPointEnd = trackPoint+trackDir*(trackDir.Dot(track->End()-trackPoint)/trackDir.Mag2());
  
  sortedTrack.trackStart=nearestPointStart;
  sortedTrack.trackEnd=nearestPointEnd;
  sortedTrack.trackDir=trackDir;
  sortedTrack.trackLength=(nearestPointEnd-nearestPointStart).Mag();

  std::vector<art::Ptr<recob::Hit>> hits=fTracksToHits[track];

  for(auto hitIter=hits.begin();hitIter!=hits.end();++hitIter){
    
    if(!fHitsToSpacePoints.count(*hitIter)) continue;
    art::Ptr<recob::SpacePoint> sp=fHitsToSpacePoints.at(*hitIter);

    TVector3 nearestPoint = trackPoint+trackDir*(trackDir.Dot(sp->XYZ()-trackPoint)/trackDir.Mag2());
    double lengthAlongTrack=(nearestPointStart-nearestPoint).Mag();
    sortedTrack.hitMap.insert(std::pair<double,art::Ptr<recob::Hit> >(lengthAlongTrack,*hitIter));
  }
}

void mvapid::MVAAlg::_Var_EValRatio(const art::Ptr<recob::Track>& track,double& eValRatio){

  std::vector< art::Ptr<recob::Hit> > hits = fTracksToHits.at(track);

  // Define the TPrincipal
  TPrincipal* principal = new TPrincipal(3,"D");
  // Define variables to hold the eigenvalues and eigenvectors
  //const TMatrixD* covar = new TMatrixD();
  //const TVectorD* meanval = new TVectorD();

  for(auto hitIter=hits.begin();hitIter!=hits.end();++hitIter){

    if(fHitsToSpacePoints.count(*hitIter)){
      principal->AddRow( fHitsToSpacePoints.at(*hitIter)->XYZ() );
    }
  }
  
  // PERFORM PCA
  principal->MakePrincipals();
  // GET EIGENVALUES AND EIGENVECTORS
  const double* evals;
  evals=principal->GetEigenValues()->GetMatrixArray();

  eValRatio=sqrt(evals[1] * evals[1] + evals[2] * evals[2]) / evals[0];

  std::cout<<"eValRatio:       "<<eValRatio<<std::endl;
}

void mvapid::MVAAlg::_Var_Shape(const mvapid::MVAAlg::SortedTrack& track,
				double& coreHaloRatio,double& concentration,
				double& conicalness){
  
  static const unsigned int conMinHits=10;
  static const double conFracRange=0.2;
  static const double MoliereRadius = 10.1;
  static const double MoliereRadiusFraction = 0.2;


  double totalCharge=0;
  double totalChargeStart=0;
  double totalChargeEnd=0;

  double chargeCore=0;
  double chargeHalo=0;
  double chargeCon=0;
  unsigned int nHits;

  //stuff for conicalness
  double chargeConStart=0;
  double chargeConEnd=0;
  unsigned int nHitsConStart=0;
  unsigned int nHitsConEnd=0;
  

  for(auto hitIter=track.hitMap.begin();hitIter!=track.hitMap.end();++hitIter){
    if(fHitsToSpacePoints.count(hitIter->second)) {
      art::Ptr<recob::SpacePoint> sp=fHitsToSpacePoints.at(hitIter->second);

      double distFromTrackFit = ((sp->XYZ() - track.trackStart).Cross(track.trackDir)).Mag();
 
      ++nHits;

      if(distFromTrackFit < MoliereRadiusFraction * MoliereRadius)
	chargeCore += hitIter->second->Integral();
      else
	chargeHalo += hitIter->second->Integral();

      totalCharge += hitIter->second->Integral();

      chargeCon += hitIter->second->Integral() / std::max(1.E-2,distFromTrackFit);
      if(hitIter->first/track.trackLength<conFracRange){
	chargeConStart+=hitIter->second->Integral() / std::max(1.E-2,distFromTrackFit);;
	++nHitsConStart;
	totalChargeStart+=hitIter->second->Integral();
      }
      else if(1.-hitIter->first/track.trackLength<conFracRange){
	chargeConEnd+=hitIter->second->Integral();
	++nHitsConEnd;
	totalChargeEnd+=hitIter->second->Integral();
      }
    }
  }
  
  coreHaloRatio=chargeHalo/TMath::Max(1E-6, chargeCore);
  concentration=chargeCon/totalCharge;
  if(nHitsConStart>=conMinHits&&nHitsConEnd>=conMinHits){
    conicalness=chargeConStart/chargeConEnd*totalChargeEnd/totalChargeStart;
  }
  else{
    conicalness=0.;
  }
  std::cout<<"coreHaloRatio:   "<<coreHaloRatio<<std::endl;
  std::cout<<"concentration:   "<<concentration<<std::endl;
  std::cout<<"conicalness:     "<<conicalness<<std::endl;
}


double mvapid::MVAAlg::CalcSegmentdEdxFrac(const mvapid::MVAAlg::SortedTrack& track,double start,double end){

  double trackLength=(track.trackEnd-track.trackStart).Mag();
  return CalcSegmentdEdxDist(track,start*trackLength,end*trackLength);
}

double mvapid::MVAAlg::CalcSegmentdEdxDist(const mvapid::MVAAlg::SortedTrack& track,double start,double end){

  art::ServiceHandle<geo::Geometry> geom;
  
  //Need to get these from geometry really
  static const TVector3 normToWiresU(0.0, 0.7071, -0.7071);
  static const TVector3 normToWiresV(0.0, 0.7071, 0.7071);
  static const TVector3 normToWiresW(0.0, 0.0, 1.0);

  //Calculate scalar product of unit fitted track vector with unit norm to wire planes
  double scalarProdTrackWires[3] = {0};
  scalarProdTrackWires[0] = track.trackDir.Unit().Dot(normToWiresU);
  scalarProdTrackWires[1] = track.trackDir.Unit().Dot(normToWiresV);
  scalarProdTrackWires[2] = track.trackDir.Unit().Dot(normToWiresW);

  double totaldEdx=0;
  unsigned int nHits=0;

  //Loop over hits again to calculate average dE/dx and shape variables
  for ( auto hitIter = track.hitMap.begin(); hitIter!=track.hitMap.end(); ++hitIter ){

    if(hitIter->first<start) continue;
    if(hitIter->first>=end) break;

    art::Ptr<recob::Hit> hit=hitIter->second;

    //Pitch to use in dEdx calculation
    double yzPitch, xComponent;
    double pitch3D;

    TVector3 dir=track.trackDir;

    yzPitch = geom->WirePitch(0,1,hit->WireID().Plane, hit->WireID().TPC) / fabs(scalarProdTrackWires[hit->WireID().Plane]);	
    xComponent = yzPitch * dir[0] / sqrt(dir[1] * dir[1] + dir[2] * dir[2]);
    pitch3D = sqrt(xComponent * xComponent + yzPitch * yzPitch);
    
    double dEdx=fCaloAlg.dEdx_AREA(*hit, pitch3D, fEventT0);
    if( dEdx < 50.){
      ++nHits;
      totaldEdx += dEdx;
    }
  }
  return nHits?totaldEdx/nHits:0;
}

int mvapid::MVAAlg::LinFit(const art::Ptr<recob::Track> track,TVector3& trackPoint,TVector3& trackDir){

  const std::vector<art::Ptr<recob::SpacePoint> >& sp = fTracksToSpacePoints.at(track);

    TGraph2D grFit(1);
    unsigned int iPt=0;
    for(auto spIter=sp.begin();spIter!=sp.end();++spIter){
      TVector3 point=(*spIter)->XYZ();
      grFit.SetPoint(iPt++,point.X(),point.Y(),point.Z());
    }

    //Lift from the ROOT line3Dfit.C tutorial
    ROOT::Fit::Fitter fitter;
    // make the functor object
    mvapid::MVAAlg::SumDistance2 sdist(&grFit);
    ROOT::Math::Functor fcn(sdist,4);

    //Initial fit parameters from track start and end...
    TVector3 trackStart=track->Vertex();
    TVector3 trackEnd=track->End();
    trackDir=(trackEnd-trackStart).Unit();

    TVector3 x0=trackStart-trackDir*(trackStart.Z()/trackDir.Z());
    TVector3 u=trackDir*(1./trackDir.Z());

    double pStart[4] = {x0.X(),u.X(),x0.Y(),u.Y()};
    fitter.SetFCN(fcn,pStart);

    bool ok = fitter.FitFCN();
    if (!ok) {
      trackPoint.SetXYZ(x0.X(),x0.Y(),0.);
      trackDir.SetXYZ(u.X(),u.Y(),1.);
      trackDir=trackDir.Unit();
      return 1;
    }
    else{
      const ROOT::Fit::FitResult & result = fitter.Result();
      const double * parFit = result.GetParams();
      trackPoint.SetXYZ(parFit[0],parFit[2],0.);
      trackDir.SetXYZ(parFit[1],parFit[3],1.);
      trackDir=trackDir.Unit();
      return 0;
    }
}
