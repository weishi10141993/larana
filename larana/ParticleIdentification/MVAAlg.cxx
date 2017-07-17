/////////////////////////////////////////////////////////////////
//  \fileMVAAlg.cxx
//  m.haigh@warwick.ac.uk
////////////////////////////////////////////////////////////////////

#include "larana/ParticleIdentification/MVAAlg.h"
#include "larcore/Geometry/Geometry.h"
//#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "TPrincipal.h"
#include "TFile.h"
//#include "TF1.h"
//#include "TF2.h"
#include "TMatrix.h"
#include "TVectorD.h"
#include "larcorealg/CoreUtils/quiet_Math_Functor.h" // remove the wrapper when ROOT header is fixed
#include <Fit/Fitter.h>
//#include "TVirtualFitter.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <cmath>

mvapid::MVAAlg::MVAAlg(fhicl::ParameterSet const& pset, const art::EDProducer* parentModule):
  fCaloAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")), 
  fParentModule(parentModule),fReader(""){
  fHitLabel=pset.get<std::string>("HitLabel");
  fTrackLabel=pset.get<std::string>("TrackLabel");
  fShowerLabel=pset.get<std::string>("ShowerLabel");
  fSpacePointLabel=pset.get<std::string>("SpacePointLabel");
  fTrackingLabel=pset.get<std::string>("TrackingLabel","");

  fCheatVertex=pset.get<bool>("CheatVertex",false);

  fReader.AddVariable("evalRatio",&fResHolder.evalRatio);
  fReader.AddVariable("coreHaloRatio",&fResHolder.coreHaloRatio);
  fReader.AddVariable("concentration",&fResHolder.concentration);
  fReader.AddVariable("conicalness",&fResHolder.conicalness);
  fReader.AddVariable("dEdxStart",&fResHolder.dEdxStart);
  fReader.AddVariable("dEdxEnd",&fResHolder.dEdxEnd);
  fReader.AddVariable("dEdxEndRatio",&fResHolder.dEdxEndRatio);

  fMVAMethods=pset.get<std::vector<std::string> >("MVAMethods");
  std::vector<std::string> weightFileBnames=pset.get<std::vector<std::string> >("WeightFiles");

  cet::search_path searchPath("FW_SEARCH_PATH");
  for(auto fileIter=weightFileBnames.begin();fileIter!=weightFileBnames.end();++fileIter){
    std::string fileWithPath;
    if(!searchPath.find_file(*fileIter,fileWithPath)){
	fWeightFiles.clear();
	fMVAMethods.clear();
	throw cet::exception("MVAPID")<<"Unable to find weight file "<<*fileIter<<" in search path "
				      <<searchPath.to_string()<<std::endl;
      }
      fWeightFiles.push_back(fileWithPath);
  }

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

int mvapid::MVAAlg::IsInActiveVol(const TVector3& pos)
{

  //DUNE 10kt workspace geometry - see https://cdcvs.fnal.gov/redmine/projects/dunetpc/wiki/DUNE_Geometries 
  //Exact coordinates of edges of active volume from email from Tyler Alion on 15-3-2016
  const double activeVolMinX = -363.376;
  const double activeVolMaxX = 363.376;
  const double activeVolMinY = -607.829;
  const double activeVolMaxY = 607.829;
  const double activeVolMinZ = -0.87625;
  const double activeVolMaxZ = 463.904;
  const double fiducialDist = 5.0;

  if(pos.X() > (activeVolMinX + fiducialDist) && pos.X() < (activeVolMaxX - fiducialDist)
     && pos.Y() > (activeVolMinY + fiducialDist) && pos.Y() < (activeVolMaxY - fiducialDist)
     && pos.Z() > (activeVolMinZ + fiducialDist) && pos.Z() < (activeVolMaxZ - fiducialDist))
    return 1;
  else
    return 0;

}

void mvapid::MVAAlg::RunPID(art::Event& evt,std::vector<anab::MVAPIDResult>& result,
			    art::Assns<recob::Track, anab::MVAPIDResult, void>& trackAssns,
			    art::Assns<recob::Shower, anab::MVAPIDResult, void>& showerAssns){

  this->PrepareEvent(evt);

  for(auto trackIter=fTracks.begin();trackIter!=fTracks.end();++trackIter){
    mvapid::MVAAlg::SortedObj sortedObj;

    std::vector<double> eVals,eVecs;
    int isStoppingReco; 
    this->RunPCA(fTracksToHits[*trackIter],eVals,eVecs);
    double evalRatio;
    if(eVals[0] < 0.0001)
      evalRatio=0.0;
    else
      evalRatio=std::sqrt(eVals[1]*eVals[1]+eVals[2]*eVals[2])/eVals[0];
    this->FitAndSortTrack(*trackIter,isStoppingReco,sortedObj);
    double coreHaloRatio,concentration,conicalness;
    this->_Var_Shape(sortedObj,coreHaloRatio,concentration,conicalness);
    double dEdxStart = CalcSegmentdEdxFrac(sortedObj,0.,0.05);
    double dEdxEnd = CalcSegmentdEdxFrac(sortedObj,0.9,1.0);
    double dEdxPenultimate = CalcSegmentdEdxFrac(sortedObj,0.8,0.9);

    /*    
    std::cout<<"coreHaloRatio:   "<<coreHaloRatio<<std::endl;
    std::cout<<"concentration:   "<<concentration<<std::endl;
    std::cout<<"conicalness:     "<<conicalness<<std::endl;
    std::cout<<"dEdxStart: "<<dEdxStart<<std::endl;
    std::cout<<"dEdxEnd: "<<dEdxEnd<<std::endl;
    std::cout<<"dEdxEndRatio: ";
    if(dEdxPenultimate < 0.1)
    std::cout<<"1.0";
    else
    std::cout<<dEdxEnd/dEdxPenultimate;
    std::cout<<std::endl;
    */

    fResHolder.isTrack=1;
    fResHolder.isStoppingReco=isStoppingReco;
    fResHolder.nSpacePoints=sortedObj.hitMap.size();
    fResHolder.trackID=(*trackIter)->ID();
    fResHolder.evalRatio=evalRatio;
    fResHolder.concentration=concentration;
    fResHolder.coreHaloRatio=coreHaloRatio;
    fResHolder.conicalness=conicalness;
    fResHolder.dEdxStart=dEdxStart;
    fResHolder.dEdxEnd=dEdxEnd;
    if(dEdxPenultimate < 0.1)
      fResHolder.dEdxEndRatio=1.0;
    else
      fResHolder.dEdxEndRatio=dEdxEnd/dEdxPenultimate;
    fResHolder.length=sortedObj.length;

    for(auto methodIter=fMVAMethods.begin();methodIter!=fMVAMethods.end();++methodIter){
      fResHolder.mvaOutput[*methodIter]=fReader.EvaluateMVA(*methodIter);
    }
    result.push_back(fResHolder);
    util::CreateAssn(*fParentModule, evt, result, *trackIter, trackAssns);
  }

  for(auto showerIter=fShowers.begin();showerIter!=fShowers.end();++showerIter){
    mvapid::MVAAlg::SortedObj sortedObj;

    std::vector<double> eVals,eVecs;
    int isStoppingReco;

    this->RunPCA(fShowersToHits[*showerIter],eVals,eVecs);

    double evalRatio;
    if(eVals[0]< 0.0001)
      evalRatio=0.0;
    else
      evalRatio=std::sqrt(eVals[1]*eVals[1]+eVals[2]*eVals[2])/eVals[0];

    //this->SortShower(*showerIter,TVector3(eVecs[0],eVecs[3],eVecs[6]),isStoppingReco,sortedObj);
    this->SortShower(*showerIter,isStoppingReco,sortedObj);

    double coreHaloRatio,concentration,conicalness;
    this->_Var_Shape(sortedObj,coreHaloRatio,concentration,conicalness);
    double dEdxStart = CalcSegmentdEdxFrac(sortedObj,0.,0.05);
    double dEdxEnd = CalcSegmentdEdxFrac(sortedObj,0.9,1.0);
    double dEdxPenultimate = CalcSegmentdEdxFrac(sortedObj,0.8,0.9);

    /*    
    std::cout<<"coreHaloRatio:   "<<coreHaloRatio<<std::endl;
    std::cout<<"concentration:   "<<concentration<<std::endl;
    std::cout<<"conicalness:     "<<conicalness<<std::endl;
    std::cout<<"dEdxStart: "<<dEdxStart<<std::endl;
    std::cout<<"dEdxEnd: "<<dEdxEnd<<std::endl;
    std::cout<<"dEdxEndRatio: ";
    if(dEdxPenultimate < 0.1)
    std::cout<<"1.0";
    else
    std::cout<<dEdxEnd/dEdxPenultimate;
    std::cout<<std::endl;
    */

    fResHolder.isTrack=0;
    fResHolder.isStoppingReco=isStoppingReco;
    fResHolder.nSpacePoints=sortedObj.hitMap.size();
    fResHolder.trackID=(*showerIter)->ID()+1000; //For the moment label showers by adding 1000 to ID

    fResHolder.evalRatio=evalRatio;
    fResHolder.concentration=concentration;
    fResHolder.coreHaloRatio=coreHaloRatio;
    fResHolder.conicalness=conicalness;
    fResHolder.dEdxStart=dEdxStart;
    fResHolder.dEdxEnd=dEdxEnd;
    if(dEdxPenultimate < 0.1)
      fResHolder.dEdxEndRatio=1.0;
    else
      fResHolder.dEdxEndRatio=dEdxEnd/dEdxPenultimate;
    fResHolder.length=sortedObj.length;

    for(auto methodIter=fMVAMethods.begin();methodIter!=fMVAMethods.end();++methodIter){
      fResHolder.mvaOutput[*methodIter]=fReader.EvaluateMVA(*methodIter);
    }
    result.push_back(fResHolder);
    util::CreateAssn(*fParentModule, evt, result, *showerIter, showerAssns);
  }

}

void mvapid::MVAAlg::PrepareEvent(const art::Event& evt){

  fHits.clear();
  fSpacePoints.clear();
  fTracks.clear();
  fShowers.clear();
  fSpacePointsToHits.clear();
  fHitsToSpacePoints.clear();
  fTracksToHits.clear();
  fTracksToSpacePoints.clear();
  fShowersToHits.clear();
  fShowersToSpacePoints.clear();

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

  art::Handle< std::vector<recob::Shower> > showersHandle;
  evt.getByLabel(fShowerLabel, showersHandle);

  for (unsigned int iShower = 0; iShower < showersHandle->size(); ++iShower){
    const art::Ptr<recob::Shower> shower(showersHandle, iShower);
    fShowers.push_back(shower);
  }

  art::Handle< std::vector<recob::SpacePoint> > spHandle;
  evt.getByLabel(fSpacePointLabel, spHandle);

  for (unsigned int iSP = 0; iSP < spHandle->size(); ++iSP){
    const art::Ptr<recob::SpacePoint> spacePoint(spHandle, iSP);
    fSpacePoints.push_back(spacePoint);
  }

  art::FindManyP<recob::Hit> findTracksToHits(fTracks,evt,fTrackLabel);
  art::FindManyP<recob::Hit> findShowersToHits(fShowers,evt,fShowerLabel);
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

  for(unsigned int iShower = 0; iShower < fShowers.size(); ++iShower){
    const art::Ptr<recob::Shower> shower=fShowers.at(iShower);
    const std::vector< art::Ptr<recob::Hit> > showerHits = findShowersToHits.at(iShower);

    for (unsigned int iHit=0; iHit<showerHits.size(); ++iHit)
    {
      const art::Ptr<recob::Hit> hit = showerHits.at(iHit);          
      fShowersToHits[shower].push_back(hit);
      if(fHitsToSpacePoints.count(hit)){
	fShowersToSpacePoints[shower].push_back(fHitsToSpacePoints.at(hit));
      }
    }
  }

  if(fCheatVertex){
    art::Handle< std::vector<simb::MCParticle> > partHandle;
    evt.getByLabel(fTrackingLabel, partHandle);

    if(partHandle->size()==0||partHandle->at(0).TrackId()!=1){
      std::cout<<"Error, ID of first track in largeant list is not 0"<<std::endl;
      exit(1);
    } 
    fVertex4Vect=partHandle->at(0).Position();
  }
}

void mvapid::MVAAlg::FitAndSortTrack(art::Ptr<recob::Track> track,int& isStoppingReco,
				     mvapid::MVAAlg::SortedObj& sortedTrack){
  
  sortedTrack.hitMap.clear();
  TVector3 trackPoint,trackDir;
  this->LinFit(track,trackPoint,trackDir);

  TVector3 nearestPointStart, nearestPointEnd;

  //For single-particle events can opt to cheat vertex from start of primary trajectory.
  //Ok since in real events it should be possible to identify the true vertex.
  if(fCheatVertex){
    if((track->End()-fVertex4Vect.Vect()).Mag()>(track->Vertex()-fVertex4Vect.Vect()).Mag()){
      nearestPointStart = trackPoint+trackDir*(trackDir.Dot(track->Vertex()-trackPoint)/trackDir.Mag2());
      nearestPointEnd = trackPoint+trackDir*(trackDir.Dot(track->End()-trackPoint)/trackDir.Mag2());
      isStoppingReco = this->IsInActiveVol(track->End());
    }
    else{
      nearestPointStart = trackPoint+trackDir*(trackDir.Dot(track->End()-trackPoint)/trackDir.Mag2());
      nearestPointEnd = trackPoint+trackDir*(trackDir.Dot(track->Vertex()-trackPoint)/trackDir.Mag2());
      isStoppingReco = this->IsInActiveVol(track->Vertex());
      trackDir*=-1.;
    }
  }
  else{
    if(track->End().Z() >= track->Vertex().Z()){ //Otherwise assume particle is forward-going for now...
      nearestPointStart = trackPoint+trackDir*(trackDir.Dot(track->Vertex()-trackPoint)/trackDir.Mag2());
      nearestPointEnd = trackPoint+trackDir*(trackDir.Dot(track->End()-trackPoint)/trackDir.Mag2());
      isStoppingReco = this->IsInActiveVol(track->End());
    }
    else{
      nearestPointStart = trackPoint+trackDir*(trackDir.Dot(track->End()-trackPoint)/trackDir.Mag2());
      nearestPointEnd = trackPoint+trackDir*(trackDir.Dot(track->Vertex()-trackPoint)/trackDir.Mag2());
      isStoppingReco = this->IsInActiveVol(track->Vertex());
    }

    if(trackDir.Z() <= 0){
      trackDir.SetX(-trackDir.X());
      trackDir.SetY(-trackDir.Y());
      trackDir.SetZ(-trackDir.Z());
    }
  }

  sortedTrack.start=nearestPointStart;
  sortedTrack.end=nearestPointEnd;
  sortedTrack.dir=trackDir;
  sortedTrack.length=(nearestPointEnd-nearestPointStart).Mag();

  std::vector<art::Ptr<recob::Hit>> hits=fTracksToHits[track];

  for(auto hitIter=hits.begin();hitIter!=hits.end();++hitIter){
    
    if(!fHitsToSpacePoints.count(*hitIter)) continue;
    art::Ptr<recob::SpacePoint> sp=fHitsToSpacePoints.at(*hitIter);

    TVector3 nearestPoint = trackPoint+trackDir*(trackDir.Dot(TVector3(sp->XYZ())-trackPoint)/trackDir.Mag2());
    double lengthAlongTrack=(nearestPointStart-nearestPoint).Mag();
    sortedTrack.hitMap.insert(std::pair<double,art::Ptr<recob::Hit> >(lengthAlongTrack,*hitIter));
  }
}

//void mvapid::MVAAlg::SortShower(art::Ptr<recob::Shower> shower,TVector3 dir,int& isStoppingReco,
//				     mvapid::MVAAlg::SortedObj& sortedShower){
void mvapid::MVAAlg::SortShower(art::Ptr<recob::Shower> shower,int& isStoppingReco,
				  mvapid::MVAAlg::SortedObj& sortedShower){
  sortedShower.hitMap.clear();
  
  std::vector<art::Ptr<recob::Hit>> hits=fShowersToHits[shower];

  TVector3 showerEnd(0, 0, 0); 
  double furthestHitFromStart = -999.9;
  for(auto hitIter=hits.begin();hitIter!=hits.end();++hitIter){

    if(!fHitsToSpacePoints.count(*hitIter)) continue;
    art::Ptr<recob::SpacePoint> sp=fHitsToSpacePoints.at(*hitIter);
    if((TVector3(sp->XYZ()) - shower->ShowerStart()).Mag() > furthestHitFromStart)
      {
	showerEnd = TVector3(sp->XYZ());
        furthestHitFromStart = (TVector3(sp->XYZ()) - shower->ShowerStart()).Mag();
      }
  }

  TVector3 showerPoint,showerDir;
  this->LinFitShower(shower,showerPoint,showerDir);

  TVector3 nearestPointStart, nearestPointEnd;

  //Ensure that shower is fitted in correct direction (assuming for now that particle moves in +z direction)

  if(fCheatVertex){
    if((showerEnd-fVertex4Vect.Vect()).Mag()>(shower->ShowerStart()-fVertex4Vect.Vect()).Mag()){
      nearestPointStart = showerPoint+showerDir*(showerDir.Dot(shower->ShowerStart()-showerPoint)/showerDir.Mag2());
      nearestPointEnd = showerPoint+showerDir*(showerDir.Dot(showerEnd-showerPoint)/showerDir.Mag2());
      isStoppingReco = this->IsInActiveVol(showerEnd);
    }
    else
      {
	nearestPointStart = showerPoint+showerDir*(showerDir.Dot(showerEnd-showerPoint)/showerDir.Mag2());
	nearestPointEnd = showerPoint+showerDir*(showerDir.Dot(shower->ShowerStart()-showerPoint)/showerDir.Mag2());
	isStoppingReco = this->IsInActiveVol(shower->ShowerStart());
	showerDir*=-1.;
      }
  }
  else{
    if(showerEnd.Z() >= shower->ShowerStart().Z()){
      nearestPointStart = showerPoint+showerDir*(showerDir.Dot(shower->ShowerStart()-showerPoint)/showerDir.Mag2());
      nearestPointEnd = showerPoint+showerDir*(showerDir.Dot(showerEnd-showerPoint)/showerDir.Mag2());
      isStoppingReco = this->IsInActiveVol(showerEnd);
    }
    else{
      nearestPointStart = showerPoint+showerDir*(showerDir.Dot(showerEnd-showerPoint)/showerDir.Mag2());
      nearestPointEnd = showerPoint+showerDir*(showerDir.Dot(shower->ShowerStart()-showerPoint)/showerDir.Mag2());
      isStoppingReco = this->IsInActiveVol(shower->ShowerStart());
    }
    
    if(showerDir.Z() <= 0){
      showerDir.SetX(-showerDir.X());
      showerDir.SetY(-showerDir.Y());
      showerDir.SetZ(-showerDir.Z());
    }  
  }

  sortedShower.start=nearestPointStart;
  sortedShower.end=nearestPointEnd;
  //sortedShower.dir=dir;
  sortedShower.dir=showerDir;
  sortedShower.length=(nearestPointEnd-nearestPointStart).Mag();

  for(auto hitIter=hits.begin();hitIter!=hits.end();++hitIter){
    
    if(!fHitsToSpacePoints.count(*hitIter)) continue;
    art::Ptr<recob::SpacePoint> sp=fHitsToSpacePoints.at(*hitIter);

    TVector3 nearestPoint = showerPoint+showerDir*(showerDir.Dot(TVector3(sp->XYZ())-showerPoint)/showerDir.Mag2());
    double lengthAlongShower=(nearestPointStart-nearestPoint).Mag();
    sortedShower.hitMap.insert(std::pair<double,art::Ptr<recob::Hit> >(lengthAlongShower,*hitIter));
  }
  
}
void mvapid::MVAAlg::RunPCA(std::vector< art::Ptr<recob::Hit> >& hits,std::vector<double>& eVals,std::vector<double>& eVecs){

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
  for(unsigned int i=0;i<3;++i){
    eVals.push_back(principal->GetEigenValues()->GetMatrixArray()[i]);
  }

  for(unsigned int i=0;i<9;++i){
    eVecs.push_back(principal->GetEigenVectors()->GetMatrixArray()[i]);
  }
}
void mvapid::MVAAlg::_Var_Shape(const mvapid::MVAAlg::SortedObj& track,
				double& coreHaloRatio,double& concentration,
				double& conicalness){
  
  static const unsigned int conMinHits=3;
  static const double minCharge = 0.1;
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

      double distFromTrackFit = ((TVector3(sp->XYZ()) - track.start).Cross(track.dir)).Mag();
 
      ++nHits;

      if(distFromTrackFit < MoliereRadiusFraction * MoliereRadius)
	chargeCore += hitIter->second->Integral();
      else
	chargeHalo += hitIter->second->Integral();

      totalCharge += hitIter->second->Integral();

      chargeCon += hitIter->second->Integral() / std::max(1.E-2,distFromTrackFit);
      if(hitIter->first/track.length<conFracRange){
	chargeConStart+=distFromTrackFit*distFromTrackFit*hitIter->second->Integral();
	++nHitsConStart;
	totalChargeStart+=hitIter->second->Integral();
      }
      else if(1.-hitIter->first/track.length<conFracRange){
	chargeConEnd+=distFromTrackFit*distFromTrackFit*hitIter->second->Integral();
	++nHitsConEnd;
	totalChargeEnd+=hitIter->second->Integral();
      }
    }
  }
  
  coreHaloRatio=chargeHalo/TMath::Max(1.0E-3, chargeCore);
  coreHaloRatio=TMath::Min(100.0, coreHaloRatio);
  concentration=chargeCon/totalCharge;
  if(nHitsConStart>=conMinHits&&nHitsConEnd>=conMinHits&&totalChargeEnd>minCharge&&sqrt(chargeConStart)>minCharge&&totalChargeStart>minCharge){
    conicalness=(sqrt(chargeConEnd)/totalChargeEnd) / (sqrt(chargeConStart)/totalChargeStart);
  }
  else{
    conicalness=1.;
  }
}

double mvapid::MVAAlg::CalcSegmentdEdxFrac(const mvapid::MVAAlg::SortedObj& track,double start,double end){

  double trackLength=(track.end-track.start).Mag();
  return CalcSegmentdEdxDist(track,start*trackLength,end*trackLength);
}

double mvapid::MVAAlg::CalcSegmentdEdxDistAtEnd(const mvapid::MVAAlg::SortedObj& track,double distAtEnd){

  double trackLength=(track.end-track.start).Mag();
  return CalcSegmentdEdxDist(track,trackLength-distAtEnd,trackLength);
}

double mvapid::MVAAlg::CalcSegmentdEdxDist(const mvapid::MVAAlg::SortedObj& track,double start,double end){

  art::ServiceHandle<geo::Geometry> geom;
  
  //Need to get these from geometry really
  static const TVector3 normToWiresU(0.0, 0.7071, -0.7071);
  static const TVector3 normToWiresV(0.0, 0.7071, 0.7071);
  static const TVector3 normToWiresW(0.0, 0.0, 1.0);

  //Calculate scalar product of unit fitted track vector with unit norm to wire planes
  double scalarProdTrackWires[3] = {0};
  scalarProdTrackWires[0] = track.dir.Unit().Dot(normToWiresU);
  scalarProdTrackWires[1] = track.dir.Unit().Dot(normToWiresV);
  scalarProdTrackWires[2] = track.dir.Unit().Dot(normToWiresW);

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

    TVector3 dir=track.dir;

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

    ROOT::Math::Functor fcn(sdist,6);

    //Initial fit parameters from track start and end...
    TVector3 trackStart=track->Vertex();
    TVector3 trackEnd=track->End();
    trackDir=(trackEnd-trackStart).Unit();

    TVector3 x0=trackStart-trackDir;
    TVector3 u=trackDir;

    double pStart[6] = {x0.X(),u.X(),x0.Y(),u.Y(),x0.Z(),u.Z()};

    fitter.SetFCN(fcn,pStart);

    bool ok = fitter.FitFCN();
    if (!ok) {
      trackPoint.SetXYZ(x0.X(),x0.Y(),x0.Z());
      trackDir.SetXYZ(u.X(),u.Y(),u.Z());
      trackDir=trackDir.Unit();
      return 1;
    }
    else{
      const ROOT::Fit::FitResult & result = fitter.Result();
      const double * parFit = result.GetParams();
      trackPoint.SetXYZ(parFit[0],parFit[2],parFit[4]);
      trackDir.SetXYZ(parFit[1],parFit[3],parFit[5]);
      trackDir=trackDir.Unit();
      return 0;
    }
}

int mvapid::MVAAlg::LinFitShower(const art::Ptr<recob::Shower> shower,TVector3& showerPoint,TVector3& showerDir){

  const std::vector<art::Ptr<recob::SpacePoint> >& sp = fShowersToSpacePoints.at(shower);

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

  ROOT::Math::Functor fcn(sdist,6);

  //Initial fit parameters from shower start and end...
  TVector3 showerStart=shower->ShowerStart();
  showerDir = shower->Direction().Unit();

  TVector3 x0=showerStart-showerDir;
  TVector3 u=showerDir;

  double pStart[6] = {x0.X(),u.X(),x0.Y(),u.Y(),x0.Z(),u.Z()};

  fitter.SetFCN(fcn,pStart);

  bool ok = fitter.FitFCN();
  if (!ok) {
    showerPoint.SetXYZ(x0.X(),x0.Y(),x0.Z());
    showerDir.SetXYZ(u.X(),u.Y(),u.Z());
    showerDir=showerDir.Unit();
    return 1;
  }
  else{
    const ROOT::Fit::FitResult & result = fitter.Result();
    const double * parFit = result.GetParams();
    showerPoint.SetXYZ(parFit[0],parFit[2],parFit[4]);
    showerDir.SetXYZ(parFit[1],parFit[3],parFit[5]);
    showerDir=showerDir.Unit();
    return 0;
  }
}

