//////////////////////////////////////////////////
//
// Calorimetry class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// ART port echurch@fnal.gov
//  This algorithm is designed to perform the calorimetric reconstruction 
//  of the 3D reconstructed tracks
////////////////////////////////////////////////////////////////////////
#ifndef CALO_H
#define CALO_H


extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "AnalysisAlg/CalorimetryAlg.h"
#include "Utilities/LArProperties.h"
#include "SimpleTypesAndConstants/PhysicalConstants.h"
#include "Utilities/DetectorProperties.h"

#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Track.h"
#include "AnalysisBase/Calorimetry.h"
#include "Utilities/AssociationUtil.h"
#include "Filters/ChannelFilter.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


///calorimetry
namespace calo {
   
  class Calorimetry : public art::EDProducer {
    
  public:
    
    explicit Calorimetry(fhicl::ParameterSet const& pset); 
    virtual ~Calorimetry();
    
    void beginJob(); 
    //    void endJob();

    void produce(art::Event& evt);

  private:
        
    double LifetimeCorrection(float time);
    double BirksCorrection(double dQdx_e);
    void   ReadCaloTree();

    bool BeginsOnBoundary(art::Ptr<recob::Track> lar_track);
    bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);

    //void GetPitch(art::Ptr<recob::Hit> hit, art::FindManyP<recob::SpacePoint> fmspts, art::Event& evt, double *xyz3d, double &pitch);
    void GetPitch(art::Ptr<recob::Hit> hit, std::vector<double> trkx, std::vector<double> trky, std::vector<double> trkz, std::vector<double> trkw, std::vector<double> trkx0, double *xyz3d, double &pitch);

    std::string fTrackModuleLabel;
    std::string fSpacePointModuleLabel;
    bool fMakeTree;
    bool fUseArea;
    CalorimetryAlg caloAlg;
	
	
    TTree *ftree;
    TH1F *fdQdx_Coll;
    TH1F *fdEdx_Coll;
    TH2F *fbirk;
    TH2F *fdEdx_Coll_vsXZangle;
    TH2F *fdEdx_vs_ResRange;
    TH2F *fdEdx_vs_ResRange_ent;
    TH2F *fdEdx_vs_ResRange_cont;
    TH2F *fdEdx_vs_ResRange_pass;
    TH2F *fdEdx_vs_ResRange_esc;
    TH2F *fdEdx_vs_ResRange_thr;
    TH1F *fKinetic_En;
    TH2F *fKinetic_En_vs_Range;
 
    int frun;           //Run 
    int fevent;         //Event
    int fntrack;        //track number
    int ftotTracks;        //track number
    double fTrkPitch;
    double fTrkPitchI;  //
    double fTrkPitchC;
    double fXStart;
    double fYStart;
    double fZStart;
    double fXEnd;
    double fYEnd;
    double fZEnd;

    int fnhits3D; 
    std::vector<double> fXHit;
    std::vector<double> fYHit;
    std::vector<double> fZHit;
    std::vector<double> fMIPs3D;   
 
    int fnhitsIND; 
    int fnspsIND;
    std::vector<int>    fwireIND;
    std::vector<double> ftimeIND;
    std::vector<double> fstimeIND;
    std::vector<double> fetimeIND;
    std::vector<double> fMIPsIND;
    std::vector<double> fdEdxIND;
    std::vector<double> fResRngIND;
	
    int fnhitsCOL;
    int fnspsCOL;
    std::vector<int>    fwireCOL;
    std::vector<double> ftimeCOL;
    std::vector<double> fstimeCOL;
    std::vector<double> fetimeCOL;
    std::vector<double> fMIPsCOL;
    std::vector<double> fdEdxCOL;
    std::vector<double> fResRngCOL;

    int fnhits;
    int fnsps;
    std::vector<int>    fwire;
    std::vector<double> ftime;
    std::vector<double> fstime;
    std::vector<double> fetime;
    std::vector<double> fMIPs;
    std::vector<double> fdQdx;
    std::vector<double> fdEdx;
    std::vector<double> fResRng;
    std::vector<double> fpitch;

    double TPCsize[3];

  protected: 
    
  
  }; // class Calorimetry

}

#endif // CALO_H

//-------------------------------------------------
calo::Calorimetry::Calorimetry(fhicl::ParameterSet const& pset)
  : fTrackModuleLabel(pset.get< std::string >("TrackModuleLabel")      ),
    fSpacePointModuleLabel (pset.get< std::string >("SpacePointModuleLabel")       ),
    fMakeTree(pset.get< bool >("MakeTree") ),
    fUseArea(pset.get< bool >("UseArea") ),
    caloAlg(pset.get< fhicl::ParameterSet >("CaloAlg"))
{
  produces< std::vector<anab::Calorimetry>              >();
  produces< art::Assns<recob::Track, anab::Calorimetry> >();
}

//-------------------------------------------------
calo::Calorimetry::~Calorimetry()
{
  
}

//-------------------------------------------------
void calo::Calorimetry::beginJob()
{
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;


  fdQdx_Coll = tfs->make<TH1F>("dQdx_Coll","dQ/dx Coll (ADC/cm)",200,0.,200.);
  fdEdx_Coll = tfs->make<TH1F>("dEdx_Coll","dEdx Coll (MeV/cm)",1000,0.,10.); 
  fbirk = tfs->make<TH2F>("Birk_Coll_correction","Birk Coll MeV/cm-v-e/cm",100,0.0,5.0e5,20,0.0,50.0);
  fdEdx_Coll_vsXZangle = tfs->make<TH2F>("dEdx_Coll_vsXZ Angle","dEdx Coll vs. XZ Angle",
					 1000, 0.0, 10.0, 200, -1.0*TMath::Pi(), TMath::Pi());
  fdEdx_vs_ResRange = tfs->make<TH2F>("dEdx_vs_ResRange","dEdx vs. Residual Range ",
				      500, 0.0, 100.0, 200, 0.0, 40.0);
  // dE/dx vs residual range plots for entering, contained, passing and 
  // escaping tracks (even if only for entering and contained tracks
  // the stopping point is know and the actual residual range can be calculated!) 
  fdEdx_vs_ResRange_ent  = tfs->make<TH2F>("dEdx_vs_ResRange_ent", "dEdx vs. Residual Range entering tracks",
					   500, 0.0, 100.0, 200, 0.0, 40.0);
  fdEdx_vs_ResRange_cont = tfs->make<TH2F>("dEdx_vs_ResRange_cont", "dEdx vs. Residual  Range contained tracks",
					   500, 0.0, 100.0, 200, 0.0, 40.0);
  fdEdx_vs_ResRange_pass = tfs->make<TH2F>("dEdx_vs_ResRange_pass", "dEdx vs. Residual Range passing tracks",
					   500, 0.0, 100.0, 200, 0.0, 40.0);
  fdEdx_vs_ResRange_esc  = tfs->make<TH2F>("dEdx_vs_ResRange_esc", "dEdx vs. Residual Range escaping tracks",
					   500, 0.0, 100.0, 200, 0.0, 40.0);
  fKinetic_En = tfs->make<TH1F>("Kinetic_En","Kinetic Energy Deposited in LAr (MeV)",500,0.,500.);
  fKinetic_En_vs_Range = tfs->make<TH2F>("Kinetic_En_vs_Range","Kinetic Energy Deposited in LAr (MeV)vs tot. range",
					 500, 0.0, 100.0, 200, 0., 500.);

  /// ROOT tree for Calorimetry 
  if(fMakeTree) {
    ftree = tfs->make<TTree>("CaloTree","CaloTree");
    ftree->Branch("run",       &frun,       "run/I");
    ftree->Branch("event",     &fevent,     "event/I");
    ftree->Branch("itrack",    &fntrack,    "itrack/I");
    ftree->Branch("ntracks",   &ftotTracks, "ntracks/I");
    ftree->Branch("TrkPitchI", &fTrkPitchI, "TrkPitchI/F");
    ftree->Branch("TrkPitchC", &fTrkPitchC, "TrkPitchC/F");
    ftree->Branch("XStart",    &fXStart,    "XStart/F");
    ftree->Branch("YStart",    &fYStart,    "YStart/F");
    ftree->Branch("ZStart",    &fZStart,    "ZStart/F");  
    ftree->Branch("XEnd",      &fXEnd,      "XEnd/F");  
    ftree->Branch("YEnd",      &fYEnd,      "YEnd/F");  
    ftree->Branch("ZEnd",      &fZEnd,      "ZEnd/F");  
    ftree->Branch("nhits3D",   &fnhits3D,   "nhits3D/I"); 
    ftree->Branch("XHit",      &fXHit);		     
    ftree->Branch("YHit",      &fYHit);		     
    ftree->Branch("ZHit",      &fZHit);                 
    ftree->Branch("nhitsIND",  &fnhitsIND, "nhitsIND/I");	 
    ftree->Branch("nspsIND",   &fnspsIND,  "nspsIND/I");	 
    ftree->Branch("wireIND",   &fwireIND);		 
    ftree->Branch("timeIND",   &ftimeIND);		 
    ftree->Branch("stimeIND",  &fstimeIND);		 
    ftree->Branch("etimeIND",  &fetimeIND);		 
    ftree->Branch("MIPsIND",   &fMIPsIND);  
    ftree->Branch("dEdxIND",   &fdEdxIND);
    ftree->Branch("ResRngIND", &fResRngIND);
    ftree->Branch("nhitsCOL",  &fnhitsCOL, "nhitsCOL/I");
    ftree->Branch("nspsCOL",   &fnspsCOL,  "nspsCOL/I");
    ftree->Branch("wireCOL",   &fwireCOL);
    ftree->Branch("timeCOL",   &ftimeCOL);
    ftree->Branch("stimeCOL",  &fstimeCOL);
    ftree->Branch("etimeCOL",  &fetimeCOL);
    ftree->Branch("MIPsCOL",   &fMIPsCOL);
    ftree->Branch("dEdxCOL",   &fdEdxCOL);
    ftree->Branch("ResRngCOL", &fResRngCOL);
  }
   
  return;
}

//------------------------------------------------------------------------------------//
void calo::Calorimetry::produce(art::Event& evt)
{ 
  art::ServiceHandle<util::LArProperties> LArProp;
  art::ServiceHandle<util::DetectorProperties> detprop;

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
  
  //  art::Handle< std::vector<recob::Hit> > hitListHandle;
  //  std::vector<art::Ptr<recob::Hit> > hitlist;
  //  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
  //    art::fill_ptr_vector(hitlist, hitListHandle);

  // Electronic calibration factor to convert from ADC to electrons
  //double fElectronsToADC = detprop->ElectronsToADC();

  // Get Geometry
  art::ServiceHandle<geo::Geometry> geom;

  //TPC dimensions  
  TPCsize[0] = (geom->DetHalfWidth())*2.;
  TPCsize[1] = (geom->DetHalfHeight())*2.;
  TPCsize[2] = (geom->DetLength());

  //look for dead wires
  filter::ChannelFilter chanFilt;

  frun = evt.id().run();
  fevent = evt.id().event();

  fntrack = 0; 
  fnhits3D = 0; fnhitsIND = 0; fnhitsCOL = 0;

  size_t nplanes = geom->Nplanes();


  //create anab::Calorimetry objects and make association with recob::Track
  std::unique_ptr< std::vector<anab::Calorimetry> > calorimetrycol(new std::vector<anab::Calorimetry>);
  std::unique_ptr< art::Assns<recob::Track, anab::Calorimetry> > assn(new art::Assns<recob::Track, anab::Calorimetry>);

  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit>        fmht(trackListHandle, evt, fTrackModuleLabel);

  for(size_t trkIter = 0; trkIter < tracklist.size(); ++trkIter){   


    // dE/dx vs Residual Range plot
    fnhits3D = fmsp.at(trkIter).size();
    std::vector<double> larStart;
    std::vector<double> larEnd;
    //put xyz coordinates at begin/end of track into vectors(?)
    tracklist[trkIter]->Extent(larStart,larEnd);
    //store track directional cosines
    double trackCosStart[3]={0.,0.,0.};
    double trackCosEnd[3]={0.,0.,0.};
    tracklist[trkIter]->Direction(trackCosStart,trackCosEnd);
        
    bool startsonboundary = BeginsOnBoundary(tracklist[trkIter]);
    bool endsonboundary   = EndsOnBoundary(tracklist[trkIter]);

    int containment = 0;
    if(!startsonboundary && !endsonboundary ) containment = 0;  // contained track
    if(!startsonboundary &&  endsonboundary ) containment = 1;  // escaping track
    if( startsonboundary && !endsonboundary ) containment = 2;  // entering track
    if( startsonboundary &&  endsonboundary ) containment = 3;  // passing track 

    bool TrackStops = (containment == 0 || containment == 2);
    //      fMIPs3D = new double[fnhits3D];

    //recover the Induction (geo::kU) hits list
    
    //int npI = 0;
    
    // Some variables for the hit
    float time;          //hit time at maximum
    float stime;         //hit start time 
    float etime;         //hit end time 
    uint32_t     channel = 0;//channel number
    unsigned int cstat   = 0;    //hit cryostat number 
    unsigned int tpc     = 0;    //hit tpc number 
    unsigned int wire    = 0;   //hit wire number 
    unsigned int plane   = 0;  //hit plane number


//    try{
//      fTrkPitchI = tracklist[trkIter]->PitchInView(geo::kU);
//    }
//    catch( cet::exception &e){
//      mf::LogWarning("Calorimetry") << "caught exception " << e << "\n setting pitch (I) to 0.";
//      fTrkPitchI = 0.;
//    }

    //get hits in each plane

    std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(trkIter);
    std::vector< std::vector< art::Ptr<recob::Hit> > > hits(nplanes);
    for (size_t ah = 0; ah< allHits.size(); ++ah){
      hits[allHits[ah]->WireID().Plane].push_back(allHits[ah]);
    }

    for (size_t ipl = 0; ipl < nplanes; ++ipl){//loop over all wire planes

      fwire.clear();
      ftime.clear();
      fstime.clear();
      fetime.clear();
      fMIPs.clear();
      fdQdx.clear();
      fdEdx.clear();
      fpitch.clear();
      fResRng.clear();

      double Kin_En = 0.;
      double Trk_Length = 0.;
      std::vector<double> vdEdx;
      std::vector<double> vresRange;
      std::vector<double> vdQdx;
      std::vector<double> deadwire; //residual range for dead wires
    
      //std::vector<int>    vhitindexC;
      //Get spacepoints associated with the hits
      art::FindManyP<recob::SpacePoint> fmspts(hits[ipl], evt, fSpacePointModuleLabel);
      //range of wire signals
      unsigned int wire0 = 100000;
      unsigned int wire1 = 0;
      double PIDA = 0;
      int nPIDA = 0;

      // determine track direction. Fill residual range array
      bool GoingDS = true;
      // find the track direction by comparing US and DS charge BB
      double USChg = 0;
      double DSChg = 0;
      // temp array holding distance betweeen space points
      std::vector<double> spdelta;
      //int nht = 0; //number of hits
      fnsps = 0; //number of space points
      std::vector<double> ChargeBeg;
      std::stack<double> ChargeEnd;     

      // find the separation between all space points
      double xx = 0.,yy = 0.,zz = 0.;

      try{
	fTrkPitch = tracklist[trkIter]->PitchInView(geom->Plane(ipl).View());
      }
      catch( cet::exception &e){
	mf::LogWarning("Calorimetry") << "caught exception " 
				      << e << "\n setting pitch (C) to "
				      << util::kBogusD;
	fTrkPitch = 0;
      }

      //save track 3d points
      std::vector<double> trkx;
      std::vector<double> trky;
      std::vector<double> trkz;
      std::vector<double> trkw;
      std::vector<double> trkx0;
      for (size_t i = 0; i<fmspts.size(); ++i){
	
	std::vector< art::Ptr<recob::SpacePoint> > sptv = fmspts.at(i);
	
	//Get hits associated with spacepts
	art::FindManyP<recob::Hit> fmhits(sptv, evt, fSpacePointModuleLabel);
	
	for (size_t j = 0; j < sptv.size(); ++j){
	  
	  std::vector< art::Ptr<recob::Hit> > hitv = fmhits.at(j);
	  
	  for (size_t k = 0; k < hitv.size(); ++k){
	    
	    if (hitv[k]->WireID().Plane != ipl) continue;

	    double t = hitv[k]->PeakTime();
	    double x = detprop->ConvertTicksToX(t, hitv[k]->WireID().Plane, hitv[k]->WireID().TPC, hitv[k]->WireID().Cryostat);
	    double w = hitv[k]->WireID().Wire;
	    trkx.push_back(sptv[j]->XYZ()[0]);
	    trky.push_back(sptv[j]->XYZ()[1]);
	    trkz.push_back(sptv[j]->XYZ()[2]);
	    trkw.push_back(w);
	    trkx0.push_back(x);
	  }
	}
      }
      for (size_t ihit = 0; ihit < hits[ipl].size(); ++ihit){//loop over all hits on each wire plane

	//std::cout<<ihit<<std::endl;

	wire = hits[ipl][ihit]->WireID().Wire;
	time = hits[ipl][ihit]->PeakTime() ;
	stime = hits[ipl][ihit]->StartTime() ;
	etime = hits[ipl][ihit]->EndTime();            
	
	double charge = hits[ipl][ihit]->Charge(true);
	if (fUseArea) charge = hits[ipl][ihit]->Charge(false);
	//get 3d coordinate and track pitch for the current hit
	//not all hits are associated with space points, the method uses neighboring spacepts to interpolate
	double xyz3d[3];
	double pitch;
	GetPitch(hits[ipl][ihit], trkx, trky, trkz, trkw, trkx0, xyz3d, pitch);
	//GetPitch(hits[ipl][ihit], fmspts, evt, xyz3d, pitch);

	if (xyz3d[2]<-100) continue; //hit not on track
	if (pitch<=0) pitch = fTrkPitch;
	if (!pitch) continue;

        if(fnsps == 0) {
          xx = xyz3d[0];
          yy = xyz3d[1];
          zz = xyz3d[2];
          spdelta.push_back(0);
        } else {
          double dx = xyz3d[0] - xx;
          double dy = xyz3d[1] - yy;
          double dz = xyz3d[2] - zz;
          spdelta.push_back(sqrt(dx*dx + dy*dy + dz*dz));
          Trk_Length += spdelta.back();
          xx = xyz3d[0];
          yy = xyz3d[1];
          zz = xyz3d[2];
        }
	
	ChargeBeg.push_back(charge);
	ChargeEnd.push(charge);

	double MIPs = charge;
	double dQdx = MIPs/pitch;
	// \todo fElectronsToADC should have the value saved in CalorimetryAlg
	//double dQdx_e = dQdx/fElectronsToADC;  // Conversion from ADC/cm to e/cm
	double dEdx = 0;
	if (fUseArea) dEdx = caloAlg.dEdx_AREA(hits[ipl][ihit], pitch);
	else dEdx = caloAlg.dEdx_AMP(hits[ipl][ihit], pitch);

	Kin_En = Kin_En + dEdx * pitch;	

	//std::cout<<dEdx<<" "<<pitch<<" "<<Kin_En<<std::endl;

	if (hits[ipl][ihit]->WireID().Wire < wire0) wire0 = hits[ipl][ihit]->WireID().Wire;
	if (hits[ipl][ihit]->WireID().Wire > wire1) wire1 = hits[ipl][ihit]->WireID().Wire;

	fMIPs.push_back(MIPs);
	fdEdx.push_back(dEdx);
	fdQdx.push_back(dQdx);
	fwire.push_back(wire);
	ftime.push_back(time);
	fstime.push_back(stime);
	fetime.push_back(etime);
	fpitch.push_back(pitch);
	++fnsps;
      }
      if (!fnsps){
	calorimetrycol->push_back(anab::Calorimetry(util::kBogusD,
						    vdEdx,
						    vdQdx,
						    vresRange,
						    deadwire,
						    util::kBogusD,
						    fpitch));
	util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);
	continue;
      }
      for (int isp = 0; isp<fnsps; ++isp){
	if (isp>3) break;
	USChg += ChargeBeg[isp];
      }
      int countsp = 0;
      while (!ChargeEnd.empty()){
	if (countsp>3) break;
	DSChg += ChargeEnd.top();
	ChargeEnd.pop();
	++countsp;
      }
      // Going DS if charge is higher at the end
      GoingDS = (DSChg > USChg);
      // determine the starting residual range and fill the array
      fResRng.resize(fnsps);
      if(GoingDS) {
        fResRng[fnsps - 1] = spdelta[fnsps - 1] / 2;
        for(int isp = fnsps - 2; isp > -1; isp--) {
          fResRng[isp] = fResRng[isp+1] + spdelta[isp+1];
        }
      } else {
        fResRng[0] = spdelta[1] / 2;
        for(int isp = 1; isp < fnsps; isp++) {
          fResRng[isp] = fResRng[isp-1] + spdelta[isp];
        }
      }
      mf::LogVerbatim CaloPrtTrk("CaloPrtTrk");
    
      CaloPrtTrk << "Calorimetry Run/Evt: "
		 << frun   << " / " 
		 << fevent << " Track #"
		 << trkIter << " Plane #"
		 << ipl;
      switch (containment) {
      case 0:
	CaloPrtTrk << ", contained.";
	break;
      case 1:
	CaloPrtTrk << ", escaping.";
	break;
      case 2:
	CaloPrtTrk << ", entering.";
	break;
      case 3:
	CaloPrtTrk << ", passing.";
	break;
      default:
	CaloPrtTrk << ", ??";
      }
      if(TrackStops) {
	if(GoingDS) {
	  CaloPrtTrk <<" Going downstream";
	} 
	else {
	  CaloPrtTrk <<" Going upstream";
	}
      }
      CaloPrtTrk<<"\n";
      mf::LogVerbatim CaloPrtHit("CaloPrtHit");
      CaloPrtHit << " pt wire  time  ResRng    MIPs   pitch   dE/dx    Ai\n";

      if (ipl==0){
	fnspsIND = fnsps;
	fTrkPitchI = fTrkPitch;
      }
      else if (ipl == nplanes-1){
	fnspsCOL = fnsps;
	fTrkPitchC = fTrkPitch;
      }
      double Ai = -1;
      for (int i = 0; i < fnsps; ++i){//loop over all 3D points
	if (ipl==0){
	  fwireIND.push_back(fwire[i]);
	  ftimeIND.push_back(ftime[i]);
	  fstimeIND.push_back(fstime[i]);
	  fetimeIND.push_back(fetime[i]);
	  fMIPsIND.push_back(fMIPs[i]);
	  fdEdxIND.push_back(fdEdx[i]);
	  fResRngIND.push_back(fResRng[i]);
	}
	else if (ipl == nplanes-1){
	  fwireCOL.push_back(fwire[i]);
	  ftimeCOL.push_back(ftime[i]);
	  fstimeCOL.push_back(fstime[i]);
	  fetimeCOL.push_back(fetime[i]);
	  fMIPsCOL.push_back(fMIPs[i]);
	  fdEdxCOL.push_back(fdEdx[i]);
	  fResRngCOL.push_back(fResRng[i]);
	}
	if (i!=0 && i!= fnsps-1){//ignore the first and last point
	  vresRange.push_back(fResRng[i]);
	  vdEdx.push_back(fdEdx[i]);
	  vdQdx.push_back(fdQdx[i]);
	  // Calculate PIDA 
	  if(TrackStops){
	    Ai = fdEdx[i] * pow(fResRng[i],0.42);
	    nPIDA++;
	    PIDA += Ai;
	  }
	}
	CaloPrtHit << std::setw(4) << i
		   <<std::setw(4)  << fwire[i]
		   << std::setw(6) << (int)ftime[i]
		   << std::setiosflags(std::ios::fixed | std::ios::showpoint)
		   << std::setprecision(2)
		   << std::setw(8) << fResRng[i]
		   << std::setprecision(1)
		   << std::setw(8) << fMIPs[i]
		   << std::setprecision(2)
		   << std::setw(8) << fpitch[i]
		   << std::setw(8) << fdEdx[i]
		   << std::setw(8) << Ai
		   << "\n";
      }//end looping over 3D points
      if(nPIDA > 0) {
	PIDA = PIDA / (double)nPIDA;
      } 
      else {
	PIDA = -1;
      }
      CaloPrtTrk << "Plane # "<< ipl
		 << "TrkPitch= "
		 << std::setprecision(2) << fTrkPitch 
		 << " nhits= "        << fnsps
		 << "\n" 
		 << std::setiosflags(std::ios::fixed | std::ios::showpoint)
		 << "Trk Length= "       << std::setprecision(1)
		 << Trk_Length           << " cm,"
		 << " KE calo= "         << std::setprecision(1)
		 << Kin_En               << " MeV,"
		 << " PIDA= "            << PIDA
		 << "\n";
      
      // look for dead wires
      for (unsigned int iw = wire0; iw<wire1+1; ++iw){
	plane = hits[ipl][0]->WireID().Plane;
	tpc   = hits[ipl][0]->WireID().TPC;
	cstat = hits[ipl][0]->WireID().Cryostat;
	channel = geom->PlaneWireToChannel(plane,iw,tpc,cstat);
	if (chanFilt.BadChannel(channel)){
	  mf::LogVerbatim("Calorimetry") << "Found dead wire at Plane = " << plane 
					 << " Wire =" << iw;
	  unsigned int closestwire = 0;
	  unsigned int endwire = 0;
	  int dwire = 100000;
	  double mindis = 100000;
	  double goodresrange = 0;
	  //hitCtr = 0;
	  for (size_t ihit = 0; ihit <hits[ipl].size(); ++ihit){
	    //	for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitsV.begin(); 
	    //	    hitIter != hitsV.end();  
	    //	    ++hitCtr, hitIter++){
	    channel = hits[ipl][ihit]->Wire()->RawDigit()->Channel();
	    if (chanFilt.BadChannel(channel)) continue;
	    // grab the space points associated with this hit
	    std::vector< art::Ptr<recob::SpacePoint> > sppv = fmspts.at(ihit);
	    if(sppv.size() < 1) continue;
	    // only use the first space point in the collection, really each hit should
	    // only map to 1 space point
	    const double xyz[3] = {sppv[0]->XYZ()[0],
				   sppv[0]->XYZ()[1],
				   sppv[0]->XYZ()[2]};
	    double dis1 = (larEnd[0]-xyz[0])*(larEnd[0]-xyz[0]) +
	      (larEnd[1]-xyz[1])*(larEnd[1]-xyz[1]) +
	      (larEnd[2]-xyz[2])*(larEnd[2]-xyz[2]);
	    if (dis1) dis1 = std::sqrt(dis1);
	    if (dis1 < mindis){
	      endwire = hits[ipl][ihit]->WireID().Wire;
	      mindis = dis1;
	    }
	    if (std::abs(wire-iw) < dwire){
	      closestwire = hits[ipl][ihit]->WireID().Wire;
	      dwire = abs(hits[ipl][ihit]->WireID().Wire-iw);
	      goodresrange = dis1;
	    }
	  }
	  if (closestwire){
	    if (iw < endwire){
	      deadwire.push_back(goodresrange+(int(closestwire)-int(iw))*fTrkPitch);
	    }
	    else{
	      deadwire.push_back(goodresrange+(int(iw)-int(closestwire))*fTrkPitch);
	    }
	  }
	}
      }
      calorimetrycol->push_back(anab::Calorimetry(Kin_En,
						  vdEdx,
						  vdQdx,
						  vresRange,
						  deadwire,
						  Trk_Length,
						  fpitch));
      util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);
      
    }//end looping over planes
  }//end looping over tracks
  
  evt.put(std::move(calorimetrycol));
  evt.put(std::move(assn));

  return;
}

//------------------------------------------------------------------------------------//
double calo::Calorimetry::LifetimeCorrection(float time){

  float t = time;

  art::ServiceHandle<util::LArProperties> LArProp;
  art::ServiceHandle<util::DetectorProperties> detprop;

  double timetick = detprop->SamplingRate()*1.e-3;    //time sample in microsec
  double presamplings = detprop->TriggerOffset();

  t -= presamplings;
  time = t * timetick;  //  (in microsec)

  double tau = LArProp->ElectronLifetime();

  double correction = exp(time/tau);
  return correction;
}

//--------------------------------------------------
bool calo::Calorimetry::BeginsOnBoundary(art::Ptr<recob::Track> lar_track)
{
  double fdBoundary = 1.5;
  std::vector<double> larStart, larEnd;
  //put xyz coordinates at begin/end of track into vectors(?)
  lar_track->Extent(larStart,larEnd);

  if(std::abs(larStart[0])            < fdBoundary ||
     std::abs(TPCsize[0]-larStart[0]) < fdBoundary ||
     std::abs(larStart[1]+TPCsize[1]/2) < fdBoundary ||
     std::abs(TPCsize[1]/2-larStart[1]) < fdBoundary ||
     std::abs(larStart[2])            < fdBoundary ||
     std::abs(TPCsize[2]-larStart[2]) < fdBoundary )   
    return true;  
  else return false;
}
    
//--------------------------------------------------
bool calo::Calorimetry::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
{
  double fdBoundary = 1.5;
  std::vector<double> larStart, larEnd;
  //put xyz coordinates at begin/end of track into vectors(?)
  lar_track->Extent(larStart, larEnd); 

  if(std::abs(larEnd[0])              < fdBoundary ||
     std::abs(larEnd[2])              < fdBoundary ||
     std::abs(TPCsize[0]-larEnd[0])   < fdBoundary ||
     std::abs(TPCsize[2]-larEnd[2])   < fdBoundary ||   
     std::abs(larEnd[1]+TPCsize[1]/2) < fdBoundary ||
     std::abs(TPCsize[1]/2-larEnd[1]) < fdBoundary )
    return true;  

  return false;
}

//------------------------------------------------------------------------------------//
void calo::Calorimetry::ReadCaloTree(){
  // This method is an example of how to read 
  // the dynamic vectors in the ROOT Tree

  ftree->Scan("run:event:TrkPitchI:TrkPitchC:XStart:nhits3D:nhitsIND:nhitsCOL");
  int nentries=(int)ftree->GetEntries();
  mf::LogVerbatim("Calorimetry") << "nentries  " << nentries;

  std::vector<double> *MIPsCOL = 0;
  TBranch *bMIPsCOL = 0;
  ftree->SetBranchAddress("MIPsCOL",&MIPsCOL,&bMIPsCOL);

  std::vector<double> *MIPsIND = 0;
  TBranch *bMIPsIND = 0;
  ftree->SetBranchAddress("MIPsIND",&MIPsIND,&bMIPsIND);
  for (int i = 0; i < nentries; ++i) {
    mf::LogVerbatim("Calorimetry") << " entry " << i;
    long int tentry = ftree->LoadTree(i);

    bMIPsCOL->GetEntry(tentry);
    mf::LogVerbatim("Calorimetry") << "# of hits COLL " << MIPsCOL->size();
    for (size_t j = 0; j < MIPsCOL->size(); ++j) {
      mf::LogVerbatim("Calorimetry") << " Coll " << MIPsCOL->at(j);
    }
    bMIPsIND->GetEntry(tentry);
    mf::LogVerbatim("Calorimetry") << "# of hits IND " << MIPsIND->size();
    for (size_t j = 0; j < MIPsIND->size(); ++j) {
      mf::LogVerbatim("Calorimetry") << " Ind " << MIPsIND->at(j);
    }
  }
}

void calo::Calorimetry::GetPitch(art::Ptr<recob::Hit> hit, std::vector<double> trkx, std::vector<double> trky, std::vector<double> trkz, std::vector<double> trkw, std::vector<double> trkx0, double *xyz3d, double &pitch){
  //Get 3d coordinates and track pitch for each hit
  //Find 5 nearest space points and determine xyz and curvature->track pitch

  // Get services
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::DetectorProperties> dp;
  
  //save distance to each spacepoint sorted by distance
  std::map<double,size_t> sptmap;
  //save the sign of distance
  std::map<size_t, int> sptsignmap;

  double wire_pitch = geom->WirePitch(0,1,0);

  double t0 = hit->PeakTime();
  double x0 = dp->ConvertTicksToX(t0, hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat);
  double w0 = hit->WireID().Wire;

  for (size_t i = 0; i<trkx.size(); ++i){
    double distance = pow((trkw[i]-w0)*wire_pitch,2)+pow(trkx0[i]-x0,2);
    if (distance>0) distance = sqrt(distance);
    sptmap.insert(std::pair<double,size_t>(distance,i));
    if (w0-trkw[i]>0) sptsignmap.insert(std::pair<size_t,int>(i,1));
    else sptsignmap.insert(std::pair<size_t,int>(i,-1));
  }

  //x,y,z vs distance
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;
  std::vector<double> vs;

  double kx = 0, ky = 0, kz = 0;

  int np = 0;
  for (auto isp = sptmap.begin(); isp!=sptmap.end(); isp++){
//    const double *xyz = new double[3];
//    xyz = isp->second->XYZ();
    double xyz[3];
    xyz[0] = trkx[isp->second];
    xyz[1] = trky[isp->second];
    xyz[2] = trkz[isp->second];
    
    double distancesign = sptsignmap[isp->second];
    //std::cout<<np<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<(*isp).first<<std::endl;
    if (np==0&&isp->first>30){//hit not on track
      xyz3d[0] = -1000;
      xyz3d[1] = -1000;
      xyz3d[2] = -1000;
      pitch = -1;
      return;
    }
    if (np<5) {
      vx.push_back(xyz[0]);
      vy.push_back(xyz[1]);
      vz.push_back(xyz[2]);
      vs.push_back(isp->first*distancesign);
    }
    else {
      break;
    }
    np++;
    //delete [] xyz;
  }
  //std::cout<<"np="<<np<<std::endl;
  if (np>=2){//at least two points
    TGraph *xs = new TGraph(np,&vs[0],&vx[0]);
    //for (int i = 0; i<np; i++) std::cout<<i<<" "<<vs[i]<<" "<<vx[i]<<" "<<vy[i]<<" "<<vz[i]<<std::endl;
    try{
      if (np>2){
	xs->Fit("pol2","Q");
      }
      else{
	xs->Fit("pol1","Q");
      }
      TF1 *pol = 0;
      if (np>2) pol = (TF1*) xs->GetFunction("pol2");
      else pol = (TF1*) xs->GetFunction("pol1");
      xyz3d[0] = pol->Eval(0);
      kx = pol->GetParameter(1);
      //std::cout<<xyz3d[0]<<" "<<kx<<std::endl;
    }
    catch(...){
      mf::LogWarning("Calorimetry::GetPitch") <<"Fitter failed";
      delete xs;
      xyz3d[0] = vx[0];
    }
    delete xs;
    TGraph *ys = new TGraph(np,&vs[0],&vy[0]);
    try{
      if (np>2){
	ys->Fit("pol2","Q");
      }
      else{
	ys->Fit("pol1","Q");
      }
      TF1 *pol = 0;
      if (np>2) pol = (TF1*) ys->GetFunction("pol2");
      else pol = (TF1*) ys->GetFunction("pol1");
      xyz3d[1] = pol->Eval(0);
      ky = pol->GetParameter(1);
      //std::cout<<xyz3d[1]<<" "<<ky<<std::endl;
    }
    catch(...){
      mf::LogWarning("Calorimetry::GetPitch") <<"Fitter failed";
      delete ys;
      xyz3d[1] = vy[0];
    }
    delete ys;
    TGraph *zs = new TGraph(np,&vs[0],&vz[0]);
    try{
      if (np>2){
	zs->Fit("pol2","Q");
      }
      else{
	zs->Fit("pol1","Q");
      }
      TF1 *pol = 0;
      if (np>2) pol = (TF1*) zs->GetFunction("pol2");
      else pol = (TF1*) zs->GetFunction("pol1");
      xyz3d[2] = pol->Eval(0);
      kz = pol->GetParameter(1);
      //std::cout<<xyz3d[2]<<" "<<kz<<std::endl;
    }
    catch(...){
      mf::LogWarning("Calorimetry::GetPitch") <<"Fitter failed";
      delete zs;
      xyz3d[2] = vz[0];
    }
    delete zs;
  }
  else if (np){
    xyz3d[0] = vx[0];
    xyz3d[1] = vy[0];
    xyz3d[2] = vz[0];
  }
  else{
    xyz3d[0] = -1000;
    xyz3d[1] = -1000;
    xyz3d[2] = -1000;
    pitch = -1;
    return;
  }
//  std::cout<<xyz3d[0]<<" "<<xyz3d[1]<<" "<<xyz3d[2]<<std::endl;
//  double dcos0[3];
//  double dcos1[3];
//  track->Direction(dcos0,dcos1);
  pitch = -1;
  if (kx*kx+ky*ky+kz*kz){
    double tot = sqrt(kx*kx+ky*ky+kz*kz);
    kx /= tot;
    ky /= tot;
    kz /= tot;
    //get pitch
    //for(unsigned int i = 0; i < geom->Nplanes(); ++i){
    //if(geom->Plane(i).View() == hit->View()){
    double wirePitch = geom->WirePitch(0,1,hit->WireID().Plane);
    double angleToVert = geom->Plane(hit->WireID().Plane).Wire(0).ThetaZ(false) - 0.5*TMath::Pi();
    double cosgamma = TMath::Abs(TMath::Sin(angleToVert)*ky+TMath::Cos(angleToVert)*kz);
    if (cosgamma>0) pitch = wirePitch/cosgamma;   
    //std::cout<<kx<<" "<<ky<<" "<<kz<<" pitch = "<<pitch<<std::endl;
	//} // end if the correct view
	//} // end loop over planes

  }
  //std::cout<<kx<<" "<<ky<<" "<<kz<<" "<<dcos0[0]<<" "<<dcos0[1]<<" "<<dcos0[2]<<std::endl;  

}


namespace calo{

  DEFINE_ART_MODULE(Calorimetry)
  
} // end namespace 

