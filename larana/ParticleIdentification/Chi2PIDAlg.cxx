////////////////////////////////////////////////////////////////////////
//
// Chi2PIDAlg class
//
// tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

//#include "RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larana/ParticleIdentification/Chi2PIDAlg.h"

// ROOT includes
#include "TFile.h"
#include "TProfile.h"
#include "TMath.h"

// Framework includes
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/GeometryUtilities.h"

//------------------------------------------------------------------------------
pid::Chi2PIDAlg::Chi2PIDAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

//------------------------------------------------------------------------------
pid::Chi2PIDAlg::~Chi2PIDAlg()
{
}

//------------------------------------------------------------------------------
std::bitset<8> pid::Chi2PIDAlg::GetBitset(geo::PlaneID planeID){

    std::bitset<8> thisBitset;

    thisBitset.set(planeID.Plane);

    return thisBitset;

}

//------------------------------------------------------------------------------
void pid::Chi2PIDAlg::reconfigure(fhicl::ParameterSet const& pset)
{
  fTemplateFile           = pset.get< std::string >("TemplateFile");
  fUseMedian              = pset.get< bool >("UseMedian");
  //fCalorimetryModuleLabel = pset.get< std::string >("CalorimetryModuleLabel");

  cet::search_path sp("FW_SEARCH_PATH");
 
  if( !sp.find_file(fTemplateFile, fROOTfile) )
    throw cet::exception("Chi2ParticleID") << "cannot find the root template file: \n" 
					   << fTemplateFile
					   << "\n bail ungracefully.\n";
  TFile *file = TFile::Open(fROOTfile.c_str());
  dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
  dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
  dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
  dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

//  std::cout<<"Chi2PIDAlg configuration:"<<std::endl;
//  std::cout<<"Template file: "<<fROOTfile<<std::endl;
//  std::cout<<"fUseMedian: "<<fUseMedian<<std::endl;

  return;
}


//------------------------------------------------------------------------------
anab::ParticleID pid::Chi2PIDAlg::DoParticleID(std::vector<art::Ptr<anab::Calorimetry>> calos){

  std::vector<anab::sParticleIDAlgScores> AlgScoresVec;
  
  for (size_t i_calo = 0; i_calo < calos.size(); i_calo++){

  art::Ptr<anab::Calorimetry> calo = calos.at(i_calo);

  int npt = 0;
  double chi2pro = 0;
  double chi2ka = 0;
  double chi2pi = 0;
  double chi2mu = 0;
  double avgdedx = 0;
  double PIDA = 0; //by Bruce Baller
  std::vector<double> vpida;
  std::vector<float> trkdedx = calo->dEdx();
  std::vector<float> trkres = calo->ResidualRange();
  std::vector<float> deadwireresrc = calo->DeadWireResRC();

  int used_trkres = 0;
  for (unsigned i = 0; i<trkdedx.size(); ++i){//hits
    //ignore the first and the last point
    if (i==0 || i==trkdedx.size()-1) continue;
    avgdedx += trkdedx[i];
    if(trkres[i] < 30) {
      PIDA += trkdedx[i]*pow(trkres[i],0.42);
      vpida.push_back(trkdedx[i]*pow(trkres[i],0.42));
      used_trkres++;
    }
    if (trkdedx[i]>1000) continue; //protect against large pulse height
    int bin = dedx_range_pro->FindBin(trkres[i]);
    if (bin>=1&&bin<=dedx_range_pro->GetNbinsX()){
      double bincpro = dedx_range_pro->GetBinContent(bin);
      if (bincpro<1e-6){//for 0 bin content, using neighboring bins
	bincpro = (dedx_range_pro->GetBinContent(bin-1)+dedx_range_pro->GetBinContent(bin+1))/2;
      }
      double bincka = dedx_range_ka->GetBinContent(bin);
      if (bincka<1e-6){
	bincka = (dedx_range_ka->GetBinContent(bin-1)+dedx_range_ka->GetBinContent(bin+1))/2;
      }
      double bincpi = dedx_range_pi->GetBinContent(bin);
      if (bincpi<1e-6){
	bincpi = (dedx_range_pi->GetBinContent(bin-1)+dedx_range_pi->GetBinContent(bin+1))/2;
      }
      double bincmu = dedx_range_mu->GetBinContent(bin);
      if (bincmu<1e-6){
	bincmu = (dedx_range_mu->GetBinContent(bin-1)+dedx_range_mu->GetBinContent(bin+1))/2;
      }
      double binepro = dedx_range_pro->GetBinError(bin);
      if (binepro<1e-6){
	binepro = (dedx_range_pro->GetBinError(bin-1)+dedx_range_pro->GetBinError(bin+1))/2;
      }
      double bineka = dedx_range_ka->GetBinError(bin);
      if (bineka<1e-6){
	bineka = (dedx_range_ka->GetBinError(bin-1)+dedx_range_ka->GetBinError(bin+1))/2;
      }
      double binepi = dedx_range_pi->GetBinError(bin);
      if (binepi<1e-6){
	binepi = (dedx_range_pi->GetBinError(bin-1)+dedx_range_pi->GetBinError(bin+1))/2;
      }
      double binemu = dedx_range_mu->GetBinError(bin);
      if (binemu<1e-6){
	binemu = (dedx_range_mu->GetBinError(bin-1)+dedx_range_mu->GetBinError(bin+1))/2;
      }
      //double errke = 0.05*trkdedx[i];   //5% KE resolution
      double errdedx = 0.04231+0.0001783*trkdedx[i]*trkdedx[i]; //resolution on dE/dx
      errdedx *= trkdedx[i];
      chi2pro += pow((trkdedx[i]-bincpro)/std::sqrt(pow(binepro,2)+pow(errdedx,2)),2);
      chi2ka += pow((trkdedx[i]-bincka)/std::sqrt(pow(bineka,2)+pow(errdedx,2)),2);
      chi2pi += pow((trkdedx[i]-bincpi)/std::sqrt(pow(binepi,2)+pow(errdedx,2)),2);
      chi2mu += pow((trkdedx[i]-bincmu)/std::sqrt(pow(binemu,2)+pow(errdedx,2)),2);
      //std::cout<<i<<" "<<trkdedx[i]<<" "<<trkres[i]<<" "<<bincpro<<std::endl;
      ++npt;
    }
  }

  anab::sParticleIDAlgScores chi2proton;
  anab::sParticleIDAlgScores chi2kaon;
  anab::sParticleIDAlgScores chi2pion;
  anab::sParticleIDAlgScores chi2muon;
  anab::sParticleIDAlgScores pida_mean;
  anab::sParticleIDAlgScores pida_median;

  //anab::ParticleID pidOut;
  if (npt){
  
    chi2proton.fAlgName = "Chi2";
    chi2proton.fVariableType = anab::kGOF;
    chi2proton.fTrackDir = anab::kForward;
    chi2proton.fAssumedPdg = 2212;
    chi2proton.fPlaneMask = GetBitset(calo->PlaneID()); 
    chi2proton.fNdf = npt;
    chi2proton.fValue = chi2pro/npt;

    chi2muon.fAlgName = "Chi2";
    chi2muon.fVariableType = anab::kGOF;
    chi2muon.fTrackDir = anab::kForward;
    chi2muon.fAssumedPdg = 13;
    chi2muon.fPlaneMask = GetBitset(calo->PlaneID()); 
    chi2muon.fNdf = npt;
    chi2muon.fValue = chi2mu/npt;

    chi2kaon.fAlgName = "Chi2";
    chi2kaon.fVariableType = anab::kGOF;
    chi2kaon.fTrackDir = anab::kForward;
    chi2kaon.fAssumedPdg = 321;
    chi2kaon.fPlaneMask = GetBitset(calo->PlaneID()); 
    chi2kaon.fNdf = npt;
    chi2kaon.fValue = chi2ka/npt;

    chi2pion.fAlgName = "Chi2";
    chi2pion.fVariableType = anab::kGOF;
    chi2pion.fTrackDir = anab::kForward;
    chi2pion.fAssumedPdg = 211;
    chi2pion.fPlaneMask = GetBitset(calo->PlaneID()); 
    chi2pion.fNdf = npt;
    chi2pion.fValue = chi2pi/npt;

    AlgScoresVec.push_back(chi2proton);
    AlgScoresVec.push_back(chi2muon);
    AlgScoresVec.push_back(chi2kaon);
    AlgScoresVec.push_back(chi2pion);
  }

  //if (trkdedx.size()) pidOut.fPIDA = PIDA/trkdedx.size();
  if(used_trkres > 0){
    if (fUseMedian){
      pida_median.fAlgName = "PIDA_median";
      pida_median.fVariableType = anab::kPIDA;
      pida_median.fTrackDir = anab::kForward;
      pida_median.fValue = TMath::Median(vpida.size(), &vpida[0]);
      pida_median.fPlaneMask = GetBitset(calo->PlaneID());
      AlgScoresVec.push_back(pida_median);
    }
    else{ // use mean
      pida_mean.fAlgName = "PIDA_mean";
      pida_mean.fVariableType = anab::kPIDA;
      pida_mean.fTrackDir = anab::kForward;
      pida_mean.fValue = PIDA/used_trkres;
      pida_mean.fPlaneMask = GetBitset(calo->PlaneID());
      AlgScoresVec.push_back(pida_mean);
    } 
  }

  }

  anab::ParticleID pidOut(AlgScoresVec);
 
  return pidOut;

}

