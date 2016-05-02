#include "TFile.h"
#include "TTree.h"
#include "TMVA/Factory.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "MVAPIDResult.h"

void BuildTree(std::string inFile,std::string outFile){
  TFile* fIn=new TFile(inFile.c_str());
  TTree* tr=(TTree*)fIn->Get("pid/MVAPID");
  std::vector<anab::MVAPIDResult>* mvares=0;
  tr->SetBranchAddress("MVAResult",&mvares);
  TFile* fOut=new TFile(outFile.c_str(),"RECREATE");
  TTree* mvaTree = new TTree("mvaTree","mvaTree");

  float evalRatio, concentration, coreHaloRatio, conicalness;
  float dEdxStart, dEdxEnd, dEdxEndRatio;
  float length;
  int isTrack, isStoppingReco;

  mvaTree->Branch("evalRatio",&evalRatio,"evalRatio/F");
  mvaTree->Branch("concentration",&concentration,"concentration/F");
  mvaTree->Branch("coreHaloRatio",&coreHaloRatio,"coreHaloRatio/F");
  mvaTree->Branch("conicalness",&conicalness,"conicalness/F");
  mvaTree->Branch("dEdxStart",&dEdxStart,"dEdxStart/F");
  mvaTree->Branch("dEdxEnd",&dEdxEnd,"dEdxEnd/F");
  mvaTree->Branch("dEdxEndRatio",&dEdxEndRatio,"dEdxEndRatio/F");
  mvaTree->Branch("length",&length,"length/F");
  mvaTree->Branch("isTrack",&isTrack,"isTrack/I");
  mvaTree->Branch("isStoppingReco",&isStoppingReco,"isStoppingReco/I");

 for(int iEntry=0;iEntry<tr->GetEntries();++iEntry){
   tr->GetEntry(iEntry);
   if(!mvares->size()) continue;

   anab::MVAPIDResult* biggestTrack=&((*mvares)[0]);
   for(unsigned int iRes=0;iRes!=mvares->size();++iRes){
     if((((*mvares)[iRes]).nSpacePoints)>(biggestTrack->nSpacePoints)){
       biggestTrack=&((*mvares)[iRes]);
     }
   }

   evalRatio=biggestTrack->evalRatio;
   concentration=biggestTrack->concentration;
   coreHaloRatio=biggestTrack->coreHaloRatio;
   conicalness=biggestTrack->conicalness;
   dEdxStart=biggestTrack->dEdxStart;
   dEdxEnd=biggestTrack->dEdxEnd;
   dEdxEndRatio=biggestTrack->dEdxEndRatio;
   length=biggestTrack->length;
   isTrack=biggestTrack->isTrack;
   isStoppingReco=biggestTrack->isStoppingReco;

   mvaTree->Fill();
 }

 fIn->Close();
 fOut->Write();
 fOut->Close();
}

void TrainMVA(std::vector<std::string> signalFiles,std::vector<std::string> backgroundFiles,std::string outputFile,std::string jobName){

  TFile* fOut = new TFile(outputFile.c_str(),"RECREATE");
  TMVA::Factory* factory = new TMVA::Factory( jobName.c_str(), fOut, "" );

  std::vector<TTree*> sigTrees;

  for(std::vector<std::string>::iterator fIter=signalFiles.begin();fIter!=signalFiles.end();++fIter){
    TFile* fIn=new TFile(fIter->c_str());
    factory->AddSignalTree((TTree*)fIn->Get("mvaTree"));
  }

  for(std::vector<std::string>::iterator fIter=backgroundFiles.begin();fIter!=backgroundFiles.end();++fIter){
    TFile* fIn=new TFile(fIter->c_str());
    factory->AddBackgroundTree((TTree*)fIn->Get("mvaTree"));
  }

   factory->AddVariable("evalRatio",'F');   
   factory->AddVariable("concentration",'F');     
   factory->AddVariable("coreHaloRatio",'F');	       
   factory->AddVariable("conicalness",'F');
   factory->AddVariable("dEdxStart",'F');
   factory->AddVariable("dEdxEnd",'F');
   factory->AddVariable("dEdxEndRatio",'F');	       

   factory->BookMethod( TMVA::Types::kTMlpANN, "ANN", "" );
   factory->BookMethod( TMVA::Types::kBDT, "BDT", "" );
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();
   fOut->Write();
   fOut->Close();
}


void PrintRes(std::string inFile){

  TFile* fIn=new TFile(inFile.c_str());
  TTree* tr=(TTree*)fIn->Get("pid/ANAB");
  std::vector<anab::MVAPIDResult>* mvares=0;
  tr->SetBranchAddress("MVAResult",&mvares);

 for(int iEntry=0;iEntry<tr->GetEntries();++iEntry){
   tr->GetEntry(iEntry);
   if(!mvares->size()) continue;

   anab::MVAPIDResult* biggestTrack=&((*mvares)[0]);
   for(unsigned int iRes=0;iRes!=mvares->size();++iRes){
     if((((*mvares)[iRes]).nSpacePoints)>(biggestTrack->nSpacePoints)){
       biggestTrack=&((*mvares)[iRes]);
     }
   }

   std::cout<<biggestTrack->mvaOutput.at(string("ANN"))<<std::endl;
 }

 fIn->Close();
}
  /*
  std::vector<std::string> methods;
  methods.push_back("TMlp_ANN");
  methods.push_back("BDT");
  //methods.push_back("Cuts");
  //methods.push_back("kNN");
  //methods.push_back("Likelihood");
  gStyle->SetOptStat("0000");

  TFile fOut(outFileName.c_str(),"RECREATE");

  TSystemDirectory dir("./","./");
  TIter nextFile(dir.GetListOfFiles());
  TSystemFile* file;
  while(file = (TSystemFile*)nextFile.Next()){
    std::string fName=file->GetName();
    std::cout<<fName<<std::endl;
    if(fName.find("mva_")!=std::string::npos){
      fOut.cd();
      std::string caName=fName.substr(0,fName.size()-5);
      TCanvas* ca=new TCanvas(caName.c_str(),caName.c_str());
      TLegend* le=new TLegend(0.2,0.1,0.4,0.3);
      TFile f(fName.c_str());
      int iMethod=0;
      for(std::vector<std::string>::iterator mIter=methods.begin();mIter!=methods.end();++mIter){
	std::string methodFolder;
	if(mIter->find("ANN")!=std::string::npos){
	  methodFolder="TMlpANN";
	}
	else{
	  methodFolder=*mIter;
	}
	std::cout<<(std::string("Method_")+methodFolder+"/"+*mIter+"/MVA_"+*mIter+"_rejBvsS").c_str()<<std::endl;
	TH1D* h=(TH1D*)f.Get((std::string("Method_")+methodFolder+"/"+*mIter+"/MVA_"+*mIter+"_rejBvsS").c_str());
	le->AddEntry(h);
	if(iMethod==0){
	  h->SetTitle(caName.c_str());
	  h->SetLineColor(kRed);
	  h->Draw();
	  h->GetYaxis().SetRangeUser(0,1.2);
	  firstMethod=false;
	}
	else{
	  h->SetLineColor(iMethod==1?kBlue:kGreen);
	  h->Draw("same");
	}
	++iMethod;
      }
      le->Draw();
      fOut.cd();
      std::cout<<caName.c_str()<<std::endl;
      ca->Write(caName.c_str());
      ca->SaveAs((std::string("./plots/")+caName+".png").c_str());
    }
  }
  fOut.Close();
}
*/
