////////////////////////////////////////////////////////////////////////
//
// \file MVAPID_module.cc
//
// m.haigh@warwick.ac.uk
//
///////////////////////////////////////////////////////////////////////

#ifndef MVAPID_H
#define MVAPID_H

// ### Generic C++ includes ###
#include <iostream>

// ### Framework includes ###
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "TTree.h"

#include "MVAAlg.h"
#include "AnalysisBase/MVAResult.h"
#include "RecoBase/Track.h"
#include "Utilities/AssociationUtil.h"

namespace mvapid {

  class MVAPID : public art::EDProducer {

  public:

    explicit MVAPID(fhicl::ParameterSet const& pset);
    virtual ~MVAPID();                               
    void beginJob();        
    void beginRun(art::Run& run); 
    void reconfigure(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);                   
   
   
 private:

    MVAAlg fAlg;
    std::vector<anab::MVAResult>* fResult;
    unsigned int fRun,fSubrun,fEvent;
    TTree* fTree;
    
  }; // class MVAPID


 
//------------------------------------------------------------------------------
  MVAPID::MVAPID(fhicl::ParameterSet const& pset): fAlg(pset,this)
{
  this->reconfigure(pset);
  produces< std::vector<anab::MVAResult> >();
  produces< art::Assns<recob::Track, anab::MVAResult, void> >();
  fResult=new std::vector<anab::MVAResult>;
}

//------------------------------------------------------------------------------
void MVAPID::reconfigure(fhicl::ParameterSet const& pset) 
{
  return;
}

//------------------------------------------------------------------------------
MVAPID::~MVAPID()
{
}

// ***************** //
void MVAPID::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  fTree =tfs->make<TTree>("MVAPID","Results");/**All-knowing tree with reconstruction information*/

  fTree->Branch("run",&fRun,"run/I");
  fTree->Branch("subrun",&fSubrun,"subrun/I");
  fTree->Branch("event",&fEvent,"event/I");
  fTree->Branch("MVAResult",&fResult);
}


void MVAPID::beginRun(art::Run&)
{
  
}

// ***************** //
void MVAPID::produce(art::Event& evt)
{ 
  std::unique_ptr<std::vector<anab::MVAResult> > result(new std::vector<anab::MVAResult>);
  std::unique_ptr< art::Assns<recob::Track, anab::MVAResult> > assns(new art::Assns<recob::Track, anab::MVAResult>);

  fRun = evt.id().run();
  fSubrun = evt.id().subRun();
  fEvent = evt.id().event();

  fAlg.RunPID(evt,*result,*assns);
  *fResult=*result;
  fTree->Fill();
  evt.put(std::move(result));
  evt.put(std::move(assns));
}
  
  DEFINE_ART_MODULE(MVAPID)
  
  
} //namespace mvapid

#endif // MVAPID_H
