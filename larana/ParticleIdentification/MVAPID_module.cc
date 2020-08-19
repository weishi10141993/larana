////////////////////////////////////////////////////////////////////////
//
// \file MVAPID_module.cc
//
// m.haigh@warwick.ac.uk
//
///////////////////////////////////////////////////////////////////////

// ### Generic C++ includes ###

// ### Framework includes ###
#include "TTree.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/fwd.h"

#include "MVAAlg.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "lardataobj/RecoBase/Track.h"

namespace mvapid {

  class MVAPID : public art::EDProducer {

  public:
    explicit MVAPID(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);

  private:
    MVAAlg fAlg;
    std::vector<anab::MVAPIDResult>* fResult;
    unsigned int fRun, fSubrun, fEvent;
    TTree* fTree;

  }; // class MVAPID

  //------------------------------------------------------------------------------
  MVAPID::MVAPID(fhicl::ParameterSet const& pset) : EDProducer{pset}, fAlg(pset)
  {
    produces<std::vector<anab::MVAPIDResult>>();
    produces<art::Assns<recob::Track, anab::MVAPIDResult, void>>();
    produces<art::Assns<recob::Shower, anab::MVAPIDResult, void>>();
    fResult = new std::vector<anab::MVAPIDResult>;
  }

  // ***************** //
  void
  MVAPID::beginJob()
  {
    art::ServiceHandle<art::TFileService const> tfs;
    fTree =
      tfs->make<TTree>("MVAPID", "Results"); /**All-knowing tree with reconstruction information*/
    fTree->Branch("run", &fRun, "run/I");
    fTree->Branch("subrun", &fSubrun, "subrun/I");
    fTree->Branch("event", &fEvent, "event/I");
    fTree->Branch("MVAResult", &fResult);
    fAlg.GetDetectorEdges();
    fAlg.GetWireNormals();
  }

  // ***************** //
  void
  MVAPID::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<anab::MVAPIDResult>> result(new std::vector<anab::MVAPIDResult>);
    std::unique_ptr<art::Assns<recob::Track, anab::MVAPIDResult>> trackAssns(
      new art::Assns<recob::Track, anab::MVAPIDResult>);
    std::unique_ptr<art::Assns<recob::Shower, anab::MVAPIDResult>> showerAssns(
      new art::Assns<recob::Shower, anab::MVAPIDResult>);
    fRun = evt.id().run();
    fSubrun = evt.id().subRun();
    fEvent = evt.id().event();
    fAlg.RunPID(evt, *result, *trackAssns, *showerAssns);
    *fResult = *result;
    fTree->Fill();
    evt.put(std::move(result));
    evt.put(std::move(trackAssns));
    evt.put(std::move(showerAssns));
  }

  DEFINE_ART_MODULE(MVAPID)

} //namespace mvapid
