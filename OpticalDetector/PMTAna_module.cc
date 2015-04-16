/**
 * \file PMTAna_module.cc
 *
 * \ingroup PMTAna
 * 
 * \brief Class definition file of PMTAna
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PMTAna

@{*/
#ifndef PMTAna_H
#define PMTAna_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

// LArSoft
//#include "OpticalDetectorData/FIFOChannel.h"
#include "RawData/OpDetWaveform.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/Geometry.h"
// STL
#include <set>
#include <vector>
#include <cmath>
#include <functional>
#include <numeric>


// ROOT
#include <TString.h>
#include <TTree.h>
// My modules
#include "AlgoThreshold.h"
#include "AlgoFixedWindow.h"
#include "PulseRecoManager.h"

namespace pmtana {

  /**
     \class PMTAna
     PMTAna module to copy LArSoft data contents into LArLight data formatted file
  */ 
  class PMTAna : public art::EDAnalyzer{

  public:

    /// Constructor
    PMTAna(const fhicl::ParameterSet&);

    /// Destructor
    virtual ~PMTAna();

    /// Function to be called before an event loop
    void beginJob(){}

    /// Function to be called per sub run
    void beginSubRun(const art::SubRun& /*srun*/){}

    /// Function to be called per event
    void analyze (const art::Event&);

  private:


    /// Function to clear event-wise variables
    void ClearEventData();

    std::string _fifo_mod_name; ///< Input FIFOChannel producer name
    TTree*      _tree;          ///< output data holder TTree

    PulseRecoManager _preco_man;
    AlgoThreshold     _th_algo;
    AlgoFixedWindow   _fw_algo;
  };
  
}

#endif//  PMTAna_H

// PMTAna.cc

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace pmtana {
  DEFINE_ART_MODULE(PMTAna)
}

namespace pmtana {

  //#######################################################################################################
  PMTAna::PMTAna(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset), 
    _preco_man(),
    _th_algo(),
    _fw_algo()
  //#######################################################################################################
  {

    // Obtain module names for input data
    _fifo_mod_name = pset.get<std::string>("fModName_FIFOChannel" );

    // Next we make storage data class objects for those data types specified in fcl files.
    art::ServiceHandle<art::TFileService>  fileService;

    // Create TTree
    _tree = fileService->make<TTree>("pmt_tree","Analysis Tree");

    //
    // Demonstration purpose ... 
    //
    _preco_man.AddRecoAlgo(&_th_algo);
    _preco_man.AddRecoAlgo(&_fw_algo);
    _preco_man.SetPedAlgo(kHEAD);

  }

  PMTAna::~PMTAna(){}

  //#######################################################################################################
  void PMTAna::ClearEventData()
  //#######################################################################################################
  {
  }


  //#######################################################################################################
  void PMTAna::analyze(const art::Event& evt) 
  //#######################################################################################################
  {

    //data_ptr->set_event(evt.id().event(), evt.run(), evt.subRun());

//    std::vector<const optdata::FIFOChannel*> pmtArray;
    std::vector<const raw::OpDetWaveform*> pmtArray;
    try{

      evt.getView(_fifo_mod_name,pmtArray);

    }catch (art::Exception const& e) {

      if (e.categoryCode() != art::errors::ProductNotFound ) throw;

    }




    for(size_t i=0; i<pmtArray.size(); ++i) {

//      const optdata::FIFOChannel* fifo_ptr(pmtArray.at(i));
      const raw::OpDetWaveform* fifo_ptr(pmtArray.at(i));

      _preco_man.RecoPulse(*fifo_ptr);

      //
      // here I add code to store reco-ed pulse w/ channel number
      // OR I may make a singleton storage manager...

      /*
      data_ptr->add_pmtfifo(fifo_ptr->ChannelNumber(),
			    fifo_ptr->Category(),
			    fifo_ptr->Frame(),
			    fifo_ptr->TimeSlice(),
			    *fifo_ptr);
			    */

      //
      //
      //
    }
    
    
  }

}
