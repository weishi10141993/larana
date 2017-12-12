
#include "larana/T0Finder/AssociationsTools/IndirectHitParticleAssns_tool.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/MCCheater/BackTracker.h"

#include <cmath>
#include <algorithm>

namespace t0
{
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
IndirectHitParticleAssns::IndirectHitParticleAssns(fhicl::ParameterSet const & pset) 
{
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("IndirectHitParticleAssns") << "Configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
IndirectHitParticleAssns::~IndirectHitParticleAssns()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void IndirectHitParticleAssns::reconfigure(fhicl::ParameterSet const & pset)
{
    fHitPartAssnsModuleLabel = pset.get<art::InputTag>("HitPartAssnsLabel");
    fHitModuleLabel          = pset.get<art::InputTag>("HitModuleLabel");
}

//----------------------------------------------------------------------------
/// Rebuild method -> rebuild the basic maps to get truth information
///
/// Arguments:
///
/// event - the art event used to extract all information
///
void IndirectHitParticleAssns::CreateHitParticleAssociations(art::Event& evt, HitParticleAssociations* hitPartAssns)
{
    // This function handles the "indirect" creation of hit<-->MCParticle associations by trying to
    // use the previously created Hit<-->MCParticle associations
    // Its pretty much a brute force approach here... but time is short!
    //
    // First step is to recover the preexisting associations

    // Get a handle for the associations...
    art::Handle<HitParticleAssociations> partHitAssnsHandle;
    
    evt.getByLabel(fHitPartAssnsModuleLabel, partHitAssnsHandle);
    
    if (!partHitAssnsHandle.isValid())
    {
        throw cet::exception("IndirectHitParticleAssns") << "===>> NO MCParticle <--> Hit associations found for run/subrun/event: " << evt.run() << "/" << evt.subRun() << "/" << evt.id().event();
    }
    
    // Look up the hits we want to process as well, since if they are not there then no point in proceeding
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    evt.getByLabel(fHitModuleLabel,hitListHandle);
    
    if(!hitListHandle.isValid()){
        throw cet::exception("IndirectHitParticleAssns") << "===>> NO Hit collection found to process for run/subrun/event: " << evt.run() << "/" << evt.subRun() << "/" << evt.id().event();
    }
    
    // Go through the associations and build out our (hopefully sparse) data structure
    using ParticleDataPair         = std::pair<const art::Ptr<simb::MCParticle>*, const anab::BackTrackerHitMatchingData*>;
    using MCParticleDataSet        = std::set<ParticleDataPair>;
    using TickToPartDataMap        = std::unordered_map<raw::TDCtick_t,   MCParticleDataSet>;
    using ChannelToTickPartDataMap = std::unordered_map<raw::ChannelID_t, TickToPartDataMap>;
    
    ChannelToTickPartDataMap chanToTickPartDataMap;
    
    // Build out the maps between hits/particles
    for(HitParticleAssociations::assn_iterator partHitItr = partHitAssnsHandle->begin(); partHitItr != partHitAssnsHandle->end(); ++partHitItr)
    {
        const art::Ptr<simb::MCParticle>&       mcParticle = partHitItr->first;
        const art::Ptr<recob::Hit>&             recoHit    = partHitItr->second;
        const anab::BackTrackerHitMatchingData* data       = &partHitAssnsHandle->data(partHitItr);
        
        TickToPartDataMap& tickToPartDataMap = chanToTickPartDataMap[recoHit->Channel()];
        
        for(raw::TDCtick_t tick = recoHit->PeakTimeMinusRMS(); tick <= recoHit->PeakTimePlusRMS(); tick++)
        {
            tickToPartDataMap[tick].insert(ParticleDataPair(&mcParticle,data));
        }
    }
    
    // Armed with the map, process the hit list
    for(size_t hitIdx = 0; hitIdx < hitListHandle->size(); hitIdx++)
    {
        art::Ptr<recob::Hit> hit(hitListHandle,hitIdx);
        
        TickToPartDataMap& tickToPartDataMap = chanToTickPartDataMap[hit->Channel()];
        
        if (tickToPartDataMap.empty())
        {
            mf::LogInfo("IndirectHitParticleAssns") << "No channel information found for hit " << hit << "\n";
            continue;
        }

        // Keep track of results
        MCParticleDataSet particleDataSet;
        
        // Go through the ticks in this hit and recover associations
        for(raw::TDCtick_t tick = hit->PeakTimeMinusRMS(); tick <= hit->PeakTimePlusRMS(); tick++)
        {
            TickToPartDataMap::iterator hitInfoItr = tickToPartDataMap.find(tick);
            
            if (hitInfoItr != tickToPartDataMap.end())
            {
                for(const auto& partData : hitInfoItr->second) particleDataSet.insert(partData);
            }
        }
        
        // Now create new associations for the hit in question
        for(const auto& partData : particleDataSet)
            hitPartAssns->addSingle(*partData.first, hit, *partData.second);
    }

    return;
}

//----------------------------------------------------------------------------
    
//DEFINE_ART_CLASS_TOOL(IndirectHitParticleAssns)
}
