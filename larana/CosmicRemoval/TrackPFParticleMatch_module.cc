////////////////////////////////////////////////////////////////////////
// Class:       TrackPFParticleMatch
// Module Type: producer
// File:        TrackPFParticleMatch_module.cc
//              This module aims to enhance cosmic ray rejection by
//              matching tracks which were fit from "other" producers
//              to PFParticles by comparing hits
//              The output will be a set of associations relating the two
//
// Generated at Mon Mar 28 19:17:00 2016 by Tracy Usher by cloning CosmicTrackTagger
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::TrackPFParticleMatch
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/geo.h"

#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"

#include "lardata/AnalysisBase/CosmicTag.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"

namespace cosmic
{

class TrackPFParticleMatch : public art::EDProducer
{
public:
    explicit TrackPFParticleMatch(fhicl::ParameterSet const & p);
    virtual ~TrackPFParticleMatch();

    void produce(art::Event & e) override;

    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void endJob() override;

private:
    std::string fPFParticleModuleLabel;
    std::string fTrackModuleLabel;
};


TrackPFParticleMatch::TrackPFParticleMatch(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    produces< art::Assns<recob::Track, recob::PFParticle>>();
}

TrackPFParticleMatch::~TrackPFParticleMatch()
{
    // Clean up dynamic memory and other resources here.
}

void TrackPFParticleMatch::produce(art::Event & evt)
{
    // Instatiate the output
    std::unique_ptr< art::Assns<recob::Track, recob::PFParticle > > trackPFParticleAssns( new art::Assns<recob::Track, recob::PFParticle>);
    
    // Recover handle for PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    evt.getByLabel( fPFParticleModuleLabel, pfParticleHandle);
    
    // Recover the clusters so we can do associations to the hits
    // In theory the clusters come from the same producer as the PFParticles
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    evt.getByLabel(fPFParticleModuleLabel, clusterHandle);
    
    if (!(pfParticleHandle.isValid() && clusterHandle.isValid()))
    {
        evt.put( std::move(trackPFParticleAssns) );
        return;
    }
    
    // Recover the handle for the tracks
    art::Handle<std::vector<recob::Track> > trackHandle;
    evt.getByLabel( fTrackModuleLabel, trackHandle);
    
    if (!trackHandle.isValid())
    {
        evt.put( std::move(trackPFParticleAssns) );
        return;
    }

    // The strategy will be to build sets of maps which will enable us to go from hit to either Track or PFParticle
    using HitToTrackMap = std::map<int, std::vector<art::Ptr<recob::Track>>>;   // Allows for hit sharing
    
    // Build out the track/hit map first
    // Recover the track/hit associations
    art::FindManyP<recob::Hit>  trackHitAssns(trackHandle, evt, fTrackModuleLabel);
    
    // Get an empty map
    HitToTrackMap hitToTrackMap;
    
    // Fill it
    for(size_t idx = 0; idx < trackHandle->size(); idx++)
    {
        art::Ptr<recob::Track> track(trackHandle, idx);
        
        std::vector<art::Ptr<recob::Hit>> trackHitVec = trackHitAssns.at(track.key());
        
        for (const auto& hitPtr : trackHitVec)
        {
            hitToTrackMap[hitPtr.key()].push_back(track);
        }
    }
    
    // Now we walk up the remaining list of associations needed to go from CR tags to hits
    // We'll need to go from tracks to PFParticles
    // From PFParticles we go to clusters
    art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, evt, fPFParticleModuleLabel);
    
    // Likewise, recover the collection of associations to hits
    art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, evt, fPFParticleModuleLabel);
    
    // We'll store this info in a rather complicated data structure...
    using TrackToPFParticleHitMap = std::map<art::Ptr<recob::Track>, std::map<art::Ptr<recob::PFParticle>,std::vector<art::Ptr<recob::Hit>>>>;
    
    TrackToPFParticleHitMap trackToPFParticleHitMap;
    
    // Now we go through the PFParticles and match hits to tracks
    for(size_t idx = 0; idx < pfParticleHandle->size(); idx++)
    {
        art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, idx);
        
        // Recover associated clusters
        std::vector<art::Ptr<recob::Cluster>> clusterVec = clusterAssns.at(pfParticle.key());
        
        // Loop through the clusters
        for(const auto& cluster : clusterVec)
        {
            // Recover associated hits
            std::vector<art::Ptr<recob::Hit>> hitVec = clusterHitAssns.at(cluster.key());
            
            // Loop through hits and look for associated track
            for(const auto& hit : hitVec)
            {
                HitToTrackMap::iterator hitToTrackItr = hitToTrackMap.find(hit.key());
                
                if (hitToTrackItr != hitToTrackMap.end())
                {
                    for(auto& track : hitToTrackItr->second)
                        trackToPFParticleHitMap[track][pfParticle].push_back(hit);
                }
            }
        }
    }
    
    // Ok, now we can create the associations
    for (auto& trackMapItr : trackToPFParticleHitMap)
    {
        std::vector<art::Ptr<recob::PFParticle>> pfParticleVec;
        
        for(auto& pfParticleMapItr : trackMapItr.second)
        {
            // We need to make sure we don't associate the case where the hits are from crossing tracks
            // which would be an illegal association.
            // Two things to check: that more than one view is matched in the hits
            // That some fraction of the matched hits are from the track
            int nHitsPerView[] = {0,0,0};
            
            for(const auto& hit : pfParticleMapItr.second)
            {
                nHitsPerView[hit->View()]++;
            }
            
            int nViewsWithHits(0);
            
            for(auto& nHits : nHitsPerView) if (nHits > 0) nViewsWithHits++;
            
            if (nViewsWithHits < 2) continue;
            
            // Get the hits associated to the track again
            std::vector<art::Ptr<recob::Hit>> trackHitVec = trackHitAssns.at(trackMapItr.first.key());
            
            // Fraction of hits from track shared with PFParticle
            float sharedHitFrac = float(pfParticleMapItr.second.size()) / float(trackHitVec.size());

            if (sharedHitFrac < 0.3) continue;
            
            // Ok, this is an association to make
            pfParticleVec.push_back(pfParticleMapItr.first);
        }
        
        util::CreateAssn(*this, evt, trackMapItr.first, pfParticleVec, *trackPFParticleAssns);
    }
    
    evt.put( std::move(trackPFParticleAssns));
    
    return;

} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////

void TrackPFParticleMatch::beginJob()
{
}

void TrackPFParticleMatch::reconfigure(fhicl::ParameterSet const & p)
{
    fPFParticleModuleLabel = p.get< std::string >("PFParticleModuleLabel");
    fTrackModuleLabel      = p.get< std::string >("TrackModuleLabel", "track");
}

void TrackPFParticleMatch::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TrackPFParticleMatch)
    
}
