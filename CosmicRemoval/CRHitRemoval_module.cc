////////////////////////////////////////////////////////////////////////
//
// Class:       CRHitRemoval
// Module Type: producer
// File:        CRHitRemoval_module.cc
//
// This module produces RecoBase/Hit objects after removing those
// deemed to be due to CR muons.
//
// Configuration parameters:
//
// CosmicProducerLabels    - a list of cosmic ray producers which should be or'ed
// FlashProducerLabel      - Flash tagger module for handling special case
// HitProducerLabel        - the producer of the recob::Hit objects
// PFParticleProducerLabel - the producer of the recob::PFParticles to consider
// TrackProducerLabel      - the producer of the recob::Track objects
// CosmicTagThresholds     - a vector of thresholds to apply to label as cosmic
//
// Created by Tracy Usher (usher@slac.stanford.edu) on September 18, 2014
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/PFParticle.h"
#include "AnalysisBase/CosmicTag.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/TimeService.h"
#include "Utilities/SimpleTimeService.h"

// Local functions.
namespace
{
    //----------------------------------------------------------------------------
    // Filter a collection of hits (set difference).
    // This function is copied from Track3DKalmanHit_module.cc
    //
    // Arguments:
    //
    // hits      - Hit collection from which hits should be removed.
    // used_hits - Hits to remove.
    //
    void FilterHits(art::PtrVector<recob::Hit>& hits, art::PtrVector<recob::Hit>& used_hits)
    {
        if(used_hits.size() > 0)
        {
            // Make sure both hit collections are sorted.
            std::stable_sort(hits.begin(), hits.end());
            std::stable_sort(used_hits.begin(), used_hits.end());

            // Do set difference operation.
            art::PtrVector<recob::Hit>::iterator it = std::set_difference(hits.begin(), hits.end(), used_hits.begin(), used_hits.end(), hits.begin());

            // Truncate hit collection.
            hits.erase(it, hits.end());
        }
    }
}

class Propagator;

class CRHitRemoval : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit CRHitRemoval(fhicl::ParameterSet const & pset);
    virtual ~CRHitRemoval();

    // Overrides.
    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

private:
    // Methods
    void removeTaggedHits(const recob::PFParticle*                            pfParticle,
                          const art::Handle<std::vector<recob::PFParticle> >& pfParticleHandle,
                          const art::FindManyP<recob::Cluster>&               partToClusAssns,
                          const art::FindManyP<recob::Hit>&                   clusToHitAssns,
                          std::set<const recob::PFParticle*>&                 taggedParticles,
                          art::PtrVector<recob::Hit>&                         hitVec);

    // Fcl parameters.
    std::vector<std::string> fCosmicProducerLabels;    ///< List of cosmic tagger producers
    std::string              fFlashProducerLabel;      ///< Name of the flash tagger module
    std::string              fHitProducerLabel;        ///< The full collection of hits
    std::string              fPFParticleProducerLabel; ///< PFParticle producer
    std::string              fTrackProducerLabel;      ///< Track producer
    
    bool                     fCorrelateToFlash;        ///< if true then correlate to flash
    std::vector<double>      fCosmicTagThresholds;     ///< Thresholds for tagging

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
    int fNumCRRejects;    ///< Number of tracks produced.
};

DEFINE_ART_MODULE(CRHitRemoval)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CRHitRemoval::CRHitRemoval(fhicl::ParameterSet const & pset) :
  fNumEvent(0),
  fNumCRRejects(0)
{
    reconfigure(pset);
    produces<std::vector<recob::Hit> >();

    // Report.
    mf::LogInfo("CRHitRemoval") << "CRHitRemoval configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
CRHitRemoval::~CRHitRemoval()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CRHitRemoval::reconfigure(fhicl::ParameterSet const & pset)
{
    fCosmicProducerLabels    = pset.get<std::vector<std::string> >("CosmicProducerLabels");
    fFlashProducerLabel      = pset.get<std::string>("FlashProducerLabel");
    fHitProducerLabel        = pset.get<std::string>("HitProducerLabel");
    fPFParticleProducerLabel = pset.get<std::string>("PFParticleProducerLabel");
    fTrackProducerLabel      = pset.get<std::string>("TrackProducerLabel");
    fCosmicTagThresholds     = pset.get<std::vector<double> >("CosmicTagThresholds");
}

//----------------------------------------------------------------------------
/// Begin job method.
void CRHitRemoval::beginJob()
{
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method. The goal is to produce a list of recob::Hit
/// objects which are a "clean" subset of all hits and which are believed to
/// be due to a neutrino interaction. It does this by considering input CosmicTag
/// objects, relating them to PFParticles/Tracks and removing the hits
/// associated to those objects which are believed to be Cosmic Rays.
///
void CRHitRemoval::produce(art::Event & evt)
{
    ++fNumEvent;
    
    // get the time service for the identification of the spill window
    util::SimpleTimeService const* time_service
      = &(*(art::ServiceHandle<util::TimeService>()));
    
    // Start by looking up the original hits
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitProducerLabel, hitHandle);
    
    // If there are no hits then there should be no output
    if (!hitHandle.isValid()) return;
    
    // If there are hits then we are going to output something so get a new
    // output hit vector
    std::unique_ptr<std::vector<recob::Hit> > outputHits(new std::vector<recob::Hit>);
    
    // And fill it with the complete original list of hits
    *outputHits = *hitHandle;
    
    // Now recover thre remaining collections of objects in the event store that we need
    // Start with tracks
    art::Handle<std::vector<recob::Track> > trackHandle;
    evt.getByLabel(fTrackProducerLabel, trackHandle);
    
    // If no tracks then no point in continuing here
    if (!trackHandle.isValid())
    {
        evt.put(std::move(outputHits));
        return;
    }
    
    // Recover the PFParticles that are responsible for making the tracks
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    evt.getByLabel(fPFParticleProducerLabel, pfParticleHandle);
    
    // Without a valid collection of PFParticles we can't do the hit removal
    if (!pfParticleHandle.isValid())
    {
        evt.put(std::move(outputHits));
        return;
    }
    
    // Recover the clusters so we can do associations to the hits
    // In theory the clusters come from the same producer as the PFParticles
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    evt.getByLabel(fPFParticleProducerLabel, clusterHandle);
    
    // If there are no clusters then something is really wrong
    if (!clusterHandle.isValid())
    {
        evt.put(std::move(outputHits));
        return;
    }
    
    // Build the list of cosmic tag producers, associations and thresholds
    std::vector<art::Handle<std::vector<anab::CosmicTag> > > cosmicHandleVec;
    std::vector<std::vector<const recob::Track*> >           cosmicToTrackVecVec;
    std::vector<double>                                      thresholdVec;
    
    for(size_t idx = 0; idx != fCosmicProducerLabels.size(); idx++)
    {
        std::string& handleLabel(fCosmicProducerLabels[idx]);
        
        art::Handle<std::vector<anab::CosmicTag> > cosmicHandle;
        evt.getByLabel(handleLabel, cosmicHandle);
        
        if (cosmicHandle.isValid())
        {
            // Look up the associations to tracks
            art::Handle< art::Assns<anab::CosmicTag, recob::Track> > cosmicToTrackHandle;
            evt.getByLabel(handleLabel, cosmicToTrackHandle);
            std::vector<const recob::Track*> cosmicToTrackVec = util::GetAssociatedVectorOneP(cosmicToTrackHandle, cosmicHandle);
            
            // Store these collections
            cosmicHandleVec.emplace_back(cosmicHandle);
            cosmicToTrackVecVec.emplace_back(cosmicToTrackVec);
            thresholdVec.push_back(fCosmicTagThresholds[idx]);
        }
    }
    
    // No cosmic tags then nothing to do here
    if (cosmicHandleVec.empty())
    {
        evt.put(std::move(outputHits));
        return;
    }
    
    // Let's do a test of the flash backing up the cosmic
    art::Handle< art::Assns<recob::Track, anab::CosmicTag> > trackToFlashHandle;
    std::vector<const anab::CosmicTag*>                      trackToFlashVec;

    evt.getByLabel(fFlashProducerLabel, trackToFlashHandle);
    if (trackToFlashHandle.isValid()) trackToFlashVec = util::GetAssociatedVectorOneP(trackToFlashHandle, trackHandle);
    
    // Now we walk up the remaining list of associations needed to go from CR tags to hits
    // We'll need to go from tracks to PFParticles
    art::Handle< art::Assns<recob::Track, recob::PFParticle> > trackToPfParticleHandle;
    evt.getByLabel(fTrackProducerLabel, trackToPfParticleHandle);
    std::vector<const recob::PFParticle*> trackToPartIdVec = util::GetAssociatedVectorOneP(trackToPfParticleHandle, trackHandle);
    
    // From PFParticles we go to clusters
    art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, evt, fPFParticleProducerLabel);
    
    // Likewise, recover the collection of associations to hits
    art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, evt, fPFParticleProducerLabel);
    
    // Container to contain the "bad" hits...
    art::PtrVector<recob::Hit> taggedHits;
    
    // No point double counting hits
    std::set<const recob::PFParticle*> taggedSet;

    // Start the identification of hits to remove. The outer loop is over the various producers of
    // the CosmicTag objects we're examininig
    for(size_t idx = 0; idx != cosmicHandleVec.size(); idx++)
    {
        // Obviously, dereference the handle and associations
        const art::Handle<std::vector<anab::CosmicTag> >& cosmicHandle(cosmicHandleVec[idx]);
        const std::vector<const recob::Track*>&           crTagTrackVec(cosmicToTrackVecVec[idx]);
        
        for(size_t crIdx = 0; crIdx != cosmicHandle->size(); crIdx++)
        {
            art::Ptr<anab::CosmicTag> cosmicTag(cosmicHandle, crIdx);
        
            // If this was tagged as a CR muon then we have work to do!
            if (cosmicTag->CosmicScore() > thresholdVec[idx])
            {
                // Recover the associated track
                const recob::Track* track(crTagTrackVec[crIdx]);
                
                // Probably a needless test but just to be sure
                if (!track) continue;
                
                // Dereference the PFParticle for what comes below
                const recob::PFParticle* pfParticle(trackToPartIdVec[track->ID()]);
                
                // Again, most likely needless
                if (!pfParticle) continue;
                
                // A cosmic ray must be a primary (by fiat)
                if (!pfParticle->IsPrimary()) continue;
                
                // Avoid double counting if more than one tagger running
                if (taggedSet.find(pfParticle) != taggedSet.end()) continue;
                
                // Remove all hits associated to this particle and its daughters
                removeTaggedHits(pfParticle, pfParticleHandle, clusterAssns, clusterHitAssns, taggedSet, taggedHits);
            }
            // If the cosmic score is non-zero then cross correlate with the Flash tagger
            else if (!trackToFlashVec.empty() && cosmicTag->CosmicScore() > 0.)
            {
                // Recover the associated track
                const recob::Track* track(crTagTrackVec[crIdx]);
                
                // Probably a needless test but just to be sure
                if (!track) continue;

                // Now recover the flash tag
                const anab::CosmicTag* flashTag(trackToFlashVec[track->ID()]);
                
                // Should be ok?
                if (!flashTag) continue;
                
                // Require over threshold
                if (flashTag->CosmicScore() > thresholdVec[idx])
                {
                    // Dereference the PFParticle for what comes below
                    const recob::PFParticle* pfParticle(trackToPartIdVec[track->ID()]);
                    
                    // Again, most likely needless
                    if (!pfParticle) continue;
                    
                    // A cosmic ray must be a primary (by fiat)
                    if (!pfParticle->IsPrimary()) continue;
                    
                    // Avoid double counting if more than one tagger running
                    if (taggedSet.find(pfParticle) != taggedSet.end()) continue;
                    
                    // Remove all hits associated to this particle and its daughters
                    removeTaggedHits(pfParticle, pfParticleHandle, clusterAssns, clusterHitAssns, taggedSet, taggedHits);
                }
            }
        }
    }
    
    // Are there any tagged hits?
    if (!taggedHits.empty())
    {
        // First order of business is to attempt to restore any hits which are shared between a tagged
        // CR PFParticle and an untagged one. We can do this by going through the PFParticles and
        // "removing" hits which are in the not tagged set.
        art::PtrVector<recob::Hit> untaggedHits;
        
        for(const auto& pfParticle : *pfParticleHandle)
        {
            if (taggedSet.find(&pfParticle) != taggedSet.end()) continue;
            
            // Recover the clusters associated to the input PFParticle
            std::vector<art::Ptr<recob::Cluster> > clusterVec = clusterAssns.at(pfParticle.Self());
            
            // Loop over the clusters and grab the associated hits
            for(const auto& cluster : clusterVec)
            {
                std::vector<art::Ptr<recob::Hit> > clusHitVec = clusterHitAssns.at(cluster->ID());
                untaggedHits.insert(untaggedHits.end(), clusHitVec.begin(), clusHitVec.end());
            }
        }
        
        // Filter out the hits we want to save
        FilterHits(taggedHits, untaggedHits);
        
        // The below is rather ugly but there is an interplay between art::Ptr's and the
        // actual pointers to objects that I might be missing and this is what I see how to do
        // First move all the original art::Ptr hits into a local art::PtrVector
        art::PtrVector<recob::Hit> originalHits;

        // Fill this one hit at a time...
        for(size_t hitIdx = 0; hitIdx != hitHandle->size(); hitIdx++)
        {
            art::Ptr<recob::Hit> hit(hitHandle, hitIdx);
            
            originalHits.push_back(hit);
        }
        
        // Remove the cosmic ray tagged hits
        FilterHits(originalHits, taggedHits);
        
        // Clear the current outputHits vector since we're going to refill...
        outputHits->clear();
        
        // Now make the new list of output hits
        for (const auto& hit : originalHits)
        {
            LOG_WARNING("CRHitRemoval")
              << "This module has experiment-specific information hard-coded!";
            // ... and yes, I want it printed this annoyingly
            
            // Kludge to remove out of time hits
            const double start_time_since_trigger
              = time_service->TPCTick2TrigTime(hit->PeakTimeMinusRMS());
            const double end_time_since_trigger
              = time_service->TPCTick2TrigTime(hit->PeakTimePlusRMS());
            
            if (end_time_since_trigger < 0.)
              continue;
            if (start_time_since_trigger > time_service->TPCClock().FramePeriod())
              continue;
            
            outputHits->emplace_back(*hit);
        }
    }
    
    // Add tracks and associations to event.
    evt.put(std::move(outputHits));
}

//----------------------------------------------------------------------------
/// Hit removal method
///
/// Arguments:
///
/// pfParticle - the top level PFParticle to have hits removed
/// pfParticleHandle - handle to the PFParticle objects
/// partToClusAssns - list of PFParticle to Cluster associations
/// clusToHitAssns - list of Cluster to Hit associations
/// hitVec - the current list of hits
///
/// This recursively called method will remove all hits associated to an input
/// PFParticle and, in addition, will call itself for all daughters of the input
/// PFParticle
///
void CRHitRemoval::removeTaggedHits(const recob::PFParticle*                            pfParticle,
                                    const art::Handle<std::vector<recob::PFParticle> >& pfParticleHandle,
                                    const art::FindManyP<recob::Cluster>&               partToClusAssns,
                                    const art::FindManyP<recob::Hit>&                   clusToHitAssns,
                                    std::set<const recob::PFParticle*>&                 taggedParticles,
                                    art::PtrVector<recob::Hit>&                         hitVec)
{
    // Recover the clusters associated to the input PFParticle
    std::vector<art::Ptr<recob::Cluster> > clusterVec = partToClusAssns.at(pfParticle->Self());
    
    // Record this PFParticle as tagged
    taggedParticles.insert(pfParticle);
    
    // Loop over the clusters and grab the associated hits
    for(const auto& cluster : clusterVec)
    {
        std::vector<art::Ptr<recob::Hit> > clusHitVec = clusToHitAssns.at(cluster->ID());
        hitVec.insert(hitVec.end(), clusHitVec.begin(), clusHitVec.end());
    }
    
    // Loop over the daughters of this particle and remove their hits as well
    for(const auto& daughterId : pfParticle->Daughters())
    {
        art::Ptr<recob::PFParticle> daughter(pfParticleHandle, daughterId);
        
        removeTaggedHits(daughter.get(), pfParticleHandle, partToClusAssns, clusToHitAssns, taggedParticles, hitVec);
    }
    
    return;
}


//----------------------------------------------------------------------------
/// End job method.
void CRHitRemoval::endJob()
{
    double aveCRPerEvent = fNumEvent > 0 ? double(fNumCRRejects) / double(fNumEvent) : 0.;
    
    mf::LogInfo("CRHitRemoval")
        << "CRHitRemoval statistics:\n"
        << "  Number of events = " << fNumEvent << "\n"
        << "  Number of Cosmic Rays found = " << fNumCRRejects
        << ", " << aveCRPerEvent << " average/event";
}
