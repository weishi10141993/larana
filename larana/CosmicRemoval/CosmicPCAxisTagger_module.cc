////////////////////////////////////////////////////////////////////////
// Class:       CosmicPCAxisTagger
// Module Type: producer
// File:        CosmicPCAxisTagger_module.cc
//              This module checks timing and TPC volume boundaries as a
//              way to tag potential cosmic rays
//              This particular module uses PFParticles as input and handles
//              special cases associated with them. Instead of tracks, it uses
//              PCAxis objects for getting the start/end points of the
//              candidate CR's
//              This module started life as CosmicTrackTagger_module, written
//              by Sarah Lockwitz, and small alterations made to handle the
//              PFParticle input
//
// Generated at Wed Sep 17 19:17:00 2014 by Tracy Usher by cloning CosmicTrackTagger
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CosmicPCAxisTagger
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

#include "Geometry/Geometry.h"
#include "Geometry/geo.h"

#include "RecoBase/PFParticle.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/PCAxis.h"
#include "RecoBase/SpacePoint.h"
#include "RecoObjects/Cluster3D.h"

#include "AnalysisBase/CosmicTag.h"
#include "RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"

#include "TVector3.h"

namespace cosmic
{
    class CosmicPCAxisTagger;
    class SpacePoint;
    class Track;
}


class cosmic::CosmicPCAxisTagger : public art::EDProducer
{
public:
    explicit CosmicPCAxisTagger(fhicl::ParameterSet const & p);
    virtual ~CosmicPCAxisTagger();

    void produce(art::Event & e) override;

    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void endJob() override;

private:
    typedef std::vector<reco::ClusterHit2D> Hit2DVector;
    
    void RecobToClusterHits(const art::PtrVector<recob::Hit>& recobHitVec,
                            Hit2DVector&                      clusterHitVec,
                            reco::Hit2DListPtr&               hit2DListPtr) const;

    std::string                           fPFParticleModuleLabel;
    std::string                           fPCAxisModuleLabel;
    
    util::DetectorProperties*             fDetector;              ///<  Pointer to the detector properties
    lar_cluster3d::PrincipalComponentsAlg fPcaAlg;                ///<  Principal Components algorithm
    
    int                                   fDetectorWidthTicks;
    float                                 fTPCXBoundary, fTPCYBoundary, fTPCZBoundary;
    float                                 fDetHalfHeight, fDetWidth, fDetLength;
};


cosmic::CosmicPCAxisTagger::CosmicPCAxisTagger(fhicl::ParameterSet const & p) :
        fPcaAlg(p.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
// :
// Initialize member data here.
{
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    produces< std::vector<anab::CosmicTag> >();
    produces< art::Assns<recob::PFParticle, anab::CosmicTag> >();
    produces< art::Assns<recob::PCAxis,     anab::CosmicTag> >();
}

cosmic::CosmicPCAxisTagger::~CosmicPCAxisTagger()
{
    // Clean up dynamic memory and other resources here.
}

void cosmic::CosmicPCAxisTagger::produce(art::Event & evt)
{
    // Instatiate the output
    std::unique_ptr< std::vector< anab::CosmicTag > >                  cosmicTagPFParticleVector(  new std::vector<anab::CosmicTag> );
    std::unique_ptr< art::Assns<recob::PCAxis,     anab::CosmicTag > > assnOutCosmicTagPCAxis(     new art::Assns<recob::PCAxis,     anab::CosmicTag>);
    std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag > > assnOutCosmicTagPFParticle( new art::Assns<recob::PFParticle, anab::CosmicTag>);
    
    // Recover handle for PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    evt.getByLabel( fPFParticleModuleLabel, pfParticleHandle);
    
    if (!pfParticleHandle.isValid())
    {
        evt.put( std::move(cosmicTagPFParticleVector)  );
        evt.put( std::move(assnOutCosmicTagPFParticle) );
        evt.put( std::move(assnOutCosmicTagPCAxis)     );
        return;
    }
    
    // We need a handle to the collection of clusters in the data store so we can
    // handle associations to hits.
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    evt.getByLabel(fPFParticleModuleLabel, clusterHandle);
    
    // Recover the handle for the PCAxes
    art::Handle<std::vector<recob::PCAxis> > pcaxisHandle;
    evt.getByLabel( fPCAxisModuleLabel, pcaxisHandle);
    
    if (!pcaxisHandle.isValid() || !clusterHandle.isValid())
    {
        evt.put( std::move(cosmicTagPFParticleVector)  );
        evt.put( std::move(assnOutCosmicTagPFParticle) );
        evt.put( std::move(assnOutCosmicTagPCAxis)     );
        return;
    }
    
    // Recover the list of associated PCA axes
    art::FindManyP<recob::PCAxis> pfPartToPCAxisAssns(pfParticleHandle, evt, fPCAxisModuleLabel);
    
    // Add the relations to recover associations cluster hits
    art::FindManyP<recob::SpacePoint> spacePointAssnVec(pfParticleHandle, evt, fPFParticleModuleLabel);
    
    // Recover the collection of associations between PFParticles and clusters, this will
    // be the mechanism by which we actually deal with clusters
    art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, evt, fPFParticleModuleLabel);
    
    // Likewise, recover the collection of associations to hits
    art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, evt, fPFParticleModuleLabel);
    
    // The outer loop is going to be over PFParticles
    for(size_t pfPartIdx = 0; pfPartIdx != pfParticleHandle->size(); pfPartIdx++)
    {
        art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, pfPartIdx);
        
        // Require the PFParticle be a primary because attached hits will be daughter particles
        //if (!pfParticle->IsPrimary()) continue;
        // Oops! Don't do this unless traversing the heirarchy!
        
        // Recover the PCAxis vector
        std::vector<art::Ptr<recob::PCAxis> > pcAxisVec = pfPartToPCAxisAssns.at(pfPartIdx);
        
        // Is there an axis associated to this PFParticle?
        if (pcAxisVec.empty()) continue;

        // *****************************************************************************************
        // For what follows below we want the "best" PCAxis object only. However, it can be that
        // there are two PCAxes for a PFParticle (depending on source) where the "first" axis will
        // be the "better" one that we want (this statement by fiat, it is defined that way in the
        // axis producer module).
        if (pcAxisVec.size() > 1 && pcAxisVec.front()->getID() > pcAxisVec.back()->getID()) std::reverse(pcAxisVec.begin(), pcAxisVec.end());
        // We need to confirm this!!
        // *****************************************************************************************
        
        // Recover the axis
        const art::Ptr<recob::PCAxis>& pcAxis = pcAxisVec.front();
        
        // Start the tagging process...
        int                    isCosmic =  0;
        anab::CosmicTagID_t    tag_id   = anab::CosmicTagID_t::kNotTagged;
        
        // There are two sections to the tagging, in the first we are going to check for hits that are
        // "out of time" and for this we only need the hit vector. If no hits are out of time then
        // we need to do a more thorough check of the positions of the hits.
        // If we find hits that are out of time we'll set the "end points" of our trajectory to
        // a scale factor past the principle eigen value
        double eigenVal0 = sqrt(pcAxis->getEigenValues()[0]);
        double maxArcLen = 3. * eigenVal0;
        
        // Recover PCA end points
        TVector3 vertexPosition(pcAxis->getAvePosition()[0],pcAxis->getAvePosition()[1],pcAxis->getAvePosition()[2]);
        TVector3 vertexDirection(pcAxis->getEigenVectors()[0][0],pcAxis->getEigenVectors()[0][1],pcAxis->getEigenVectors()[0][2]);
        
        TVector3 pcAxisStart = vertexPosition - maxArcLen * vertexDirection;
        TVector3 pcAxisEnd   = vertexPosition + maxArcLen * vertexDirection;
        
        // "Track" end points in easily readable form
        float trackEndPt1_X = pcAxisStart [0];
        float trackEndPt1_Y = pcAxisStart [1];
        float trackEndPt1_Z = pcAxisStart [2];
        float trackEndPt2_X = pcAxisEnd[0];
        float trackEndPt2_Y = pcAxisEnd[1];
        float trackEndPt2_Z = pcAxisEnd[2];
        
        // Now we get the 2D clusters associated to this PFParticle
        std::vector<art::Ptr<recob::Cluster> > clusterVec = clusterAssns.at(pfParticle.key());
        
        bool dumpMe(false);
        
        // Once we have the clusters then we can loop over them to find the associated hits
        for(const auto& cluster : clusterVec)
        {
            // Recover the 2D hits associated to a given cluster
            std::vector<art::Ptr<recob::Hit> > hitVec = clusterHitAssns.at(cluster->ID());
            
            // Once we have the hits the first thing we should do is to check if any are "out of time"
            // If there are out of time hits then we are going to reject the cluster so no need to do
            // any further processing.
            /////////////////////////////////////
            // Check that all hits on particle are "in time"
            /////////////////////////////////////
            for(const auto& hit : hitVec)
            {
                if (dumpMe)
                {
                    std::cout << "***>> Hit key: " << hit.key() << ", peak - RMS: " << hit->PeakTimeMinusRMS() << ", peak + RMS: " << hit->PeakTimePlusRMS() << ", det width: " << fDetectorWidthTicks << std::endl;
                }
                //if( hit->StartTick() < fDetectorWidthTicks || hit->EndTick() > 2.*fDetectorWidthTicks)
                if( hit->PeakTimeMinusRMS() < fDetectorWidthTicks || hit->PeakTimePlusRMS() > 2.*fDetectorWidthTicks)
                {
                    isCosmic = 1;
                    tag_id   = anab::CosmicTagID_t::kOutsideDrift_Partial;
                    break;     // If one hit is out of time it must be a cosmic ray
                }
            }
        }
        
        // Recover the space points associated to this PFParticle.
        std::vector<art::Ptr<recob::SpacePoint> > spacePointVec = spacePointAssnVec.at(pfParticle.key());
        
        /////////////////////////////////
        // Now check the TPC boundaries:
        /////////////////////////////////
        if(isCosmic==0 && !spacePointVec.empty())
        {
            // Do a check on the transverse components of the PCA axes to make sure we are looking at long straight
            // tracks and not the kind of events we might want to keep
            double transRMS  = sqrt(std::pow(pcAxis->getEigenValues()[1],2) + std::pow(pcAxis->getEigenValues()[1],2));
            
            //if (eigenVal0 > 10. && transRMS < 10. && transRMS / eigenVal0 < 0.1)
            if (eigenVal0 > 0. && transRMS > 0.)
            {
                // The idea is to find the maximum extents of this PFParticle using the PCA axis which we
                // can then use to determine proximity to a TPC boundary.
                // We implement this by recovering the 3D Space Points and then make a pass through them to
                // find the space points at the extremes of the distance along the principle axis.
                // We'll loop through the space points looking for those which have the largest arc lengths along
                // the principle axis. Set up to do that
                double arcLengthToFirstHit(9999.);
                double arcLengthToLastHit(-9999.);
                
                for(const auto spacePoint : spacePointVec)
                {
                    TVector3 spacePointPos(spacePoint->XYZ()[0],spacePoint->XYZ()[1],spacePoint->XYZ()[2]);
                    TVector3 deltaPos    = spacePointPos - vertexPosition;
                    double   arcLenToHit = deltaPos.Dot(vertexDirection);
                    
                    if (arcLenToHit < arcLengthToFirstHit)
                    {
                        arcLengthToFirstHit = arcLenToHit;
                        pcAxisStart         = spacePointPos;
                    }
                    
                    if (arcLenToHit > arcLengthToLastHit)
                    {
                        arcLengthToLastHit = arcLenToHit;
                        pcAxisEnd          = spacePointPos;
                    }
                }
                
                // "Track" end points in easily readable form
                trackEndPt1_X = pcAxisStart [0];
                trackEndPt1_Y = pcAxisStart [1];
                trackEndPt1_Z = pcAxisStart [2];
                trackEndPt2_X = pcAxisEnd[0];
                trackEndPt2_Y = pcAxisEnd[1];
                trackEndPt2_Z = pcAxisEnd[2];
                
                // In below we check entry and exit points. Note that a special case of a particle entering
                // and exiting the same surface is considered to be running parallel to the surface and NOT
                // entering and exiting.
                // Also, in what follows we make no assumptions on which end point is the "start" or
                // "end" of the track being considered.
                int nBdX(0);
                int nBdY(0);
                int nBdZ(0);
            
                // Check x extents - note that uboone coordinaes system has x=0 at edge
                // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
                // Also note that, in theory, any cosmic ray entering or exiting the X surfaces will have presumably
                // been removed already by the checking of "out of time" hits... but this will at least label
                // neutrino interaction tracks which exit through the X surfaces of the TPC
                //if (fDetWidth - trackEndPt1_X < fTPCXBoundary || fDetWidth - trackEndPt2_X < fTPCXBoundary) nBdX++;
                //if (            trackEndPt1_X < fTPCXBoundary ||             trackEndPt2_X < fTPCXBoundary) nBdX++;
            
                // Check y extents (note coordinate system change)
                // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
                if (fDetHalfHeight - trackEndPt1_Y < fTPCYBoundary || fDetHalfHeight - trackEndPt2_Y < fTPCYBoundary) nBdY++;  // one end of track exits out top
                if (fDetHalfHeight + trackEndPt1_Y < fTPCYBoundary || fDetHalfHeight + trackEndPt2_Y < fTPCYBoundary) nBdY++;  // one end of track exist out bottom
            
                // Check z extents
                // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
                if (fDetLength - trackEndPt1_Z < fTPCZBoundary || fDetLength - trackEndPt2_Z < fTPCZBoundary ) nBdZ++;
                if (             trackEndPt1_Z < fTPCZBoundary ||              trackEndPt2_Z < fTPCZBoundary ) nBdZ++;
                
                // This should check for the case of a track which is both entering and exiting
                // but we consider entering and exiting the z boundaries to be a special case (should it be?)
                if( (nBdX + nBdY + nBdZ) > 1 && nBdZ != 2)
                {
                    isCosmic = 2;
                    if      (nBdX > 1)             tag_id = anab::CosmicTagID_t::kGeometry_XX;
                    else if (nBdY > 1)             tag_id = anab::CosmicTagID_t::kGeometry_YY;
                    else if (nBdX > 0 && nBdY > 0) tag_id = anab::CosmicTagID_t::kGeometry_XY;
                    else if (nBdX > 0 && nBdZ > 0) tag_id = anab::CosmicTagID_t::kGeometry_XZ;
                    else                           tag_id = anab::CosmicTagID_t::kGeometry_YZ;
                }
                // This is the special case of track which appears to enter/exit z boundaries
                else if ( nBdZ > 1)
                {
                    isCosmic = 3;
                    tag_id = anab::CosmicTagID_t::kGeometry_ZZ;
                }
                // This looks for track which enters/exits a boundary but has other endpoint in TPC
                else if ( (nBdX + nBdY + nBdZ) == 1)
                {
                    isCosmic = 4 ;
                    if      (nBdX == 1) tag_id = anab::CosmicTagID_t::kGeometry_X;
                    else if (nBdY == 1) tag_id = anab::CosmicTagID_t::kGeometry_Y;
                    else if (nBdZ == 1) tag_id = anab::CosmicTagID_t::kGeometry_Z;
                }
            }
        }
        
        std::vector<float> endPt1;
        std::vector<float> endPt2;
        endPt1.push_back( trackEndPt1_X );
        endPt1.push_back( trackEndPt1_Y );
        endPt1.push_back( trackEndPt1_Z );
        endPt2.push_back( trackEndPt2_X );
        endPt2.push_back( trackEndPt2_Y );
        endPt2.push_back( trackEndPt2_Z );
        
        float cosmicScore = isCosmic > 0 ? 1. : 0.;

        // Handle special cases
        if      (isCosmic == 3) cosmicScore = 0.4;   // Enter/Exit at opposite Z boundaries
        else if (isCosmic == 4) cosmicScore = 0.5;   // Enter or Exit but not both
        
        // Create the tag object for this PFParticle and make the corresponding association
        cosmicTagPFParticleVector->emplace_back( endPt1, endPt2, cosmicScore, tag_id);
        
        util::CreateAssn(*this, evt, *cosmicTagPFParticleVector, pfParticle, *assnOutCosmicTagPFParticle );

        // Loop through the tracks resulting from this PFParticle and mark them
        for(const auto& axis : pcAxisVec)
        {
            util::CreateAssn(*this, evt, *cosmicTagPFParticleVector, axis, *assnOutCosmicTagPCAxis );
        }
    }
    
    evt.put( std::move(cosmicTagPFParticleVector)  );
    evt.put( std::move(assnOutCosmicTagPFParticle) );
    evt.put( std::move(assnOutCosmicTagPCAxis)     );
    
    return;

} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////

void cosmic::CosmicPCAxisTagger::beginJob()
{
}

void cosmic::CosmicPCAxisTagger::reconfigure(fhicl::ParameterSet const & p)
{
    // Implementation of optional member function here.
  
    ////////  fSptalg  = new cosmic::SpacePointAlg(p.get<fhicl::ParameterSet>("SpacePointAlg"));
    art::ServiceHandle<util::DetectorProperties> detp;
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<geo::Geometry> geo;
    
    fDetector      = &*detp;

    fDetHalfHeight = geo->DetHalfHeight();
    fDetWidth      = 2.*geo->DetHalfWidth();
    fDetLength     = geo->DetLength();

    float fSamplingRate = detp->SamplingRate();

    fPFParticleModuleLabel = p.get<std::string >("PFParticleModuleLabel");
    fPCAxisModuleLabel     = p.get< std::string >("PCAxisModuleLabel");
    
    fPcaAlg.reconfigure(p.get<fhicl::ParameterSet>("PrincipalComponentsAlg"));

    fTPCXBoundary = p.get< float >("TPCXBoundary", 5);
    fTPCYBoundary = p.get< float >("TPCYBoundary", 5);
    fTPCZBoundary = p.get< float >("TPCZBoundary", 5);

    const double driftVelocity = larp->DriftVelocity( larp->Efield(), larp->Temperature() ); // cm/us

    //std::cerr << "Drift velocity is " << driftVelocity << " cm/us.  Sampling rate is: "<< fSamplingRate << " detector width: " <<  2*geo->DetHalfWidth() << std::endl;
    fDetectorWidthTicks = 2*geo->DetHalfWidth()/(driftVelocity*fSamplingRate/1000); // ~3200 for uB
}

void cosmic::CosmicPCAxisTagger::endJob() {
  // Implementation of optional member function here.
}

//------------------------------------------------------------------------------------------------------------------------------------------
void cosmic::CosmicPCAxisTagger::RecobToClusterHits(const art::PtrVector<recob::Hit>& recobHitVec,
                                                    Hit2DVector&                      clusterHitVec,
                                                    reco::Hit2DListPtr&               hit2DListPtr) const
{
    /**
     *  @brief Takes the input vector of recob::Hits and returns a vector of Cluster2D hits
     */

    
    // We'll need the offsets for each plane
    std::map<geo::View_t, double> viewOffsetMap;
    
    viewOffsetMap[geo::kU] = fDetector->GetXTicksOffset(geo::kU, 0, 0)-fDetector->TriggerOffset();
    viewOffsetMap[geo::kV] = fDetector->GetXTicksOffset(geo::kV, 0, 0)-fDetector->TriggerOffset();
    viewOffsetMap[geo::kW] = fDetector->GetXTicksOffset(geo::kW, 0, 0)-fDetector->TriggerOffset();
    
    // Reserve memory for the hit vector - this is the container of the actual hits
    clusterHitVec.reserve(recobHitVec.size());
    
    // Reserve memory for the list of pointers to these hits - what we'll operate on
    
    // Cycle through the recob hits to build ClusterHit2D objects and insert
    // them into the map
    for(const auto& recobHit : recobHitVec)
    {
        const geo::WireID& hitWireID(recobHit->WireID());
        
        double hitPeakTime(recobHit->PeakTime() - viewOffsetMap[recobHit->View()]);
        double xPosition(fDetector->ConvertTicksToX(recobHit->PeakTime(), hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat));
        
        clusterHitVec.emplace_back(reco::ClusterHit2D(0, 0., 0., xPosition, hitPeakTime, *recobHit));
        hit2DListPtr.emplace_back(&clusterHitVec.back());
    }
    
    return;
}

DEFINE_ART_MODULE(cosmic::CosmicPCAxisTagger)
