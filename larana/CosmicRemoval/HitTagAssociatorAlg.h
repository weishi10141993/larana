#ifndef HITTAGASSOCIATIORALG_H
#define HITTAGASSOCIATIORALG_H
/*!
 * Title:   Hit <--> Cosmic Tag Associator Alg
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Algorithm that will provide associations of Hits to 
 *              cosmic tags, where both of those are associated to some
 *              intermediate object (like a track or cluster)
 * Input:       Assn<recob::Hit,???> and Assn<???,anab::CosmicTag>
 * Output:      Assn<recob::Hit,anab::CosmicTag> 
*/
#include <iostream>

#include "fhiclcpp/ParameterSet.h"

#include "lardata/RecoBase/Hit.h"
#include "lardata/AnalysisBase/CosmicTag.h"

namespace cosmic{
  class HitTagAssociatorAlg;
}

class cosmic::HitTagAssociatorAlg{
 public:
  HitTagAssociatorAlg(fhicl::ParameterSet const& p);
  void reconfigure(fhicl::ParameterSet const& p);

  //possiblity of multiple tags per bridge object
  void MakeHitTagAssociations(std::vector< std::vector<size_t> > const& bridges_per_hit,
			      std::vector< std::vector<size_t> > const& tags_per_bridges,
			      std::vector< std::vector<size_t> >& tags_per_hit);

  //exactly one tag per bridge object
  void MakeHitTagAssociations(std::vector< std::vector<size_t> > const& bridges_per_hit,
			      std::vector<size_t> const& tag_per_bridge,
			      std::vector< std::vector<size_t> >& tags_per_hit);

 private:

  //anything need to be private?

};

#endif
