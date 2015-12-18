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

#include "HitTagAssociatorAlg.h"
#include <limits>

cosmic::HitTagAssociatorAlg::HitTagAssociatorAlg(fhicl::ParameterSet const& p) 
{
  this->reconfigure(p);
}

void cosmic::HitTagAssociatorAlg::reconfigure(fhicl::ParameterSet const& p){}

void cosmic::HitTagAssociatorAlg::MakeHitTagAssociations(std::vector< std::vector<size_t> > const& bridges_per_hit,
							 std::vector< std::vector<size_t> > const& tags_per_bridges,
							 std::vector< std::vector<size_t> >& tags_per_hit)
{

  const size_t N_HITS = bridges_per_hit.size();
  tags_per_hit.clear();
  tags_per_hit.resize(N_HITS);

  for(size_t i_hit=0; i_hit<N_HITS; i_hit++){
    for(size_t i_bridge=0; i_bridge<bridges_per_hit[i_hit].size(); i_bridge++){
      tags_per_hit[i_hit].insert(tags_per_hit[i_hit].end(),
				 tags_per_bridges[i_bridge].begin(),
				 tags_per_bridges[i_bridge].end());  
    }
  }

}

void cosmic::HitTagAssociatorAlg::MakeHitTagAssociations(std::vector< std::vector<size_t> > const& bridges_per_hit,
							 std::vector<size_t> const& tag_per_bridge,
							 std::vector< std::vector<size_t> >& tags_per_hit)
{

  const size_t N_HITS = bridges_per_hit.size();

  tags_per_hit.clear();
  tags_per_hit.resize(N_HITS);

  for(size_t i_hit=0; i_hit<N_HITS; i_hit++){

    for(size_t i_bridge=0; i_bridge<bridges_per_hit[i_hit].size(); i_bridge++){

      if(i_bridge >= tag_per_bridge.size()) continue;

      if(tag_per_bridge[i_bridge]==std::numeric_limits<size_t>::max()) continue;

      tags_per_hit[i_hit].push_back(tag_per_bridge[i_bridge]); 
    }//end loop over bridges

  }//end loop over hits

}
