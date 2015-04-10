#ifndef OPFLASHANAALG_H
#define OPFLASHANAALG_H

/*!
 * Title:   OpFlashAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Alg that stores flash information in a handy tree.
 * 
*/

#include <vector>
#include <memory>

#include "RecoBase/OpFlash.h"

#include "TTree.h"

namespace opdet{

  class OpFlashAnaAlg{

  public:

    OpFlashAnaAlg(){}
    void SetOpFlashTree(TTree*);

    void FillOpFlashes(const std::vector<recob::OpFlash>&);
    
  private:

    TTree* fOpFlashTree;
    std::unique_ptr<recob::OpFlash> fOpFlashDataPtr;

  };
  
}


#endif
