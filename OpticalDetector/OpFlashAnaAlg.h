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
#include "TH1.h"

namespace opdet{

  class OpFlashAnaAlg{

  public:

    OpFlashAnaAlg(){}
    void SetOpFlashTree(TTree*);
    void SetFlashTimeHist(TH1F*);
    
    void FillOpFlashes(const std::vector<recob::OpFlash>&);
    
  private:

    TTree* fOpFlashTree;
    std::unique_ptr<recob::OpFlash> fOpFlashDataPtr;

    TH1F* fFlashTimeHist;
  };
  
}


#endif
