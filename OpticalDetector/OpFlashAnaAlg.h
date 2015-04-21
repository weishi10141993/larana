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

#include "TTree.h"
#include "TH1.h"

#include "RecoBase/OpHit.h"
#include "RecoBase/OpFlash.h"

namespace opdet{

  class OpFlashAnaAlg{

  public:

    OpFlashAnaAlg(){ fMakeOpHitHist=false; }
    void SetOpFlashTree(TTree*,bool makeOpHitHist=true);
    void SetOpHitTree(TTree*);
    
    void FillOpFlashes(const std::vector<recob::OpFlash>&);
    void FillOpHits(const std::vector<recob::OpHit>&);
    
  private:

    std::unique_ptr<recob::OpFlash> fOpFlashDataPtr;
    std::unique_ptr<recob::OpHit> fOpHitDataPtr;

    bool fMakeOpHitHist;
    std::unique_ptr<TH1D> fOpFlashHitHistPtr;

    TTree* fOpFlashTree;
    void   FillOpFlashTree(const std::vector<recob::OpFlash>&);
    
    TTree* fOpHitTree;
    void   FillOpHitTree(const std::vector<recob::OpHit>&);

  };
  
}


#endif
