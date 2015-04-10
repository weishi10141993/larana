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

#include "RecoBase/OpFlash.h"

#include "TTree.h"

namespace opdet{

  class OpFlashAnaAlg{

  public:

    OpFlashAnaAlg(){}
    void SetOutputTree(TTree*);
    
  private:

    TTree* fTree;

    //timing info
    double fTime;
    double fTimeWidth;
    double fAbsTime;
    unsigned int fFrame;

    //position info
    float fYCenter;
    float fYWidth;
    float fZCenter;
    float fZWidth;
    
    //PE info
    float fFastToTotal;
    
  };
  
}


#endif
