#ifndef FLASHHYPOTHESISCALCULATOR_H
#define FLASHHYPOTHESISCALCULATOR_H

/*!
 * Title:   FlashHypothesis Calculator Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Simple class for calculating flash hypotheses
*/

#include<vector>
#include<iostream>

#include "FlashHypothesis.h"
#include "TVector3.h"

namespace opdet{

  class FlashHypothesisCalculator{

  public:

    FlashHypothesisCalculator(){}

    std::vector<double> SegmentMidpoint(const TVector3 &pt1, const TVector3 &pt2, float XOffset=0);
    void FillFlashHypothesis(const float& yield,
			     const float& dEdx,
			     const TVector3& pt1,
			     const TVector3& pt2,
			     const std::vector<float>& qe_vector,
			     const std::vector<float>& vis_vector,
			     FlashHypothesis& hyp);

  };

}

#endif
