#ifndef FLASHHYPOTHESISCALCULATOR_CXX
#define FLASHHYPOTHESISCALCULATOR_CXX

#include "FlashHypothesisCalculator.h"

std::vector<double> opdet::FlashHypothesisCalculator::SegmentMidpoint(TVector3 const& pt1, TVector3 const& pt2, float XOffset)
{
  std::vector<double> xyz_segment(3);
  xyz_segment[0] = 0.5*(pt2.x()+pt1.x()) + XOffset;
  xyz_segment[1] = 0.5*(pt2.y()+pt1.y());
  xyz_segment[2] = 0.5*(pt2.z()+pt1.z());
  return xyz_segment;
}

void opdet::FlashHypothesisCalculator::FillFlashHypotheses(const float& yield,
							   const float& prompt_frac,
							   const float& dEdx,
							   const TVector3& pt1,
							   const TVector3& pt2,
							   const std::vector<float>& qe_vector,
							   const std::vector<float>& vis_vector,
							   FlashHypothesis& prompt_hyp,
							   FlashHypothesis& late_hyp)
{

  if(qe_vector.size()!=prompt_hyp.GetVectorSize() ||
     vis_vector.size()!=prompt_hyp.GetVectorSize() || 
     late_hyp.GetVectorSize()!=prompt_hyp.GetVectorSize())
    throw std::runtime_error("ERROR in FlashHypothesisCalculator: vector sizes not equal!");

  if(prompt_frac<0 || prompt_frac>1)
    throw std::runtime_error("ERROR in FlashHypothesisCalculator: prompt fraction not between zero and 1.");

  const float total_yield = yield*dEdx*(pt2-pt1).Mag();
  const float prompt_yield = total_yield*prompt_frac;
  const float late_yield = total_yield-prompt_yield;

  for(size_t i_chan=0; i_chan<prompt_hyp.GetVectorSize(); i_chan++){
    prompt_hyp.SetHypothesisAndError(i_chan,prompt_yield*vis_vector[i_chan]*qe_vector[i_chan]);
    late_hyp.SetHypothesisAndError(i_chan,late_yield*vis_vector[i_chan]*qe_vector[i_chan]);
  }

}

#endif
