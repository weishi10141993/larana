#ifndef FLASHHYPOTHESISALG_H
#define FLASHHYPOTHESISALG_H

/*!
 * Title:   FlashHypothesis Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Algorithm that produces a flash hypothesis for a trajectory.
 * Input:       Trajectory (std::vector<TVector3> objects)
 * Output:      FlashHypotheses
*/

#include <iostream>
#include <numeric>

#include "fhiclcpp/ParameterSet.h"

#include "RecoBase/Track.h"
#include "SimulationBase/MCParticle.h"

#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "Utilities/LArProperties.h"
#include "OpticalDetector/OpDigiProperties.h"

#include "TVector3.h"

namespace opdet{

class FlashHypothesis{

 public:
  FlashHypothesis(){}
  FlashHypothesis(size_t s):
   _NPEs_Vector(std::vector<float>(s,0.0)),_NPEs_ErrorVector(std::vector<float>(s,0.0)){}
  FlashHypothesis(std::vector<float> const& vector,
		 std::vector<float> const& vector_error):
   _NPEs_Vector(vector),_NPEs_ErrorVector(vector_error)
  {
    if( vector.size()!=vector_error.size())
      throw "ERROR in FlashHypothesisConstructor: Vector sizes not equal";
  }

  std::vector<float> const& GetHypothesisVector() const { return _NPEs_Vector; }
  std::vector<float> const& GetHypothesisErrorVector() const { return _NPEs_ErrorVector; }
  void SetHypothesisVector( std::vector<float> v ) { _NPEs_Vector=v; }
  void SetHypothesisErrorVector( std::vector<float> v ) { _NPEs_ErrorVector = v; }

  float const& GetHypothesis(size_t i_opdet) const { return _NPEs_Vector.at(i_opdet); }
  float const& GetHypothesisError(size_t i_opdet) const { return _NPEs_ErrorVector.at(i_opdet); }
  void SetHypothesis( size_t i_opdet, float pe ) { _NPEs_Vector.at(i_opdet)=pe; }
  void SetHypothesisError( size_t i_opdet, float err ) { _NPEs_ErrorVector.at(i_opdet) = err; }

  void SetHypothesisAndError( size_t i_opdet, float pe )
  {  SetHypothesis(i_opdet,pe); SetHypothesisError(i_opdet,std::sqrt(pe)); }

  float GetTotalPEs() const
  { return std::accumulate(_NPEs_Vector.begin(),_NPEs_Vector.end(),0.0); }
  float GetTotalPEsError() const
  { return std::sqrt( std::inner_product(_NPEs_ErrorVector.begin(),_NPEs_Vector.end(),_NPEs_Vector.begin(),0.0) ); }

  size_t GetVectorSize() const { return _NPEs_Vector.size(); }
  
  void Normalize(float const& totalPE_target);

  FlashHypothesis operator+(const FlashHypothesis& fh){

    if( _NPEs_Vector.size() != fh.GetVectorSize() )
      throw "ERROR in FlashHypothesisAddition: Cannot add hypothesis of different size";
    
    FlashHypothesis flashhyp(_NPEs_Vector.size());
    for(size_t i=0; i<_NPEs_Vector.size(); i++){
      flashhyp._NPEs_Vector[i] = _NPEs_Vector[i] + fh._NPEs_Vector[i];
      flashhyp._NPEs_ErrorVector[i] =
	std::sqrt(this->_NPEs_ErrorVector[i]*this->_NPEs_ErrorVector[i] +
		  fh._NPEs_ErrorVector[i]*fh._NPEs_ErrorVector[i]);
    }
    return flashhyp;
  }

 private:
  std::vector<float> _NPEs_Vector;
  std::vector<float> _NPEs_ErrorVector;
};

class FlashHypothesisAlg{

 public:
  FlashHypothesisAlg() {}

  void CreateFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
					float dEdx,
					geo::Geometry const& geom,
					phot::PhotonVisibilityService const& pvs,
					util::LArProperties const& larp,
					opdet::OpDigiProperties const& opdigip,
					float XOffset,
					FlashHypothesis &prompt_flash_hyp,
					FlashHypothesis &late_flash_hyp);
  
  void AddFlashHypothesesFromSegment(TVector3 const& pt1, TVector3 const& pt2, 
				     float dEdx,
				     geo::Geometry const& geom,
				     phot::PhotonVisibilityService const& pvs,
				     util::LArProperties const& larp,
				     opdet::OpDigiProperties const& opdigip,
				     float XOffset,
				     FlashHypothesis &prompt_flash_hyp,
				     FlashHypothesis &late_flash_hyp);
  
  private:
    
};

 class FlashHypothesisCollection{
  
 public:
  
  FlashHypothesisCollection(recob::Track const& track, 
			    std::vector<float> dEdxVector,
			    geo::Geometry const& geom,
			    phot::PhotonVisibilityService const& pvs,
			    util::LArProperties const& larp,
			    opdet::OpDigiProperties const& opdigip,
			    float XOffset=0);
  
  FlashHypothesisCollection(std::vector<TVector3> const& trajVector, 
			    std::vector<float> dEdxVector,
			    geo::Geometry const& geom,
			    phot::PhotonVisibilityService const& pvs,
			    util::LArProperties const& larp,
			    opdet::OpDigiProperties const& opdigip,
			    float XOffset=0);
  
  FlashHypothesisCollection(TVector3 const& pt1, TVector3 const& pt2, 
			    float dEdx,
			    geo::Geometry const& geom,
			    phot::PhotonVisibilityService const& pvs,
			    util::LArProperties const& larp,
			    opdet::OpDigiProperties const& opdigip,
			    float XOffset=0);
  
 FlashHypothesisCollection(FlashHypothesis const& prompt,
			   FlashHypothesis const& late):
  _prompt_hypothesis(prompt), _late_hypothesis(late)
  {
    _total_hypothesis = _prompt_hypothesis + _late_hypothesis;
  }
  
  void Initialize(size_t s)
  {
    _prompt_hypothesis = FlashHypothesis(s);
    _late_hypothesis = FlashHypothesis(s);
    _total_hypothesis = _prompt_hypothesis + _late_hypothesis;
  }

  void Normalize(float const& totalPE, util::LArProperties const& larp);
  
  FlashHypothesis const& GetPromptHypothesis() { return _prompt_hypothesis; }
  FlashHypothesis const& GetLateHypothesis() { return _late_hypothesis; }
  FlashHypothesis const& GetTotalHypothesis() { return _total_hypothesis; }

  float GetPromptPEs() const { return _prompt_hypothesis.GetTotalPEs(); }
  float GetPromptPEsError() const { return _prompt_hypothesis.GetTotalPEsError(); }
  float GetLatePEs() const { return _late_hypothesis.GetTotalPEs(); }
  float GetLatePEsError() const { return _late_hypothesis.GetTotalPEsError(); }
  float GetTotalPEs() const { return _total_hypothesis.GetTotalPEs(); }
  float GetTotalPEsError() const { return _total_hypothesis.GetTotalPEsError(); }
  
 private:
  FlashHypothesis _prompt_hypothesis;
  FlashHypothesis _late_hypothesis;
  FlashHypothesis _total_hypothesis;

  FlashHypothesisAlg _alg;
};


}

#endif
