#ifndef FLASHHYPOTHESIS_H
#define FLASHHYPOTHESIS_H

/*!
 * Title:   FlashHypothesis Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Class that represents a flash hypothesis (PEs per opdet) and its error
*/

#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace opdet{

class FlashHypothesis{

 public:
  FlashHypothesis(){}
  FlashHypothesis(size_t s):
   _NPEs_Vector(std::vector<float>(s,0.0)),_NPEs_ErrorVector(std::vector<float>(s,0.0)){}
  FlashHypothesis(std::vector<float> const& vector,
		  std::vector<float> const& vector_error=std::vector<float>())
  {
    SetHypothesisVectorAndErrorVector(vector,vector_error);
  }

  std::vector<float> const& GetHypothesisVector() const { return _NPEs_Vector; }
  std::vector<float> const& GetHypothesisErrorVector() const { return _NPEs_ErrorVector; }
  void SetHypothesisVector( std::vector<float> v ) { _NPEs_Vector=v; _NPEs_ErrorVector.resize(v.size()); }
  void SetHypothesisErrorVector( std::vector<float> v ) { _NPEs_ErrorVector = v; _NPEs_Vector.resize(v.size()); }
  void SetHypothesisVectorAndErrorVector( std::vector<float> v , std::vector<float> err=std::vector<float>(0)) 
  { 
    if(err.size()!=0 && err.size()!=v.size())
      throw std::runtime_error("ERROR in FlashHypothesisVectorSetter: Vector sizes not equal");

    _NPEs_Vector = v;
    if(err.size()==0){
      _NPEs_Vector = v; 
      _NPEs_ErrorVector.resize(v.size()); 
      for(size_t i=0; i<_NPEs_Vector.size(); i++)
	_NPEs_ErrorVector[i] = std::sqrt(_NPEs_Vector[i]);
    }
    else
      _NPEs_ErrorVector = err;

  }

  float const& GetHypothesis(size_t i_opdet) const { return _NPEs_Vector.at(i_opdet); }
  float const& GetHypothesisError(size_t i_opdet) const { return _NPEs_ErrorVector.at(i_opdet); }
  void SetHypothesis( size_t i_opdet, float pe ) { _NPEs_Vector.at(i_opdet)=pe; }
  void SetHypothesisError( size_t i_opdet, float err ) { _NPEs_ErrorVector.at(i_opdet) = err; }

  void SetHypothesisAndError( size_t i_opdet, float pe , float err=-999 )
  {  
    SetHypothesis(i_opdet,pe); 
    if(err>0) SetHypothesisError(i_opdet,err); 
    else SetHypothesisError(i_opdet,std::sqrt(pe)); 
  }

  float GetTotalPEs() const
  { return std::accumulate(_NPEs_Vector.begin(),_NPEs_Vector.end(),0.0); }
  float GetTotalPEsError() const
  { return std::sqrt( std::inner_product(_NPEs_ErrorVector.begin(),_NPEs_ErrorVector.end(),_NPEs_ErrorVector.begin(),0.0) ); }

  size_t GetVectorSize() const { return _NPEs_Vector.size(); }
  
  void Normalize(float const& totalPE_target){

    if( GetTotalPEs() < std::numeric_limits<float>::epsilon() ) return;
    
    const float PE_ratio = totalPE_target/GetTotalPEs();
    for(size_t i_opdet=0; i_opdet<_NPEs_Vector.size(); i_opdet++){
      _NPEs_Vector[i_opdet] *= PE_ratio;
      _NPEs_ErrorVector[i_opdet] = std::sqrt(_NPEs_Vector[i_opdet]);
    }
    
  }

  FlashHypothesis operator+(const FlashHypothesis& fh){

    if( _NPEs_Vector.size() != fh.GetVectorSize() )
      throw std::runtime_error("ERROR in FlashHypothesisAddition: Cannot add hypothesis of different size");
    
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

}

#endif
