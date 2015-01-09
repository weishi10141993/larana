/*!
 * Title:   FlashHypothesis Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Class that represents a flash hypothesis (PEs per opdet) and its error
*/

#include <iostream>
#include <numeric>


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

}
