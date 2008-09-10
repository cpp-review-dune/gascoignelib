#ifndef  __NonstationaryAlgorithm_h
#define  __NonstationaryAlgorithm_h

#include  "multilevelalgorithm.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class NonstationaryAlgorithm : public MultiLevelAlgorithm
{
 private:

  double  dt, time, theta;

 protected:

  void InitSolution(const std::string& initial, Gascoigne::VectorInterface& u) const;
  void TimeInfoBroadcast();

public:

  NonstationaryAlgorithm() :  MultiLevelAlgorithm() {}
  virtual ~NonstationaryAlgorithm() {}

  virtual void BasicInit(const Gascoigne::ParamFile* paramfile, 
			 const Gascoigne::NumericInterface* NI,
			 const Gascoigne::ProblemContainer* PC);

  void ImplicitEuler(const std::string& problemlabel);
  void CrankNicholson(const std::string& problemlabel);
  void ThetaScheme(const std::string& problemlabel);
};
}

/*-----------------------------------------*/

#endif
