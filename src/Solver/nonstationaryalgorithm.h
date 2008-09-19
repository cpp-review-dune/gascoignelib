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
 protected:

  double  dt, time, theta;

  void InitSolution(const std::string& initial, Gascoigne::VectorInterface& u) const;
  void TimeInfoBroadcast();

public:

  NonstationaryAlgorithm() :  MultiLevelAlgorithm() {}
  virtual ~NonstationaryAlgorithm() {}

  virtual void BasicInit(const Gascoigne::ParamFile* paramfile, MultiLevelSolver* MLS,
			 const Gascoigne::NumericInterface* NI,
			 const Gascoigne::ProblemContainer* PC);

  void ImplicitEuler            (const std::string&);
  void CrankNicholson           (const std::string&);
  void ThetaScheme              (const std::string&);
  void FractionalStepThetaScheme(const std::string&);
};
}

/*-----------------------------------------*/

#endif
