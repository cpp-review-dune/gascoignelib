#ifndef  __DwrAlgorithm_h
#define  __DwrAlgorithm_h

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

class DwrAlgorithm : public MultiLevelAlgorithm
{
 protected:

  DiscretizationInterface* CreateOtherDiscretization() const;

  void PrimalResidualsHigher(VectorInterface& f, const VectorInterface& u);
  void DualResidualsHigher  (VectorInterface& f, const VectorInterface& u, const VectorInterface& z);

public:

  DwrAlgorithm() :  MultiLevelAlgorithm()  {}
  virtual ~DwrAlgorithm() {}
  void AdaptiveLoop(const std::string& problemlabel, const std::string& duallabel, Functional& J);
};

}

/*-----------------------------------------*/

#endif
