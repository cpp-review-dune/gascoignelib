#ifndef  __MgInterpolatorInterface_h
#define  __MgInterpolatorInterface_h

#include  "gascoigne.h"

using namespace Gascoigne;

/*--------------------------------------------------------*/

class MgInterpolatorInterface
{
  
public:
  
  MgInterpolatorInterface() {}
  virtual ~MgInterpolatorInterface() {}

  virtual void restrict_zero   (GlobalVector&, const GlobalVector&) const=0;
  virtual void prolongate_add  (GlobalVector&, const GlobalVector&) const=0;
  virtual void SolutionTransfer(GlobalVector&, const GlobalVector&) const=0;
  virtual void SolutionTransferUp(GlobalVector& ul, const GlobalVector& uL) const {
    prolongate_add(ul,uL); 
  }
  virtual void Pi(GlobalVector& u) const {assert(0);}
};

#endif
