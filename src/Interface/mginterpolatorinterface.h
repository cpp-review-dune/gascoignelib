#ifndef  __MgInterpolatorInterface_h
#define  __MgInterpolatorInterface_h

#include  "gascoigne.h"

/*--------------------------------------------------------*/

class MgInterpolatorInterface
{
public:
  
  MgInterpolatorInterface() {}
  virtual ~MgInterpolatorInterface() {}

  virtual void restrict_zero   (Gascoigne::GlobalVector&, const Gascoigne::GlobalVector&) const=0;
  virtual void prolongate_add  (Gascoigne::GlobalVector&, const Gascoigne::GlobalVector&) const=0;
  virtual void SolutionTransfer(Gascoigne::GlobalVector&, const Gascoigne::GlobalVector&) const=0;
  virtual void SolutionTransferUp(Gascoigne::GlobalVector& ul, const Gascoigne::GlobalVector& uL) const {
    prolongate_add(ul,uL); 
  }
  virtual void Pi(Gascoigne::GlobalVector& u) const {assert(0);}
};

#endif
