#ifndef  __MgInterpolatorInterface_h
#define  __MgInterpolatorInterface_h

#include  "compvector.h"
#include  "meshtransferinterface.h"

/*--------------------------------------------------------*/

class MgInterpolatorInterface
{
protected:
  
  typedef  CompVector<double>   Vector;
  
public:
  
  virtual ~MgInterpolatorInterface() {}
  virtual void restrict_zero   (Vector&, const Vector&) const=0;
  virtual void prolongate_add  (Vector&, const Vector&) const=0;
  virtual void SolutionTransfer(Vector&, const Vector&) const=0;

  virtual void Pi    (Vector& u) const {assert(0);}
};

#endif
