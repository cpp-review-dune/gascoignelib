#ifndef __visucompvector_h
#define __visucompvector_h

#include  "compvector.h"
#include  "visudata.h"

/***************************************************************/

namespace Gascoigne
{
class VisuCompVector : public VisuData
{
 protected:

  typedef DoubleVector::iterator         pointer;
  typedef DoubleVector::const_iterator   const_pointer;
  
  const DoubleVector *uR, *zR;
  int           ncomp;

 public:
  
  VisuCompVector(const DoubleVector& u) : uR(&u), zR(0)
    {
      ncomp = u.ncomp();
    }
  VisuCompVector(const DoubleVector& u, const DoubleVector& z) : uR(&u), zR(&z) 
    {
      ncomp = u.ncomp()+z.ncomp();
    }

  int    visucomp()            const { return ncomp;}
  int    visun()               const { return uR->n();}
  double visudata(int i,int c) const 
    { 
      if ( c < uR->ncomp() ) return (*uR)(i,c);
      return (*zR)(i,c-uR->ncomp());
    }
};
}

/***************************************************************/

#endif
