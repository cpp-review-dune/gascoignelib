#ifndef __visucompvector_h
#define __visucompvector_h

#include  "compvector.h"
#include  "visudata.h"

/***************************************************************/

class VisuCompVector : public VisuData
{
 protected:

  typedef CompVector<double>       Vector;
  typedef Vector::iterator         pointer;
  typedef Vector::const_iterator   const_pointer;
  
  const Vector *uR, *zR;
  int           ncomp;

 public:
  
  VisuCompVector(const Vector& u) : uR(&u), zR(0)
    {
      ncomp = u.ncomp();
    }
  VisuCompVector(const Vector& u, const Vector& z) : uR(&u), zR(&z) 
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

/***************************************************************/

#endif
