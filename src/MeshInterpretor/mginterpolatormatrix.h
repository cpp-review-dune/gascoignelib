#ifndef  __MgInterpolatorMatrix_h
#define  __MgInterpolatorMatrix_h

#include  "mginterpolatorinterface.h"
#include  "columnstencil.h"
#include  "gascoigne.h"


/*-----------------------------------------*/


class MgInterpolatorMatrix : public virtual MgInterpolatorInterface
{
private:

  ColumnStencil  ST;
  nvector<double>      val;

public:


  MgInterpolatorMatrix() : MgInterpolatorInterface() {}

  ColumnStencil& GetStencil() {return  ST;}
  const ColumnStencil& GetStencil() const {return  ST;}

  nvector<double>& GetAlpha() {return val;}
  double Alpha(int pos) const {return val[pos];}

  void restrict_zero   (Gascoigne::GlobalVector& uL, const Gascoigne::GlobalVector& ul) const;
  void prolongate_add  (Gascoigne::GlobalVector& ul, const Gascoigne::GlobalVector& uL) const;
  void SolutionTransfer(Gascoigne::GlobalVector& uL, const Gascoigne::GlobalVector& ul) const;

};


#endif
