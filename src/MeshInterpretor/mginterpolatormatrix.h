#ifndef  __MgInterpolatorMatrix_h
#define  __MgInterpolatorMatrix_h

#include  "mginterpolatorinterface.h"
#include  "columnstencil.h"
#include  "gascoigne.h"

using namespace Gascoigne;


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

  void restrict_zero   (GlobalVector& uL, const GlobalVector& ul) const;
  void prolongate_add  (GlobalVector& ul, const GlobalVector& uL) const;

  void SolutionTransfer(GlobalVector& uL, const GlobalVector& ul) const;

};


#endif