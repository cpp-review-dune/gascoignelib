#ifndef  __MgInterpolator_h
#define  __MgInterpolator_h


/////////////////////////////////////////////
////
////@brief
////  ... comments MgInterpolator

////
////
/////////////////////////////////////////////

#include  "gascoigne.h"
#include  "columnstencil.h"
#include  "meshtransferinterface.h"
#include  "mginterpolatorinterface.h"

using namespace Gascoigne;

class MgInterpolator : public MgInterpolatorInterface
{
protected:

  ColumnStencil  ST;
  nvector<double>      val;

public:

//
////  Con(De)structor 
//

  MgInterpolator() {}
  virtual ~MgInterpolator() {}

  virtual void restrict_zero   (GlobalVector& uL, const GlobalVector& ul) const;
  virtual void prolongate_add  (GlobalVector& ul, const GlobalVector& uL) const;
  virtual void SolutionTransfer(GlobalVector& uL, const GlobalVector& ul) const;
  virtual void SolutionTransferUp(GlobalVector& ul, const GlobalVector& uL) const {
    prolongate_add(ul,uL); 
  }

  ColumnStencil& GetStencil() {return  ST;}
  const ColumnStencil& GetStencil() const {return  ST;}

  nvector<double>& GetAlpha() {return val;}

  double Alpha(int pos) const {return val[pos];}

  void init(const MeshTransferInterface* MT) {assert(0);}

};


#endif
