#ifndef  __Q22d_h
#define  __Q22d_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Q22d

////
////
/////////////////////////////////////////////

#include  "q2.h"

namespace Gascoigne
{

class Q22d : public Q2
{
protected:

  nmatrix<double> GetLocalInterpolationWeights(int iq) const;

public:

//
////  Con(De)structor 
//

  Q22d();
  ~Q22d();

  std::string GetName() const {return "Q22d";}
  
  void BasicInit(const Gascoigne::ParamFile* paramfile);

  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* GMT);
};

}
#endif
