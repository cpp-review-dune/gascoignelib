#ifndef  __Q23d_h
#define  __Q23d_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Q23d

////
////
/////////////////////////////////////////////

#include  "q2.h"

namespace Gascoigne
{

class Q23d : public Q2
{
protected:

  nmatrix<double> GetLocalInterpolationWeights(int iq) const;

public:

//
////  Con(De)structor 
//

  Q23d();
  ~Q23d();

  std::string GetName() const {return "Q23d";}
  
  void BasicInit(const Gascoigne::ParamFile* paramfile);

  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* GMT);
};

}

#endif
