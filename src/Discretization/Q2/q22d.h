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

  int GetPatchNumber(const Vertex2d& p0, Vertex2d& p) const;
  void VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq) const;

public:

//
////  Con(De)structor 
//

  Q22d();
  ~Q22d();

  std::string GetName() const {return "Q22d";}
  
  void BasicInit(const ParamFile* paramfile);

  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* GMT);

  nmatrix<double> GetLocalInterpolationWeights(int iq) const;
};
}

#endif
