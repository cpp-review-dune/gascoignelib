#ifndef  __Q1_h
#define  __Q1_h


/////////////////////////////////////////////
////
////@brief
////  ... comments Q1

////
////
/////////////////////////////////////////////

#include  "cellmeshinterpretor.h"
#include  "hnstructureinterface.h"

namespace Gascoigne
{
class Q1 : public CellMeshInterpretor
{
protected:

  HNStructureInterface*    HN;

  IntVector GetLocalIndices(int iq) const {
    return GetMesh()->IndicesOfCell(iq);
  }
  void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;
  
  virtual HNStructureInterface* NewHNStructure()=0;

public:

//
////  Con(De)structor 
//

  Q1();
  ~Q1();

  int n() const {
    return GetMesh()->nnodes();
  }
  //  const HNStructureQ1* GetHNStructure() const { return HN;}

  void ReInit   (const MeshInterface* MP);
  void StrongDirichletMatrix       (MatrixInterface& A, int col, const std::vector<int>& comp) const;
  void StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const std::vector<int>& comp) const;
  void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const;
  void InterpolateSolution(GlobalVector& u, const GlobalVector& uold)const;
  virtual void InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const {}
  void HNAverage   (GlobalVector& x) const;
  void HNDistribute(GlobalVector& x) const;
  void HNZero      (GlobalVector& x) const;
  bool HNZeroCheck (const GlobalVector& x) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const;
  void MassMatrix(MatrixInterface& A) const;
  void Structure(SparseStructureInterface* SI) const;
  void InitFilter(DoubleVector& F) const;
};
}


#endif
