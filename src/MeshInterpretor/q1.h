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

class Q1 : public CellMeshInterpretor
{
protected:

  HNStructureInterface*    HN;

  nvector<int> GetLocalIndices(int iq) const;
  void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;
  
  virtual HNStructureInterface* NewHNStructure()=0;

public:

//
////  Con(De)structor 
//

  Q1();
  ~Q1();

  int n() const;
  //  const HNStructureQ1* GetHNStructure() const { return HN;}

  void ReInit   (const MeshInterface* MP);
  void StrongDirichletMatrix(MatrixInterface& A, int col, const std::vector<int>& comp) const;
  void StrongDirichletVectorZero(Gascoigne::GlobalVector& u, int col, const std::vector<int>& comp) const;
  void InterpolateSolution(Gascoigne::GlobalVector& u, const Gascoigne::GlobalVector& uold)const;
  void HNAverage   (Gascoigne::GlobalVector& x) const;
  void HNDistribute(Gascoigne::GlobalVector& x) const;
  void HNZero      (Gascoigne::GlobalVector& x) const;
  bool HNZeroCheck (const Gascoigne::GlobalVector& x) const;
  void Matrix(MatrixInterface& A, const Gascoigne::GlobalVector& u, const Equation& EQ, double d) const;
  void MassMatrix(MatrixInterface& A) const;
  void Structure(SparseStructureInterface* SI) const;
};


#endif
