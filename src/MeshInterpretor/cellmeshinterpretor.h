#ifndef  __CellMeshInterpretor_h
#define  __CellMeshInterpretor_h

#include  "basicmeshinterpretor.h"

#include  "hierarchicalmesh.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments CellMeshInterpretor
///
///
/////////////////////////////////////////////

namespace Gascoigne
{
class CellMeshInterpretor : public BasicMeshInterpretor
{
protected:

  FemInterface*         __FEM;
  IntegratorInterface*  __INT;

  const FemInterface* GetFem() const {assert(__FEM); return __FEM;}
  const IntegratorInterface* GetIntegrator() const {assert(__INT); return __INT;}
  IntegratorInterface*& GetIntegratorPointer() {return __INT;}
  FemInterface*& GetFemPointer() {return __FEM;}

  void Transformation(FemInterface::Matrix& T, int iq) const;

  int RhsPoint(GlobalVector& f, const Vertex2d& p0, int comp, double d) const;
  int RhsPoint(GlobalVector& f, const Vertex3d& p0, int comp, double d) const;

  double ComputePointValue(const GlobalVector& u, const Vertex2d& p0,int comp) const;
  double ComputePointValue(const GlobalVector& u, const Vertex3d& p0,int comp) const; 

  virtual nmatrix<double> GetLocalInterpolationWeights() const { assert(0); return nmatrix<double>();}
  int GetCellNumber(const Vertex2d& p0, nvector<Vertex2d>& p) const;
  void VertexTransformation(const nvector<Vertex2d>& p, 
			    const Vertex2d& p0, Vertex2d& tp) const;

  /////

  void Transformation_HM(FemInterface::Matrix& T, const HierarchicalMesh* HM, int iq) const;
  void GlobalToLocal_HM(LocalVector& U, const GlobalVector& u, const HierarchicalMesh* HM, int iq) const;
  void swapIndices(IntVector& indices) const;


public:

  //
  ////  Constructor 
  //

  CellMeshInterpretor() : BasicMeshInterpretor(), __FEM(NULL), __INT(NULL) {}
  ~CellMeshInterpretor(){
    if(__FEM) {delete __FEM; __FEM=NULL;}
    if(__INT) {delete __INT; __INT=NULL;}
  }

  std::string GetName() const {return "CellMeshInterpretor";}

  void Structure(SparseStructureInterface* S) const;

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double) const;
  void MassMatrix(MatrixInterface& M) const;

  void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;

  void Rhs(GlobalVector& f, const RightHandSideData& RHS, double s) const;
  void DiracRhs(GlobalVector& f, const RightHandSideData& RHS, double s) const;
  int RhsPoint(GlobalVector& f, const Functional* F) const;
  void RhsNeumann(GlobalVector& f, const IntSet& Colors,  const NeumannData& NRHS, double s) const;

  void InitFilter(DoubleVector&) const;

  // Functionals
  double ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const;
  double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;

  double ComputeNewPointFunctional(const GlobalVector& u, const NewPointFunctional* FP) const;

};
}

#endif
