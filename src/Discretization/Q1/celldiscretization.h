#ifndef  __CellDiscretization_h
#define  __CellDiscretization_h

#include  "basicdiscretization.h"
#include  "hierarchicalmesh.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments CellDiscretization
///
///
/////////////////////////////////////////////

class CellDiscretization : public BasicDiscretization
{
protected:

  FemInterface*         __FEM;
  IntegratorInterface*  __INT;

  const FemInterface* GetFem() const {assert(__FEM); return __FEM;}
  const IntegratorInterface* GetIntegrator() const {assert(__INT); return __INT;}
  IntegratorInterface*& GetIntegratorPointer() {return __INT;}
  FemInterface*& GetFemPointer() {return __FEM;}

  virtual void Transformation(FemInterface::Matrix& T, int iq) const;

  double ComputePointValue(const GlobalVector& u, const Vertex2d& p0,int comp) const;
  double ComputePointValue(const GlobalVector& u, const Vertex3d& p0,int comp) const; 

  virtual int GetCellNumber(const Vertex2d& p0, Vertex2d& p) const { assert(0); }
  virtual int GetCellNumber(const Vertex3d& p0, Vertex3d& p) const { assert(0); }

  /////

  void Transformation_HM(FemInterface::Matrix& T, const HierarchicalMesh* HM, int iq) const;
  void GlobalToLocal_HM(LocalVector& U, const GlobalVector& u, const HierarchicalMesh* HM, int iq) const;
  void swapIndices(IntVector& indices) const;
  
  virtual void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex2d& p0,int i,double s) const;
  virtual void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex3d& p0,int i,double s) const;


public:

  //
  ////  Constructor 
  //

  CellDiscretization() : BasicDiscretization(), __FEM(NULL), __INT(NULL) {}
  ~CellDiscretization(){
    if(__FEM) {delete __FEM; __FEM=NULL;}
    if(__INT) {delete __INT; __INT=NULL;}
  }

  std::string GetName() const {return "CellDiscretization";}

  void Structure(SparseStructureInterface* S) const;

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double) const;
  void BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
  void MassMatrix(MatrixInterface& M) const;

  void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;

  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
  void DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const;

  void BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, double s) const;

  void InitFilter(DoubleVector&) const;

  // Functionals
  double ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const;
  double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;

  double ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const;

  virtual nmatrix<double> GetLocalInterpolationWeights() const { assert(0); return nmatrix<double>();}

};
}

#endif
