#ifndef  __PatchDiscretization_h
#define  __PatchDiscretization_h

#include  "basicdiscretization.h"
#include  "patchmesh.h"
#include  "gascoignemesh.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments PatchDiscretization

///
///
/////////////////////////////////////////////

class PatchDiscretization : public BasicDiscretization
{
protected:

  FemInterface*         __FEM;
  IntegratorInterface*  __INT;

  const FemInterface* GetFem() const {return __FEM;}
  const IntegratorInterface* GetIntegrator() const {return __INT;}
  IntegratorInterface*& GetIntegratorPointer() {return __INT;}
  FemInterface*& GetFemPointer() {return __FEM;}

  virtual void Transformation(FemInterface::Matrix& T, int iq) const;
  virtual double compute_element_mean_matrix(int iq, EntryMatrix& E) const;
  virtual nmatrix<double> GetLocalInterpolationWeights(int iq) const { assert(0);}

  virtual int GetPatchNumber(const Vertex2d& p0, Vertex2d& p) const { assert(0); }
  virtual int GetPatchNumber(const Vertex3d& p0, Vertex3d& p) const { assert(0); }
  
  const PatchMesh* GetPatchMesh() const {
    const PatchMesh* MP = dynamic_cast<const PatchMesh*>(GetMesh());
    assert(MP);
    return MP;
  }

  const GascoigneMesh*  GetGascoigneMesh() const {
    const GascoigneMesh* MP = dynamic_cast<const GascoigneMesh*>(GetMesh());
    assert(MP);
    return MP;
  }
  
  double ComputePointValue(const GlobalVector& u, const Vertex2d& p0,int comp) const;
  double ComputePointValue(const GlobalVector& u, const Vertex3d& p0,int comp) const; 

public:

//
///  Constructor 
//
  PatchDiscretization() : BasicDiscretization(), __FEM(NULL), __INT(NULL) {}
  std::string GetName() const {return "PatchDiscretization";}

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
  void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex2d& p0,int i,double s) const;
  void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex3d& p0,int i,double s) const;
  
  void BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, double s) const;

  void InitFilter(nvector<double>&) const;

  // Functionals
  double ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const;
  double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;
  double ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const;
};
}

#endif
