#ifndef  __PatchMeshInterpretor_h
#define  __PatchMeshInterpretor_h



/////////////////////////////////////////////
///
///@brief
///  ... comments PatchMeshInterpretor

///
///
/////////////////////////////////////////////


#include  "basicmeshinterpretor.h"
#include  "patchmesh.h"


namespace Gascoigne
{
class PatchMeshInterpretor : public BasicMeshInterpretor
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

public:

//
///  Constructor 
//
  PatchMeshInterpretor() : BasicMeshInterpretor(), __FEM(NULL), __INT(NULL) {}
  std::string GetName() const {return "PatchMeshInterpretor";}

  void Structure(SparseStructureInterface* S) const;

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double) const;
  void MassMatrix(MatrixInterface& M) const;

  void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;

  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
  void DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const;
  void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex2d& p0,int i,double s) const;
  void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex3d& p0,int i,double s) const;
  
  void RhsNeumann(GlobalVector& f, const IntSet& Colors,  const NeumannData& NRHS, double s) const;

  void InitFilter(nvector<double>&) const;

  // Functionals
  double ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const;
  double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;
};
}

#endif
