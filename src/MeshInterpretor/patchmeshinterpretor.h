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


class PatchMeshInterpretor : public BasicMeshInterpretor
{
protected:

  virtual void Transformation(FemInterface::Matrix& T, int iq) const;
  virtual double compute_element_mean_matrix(int iq, EntryMatrix& E) const;
  int RhsPoint(Gascoigne::GlobalVector& f, const Vertex2d& p0, int comp, double d) const;
  int RhsPoint(Gascoigne::GlobalVector& f, const Vertex3d& p0, int comp, double d) const;
  virtual nmatrix<double> GetLocalInterpolationWeights(int iq) const { assert(0);}

  const PatchMesh* GetPatchMesh() const {
    const PatchMesh* MP = dynamic_cast<const PatchMesh*>(GetMesh());
    assert(MP);
    return MP;
  }

public:

//
///  Constructor 
//
  PatchMeshInterpretor() : BasicMeshInterpretor() {}
  std::string GetName() const {return "PatchMeshInterpretor";}

  void Structure(SparseStructureInterface* S) const;

  void Form(Gascoigne::GlobalVector& f, const Gascoigne::GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const Gascoigne::GlobalVector& u, const Equation& EQ, double) const;
  void MassMatrix(MatrixInterface& M) const;

  void ComputeError(const Gascoigne::GlobalVector& u, Gascoigne::LocalVector& err, const ExactSolution* ES) const;

  void Rhs(Gascoigne::GlobalVector& f, const RightHandSideData& RHS, double s) const;
  void DiracRhs(Gascoigne::GlobalVector& f, const RightHandSideData& RHS, double s) const;
  int RhsPoint(Gascoigne::GlobalVector& f, const std::vector<Vertex2d>& p0, int comp, const nvector<double>& d) const;
  int RhsPoint(Gascoigne::GlobalVector& f, const std::vector<Vertex3d>& p0, int comp, const nvector<double>& d) const;
  void RhsNeumann(Gascoigne::GlobalVector& f, const Equation& EQ, const Gascoigne::IntSet& Colors,  const NeumannData& NRHS, double s) const;

  double PressureFilter(nvector<double>&) const;

  // Functionals
  double ComputeBoundaryFunctional(const Gascoigne::GlobalVector& u, const BoundaryFunctional& BF) const;
  double ComputeDomainFunctional(const Gascoigne::GlobalVector& u, const DomainFunctional& F) const;
};

#endif
