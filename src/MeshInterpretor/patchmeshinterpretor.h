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
  int RhsPoint(GlobalVector& f, const Vertex2d& p0, int comp, double d) const;
  int RhsPoint(GlobalVector& f, const Vertex3d& p0, int comp, double d) const;
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

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double) const;
  void MassMatrix(MatrixInterface& M) const;

  void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;

  void Rhs(GlobalVector& f, const RightHandSideData& RHS, double s) const;
  void DiracRhs(GlobalVector& f, const RightHandSideData& RHS, double s) const;
  int RhsPoint(GlobalVector& f, const std::vector<Vertex2d>& p0, int comp, const nvector<double>& d) const;
  int RhsPoint(GlobalVector& f, const std::vector<Vertex3d>& p0, int comp, const nvector<double>& d) const;
  void RhsNeumann(GlobalVector& f, const Equation& EQ, const IntSet& Colors,  const NeumannData& NRHS, double s) const;

  double PressureFilter(nvector<double>&) const;

  // Functionals
  double ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const;
  double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;
};

#endif
