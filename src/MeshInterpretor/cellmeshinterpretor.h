#ifndef  __CellMeshInterpretor_h
#define  __CellMeshInterpretor_h

#include  "basicmeshinterpretor.h"

/////////////////////////////////////////////
///
///@brief
///  ... comments CellMeshInterpretor
///
///
/////////////////////////////////////////////

class CellMeshInterpretor : public BasicMeshInterpretor
{
protected:

  void Transformation(FemInterface::Matrix& T, int iq) const;

  virtual double compute_element_mean_matrix(int iq, EntryMatrix& E) const;
  int RhsPoint(GlobalVector& f, const Vertex2d& p0, int comp, double d) const;
  int RhsPoint(GlobalVector& f, const Vertex3d& p0, int comp, double d) const;
  virtual nmatrix<double> GetLocalInterpolationWeights() const { assert(0);}
  int GetCellNumber(const Vertex2d& p0, nvector<Vertex2d>& p) const;
  void VertexTransformation(const nvector<Vertex2d>& p, 
			    const Vertex2d& p0, Vertex2d& tp) const;

public:

  //
  ////  Constructor 
  //

  CellMeshInterpretor() : BasicMeshInterpretor() {}
  ~CellMeshInterpretor() {}

  std::string GetName() const {return "CellMeshInterpretor";}

  void Structure(SparseStructureInterface* S) const;

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double) const;
  void MassMatrix(MatrixInterface& M) const;

  void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;

  void Rhs(GlobalVector& f, const RightHandSideData& RHS, double s) const;
  void DiracRhs(GlobalVector& f, const RightHandSideData& RHS, double s) const;
  int RhsPoint(GlobalVector& f, const Functional* F) const;
  void RhsNeumann(GlobalVector& f, const Equation& EQ, const IntSet& Colors,  const NeumannData& NRHS, double s) const;

  double PressureFilter(nvector<double>&) const;

  // Functionals
  double ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const;
  double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;
};


#endif
