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

  int RhsPoint(Gascoigne::GlobalVector& f, const Vertex2d& p0, int comp, double d) const;
  int RhsPoint(Gascoigne::GlobalVector& f, const Vertex3d& p0, int comp, double d) const;
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

  void Form(Gascoigne::GlobalVector& f, const Gascoigne::GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const Gascoigne::GlobalVector& u, const Equation& EQ, double) const;
  void MassMatrix(MatrixInterface& M) const;

  void ComputeError(const Gascoigne::GlobalVector& u, Gascoigne::LocalVector& err, const ExactSolution* ES) const;

  void Rhs(Gascoigne::GlobalVector& f, const RightHandSideData& RHS, double s) const;
  void DiracRhs(Gascoigne::GlobalVector& f, const RightHandSideData& RHS, double s) const;
  int RhsPoint(Gascoigne::GlobalVector& f, const Functional* F) const;
  void RhsNeumann(Gascoigne::GlobalVector& f, const Equation& EQ, const Gascoigne::IntSet& Colors,  const NeumannData& NRHS, double s) const;

  double PressureFilter(nvector<double>&) const;

  // Functionals
  double ComputeBoundaryFunctional(const Gascoigne::GlobalVector& u, const BoundaryFunctional& BF) const;
  double ComputeDomainFunctional(const Gascoigne::GlobalVector& u, const DomainFunctional& F) const;
};


#endif
