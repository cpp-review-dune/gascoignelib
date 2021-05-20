/*----------------------------   dg.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dg_H
#define __dg_H
/*----------------------------   dg.h     ---------------------------*/

#include "dgdofhandler.h"
#include "dgequation.h"
#include "dgintegrator.h"
#include "discretizationinterface.h"
#include "finiteelement.h"
#include "gascoignemesh.h"
#include "nvector.h"
#include "transformation2d.h"

namespace Gascoigne {
template<class BASE>
class DG : public DiscretizationInterface
{
private:
  DGDofHandler<BASE> _dofhandler;
  const GascoigneMesh* _mesh;

  Transformation2d<BASEQ12D> _trans;
  FiniteElement<2, 1, Transformation2d<BASEQ12D>, BASE> _fe;
  FiniteElement<2, 1, Transformation2d<BASEQ12D>, BASE> _feslave;

  DGIntegrator<2> _integrator;

  mutable DataContainer _datacontainer;

  mutable LocalData __QN_master, __QN_slave;

protected:
public:
  // Init
  DG();
  void BasicInit(const ParamFile* pf);
  void ReInit(const GascoigneMesh* M);

  void GlobalToLocalData(LocalData& QN, const DataContainer& DC, int iq) const;

  // Handling of Vectors & Data
  const DataContainer& GetDataContainer() const { abort(); }
  void SetDataContainer(const DataContainer& q) { abort(); }
  void AddNodeVector(const std::string& name, const GlobalVector* q) const
  {
    _datacontainer.AddNodeVector(name, q);
  };
  void DeleteNodeVector(const std::string& name) const
  {
    _datacontainer.DeleteNodeVector(name);
  };
  void AddCellVector(const std::string& name, const GlobalVector* q) const
  {
    assert(0);
  }
  void DeleteCellVector(const std::string& name) const { assert(0); }
  void AddParameterVector(const std::string& name,
                          const GlobalParameterVector* q) const
  {
    assert(0);
  }
  void DeleteParameterVector(const std::string& name) const { assert(0); }
  void GlobalToGlobalData(LocalParameterData& QP) const;

  // Info
  std::string GetName() const { return "DG"; }
  int ndofs() const // number of unknowns
  {
    return _dofhandler.ndofs();
  };
  int nelements() const { return _dofhandler.nelements(); };
  int n_withouthanging() const
  {
    std::cerr << "DG::n_withouthanging() not implemented" << std::endl;
    abort();
  }

  // Dof-Handling

  void Transformation(FemInterface::Matrix& T, int iq) const;

  void Structure(SparseStructureInterface* S) const
  {
    _dofhandler.Structure(S);
  }

  // Integration functions
  void EdgeForm(GlobalVector& f,
                const GlobalVector& u,
                const DGEquation& EQ,
                double d) const;
  void EdgeRhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
  void EdgeMatrix(MatrixInterface& A,
                  const GlobalVector& u,
                  const DGEquation& EQ,
                  double) const;

  void Form(GlobalVector& f,
            const GlobalVector& u,
            const Equation& EQ,
            double d) const;
  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
  void Matrix(MatrixInterface& A,
              const GlobalVector& u,
              const Equation& EQ,
              double) const;

  void AdjointForm(GlobalVector& f,
                   const GlobalVector& u,
                   const Equation& EQ,
                   double d) const
  {
    std::cerr << "\"DiscretizationInterface::AdjointForm\" not written!"
              << std::endl;
    abort();
  }
  void BoundaryForm(GlobalVector& f,
                    const GlobalVector& u,
                    const IntSet& Colors,
                    const BoundaryEquation& BE,
                    double d) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryForm\" not written!"
              << std::endl;
    abort();
  }
  void BoundaryMatrix(MatrixInterface& A,
                      const GlobalVector& u,
                      const IntSet& Colors,
                      const BoundaryEquation& BE,
                      double d) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryMatrix\" not written!"
              << std::endl;
    abort();
  }
  void MassMatrix(MatrixInterface& M) const
  {
    std::cerr << "\"DiscretizationInterface::MassMatrix\" not written!"
              << std::endl;
    abort();
  }

  void BoundaryMassMatrix(MatrixInterface& A, const IntSet& Colors) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryMassMatrix\" not written!"
              << std::endl;
    abort();
  }
  void MassForm(GlobalVector& f,
                const GlobalVector& u,
                const TimePattern& TP,
                double s) const
  {
    std::cerr << "\"DiscretizationInterface::MassForm\" not written!"
              << std::endl;
    abort();
  }
  void DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const
  {
    std::cerr << "\"DiscretizationInterface::DiracRhs\" not written!"
              << std::endl;
    abort();
  }
  void BoundaryRhs(GlobalVector& f,
                   const IntSet& Colors,
                   const BoundaryRightHandSide& NRHS,
                   double s) const
  {
    std::cerr << "\"DiscretizationInterface::BoundaryRhs\" not written!"
              << std::endl;
    abort();
  }
  void StabForm(GlobalVector& f,
                const GlobalVector& u,
                const Equation& EQ,
                double d) const
  {
    std::cerr << "\"DiscretizationInterface::StabForm\" not written!"
              << std::endl;
    abort();
  }

  // Hanging Nodes
  void HNAverage(GlobalVector& x) const {}
  void HNDistribute(GlobalVector& x) const {}
  void HNZero(GlobalVector& x) const {}
  bool HNZeroCheck(const GlobalVector& x) const { return false; }
  void HNAverageData() const {}
  void HNZeroData() const {}

  // Interpolation?
  void Interpolate(GlobalVector& u, const DomainInitialCondition& U) const
  {
    std::cerr << "\"DiscretizationInterface::Interpolate\" not written!"
              << std::endl;
    abort();
  }
  void InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateSolution\" not written!"
              << std::endl;
  }
  void InterpolateDirac(GlobalVector& u, const GlobalVector& uold) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateDirac\" not written!"
              << std::endl;
    abort();
  }
  void ConstructInterpolator(MgInterpolatorInterface* I,
                             const MeshTransferInterface* MT)
  {
    std::cerr
      << "\"DiscretizationInterface::ConstructInterpolator\" not written!"
      << std::endl;
  }

  /// ??
  void InitFilter(DoubleVector&) const
  {
    std::cerr << "\"DiscretizationInterface::InitFilter\" not written!"
              << std::endl;
  }

  // Functional & Error
  void ComputeError(const GlobalVector& u,
                    LocalVector& err,
                    const ExactSolution* ES) const
  {
    std::cerr << "\"DiscretizationInterface::ComputeError\" not written!"
              << std::endl;
    abort();
  }
  void AssembleError(GlobalVector& eta,
                     const GlobalVector& u,
                     LocalVector& err,
                     const ExactSolution* ES) const
  {
    std::cerr << "\"DiscretizationInterface::AssembleError\" not written!"
              << std::endl;
    abort();
  }
  double ComputeBoundaryFunctional(const GlobalVector& u,
                                   const IntSet& Colors,
                                   const BoundaryFunctional& BF) const
  {
    std::cerr << "\"DiscretizationInterface::ComputeBoundaryFunctional\" not "
                 "written!"
              << std::endl;
    abort();
  }
  double ComputeDomainFunctional(const GlobalVector& u,
                                 const DomainFunctional& F) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputeDomainFunctional\" not written!"
      << std::endl;
    abort();
  }
  double ComputeErrorDomainFunctional(const GlobalVector& u,
                                      const DomainFunctional& F) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputeDomainFunctional\" not written!"
      << std::endl;
    abort();
  }

  double ComputePointFunctional(const GlobalVector& u,
                                const PointFunctional& FP) const
  {
    std::cerr
      << "\"DiscretizationInterface::ComputePointFunctional\" not written!"
      << std::endl;
    abort();
  }

  void EvaluateCellRightHandSide(GlobalVector& f,
                                 const DomainRightHandSide& CF,
                                 double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::EvaluateCellRighthandside\" not "
                 "written!"
              << std::endl;
    abort();
  }

  void EvaluateBoundaryCellRightHandSide(GlobalVector& f,
                                         const IntSet& Colors,
                                         const BoundaryRightHandSide& CF,
                                         double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::"
                 "EvaluateBoundaryCellRighthandside\" not written!"
              << std::endl;
    abort();
  }

  void EvaluateParameterRightHandSide(GlobalVector& f,
                                      const DomainRightHandSide& CF,
                                      double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::"
                 "EvaluateParameterRighthandside\" not written!"
              << std::endl;
    abort();
  }

  void EvaluateBoundaryParameterRightHandSide(GlobalVector& f,
                                              const IntSet& Colors,
                                              const BoundaryRightHandSide& CF,
                                              double d = 1.) const
  {
    std::cerr << "\"DiscretizationInterface::"
                 "EvaluateBoundaryParameterRighthandside\" not written!"
              << std::endl;
    abort();
  }

  void InterpolateDomainFunction(GlobalVector& f,
                                 const DomainFunction& DF) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateDomainFunction\" not "
                 "written!"
              << std::endl;
    abort();
  }

  void InterpolateCellDomainFunction(GlobalVector& f,
                                     const DomainFunction& DF) const
  {
    std::cerr << "\"DiscretizationInterface::InterpolateCellDomainFunction\" "
                 "not written!"
              << std::endl;
    abort();
  }
  void GetVolumes(DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetVolumes\" not written!"
              << std::endl;
    abort();
  }

  void GetAreas(DoubleVector& a, const IntSet& Colors) const
  {
    std::cerr << "\"DiscretizationInterface::GetAreas\" not written!"
              << std::endl;
    abort();
  }

  void GetMassDiag(DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetMassDiag\" not written!"
              << std::endl;
    abort();
  }

  void GetBoundaryMassDiag(DoubleVector& a) const
  {
    std::cerr << "\"DiscretizationInterface::GetBoundaryMassDiag\" not written!"
              << std::endl;
    abort();
  }

  void RhsCurve(GlobalVector& F, const Curve& C, int comp, int N) const
  {
    std::cerr << "\"DiscretizationInterface::RhsCurve\" not written!"
              << std::endl;
    abort();
  }

  // This functions are not used in DG methods
  void StrongDirichletMatrix(MatrixInterface& A,
                             int col,
                             const std::vector<int>& comp) const
  {}
  void StrongDirichletMatrixOnlyRow(MatrixInterface& A,
                                    int col,
                                    const std::vector<int>& comp) const
  {}
  void StrongDirichletVector(GlobalVector& u,
                             const DirichletData& BF,
                             int col,
                             const std::vector<int>& comp,
                             double d = 1.) const
  {}
  void StrongDirichletVectorZero(GlobalVector& u,
                                 int col,
                                 const std::vector<int>& comp) const
  {}
  void StrongPeriodicVector(GlobalVector& u,
                            const PeriodicData& BF,
                            int col,
                            const std::vector<int>& comp,
                            double d = 1.) const
  {}
};
} // namespace Gascoigne

/*----------------------------   dg.h     ---------------------------*/
/* end of #ifndef __dg_H */
#endif
/*----------------------------   dg.h     ---------------------------*/
