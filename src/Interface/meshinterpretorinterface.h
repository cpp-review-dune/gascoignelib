#ifndef  __MeshInterpretorInterface_h
#define  __MeshInterpretorInterface_h


#include  <string>
#include  "meshinterface.h"
#include  "gascoigne.h"
#include  "equation.h"
#include  "matrixinterface.h"
#include  "dirichletdata.h"
#include  "boundaryrighthandside.h"
#include  "boundaryequation.h"
#include  "exactsolution.h"
#include  "boundaryfunctional.h"
#include  "domainfunctional.h"
#include  "mginterpolatorinterface.h"
#include  "meshtransferinterface.h"
#include  "globaldata.h"
#include  "paramfile.h"
#include  "pointfunctional.h"
#include  "domainrighthandside.h"
#include  "diracrighthandside.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments MeshInterpretorInterface

  ///
  ///
  /////////////////////////////////////////////

  class MeshInterpretorInterface
  {
    private:
      mutable GlobalData __q;

    protected:
      const GlobalData& GetGlobalData() const {return __q;}

    public:
      MeshInterpretorInterface() {}
      virtual ~MeshInterpretorInterface() {}

      //
      //// Functions called from the Solver
      //
      virtual std::string GetName() const=0;

      virtual void AddNodeVector(const std::string& name, const GlobalVector* q) const {
        __q.AddNodeVector(name,q);
      }
      virtual void DeleteNodeVector(const std::string& name) const {
        __q.DeleteNodeVector(name);
      }
      virtual void AddCellVector(const std::string& name, const GlobalCellVector* q) const {
        __q.AddCellVector(name,q);
      }
      virtual void DeleteCellVector(const std::string& name) const {
        __q.DeleteCellVector(name);
      }
      virtual void AddParameterVector(const std::string& name, const GlobalParameterVector* q) const {
        __q.AddParameterVector(name,q);
      }
      virtual void DeleteParameterVector(const std::string& name) const {
        __q.DeleteParameterVector(name);
      }

      virtual void BasicInit(const ParamFile* pf)=0;
      virtual void ReInit   (const MeshInterface* M)=0;

      virtual int n() const=0;
      virtual int n_withouthanging()const {
        return n();
      }

      virtual void Structure(SparseStructureInterface* S) const=0;
      virtual void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const=0;
      virtual void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const=0;
      virtual void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, 
          const BoundaryEquation& BE, double d) const {
        std::cerr << "\"MeshInterpretorInterface::BoundaryForm\" not written!" << std::endl;
        abort();
      }
      virtual void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double) const=0;
      virtual void BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, 
          double d) const {
        std::cerr << "\"MeshInterpretorInterface::BoundaryMatrix\" not written!" << std::endl;
        abort();
      }
      virtual void MassMatrix(MatrixInterface& M) const=0;
      virtual void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const=0;
      virtual void DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const=0;
      virtual void BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, double s) const=0;

      virtual void HNAverage   (GlobalVector& x) const {}
      virtual void HNDistribute(GlobalVector& x) const {}
      virtual void HNZero      (GlobalVector& x) const {}
      virtual bool HNZeroCheck (const GlobalVector& x) const {
        return false;
      }
      virtual void HNAverageData() const {}
      virtual void HNZeroData   () const {}
      virtual void Interpolate(GlobalVector& u, const DomainInitialCondition& U) const {
        std::cerr << "\"MeshInterpretorInterface::Interpolate\" not written!" << std::endl;
        abort();
      }
      virtual void InterpolateSolution(GlobalVector& u, const GlobalVector& uold)const {
        std::cerr << "\"MeshInterpretorInterface::InterpolateSolution\" not written!" << std::endl;
        abort();
      }
      virtual void StrongDirichletMatrix(MatrixInterface& A, int col, const std::vector<int>& comp) const {
        std::cerr << "\"MeshInterpretorInterface::StrongDirichletmatrix\" not written!" << std::endl;
        abort();
      }
      virtual void StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const std::vector<int>& comp) const {
        std::cerr << "\"MeshInterpretorInterface::StrongDirichletMatrixOnlyRow\" not written!" << std::endl;
        abort();
      }
      virtual void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const {
        std::cerr << "\"MeshInterpretorInterface::StronDirichletVector\" not written!" << std::endl;
        abort();
      }
      virtual void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const {
        std::cerr << "\"MeshInterpretorInterface::StrongDirichletVectorZero\" not written!" << std::endl;
        abort();
      }

      virtual void InitFilter(DoubleVector&) const=0;

      virtual void StabForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const {
        std::cerr << "\"MeshInterpretorInterface::StabForm\" not written!" << std::endl;
        abort();
      }

      // Functionals
      virtual void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const=0;
      virtual double ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const=0;
      virtual double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const=0;
      virtual double ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const=0;

      virtual void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT) {
        std::cerr << "\"MeshInterpretorInterface::ConstructInterpolator\" not written!" << std::endl;
        abort();
      }
  };
}

#endif
