#ifndef  __DiscretizationInterface_h
#define  __DiscretizationInterface_h


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
#include  "domainfunction.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments DiscretizationInterface

  ///
  ///
  /////////////////////////////////////////////

  class DiscretizationInterface
  {
    private:

    protected:

    public:
      DiscretizationInterface() {}
      virtual ~DiscretizationInterface() {}

      virtual const GlobalData& GetGlobalData() const=0;
      virtual void SetGlobalData(const GlobalData& q) const=0;

      //
      //// Functions called from the Solver
      //
      virtual std::string GetName() const=0;

      virtual void AddNodeVector(const std::string& name, const GlobalVector* q) const=0;
      virtual void DeleteNodeVector(const std::string& name) const=0;

      virtual void AddCellVector(const std::string& name, const GlobalCellVector* q) const=0;
      virtual void DeleteCellVector(const std::string& name) const=0;

      virtual void AddParameterVector(const std::string& name, const GlobalParameterVector* q) const=0;
      virtual void DeleteParameterVector(const std::string& name) const=0;

      virtual void BasicInit(const ParamFile* pf)=0;
      virtual void ReInit   (const MeshInterface* M)=0;

      virtual int n() const=0;
      virtual int n_withouthanging()const {
        return n();
      }

      virtual void Structure(SparseStructureInterface* S) const=0;
      virtual void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const=0;
      virtual void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const=0;
      virtual void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double) const=0;

      virtual void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const {
        std::cerr << "\"DiscretizationInterface::AdjointForm\" not written!" << std::endl;
        abort();
			}
      virtual void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, 
          const BoundaryEquation& BE, double d) const {
        std::cerr << "\"DiscretizationInterface::BoundaryForm\" not written!" << std::endl;
        abort();
      }
      virtual void BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, 
          const BoundaryEquation& BE, double d) const {
        std::cerr << "\"DiscretizationInterface::BoundaryMatrix\" not written!" << std::endl;
        abort();
      }
      virtual void MassMatrix(MatrixInterface& M) const {
        std::cerr << "\"DiscretizationInterface::MassMatrix\" not written!" << std::endl;
        abort();
      }
      virtual void DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const {
        std::cerr << "\"DiscretizationInterface::DiracRhs\" not written!" << std::endl;
        abort();
      }		
      virtual void BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, 
          double s) const{
        std::cerr << "\"DiscretizationInterface::BoundaryRhs\" not written!" << std::endl;
        abort();
      }
      virtual void HNAverage   (GlobalVector& x) const {}
      virtual void HNDistribute(GlobalVector& x) const {}
      virtual void HNZero      (GlobalVector& x) const {}
      virtual bool HNZeroCheck (const GlobalVector& x) const {
        return false;
      }
      virtual void HNAverageData() const {}
      virtual void HNZeroData   () const {}
      virtual void Interpolate(GlobalVector& u, const DomainInitialCondition& U) const {
        std::cerr << "\"DiscretizationInterface::Interpolate\" not written!" << std::endl;
        abort();
      }
      virtual void InterpolateSolution(GlobalVector& u, const GlobalVector& uold)const {
        std::cerr << "\"DiscretizationInterface::InterpolateSolution\" not written!" << std::endl;
        abort();
      }
      virtual void StrongDirichletMatrix(MatrixInterface& A, int col, const std::vector<int>& comp) const {
        std::cerr << "\"DiscretizationInterface::StrongDirichletmatrix\" not written!" << std::endl;
        abort();
      }
      virtual void StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const std::vector<int>& comp) const {
        std::cerr << "\"DiscretizationInterface::StrongDirichletMatrixOnlyRow\" not written!" << std::endl;
        abort();
      }
      virtual void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp, double d=1.) const {
        std::cerr << "\"DiscretizationInterface::StronDirichletVector\" not written!" << std::endl;
        abort();
      }
      virtual void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const {
        std::cerr << "\"DiscretizationInterface::StrongDirichletVectorZero\" not written!" << std::endl;
        abort();
      }

      virtual void InitFilter(DoubleVector&) const {
        std::cerr << "\"DiscretizationInterface::InitFilter\" not written!" << std::endl;
        abort();
      }
			
      virtual void StabForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const {
        std::cerr << "\"DiscretizationInterface::StabForm\" not written!" << std::endl;
        abort();
      }

      // Functionals
      virtual void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const{
        std::cerr << "\"DiscretizationInterface::ComputeError\" not written!" << std::endl;
        abort();
      }
      virtual double ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const{
        std::cerr << "\"DiscretizationInterface::ComputeBoundaryFunctional\" not written!" << std::endl;
        abort();
      }
      virtual double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const{
        std::cerr << "\"DiscretizationInterface::ComputeDomainFunctional\" not written!" << std::endl;
        abort();
      }
			
      virtual double ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const{
        std::cerr << "\"DiscretizationInterface::ComputePointFunctional\" not written!" << std::endl;
        abort();
      }

      virtual void EvaluateCellRightHandSide(GlobalCellVector& f, const DomainRightHandSide& CF, double d = 1.) const{
        std::cerr << "\"DiscretizationInterface::EvaluateCellRighthandside\" not written!" << std::endl;
        abort();
      }

      virtual void InterpolateDomainFunction(GlobalVector& f, const DomainFunction& DF) const{
        std::cerr << "\"DiscretizationInterface::InterpolateDomainFunction\" not written!" << std::endl;
        abort();
      }

      virtual void InterpolateCellDomainFunction(GlobalCellVector& f, const DomainFunction& DF) const{
        std::cerr << "\"DiscretizationInterface::InterpolateCellDomainFunction\" not written!" << std::endl;
        abort();
      }
			
      virtual void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT) {
        std::cerr << "\"DiscretizationInterface::ConstructInterpolator\" not written!" << std::endl;
        abort();
      }

      virtual void GetVolumes(DoubleVector& a) const {
	std::cerr << "\"DiscretizationInterface::GetVolumes\" not written!" << std::endl;
        abort();
      }
  };
}

#endif
