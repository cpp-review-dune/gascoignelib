#ifndef  __SolverInterface_h
#define  __SolverInterface_h

#include  "gascoigne.h"
#include  "matrixinterface.h"
#include  "iluinterface.h"

#include  "multigridmeshinterface.h"
#include  "meshinterpretorinterface.h"
#include  "functional.h"

#include  "problemdescriptorinterface.h"
#include  "meshinterface.h"
#include  "mginterpolatorinterface.h"
#include  "meshtransferinterface.h"
#include  "basicghostvector.h"
#include  "paramfile.h"


/*---------------------------------------------------------*/

namespace Gascoigne
{
  
  //////////////////////////////////////////////
  ///
  ///@brief
  /// Interface class for Solver

  ///  Some porporties
  ///  - lives on one level of the hierarchy
  ///  - stores the matrices
  ///  - holds memory for vectors
  ///  - provides e.g. nonlinear and linear residuals of the equations
  ///  - calls class MeshInterpretor
  ///
  //////////////////////////////////////////////
  
  class SolverInterface
  {
    private: 

    protected:

    public:
      SolverInterface() {}
      virtual ~SolverInterface(){}

      virtual std::string GetName() const=0;

      virtual void BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP)=0;
  
      virtual void SetProblem(const ProblemDescriptorInterface& PDX)=0;
      virtual const ProblemDescriptorInterface* GetProblemDescriptor() const=0;

//      virtual void SetState(const std::string& s)=0;

      virtual bool DirectSolver() const=0;
  
      virtual void NewMesh(int l, const MeshInterface* MP)=0;
      virtual const MeshInterface* GetMesh() const=0;

      virtual void RegisterMatrix()=0;
      virtual void ReInitVector()=0;
      virtual void ReInitMatrix()=0;

      virtual void OutputSettings() const=0;
      virtual void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)=0;

      virtual void VisuGrid(const std::string& name, int i) const=0;

      virtual void AddNodeVector(const std::string&, const GlobalVector* q)=0;
      virtual void AddCellVector(const std::string&, const GlobalCellVector* q)=0;
      virtual void AddParameterVector(const std::string&, const GlobalParameterVector* q)=0;
      virtual void DeleteNodeVector(const std::string&)=0;
      virtual void DeleteCellVector(const std::string&)=0;
      virtual void DeleteParameterVector(const std::string&)=0;

      virtual const MeshInterpretorInterface* GetMeshInterpretor() const=0;
      virtual       MeshInterpretorInterface* GetMeshInterpretor()      =0;

      //
      /// vector - manamgement
      //
      virtual void ResizeVector(GlobalVector* x, std::string type) const=0;
      virtual void RegisterVector(const BasicGhostVector& g)=0;
      virtual GlobalVector& GetGV(BasicGhostVector& u) const=0;
      virtual const GlobalVector& GetGV(const BasicGhostVector& u) const=0;

      //
      /// vector - hanging nodes
      //
      virtual void HNAverage   (const BasicGhostVector& x) const=0;
      virtual void HNZero      (const BasicGhostVector& x) const=0;
      virtual void HNDistribute(BasicGhostVector& x) const=0;

      virtual void HNAverage   (const GlobalVector& x) const=0;
      virtual void HNZero      (const GlobalVector& x) const=0;
      virtual bool HNZeroCheck(const GlobalVector& x) const=0;
      virtual void HNDistribute(GlobalVector& x) const=0;
      virtual void HNAverageData() const=0;
      virtual void HNZeroData() const=0;

      //
      /// vector - io
      //
      virtual void Visu(const std::string& name, const BasicGhostVector& u, int i) const=0;
      virtual void Write(const BasicGhostVector& u, const std::string& filename) const=0;
      virtual void Read(BasicGhostVector& u, const std::string& filename) const=0;

      //
      /// vector - interpolation
      //
      virtual void InterpolateSolution(BasicGhostVector& u, const GlobalVector& uold) const=0;

      //
      /// vector - rhs (integration)
      //
      virtual void Rhs(BasicGhostVector& f, double d=1.) const=0;
      virtual void TimeRhsOperator(BasicGhostVector& f, const BasicGhostVector& u) const {
        std::cerr << "\"SolverInterface::TimeRhsOperator\" not written!" << std::endl;
        abort();
      }
      virtual void TimeRhs(int k, BasicGhostVector& f) const {
        std::cerr << "\"SolverInterface::TimeRhs\" not written!" << std::endl;
        abort();
      }
      virtual void IC(BasicGhostVector& f, double d=1.) const {
        std::cerr << "\"SolverInterface::IC\" not written!" << std::endl;
        abort();
      }
      virtual void L2Projection(BasicGhostVector& u) {
        std::cerr << "\"SolverInterface::L2Projection\" not written!" << std::endl;
        abort();
      }

      //
      /// vector - residual (integration)
      //
      virtual void Form(BasicGhostVector& y, const BasicGhostVector& x, double d) const=0;

      //
      /// vector - boundary condition
      //
      virtual void SetBoundaryVector(BasicGhostVector& f) const=0;
      virtual void SetBoundaryVectorZero(BasicGhostVector& f) const=0;
      virtual void SetBoundaryVectorStrong(BasicGhostVector& f, const BoundaryManager& BM, const DirichletData& DD) const=0;
  
      //
      /// vector - linear algebra
      //
      virtual double NewtonNorm(const BasicGhostVector& u) const=0;

      virtual void residualgmres(BasicGhostVector& y, const BasicGhostVector& x, const BasicGhostVector& b) const=0;
      virtual void MatrixResidual(BasicGhostVector& y, const BasicGhostVector& x, const BasicGhostVector& b) const=0;
      virtual void vmult(BasicGhostVector& y, const BasicGhostVector& x, double d) const=0;
      virtual void vmulteq(BasicGhostVector& y, const BasicGhostVector& x) const=0;
      virtual void smooth_pre(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const=0;
      virtual void smooth_exact(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const=0;
      virtual void smooth_post(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const=0;

      //
      /// vector - additional
      //
      virtual void SubtractMean(BasicGhostVector& x) const=0;
      virtual void SubtractMeanAlgebraic(BasicGhostVector& x) const=0;

      //
      /// vector - matrix
      //
      virtual void AssembleMatrix(const BasicGhostVector& u, double d=1.)=0;
      virtual void DirichletMatrix() const=0;
      virtual void MatrixZero() const=0;
      virtual void ComputeIlu(const BasicGhostVector& u) const=0;
      virtual void ComputeIlu() const=0;
      virtual void AssembleDualMatrix(const BasicGhostVector& gu, double d)=0;

      //
      /// vector - "postprocessing"
      //
      virtual void ComputeError(const BasicGhostVector& u, GlobalVector& err) const=0;
      virtual double ComputeFunctional(BasicGhostVector& f, const BasicGhostVector& u, const Functional* FP) const=0;

      //
      /// vector - initialize
      //
      virtual void BoundaryInit(BasicGhostVector& u) const=0;
      virtual void SolutionInit(BasicGhostVector& u) const=0;
};
}
#endif

