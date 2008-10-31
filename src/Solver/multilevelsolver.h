#ifndef  __MultiLevelSolver_h
#define  __MultiLevelSolver_h

#include  "numericinterface.h"
#include  "mginterpolatorinterface.h"
#include  "problemcontainer.h"
#include  "problemdescriptorinterface.h"
#include  "paramfile.h"
#include  "functionalcontainer.h"
#include  "solverinterface.h"
#include  "meshagentinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear MultilevelSolver

/// - stores MultiGridMeshInterace
/// - stores array of MGInterpolator
/// - stores array of SolverInterface
///
//////////////////////////////////////////////

class MultiLevelSolver
{
  protected :

  std::vector<SolverInterface*>           _SP;
  std::vector<MgInterpolatorInterface*>   _Interpolator;

  const MeshAgentInterface*               _MAP;
  const ParamFile*                        _paramfile;
  const ProblemContainer*                 _PC;
  const ProblemDescriptorInterface*       _PD;
  const NumericInterface*                 _NI;

  int ComputeLevel;

  const MeshAgentInterface* GetMeshAgent()    const { return _MAP;}

  void NewMgInterpolator();
  void SolverNewMesh();
  void ReInitMatrix();
  void RegisterMatrix();
  virtual void NewSolvers();
  //virtual SolverInterface* NewSolver(int solverlevel);
  void Transfer(int high, int low, VectorInterface& u) const;
  void SolutionTransfer(VectorInterface& u) const;

 public:

  // Constructor

  MultiLevelSolver();
  ~MultiLevelSolver();

  std::string GetName() const {return "MultiLevelSolver";}
  int nlevels() const { assert(GetMeshAgent()); return GetMeshAgent()->nlevels();}

  const ProblemContainer* GetProblemContainer() const { return _PC;};

  int FinestLevel()                       const { return nlevels()-1;}

        SolverInterface* GetSolver(int l)       {assert(l<_SP.size()); return _SP[l];}
  const SolverInterface* GetSolver(int l) const {assert(l<_SP.size()); return _SP[l];}
        SolverInterface* GetSolver()            {assert(_SP.size()==nlevels()); 
	                                         return _SP[FinestLevel()];}
  const SolverInterface* GetSolver() const      {assert(_SP.size()==nlevels()); 
                                                 return _SP[FinestLevel()];}
  
        std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers()       { return _Interpolator; }
  const std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers() const { return _Interpolator; }

  void BasicInit(const NumericInterface* NI, const MeshAgentInterface*, const ParamFile* paramfile,
		 const ProblemContainer* PC);

  void ReInitVector(VectorInterface& v);
  void DeleteVector(VectorInterface& g);
  void SetProblem(const std::string& label);
  void ReInit(const std::string& problemlabel);
  void AssembleMatrix    (VectorInterface& u);
  void AssembleDualMatrix(VectorInterface& u);
  void ComputeIlu();
  void ComputeIlu(VectorInterface& u);
  const DoubleVector ComputeFunctionals(VectorInterface& f, const VectorInterface& u,
					FunctionalContainer* FC) const;
};
}

#endif
