#ifndef __MeshInterpolator_h
#define __MeshInterpolator_h

#include "hierarchicalmesh.h"
#include "meshagent.h"
#include "problemdescriptorinterface.h"
#include "solverinterface.h"
#include <cassert>
#include <string>

/**********************************************************/

namespace Gascoigne
{
class MeshInterpolator
{
  protected:
    BasicGhostVector                  _help;
    HierarchicalMesh                 *_Old,*_New;
    IntVector                         _ToBeRef,_ToBeRefNew;
    MeshAgent                        *_MA;
    MeshAgentInterface               *_OMA;
    const ProblemDescriptorInterface *_PD;
    SolverInterface                  *_OSI,*_SI;
    int                               _dim;
    nmatrix<double>                   _weights;
    std::set<int>                     _BaseCells;
    std::string                       _name;
    std::vector<GlobalVector*>        _VecOld,_VecNew;
    std::vector<int>                  _NewNodeNumber;

    void CheckCell(int oldNumber, int newNumber);
    void Coarsen(GlobalVector& f, int newNumber);
    void Distribute(GlobalVector& f, int oldNumber, int newNumber);
    void RefineAndInterpolate(HierarchicalMesh* Mesh, std::vector<GlobalVector*>& u, const IntVector& refine, std::vector<bool>& done);

          MeshAgent*& GetMeshAgentPointer()       { return _MA; }
          MeshAgent*  GetMeshAgent()              { assert(_MA); return _MA; }
    const MeshAgent*  GetMeshAgent()        const { assert(_MA); return _MA; }

          MeshAgentInterface*& GetOriginalMeshAgentPointer()       { return _OMA; }
          MeshAgent*           GetOriginalMeshAgent()              { assert(dynamic_cast<MeshAgent*>(_OMA)); return dynamic_cast<MeshAgent*>(_OMA); }
    const MeshAgent*           GetOriginalMeshAgent()        const { assert(dynamic_cast<const MeshAgent*>(_OMA)); return dynamic_cast<const MeshAgent*>(_OMA); }

    const ProblemDescriptorInterface*& GetProblemPointer()       { return _PD; }
    const ProblemDescriptorInterface*  GetProblem()        const { assert(_PD); return _PD; }

          SolverInterface*& GetOriginalSolverPointer()       { return _OSI; }
          SolverInterface*  GetOriginalSolver()              { assert(_OSI); return _OSI; }
    const SolverInterface*  GetOriginalSolver()        const { assert(_OSI); return _OSI; }

          SolverInterface*& GetSolverPointer()       { return _SI; }
          SolverInterface*  GetSolver()              { assert(_SI); return _SI; }
    const SolverInterface*  GetSolver()        const { assert(_SI); return _SI; }

  public:
    MeshInterpolator();
    virtual ~MeshInterpolator();

    virtual void AddVectorOld(GlobalVector* u);
    virtual void AddVectorNew(GlobalVector* u);
    virtual void BasicInit(SolverInterface* SI, MeshAgentInterface* MA, const std::string& name, const ProblemDescriptorInterface* PD=NULL);
    virtual void RhsForProjection(BasicGhostVector& gf);
};
}

/**********************************************************/

#endif
