#ifndef __MeshInterpolator_h
#define __MeshInterpolator_h

#include "hierarchicalmesh.h"
#include "meshagent.h"
#include "discretizationinterface.h"
#include <cassert>
#include <string>

/**********************************************************/

namespace Gascoigne
{

class MeshInterpolator
{
  protected:

    HierarchicalMesh                 *_Old,*_New;
    IntVector                         _ToBeRef,_ToBeRefNew;
    MeshAgent                        *_MA;
    MeshAgentInterface               *_OMA;
    DiscretizationInterface          *_ODI,*_DI;
    nmatrix<double>                   _weights;
    std::set<int>                     _BaseCells;
    std::vector<GlobalVector>         _VecInt,_VecOld,_VecNew;
    std::vector<int>                  _NewNodeNumber;

    void CheckCell (int oldNumber, int newNumber);
    void Coarsen   (int newNumber);
    void Distribute(int oldNumber, int newNumber);
    void RefineAndInterpolate(HierarchicalMesh* Mesh, std::vector<GlobalVector>& u, const IntVector& refine, std::vector<bool>& done);

          MeshAgent*  GetMeshAgent()              { assert(_MA); return _MA; }
    const MeshAgent*  GetMeshAgent()        const { assert(_MA); return _MA; }

          MeshAgent*           GetOriginalMeshAgent()           { return dynamic_cast<MeshAgent*>(_OMA); }
    const MeshAgent*           GetOriginalMeshAgent()     const { return dynamic_cast<const MeshAgent*>(_OMA); }

          DiscretizationInterface*  GetOriginalDiscretization()       { assert(_ODI); return _ODI; }
    const DiscretizationInterface*  GetOriginalDiscretization() const { assert(_ODI); return _ODI; }

          DiscretizationInterface*  GetDiscretization()              { assert(_DI); return _DI; }
    const DiscretizationInterface*  GetDiscretization()        const { assert(_DI); return _DI; }

  public:

    MeshInterpolator();
    virtual ~MeshInterpolator();

    virtual void BasicInit(DiscretizationInterface* DI, MeshAgentInterface* MA, 
			   const std::string& name);
    virtual void RhsForProjection(GlobalVector& gf, const GlobalVector& u);

    virtual void AddVectorIntermediate(const GlobalVector& u);
    virtual void AddVectorOld         (const GlobalVector& u);
    virtual void AddVectorNew         (const GlobalVector& u);
};
}

/**********************************************************/

#endif
