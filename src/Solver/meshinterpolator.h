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

    HierarchicalMesh                          *_Old,*_New;
    IntSet                                     _BaseCells,_ToBeRef,_ToBeRefNew;
    IntVector                                  _NewNodeNumber,_NewCellNumber;
    MeshAgent                                 *_MA, *_OMA;
    DiscretizationInterface                   *_ODI,*_DI;
    nmatrix<double>                            _wq1;
    std::vector<nmatrix<double> >              _wq2;
    std::vector<std::pair<int,int> >           _iq2;
    std::vector<std::pair<GlobalVector,int> >  _VecInt,_VecOld,_VecNew;

    void CheckCell (int oldNumber, int newNumber);
    void Coarsen   (int newNumber);
    void Distribute(int oldNumber, int newNumber);
    void InitIndizes(int dim);
    void InitInterpolationWeights(int dim);
    void RefineAndInterpolate(HierarchicalMesh* Mesh, std::vector<std::pair<GlobalVector,int> >& u, const IntSet& refine, std::vector<std::vector<bool> >& done);

          MeshAgent* GetMeshAgent()       { assert(_MA); return _MA; }
    const MeshAgent* GetMeshAgent() const { assert(_MA); return _MA; }

          MeshAgent* GetOriginalMeshAgent()       { assert(_OMA); return _OMA; }
    const MeshAgent* GetOriginalMeshAgent() const { assert(_OMA); return _OMA; }

          DiscretizationInterface* GetOriginalDiscretization()       { assert(_ODI); return _ODI; }
    const DiscretizationInterface* GetOriginalDiscretization() const { assert(_ODI); return _ODI; }

          DiscretizationInterface* GetDiscretization()       { assert(_DI); return _DI; }
    const DiscretizationInterface* GetDiscretization() const { assert(_DI); return _DI; }

  public:

    MeshInterpolator();
    virtual ~MeshInterpolator();

    virtual void BasicInit(DiscretizationInterface* DI, MeshAgentInterface* MA, 
			   const std::string& name);
    virtual void InterpolateCellVector(GlobalCellVector& uNew, const GlobalCellVector& uOld);
    virtual void RhsForProjection(GlobalVector& gf, const GlobalVector& u);

    virtual void AddVectorIntermediate(const GlobalVector& u, int order);
    virtual void AddVectorOld         (const GlobalVector& u, int order);
    virtual void AddVectorNew         (const GlobalVector& u, int order);

    virtual void AddCellVectorOld(const GlobalCellVector& u);
    virtual void AddCellVectorNew(const GlobalCellVector& u);
};
}

/**********************************************************/

#endif
