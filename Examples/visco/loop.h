/*----------------------------   loop.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/

#include "stdloop.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"
#include  "usefullfunctionsbd.h"
#include "vertex.h"
#include  "gascoignemesh2d.h"
#include  "stdsolver.h"
#include  "stdmultilevelsolver.h"


#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HASHMAP   std::tr1::unordered_map
#else
#include  <ext/hash_map>
#define HASHMAP  __gnu_cxx::hash_map
#endif

extern double BOUNDARY;


namespace Gascoigne
{

  class MySolver : public StdSolver
  {
  public:
    
    void NewMesh(const MeshInterface* mp)
    {
      BOUNDARY = true;
      StdSolver::NewMesh(mp);
    }
    void SetBoundaryVectorZero(VectorInterface& gf) const
    {
      if (!BOUNDARY) return;
      StdSolver::SetBoundaryVectorZero(gf);
    }
    
    void SetBoundaryVector(VectorInterface& gf) const
    {
      if (!BOUNDARY) return;
      StdSolver::SetBoundaryVector(gf);
    }
    
    void SetPeriodicVector(VectorInterface& gf)  const
    {
      if (!BOUNDARY) return;
      StdSolver::SetPeriodicVector(gf);
    }
    
    void SetBoundaryVectorStrong(VectorInterface& gf, const BoundaryManager& BM, const DirichletData& DD, double d)  const
    {
      if (!BOUNDARY) return;
      StdSolver::SetBoundaryVectorStrong(gf, BM, DD, d);
    }
    
    void SetPeriodicVectorStrong(VectorInterface& gf, const BoundaryManager& BM, const PeriodicData& PD, double d) const
    {
      if (!BOUNDARY) return;
      StdSolver::SetPeriodicVectorStrong(gf, BM, PD, d);
    }
    
    void DirichletMatrix() const
    {
      if (!BOUNDARY) return;
      StdSolver::DirichletMatrix();
    }
    
    void PeriodicMatrix() const
    {
      if (!BOUNDARY) return;
      StdSolver::PeriodicMatrix();
    }
    
    void SubtractMean(VectorInterface& gx) const
    {
      if (!BOUNDARY) return;
      StdSolver::SubtractMean(gx);
    }
    

    void SubtractMeanAlgebraic(VectorInterface& gx) const
    {
      if (!BOUNDARY) return;
      StdSolver::SubtractMeanAlgebraic(gx);
    }
    
  };
  
  class MyMLS : public StdMultiLevelSolver
  {
  public:
    SolverInterface* NewSolver(int solverlevel)
    {
      return new MySolver;
    }
    

  };
  

  

  template<int DIM>
    class Loop : public StdLoop
    {
      
    public:
      
      void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,const FunctionalContainer* FC)
      {
	GetMultiLevelSolverPointer() = new MyMLS;
	
	StdLoop::BasicInit(paramfile,PC,FC);
      }
      
      void run(const std::string& problemlabel);
    };
  
}


/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
