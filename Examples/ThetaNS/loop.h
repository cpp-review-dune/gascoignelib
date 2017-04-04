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
#include "solvers.h"

#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HASHMAP   std::tr1::unordered_map
#else
#include  <ext/hash_map>
#define HASHMAP  __gnu_cxx::hash_map
#endif



namespace Gascoigne
{
  


  class Cyl : public BoundaryFunction<3>
  {
    double squareradius;
    Vertex2d center;
  public :
  
    std::string GetName() const { return "Cyl";}

    void BasicInit(Vertex2d c, double r) 
    {
      center = c; 
      squareradius = r;
    }
  
    double operator()(const Vertex3d& c) const 
    {
      double r = - squareradius;
      for (int i=0; i<2; i++)
	{
	  double dx = c[i]-center[i];
	  r += dx * dx;
	}
      return r;
    }
  };

  class Cir : public BoundaryFunction<2>
  {
    double squareradius;
    Vertex2d center;
  public :
  
    std::string GetName() const { return "Cir";}

    void BasicInit(Vertex2d c, double r) 
    {
      center = c; 
      squareradius = r;
    }
  
    double operator()(const Vertex2d& c) const 
    {
      double r = - squareradius;
      for (int i=0; i<2; i++)
	{
	  double dx = c[i]-center[i];
	  r += dx * dx;
	}
      return r;
    }
  };

  /*---------------------------------------------------*/


  class MA2d : public MeshAgent
  {
  protected:
  
    Cir RK;
  
  public:
  
    MA2d() : MeshAgent()
    {
      double r2 = 0.05*0.05;
      Vertex2d v(0.2,0.2);
      RK.BasicInit(v,r2);
      AddShape(80,&RK);
    }
  };


  class Loop : public StdLoop
  {
    vector<GlobalVector> __U,__Z;

    nvector<double> __dt, __theta, __start, __stop,__time;
    public:
    
    void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,const FunctionalContainer* FC)
    {
      GetMeshAgentPointer() = new MA2d;
      GetMultiLevelSolverPointer() = new MLS;
      StdLoop::BasicInit(paramfile,PC,FC);
    }
    void AdaptFSTMesh(nvector<double>& fst);
    void InitUniformTimeMesh(double start, double stop, double dt,double theta);
    void InitUniformFSTMesh(double start, double stop, double dt);

    string Solve(VectorInterface& u, VectorInterface& f, string name);
    
    void AdjointLoop(VectorInterface& z, VectorInterface& uu, VectorInterface& f, VectorInterface& zold);
    
    void PrimalLoop(VectorInterface& u, VectorInterface& f, VectorInterface& old);
    nvector<double> GetEndFunctional(const GlobalVector& U);
    nvector<double> GetMeanFunctional();
    void WriteFunctionals();
    
    
    void extrapolate(const vector<nvector<double> >& j);
    void run(const std::string& problemlabel);
    void run_fst_adaptive(const std::string& problemlabel);
    void run_exact(const std::string& problemlabel);
    void run_exact_fst(const std::string& problemlabel);


    ///////// ESTIMATE

    void InitPrimal(GlobalVector& U, int m, int q);
    void InitPrimalTime(GlobalVector& dtU, int m);
    void InitDual(GlobalVector& Z, int m, int q);
    void InitDualHigher(GlobalVector& Z, int m, int q);
    void InitPrimalHigher(GlobalVector& U, int m, int q);
    void InitDualHigherFST(GlobalVector& Z, int m, int q);
    void InitPrimalHigherFST(GlobalVector& U, int m, int q);
    void InitPrimalHigherEnd(GlobalVector& U);
    void InitPrimalTimeHigher(GlobalVector& U, int m, int q);
    void InitPrimalTimeHigherFST(GlobalVector& U, int m, int q);

    void EstimateQuadrature(nvector<double>& estA, VectorInterface& u, VectorInterface& f);
    void EstimatePrimal(nvector<double>& estP, VectorInterface& u, VectorInterface& f);
    void EstimateDual(nvector<double>& estD, VectorInterface& u, VectorInterface& f);
    

  };
  
}


/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
