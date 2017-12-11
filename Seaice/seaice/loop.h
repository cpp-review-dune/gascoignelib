/*----------------------------   loop.h     ---------------------------*/
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/


#include "stdloop.h"
#include  "periodicmapping.h"
#include "meshagent.h"
#include "gascoignemesh2d.h"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"
#include  "simplematrix.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"


namespace Gascoigne
{

  class MyPeriodicMapping : public virtual PeriodicMapping
  {
    public:
    std::string GetName() const { return "Periodic Data"; }
    
    void transformCoords(Vertex2d& w, const Vertex2d& v) const
    {
      w.x()=v.x()-40;
      w.y()=v.y();
    }
    void transformCoords(Vertex3d& w, const Vertex3d& v) const
    {
      abort();
    }

  };
  

  
  
  class Insel : public BoundaryFunction<2>
  {
  protected:
    
    typedef Vertex<2>  Vector;
    
  public :
    
    std::string GetName() const { return "Insel"; }
    virtual double operator()(const Vector& c) const
    {
      double x = c.x()-2.5;
      double y = c.y()-2.5;
      return pow(0.5,8.0)-pow(x,8.0)-pow(y,8.0);
    }
    
    /* virtual void newton(Vector& x) const */
    /* { */
    /*   if (x.x()<-3) x.y() = 10./9. * (x.x()+11); */
    /*   else if (x.x()<-1) x.y() = 175./18.  - 5./9. * x.x() - 5./18. * x.x()*x.x(); */
    /*   else if (x.x()<1) x.y() = 10.0; */
    /*   else if (x.x()<3) x.y() = 175./18.  + 5./9. * x.x() - 5./18. * x.x()*x.x(); */
    /*   else x.y() = -10./9.* (x.x()-11.);	 */
    /* } */

  };


  class Loop : public StdLoop
  {
    SimpleMatrix    MM;
    nvector<double> LMM;
   
   MyPeriodicMapping mymap;
   Insel insel;
  public:
   void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC)
   {
     GetMeshAgentPointer() = new MeshAgent();
     
     GetMeshAgent()->AddShape(5,&insel);
     GetMeshAgent()->AddShape(6,&insel);
     
     StdLoop::BasicInit(paramfile, PC, FC);
   }
   
   string PrecondCGMass(GlobalVector& u, GlobalVector& f, const TimePattern& TP, double s);
   void run(const std::string& problemlabel);
   void runwater(const std::string& problemlabel);
   void AssembleMassMatrix();
   void SolveDIV(VectorInterface& div,VectorInterface& vel,VectorInterface& f);
   void SolveTransport(VectorInterface& h, VectorInterface& hl, VectorInterface& hh, 
		       VectorInterface& div, VectorInterface& vel, VectorInterface& f);
  
  };
  
}



/*----------------------------   loop.h     ---------------------------*/
#endif
/*----------------------------   loop.h     ---------------------------*/
