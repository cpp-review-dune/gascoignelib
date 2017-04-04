/*----------------------------   aleloop.h     ---------------------------*/
/*      $Id: aleloop.h,v 1.4 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __aleloop_H
#define __aleloop_H
/*----------------------------   aleloop.h     ---------------------------*/

#include "stdloop.h"
#include "meshagent.h"
#include "alemultilevelsolver.h"

namespace Gascoigne
{

  class Circle : public BoundaryFunction<2>
  {
  public:
    std::string GetName() const {return "tut";}
    
    double operator()(const Vector& c) const
    {
      return  pow(c.x()-0.2,2.0) + pow(c.y()-0.2,2.0) - 0.05*0.05;
    }

  };
  class Wave : public BoundaryFunction<2>
  {
    double __y0,__a0,__x0;
    
  public:
    void BasicInit(double y0,double a0, double x0)
    { __y0 = y0; __a0 = a0; __x0 = x0; }
    
    std::string GetName() const {return "tut";}
    
    double operator()(const Vector& c) const
    { return  c.y() - (__y0 - __a0 * sin(M_PI * c.x() / __x0)); }
  };
  
  
  class AleMeshAgent : public MeshAgent
  {
    Circle circle;
    Wave w2,w3;
    
  public:
  AleMeshAgent() : MeshAgent()
      {
	/* w2.BasicInit(0.0,  0.03, 1.64/2.0); */
	/* w3.BasicInit(0.41, 0.03, 1.64/3.0); */
	/* AddShape(2,&w2); */
	/* AddShape(3,&w3); */
	
	AddShape(80, &circle);
	AddShape(81, &circle);
      }
  };
  
  

  class AleLoop : public StdLoop
  {

    
  public:
    void BasicInit(const ParamFile* paramfile,
		 const ProblemContainer* PC,
		 const FunctionalContainer* FC=NULL)
    {
      GetMeshAgentPointer() = new AleMeshAgent;
      GetMultiLevelSolverPointer() = new AleMultiLevelSolver;
      StdLoop::BasicInit(paramfile,PC,FC);
    }

    
    void SolveDualProblem(const std::string& duallabel,
			  VectorInterface& z, VectorInterface& u, VectorInterface& f);
    
    double Adjoint(const std::string& problemlabel,
		   const std::string& label, DoubleVector& eta,
		   const FunctionalContainer* FC,
		   VectorInterface& z, VectorInterface& u, VectorInterface& f);
    
    void AdaptMesh(const DoubleVector& eta) ;

    void run(const std::string& problemlabel);

  };
  
}



/*----------------------------   aleloop.h     ---------------------------*/
/* end of #ifndef __aleloop_H */
#endif
/*----------------------------   aleloop.h     ---------------------------*/
