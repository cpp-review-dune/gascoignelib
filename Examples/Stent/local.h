
#ifndef  __local_h
#define  __local_h

#include  "fsi.h"
#include  "dirichletdata.h"
#include  "problemdescriptorbase.h"
#include  "domainrighthandside.h"
#include  "componentinformationbase.h"
#include <stdio.h>
#include <stdlib.h>


using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */
extern double THETA,TIME,DT;

class MyDD : public DirichletData
{
 protected:

  
 public:
  MyDD(const ParamFile* pf)
    {
      DataFormatHandler DFH;
      FileScanner FS(DFH, pf, "Equation");
    }
  
  std::string GetName() const {return "MyDD";}
  
  void operator()(DoubleVector& b, const Vertex3d& v, int color) const
  {
    b.zero();

    double x = v.x();  double y = v.y(); double z = v.z();
    double fact_time=0;
      if (TIME >= 0.0 && TIME < 0.2)
        {
	  
	  fact_time =  (-11292 * TIME*TIME+ 2493.4 * TIME);
	    
//           b[1] =  (-9291.6667 * TIME * TIME
//                + 1893.3333 * TIME
//                + 40)
//             *__EINSTROMGESCHWINDIGKEIT* v.y() * (1.61-v.y()) *4/1.61/1.61;
        }
          else if (TIME >= 0.9 && TIME < 1.1)
        {
          fact_time =  (-9291.6667 * (TIME-0.9) * (TIME-0.9)
               + 1893.3333 * (TIME-0.9)
               + 40);
        }
          else if (TIME >= 1.8 && TIME < 2.0)
        {
          fact_time=  (-9291.6667 * (TIME-1.8) * (TIME-1.8)
               + 1893.3333 * (TIME-1.8)
               + 40);
        }
          else if (TIME >= 2.7 && TIME < 2.9)
        {
          fact_time =  (-9291.6667 * (TIME-2.7) * (TIME-2.7)
               + 1893.3333 * (TIME-2.7)
               + 40);
        }
          else if (TIME >= 0.2 && TIME < 0.25)
        {
          fact_time = (120 * (TIME - 0.2) + 47);
        }
          else if (TIME >= 1.1 && TIME < 1.15)
        {
          fact_time =  (120 * (TIME - 1.1) + 47);
        }
          else if (TIME >= 2.0 && TIME < 2.05)
        {
          fact_time=  (120 * (TIME - 2.0) + 47);
        }
          else if (TIME >= 2.9 && TIME < 2.95)
        {
          fact_time =  (120 * (TIME - 2.9) + 47);
        }
          else if (TIME >= 0.25 && TIME < 0.3)
        {
          fact_time = (-220 * (TIME - 0.25) + 53);
        }
          else if (TIME >= 1.15 && TIME < 1.2)
        {
          fact_time = (-220 * (TIME - 1.15) + 53);
        }
          else if (TIME >= 2.05 && TIME < 2.1)
        {
          fact_time =  (-220 * (TIME - 2.05) + 53);
        }
          else if (TIME >= 2.95 && TIME < 3.0)
        {
          fact_time = (-220 * (TIME - 2.95) + 53);
        }
          else if (TIME >= 0.3 && TIME < 0.37)
        {
          fact_time = (114.2857142 * (TIME - 0.3) + 42);
        }
          else if (TIME >= 1.2 && TIME < 1.27)
        {
          fact_time = (114.2857142 * (TIME - 1.2) + 42);
        }
          else if (TIME >= 2.1 && TIME < 2.17)
        {
          fact_time = (114.2857142 * (TIME - 2.1) + 42);
        }
          else if (TIME >= 3.0 && TIME < 3.07)
        {
           fact_time = (114.2857142 * (TIME - 3.0) + 42);
        }
          else if (TIME >= 0.37 && TIME < 0.9)
        {
           fact_time = (50*exp(-0.4210255685 * (TIME - 0.37)));
        }
          else if (TIME >= 1.27 && TIME < 1.8)
        {
           fact_time = (50*exp(-0.4210255685 * (TIME - 1.27)));
        }
          else if (TIME >= 2.17 && TIME < 2.7)
        {
           fact_time = (50*exp(-0.4210255685 * (TIME - 2.17)));
        }
          else if (TIME >= 3.07 )
        {
           fact_time = (50*exp(-0.4210255685 * (TIME - 3.07)));
        } 



    if (color==7)
      {
	b[1] = 1.0e-2*-499.5459e-003*fact_time * ( (x-(5.05))*(x-(5.05))+(y-(-5.15))*(y-(-5.15))+(z-(4.3))*(z-(4.3))-0.3*0.3)/(-0.3*0.3);
 	b[2] = 1.0e-2*500.6380e-003 *fact_time * ((x-(5.05))*(x-(5.05))+(y-(-5.15))*(y-(-5.15))+(z-(4.3))*(z-(4.3))-0.3*0.3)/(-0.3*0.3);
	b[3] = 1.0e-2*-706.9763e-003*fact_time * ((x-(5.05))*(x-(5.05))+(y-(-5.15))*(y-(-5.15))+(z-(4.3))*(z-(4.3))-0.3*0.3)/(-0.3*0.3);
      }
    if (color==2)
      {
	double high = -0.1;
	b[3] = high * ((x-(-1))*(x-(-1))+(y-(0))*(y-(0))-0.3*0.3)/(-0.3*0.3);
      }
  }
};

// -----------------------------------------

class PD : public ProblemDescriptorBase
{
 public:
    
  std::string GetName() const {return "fsi";}
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new FSI(GetParamFile());
    GetBoundaryEquationPointer() = new FSI(GetParamFile());

    GetDirichletDataPointer() = new MyDD(GetParamFile());
    
    ProblemDescriptorBase::BasicInit(pf);

  }
};


#endif
