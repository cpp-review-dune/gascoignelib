#/*----------------------------   problem.h     ---------------------------*/
/*      $Id: problem.h,v 1.2 2009/10/26 09:07:16 richter Exp $                 */
#ifndef __problem_H
#define __problem_H
/*----------------------------   problem.h     ---------------------------*/


#include  "problemdescriptorbase.h"
#include  "myequation.h"
#include  "water.h"
#include  "dirichletdata.h"
#include  "domainrighthandside.h"
#include  "periodicdata.h"
#include  "stress.h"
#include<fstream>
#include  "zerodirichletdata.h"
#include "componentinformationbase.h"
extern double TIME, DT,DDD;

using namespace std;

namespace Gascoigne
{

/*---------------------------------------------------*/
  
 /* class MyPeriodic : public virtual PeriodicData
  {
    public:
    //std::string GetName() const { return "Periodic Data"; }
    
    virtual std::set<int> preferred_colors()const 
    {
      return std::set<int>();
    }
  };
*/

// class for specifying Dirichlet data
  class SeaIceDirichletData : public DirichletData
  {
  public:
    std::string GetName() const { return "SeaIce Dirichlet Data"; }
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();
      // double bb = v.x()*v.y()*(1.28-v.x())*(1.28-v.y());
      // b[0] = 0.0005*bb;
      
      //  b[1] = 0.0005*bb;
    }
 };

  // class for specifying Dirichlet data
  class TransportDirichletData : public DirichletData
  {
  public:
    std::string GetName() const { return "Transport Dirichlet Data"; }
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();
      // Solver
      //  b[1] = 0.625*v.x()+0.2;
      //b[0] = 1.40625*v.x()+0.2;
      b[0]=0.5;
      //  if ( v.y()>0.78125*v.x())
      //     if (v.x()>0.64)
      //{
      //   b[0]=1;
      //  }
    
      //b[1] = 0.78125*v.x();
      // b[0] = 1.5625*v.x();
      // b[0] = 1.0;
      
	  b[1] = 1.0;
    }
 };

  
 


  
  /*---------------------------------------------------*/

 class MyRHSTransport : public virtual DomainRightHandSide
  {
    mutable FemFunction* oldh;
    
    mutable FemFunction*  oldu;

    void SetFemData(FemData& q) const 
    { 
	assert(q.find("oldh")!=q.end()); 
	oldh = &q["oldh"];

	assert(q.find("oldu")!=q.end()); 
	oldu = &q["oldu"];
    }

  public:
    
    MyRHSTransport(const ParamFile* pf) : DomainRightHandSide()
    {

      DataFormatHandler DFH;
      FileScanner FS(DFH,pf,"Equation");
      FS.NoComplain();
	
	

    }
    


    std::string GetName() const {return "RHST";}    
    int GetNcomp() const {return 2; }
    
    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    
 {
// Gleichung fuer h

    double div=((*oldu)[0].x()+(*oldu)[1].y());
    double d_t1= (*oldu)[0].m()*(*oldh)[0].x()+(*oldu)[1].m()*(*oldh)[0].y() ; 
    double d_t2= (*oldh)[0].m() * div;
    double d_t= d_t1+d_t2;

    b[0] += -  d_t*N.m();

    // Gleichung für A
    

    double d_At1= (*oldu)[0].m()*(*oldh)[1].x()+(*oldu)[1].m()*(*oldh)[1].y() ; 
    double d_At2= (*oldh)[1].m() * div;
    double d_At= d_At1+d_At2;
    
    b[1] +=  - d_At*N.m();
    }
  };

class MyDD : virtual public DirichletData
  {
  public:
    std::string GetName() const {return "DD";}  
  
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const 
    {
      b.zero();
    }
    void operator()(DoubleVector& b, const Vertex3d& v, int col) const 
    {
      b.zero();
    }
  };





  // class for specifying problem's right hand side
  class MyRightHandSide : public DomainRightHandSide
  {

  public:

    void SetFemData(FemData& q) const 
    {
    }

    int GetNcomp()        const { return 4; }
    std::string GetName() const { return "Rhs"; }
    
    double rho,rhoa, Lref,Tref,Cda,theta_a, windx,windy;

    MyRightHandSide(const ParamFile* paramfile) : DomainRightHandSide()
      {
	DataFormatHandler DFH;
	DFH.insert("rho", &rho, 0.);
	DFH.insert("rhoa", &rhoa, 0.);
	DFH.insert("Lref",&Lref,0.0);
	DFH.insert("Tref",&Tref,0.0);
	DFH.insert("Cda",&Cda,0.0);
	DFH.insert("theta_a",&theta_a,0.0);
	DFH.insert("windx",&windx,0.0);
	DFH.insert("windy",&windy,0.0);
	FileScanner FS(DFH);
	FS.NoComplain();
	FS.readfile(paramfile, "Equation");
	assert(rho>0);
	assert(rhoa>0);
	assert(Lref>0);
	assert(Tref>0);

      }

    double operator()(int comp, const Vertex2d& v) const
    {
      




      //double Cda   = 2e-3;

   //   double vw_x  = windx  * Tref / Lref;
    // double vw_y  = windy  * Tref / Lref;
      
           double SZ = rhoa/rho * Cda * Lref;
      
           double time = TIME * Tref;
     double X = v.x()*Lref;
      double Y = v.y()*Lref;
      double Lx = 1.28e6;
      double Ly = 1.28e6;
      double THETA = 4.0 * 24.0 * 60.0 * 60.0;

      double U_x = 5.0 +(sin(2*M_PI*TIME*Tref/(4*24*60*60))- 3.0) * sin(2.0 * M_PI * X / Lx) * sin(M_PI*Y/Ly);
      
      double U_y = 5.0 +(sin(2*M_PI*TIME*Tref/(4*24*60*60))- 3.0) * sin(2.0 * M_PI * Y/ Ly) * sin(M_PI*X/Lx);
 
    // double U_x=0;
    // double U_y=0;
      
    double ux = U_x/Lref * Tref;
    double uy = U_y/Lref * Tref;
   
  
 // anfang im September
    // double zeit = (0.5+0.5*cos(2.0*M_PI* (TIME*Tref + 3.0*30.0*24.*60*60 )/(365*24*60*60)));
    // double ux=0.0008*zeit;

  //   double uy=(1-zeit)*(-0.0014);
       double vw    = sqrt(ux*ux + uy*uy);
    
   
      
      double time_factor = 0.5*(1.0-cos(TIME * M_PI/10.0));
      //double time_factor=1.0;
      if (TIME>10) time_factor = 1.0;
      
      
      if (comp==0)
	{
	 double air_x =SZ * vw * (ux * cos(theta_a) - uy * sin(theta_a));
	 // double air_x=0; 
	  return time_factor * air_x;
	}
      else if (comp==1)
	{
	   double air_y  =SZ * vw * (uy * cos(theta_a) + ux * sin(theta_a));
	  //double air_y=0;  
	  return time_factor * air_y;
	}
      return 0;
    }
  };
  

  /*--------------------------------------------------*/
  
  
  class MyCI : public ComponentInformationBase 
  {
    
    public:
      
      virtual void      GetScalarName   (int i, std::string& s_name) 
      {
	if (i==0) s_name = "vx";
	if (i==1) s_name = "vy";
      }
      virtual const int GetNVectors     () const { return 1; }
      virtual void      GetVectorName   (int i, std::string& s_name) const { s_name = "V";}
      virtual void      GetVectorIndices(int i, fixarray<3,int>& fa) const
      {
	fa[0]=0; fa[1]=1; fa[2]=-1;
      }
  };

  
  class MyCIWater : public ComponentInformationBase 
  {
    
    public:
      
      virtual void      GetScalarName   (int i, std::string& s_name) 
      {
	if (i==0) s_name = "p";
	if (i==1) s_name = "vx";
	if (i==2) s_name = "vy";
      }
      virtual const int GetNVectors     () const { return 1; }
      virtual void      GetVectorName   (int i, std::string& s_name) const { s_name = "V";}
      virtual void      GetVectorIndices(int i, fixarray<3,int>& fa) const
      {
	fa[0]=1; fa[1]=2; fa[2]=-1;
      }
  };

  
  class MyCITr : public ComponentInformationBase 
  {
    
    public:
      
      virtual void      GetScalarName   (int i, std::string& s_name) 
      {
	if (i==0) s_name = "h";
	if (i==1) s_name = "A";
      }
      virtual const int GetNVectors     () const { return 0; }
  };

    class WaterDirichletData : public DirichletData
  {
  public:
    std::string GetName() const { return "Dirichlet Data"; }
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();
    }
  };
  
  

  // main class for defining the problem to solve
  class SeaIceProblem : public ProblemDescriptorBase
    {
    public:
  
      std::string GetName() const {return "SeaIce Problem";}
      void BasicInit(const ParamFile* pf)
      {
	// equation to solve
	GetEquationPointer()      = new MyEquation(pf);
	GetBoundaryEquationPointer()      = new MyEquation(pf);

	// definition for right hand side
	GetRightHandSidePointer() = new MyRightHandSide(pf);

	// definition of dirichlet boundary data
	//	GetDirichletDataPointer() = new ZeroDirichletData;
	GetDirichletDataPointer() = new SeaIceDirichletData();
	GetComponentInformationPointer() = new MyCI;
	// 
	ProblemDescriptorBase::BasicInit(pf);
      }
    };
  



  // --------------------------------------------------
  // --------------------------------------------------
  // --------------------------------------------------




  
  // class for specifying Dirichlet data
  class OtherDirichletData : public DirichletData
  {
  public:

    std::string GetName() const { return "Transport Dirichlet Data"; }
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();

    }

  };

  // main class for defining the problem to solve
  class TransportProblem : public ProblemDescriptorBase
    {
    public:
  
      std::string GetName() const {return "Transport Problem";}
      void BasicInit(const ParamFile* pf)
      {
	// equation to solve
	GetEquationPointer()      = new TransportEquation(pf);
        GetBoundaryEquationPointer()      = new TransportEquation(pf);
        GetDirichletDataPointer() = new TransportDirichletData();

	ProblemDescriptorBase::BasicInit(pf);
	//GetPeriodicDataPointer() = new MyPeriodic;
      }
    };
  // main class for defining the problem to solve
  class OtherProblem : public ProblemDescriptorBase
    {
    public:
  
      std::string GetName() const {return "Other Problem";}
      void BasicInit(const ParamFile* pf)
      {
	// equation to solve
	GetEquationPointer()      = new OtherEquation(pf);
	//		GetDirichletDataPointer() = new OtherDirichletData();
	   GetDirichletDataPointer() = new TransportDirichletData();
	ProblemDescriptorBase::BasicInit(pf);
	//GetPeriodicDataPointer() = new MyPeriodic;
      }
    };

  
  class WaterProblem : public ProblemDescriptorBase
    {
    public:
  
      std::string GetName() const {return "Water Problem";}
      void BasicInit(const ParamFile* pf)
      {
    	// equation to solve
    	GetEquationPointer()      = new Water(pf);
    	GetBoundaryEquationPointer()      = new Water(pf);

    	// definition for right hand side
    	//GetRightHandSidePointer() = new WaterRightHandSide(pf);

    	// definition of dirichlet boundary data
    	GetDirichletDataPointer() = new WaterDirichletData;

	GetComponentInformationPointer() = new MyCIWater;

    	ProblemDescriptorBase::BasicInit(pf);
      }
    };

  class MyBM : public BoundaryManager
  {
    IntSet _nocolors;
  public:
    const IntSet& GetDirichletDataColors        () const
    { return _nocolors; }
  };

class TransportTest : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "test";}
    void BasicInit(const ParamFile* pf) 
    {
      GetParamFilePointer() = pf;

      // Hier ist die Cij-Matrix definiert
      GetEquationPointer() = new MyEq(GetParamFile());
      // Die hat keine Dirichlet-Daten, daher neuer Boundary-Manager ohne Randfarben
      GetBoundaryManagerPointer() = new MyBM();

      GetRightHandSidePointer() = new MyRHSTransport(GetParamFile());

      GetDirichletDataPointer() = new TransportDirichletData;
      //GetInitialConditionPointer() = new MyInitial;
      ProblemDescriptorBase::BasicInit(pf);
    }
  };
  
}

  /*----------------------------   problem.h     ---------------------------*/
  /* end of #ifndef __problem_H */
#endif
  /*----------------------------   problem.h     ---------------------------*/
  
