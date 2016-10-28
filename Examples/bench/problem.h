#/*----------------------------   problem.h     ---------------------------*/
/*      $Id: problem.h,v 1.2 2009/10/26 09:07:16 richter Exp $                 */
#ifndef __problem_H
#define __problem_H
/*----------------------------   problem.h     ---------------------------*/


#include  "problemdescriptorbase.h"
#include  "myequation.h"
#include  "dirichletdata.h"
#include  "domainrighthandside.h"
#include  "periodicdata.h"
#include<fstream>
#include  "zerodirichletdata.h"
#include "componentinformationbase.h"

extern double TIME, DT,DTSUB,DDD;

using namespace std;

namespace Gascoigne
{
  double WindX(double x, double y, double t)
  {
    double tP=t; while (tP>8) { tP-=8.0; }
    int p = ((static_cast<int>(t/8)))%2;
    double vmax = 15.0; // maximale windgeschwindigkeit in m/s
    double mx,my;
    double alpha = M_PI/2.0; // 90 grad ist ohne Konvergenz oder Divergenz
    // windstarke
    double ws = tanh(tP*(8.0-tP)/2.0);  
    //   double ws = 1.0;
    if ((p%2) == 0)
      {
	mx = 50+800/16.0*tP;
	my = 50+800/16.0*tP;
	alpha -=  M_PI/2.0/5; // 18 Grad Konvergenz
      }
    else
      { ws=-ws;
	mx = 450-800/16.0*tP;
	my = 450-800/16.0*tP;
	alpha -=  M_PI/2.0/10; // 9 Grad Divergenz 10
      }

    double wx = cos(alpha)*(x-mx) + sin(alpha)*(y-my);
    double wy = -sin(alpha)*(x-mx) + cos(alpha)*(y-my);
    double r = sqrt((mx-x)*(mx-x)+(my-y)*(my-y));
    double s = 1.0/50.0*exp(-0.01*r);

    return -wx*s*ws*vmax;
  }
  double WindY(double x, double y, double t)
  {
    double tP=t; while (tP>8) { tP-=8.0; }
    int p = ((static_cast<int>(t/8)))%2;
    double vmax = 15.0; // maximale windgeschwindigkeit
    double mx,my;
    double alpha = M_PI/2.0;  // 90 grad ist ohne Konvergenz oder Divergenz
    // windstarke
    //  double ws = 1.0;
    double ws=tanh(tP*(8.0-tP)/2.0);

    if ((p%2) == 0)
      {
	mx = 50+800/16.0*tP;
	my = 50+800/16.0*tP;
	alpha -=  M_PI/2.0/5; // 18 Grad Konvergenz
      }
    else
      {ws=-ws;
	mx = 450-800/16.0*tP;
	my = 450-800/16.0*tP;
	alpha -=  M_PI/2.0/10; // 9 Grad Divergenz
      }

    double wx = cos(alpha)*(x-mx) + sin(alpha)*(y-my);
    double wy = -sin(alpha)*(x-mx) + cos(alpha)*(y-my);
    double r = sqrt((mx-x)*(mx-x)+(my-y)*(my-y));
    double s = 1.0/50.0*exp(-0.01*r);

    return -wy * s * ws*vmax;
  }

  

  // class for specifying Dirichlet data
  class SeaIceDirichletData : public DirichletData
  {
  public:
    std::string GetName() const { return "SeaIce Dirichlet Data"; }
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();
    }
  };

  // Startwerte fuer das Transportproblem
  class TGDirichletData : public DirichletData
  {
  public:
    std::string GetName() const { return "TG Dirichlet Data"; }
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();
    
      // b[0]=0.3 + 0.01*static_cast<double>(rand()%20001-10000)/10000.0; 
      b[0] =0.005*sin(500*v.x())+0.005*sin(500*v.y())+0.3;
      b[1]=1.0;
    }
  };

  
 


  
  /*---------------------------------------------------*/

  class TGRhs : public virtual DomainRightHandSide
  {
    mutable FemFunction *oldh, *vel, *div;

    void SetFemData(FemData& q) const 
    { 
      assert(q.find("oldh")!=q.end()); 
      oldh = &q["oldh"];

      assert(q.find("vel")!=q.end()); 
      vel = &q["vel"];
    }

  public:
    
  TGRhs(const ParamFile* pf) : DomainRightHandSide()
      {

	DataFormatHandler DFH;
	FileScanner FS(DFH,pf,"Equation");
	FS.NoComplain();
      }
    
    std::string GetName() const {return "TG Rhs";}    
    int GetNcomp() const {return 2; }
    
    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    {
      // Gleichung fuer h
      b[0] += (*oldh)[0].m() * N.m() 
	- DTSUB *( (*vel)[0].x() + (*vel)[1].y()) * (*oldh)[0].m() * N.m() // eigentlich muesste hier (div_ h) eine Variable sein
	- DTSUB * ((*vel)[0].m() * (*oldh)[0].x() + (*vel)[1].m()*(*oldh)[0].y()) * N.m() 
	- 0.5 * DTSUB * DTSUB * ((*vel)[0].m() * (*oldh)[0].x() + (*vel)[1].m()*(*oldh)[0].y()) 
	* ((*vel)[0].m() * N.x() + (*vel)[1].m()*N.y());
      
      // Gleichung fuer A
      b[1] += (*oldh)[1].m() * N.m() 
	- DTSUB * ((*vel)[0].x() + (*vel)[1].y()) * (*oldh)[1].m() * N.m() // eigentlich muesste hier (div_ h) eine Variable sein
	- DTSUB * ((*vel)[0].m() * (*oldh)[1].x() + (*vel)[1].m()*(*oldh)[1].y()) * N.m() 
	- 0.5 * DTSUB * DTSUB * ((*vel)[0].m() * (*oldh)[1].x() + (*vel)[1].m()*(*oldh)[1].y()) 
	* ((*vel)[0].m() * N.x() + (*vel)[1].m()*N.y());
    }
  };


  
  /*---------------------------------------------------*/

  class DivRhs : public virtual DomainRightHandSide
  {
    mutable FemFunction *vel;

    void SetFemData(FemData& q) const 
    { 
	
    }

  public:
    
  DivRhs(const ParamFile* pf) : DomainRightHandSide()
      {
	
	DataFormatHandler DFH;
	FileScanner FS(DFH,pf,"Equation");
	FS.NoComplain();
      }
    
    std::string GetName() const {return "Div Rhs";}    
    int GetNcomp() const {return 2; }
    
    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    {
      
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
      //   double Lx = 1.28e6;
      // double Ly = 1.28e6;
      double THETA = 4.0 * 24.0 * 60.0 * 60.0;
      
      double U_x = WindX(X/1000, Y/1000, TIME*Tref/(60.0*60.0*24.0)); // umrechnen in km 
      double U_y = WindY(X/1000, Y/1000, TIME*Tref/(60.0*60.0*24.0));// umrechnen in km 
 
      double ux = U_x/Lref * Tref;
      double uy = U_y/Lref * Tref;
      
  
      // anfang im September
      // double zeit = (0.5+0.5*cos(2.0*M_PI* (TIME*Tref + 3.0*30.0*24.*60*60 )/(365*24*60*60)));
      // double ux=0.0008*zeit;

      //   double uy=(1-zeit)*(-0.0014);
      double vw    = sqrt(ux*ux + uy*uy);
    

      
      if (comp==0)
	return SZ * vw * (ux * cos(theta_a) - uy * sin(theta_a));
      else if (comp==1)
	return SZ * vw * (uy * cos(theta_a) + ux * sin(theta_a));
      else
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
  class TGProblem : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "TG Problem";}
    void BasicInit(const ParamFile* pf)
    {
      GetRightHandSidePointer() = new TGRhs(pf);
      // equation to solve
      GetEquationPointer()      = new TransportEquation(pf);// wird nicht verwendet
      GetDirichletDataPointer() = new TGDirichletData();
	
      ProblemDescriptorBase::BasicInit(pf);

    }
  };

  class OtherProblem : public ProblemDescriptorBase
  {
  public:
  
    std::string GetName() const {return "Other Problem";}
    void BasicInit(const ParamFile* pf)
    {
      // equation to solve
      GetEquationPointer()      = new OtherEquation(pf);
      GetDirichletDataPointer() = new OtherDirichletData();

      ProblemDescriptorBase::BasicInit(pf);
      //GetPeriodicDataPointer() = new MyPeriodic;
    }
  };


  
  class MyBM : public BoundaryManager
  {
    IntSet _nocolors;
  public:
    const IntSet& GetDirichletDataColors        () const
    { return _nocolors; }
  };

}

/*----------------------------   problem.h     ---------------------------*/
/* end of #ifndef __problem_H */
#endif
/*----------------------------   problem.h     ---------------------------*/
  
