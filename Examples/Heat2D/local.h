#ifndef  __local_h
#define  __local_h

#include  "initialcondition.h"
#include  "equation.h"
#include  "filescanner.h"
#include  "problemdescriptorbase.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

class LocalEquation : public Equation
{
 private:
  mutable double _visc;
  mutable double _us, _vs;
  mutable double _h, _k, _r;
 public:
  LocalEquation(const Gascoigne::ParamFile* paramfile) {
    DataFormatHandler DFH;
    DFH.insert("visc",&_visc,1.);
    DFH.insert("h",&_h,0.4);
    DFH.insert("k",&_k,2.);
    DFH.insert("r",&_r,0.3);
    FileScanner FS(DFH,paramfile,"Equation");
    
    _us = (_r*_h)/(1.-_r);
    _vs = (1.-_us)*(_h+_us);
    cerr << "****** us, vs = **** " << _us << " " << _vs << endl; 
  }

  double GetUs() const {return _us;}
  double GetVs() const {return _vs;}

  std::string GetName() const { return "Local";}
  int    ncomp      () const {return 2;}
  void SetTimePattern(TimePattern& P) const {
    P.reservesize(ncomp(),ncomp(),0.);
    P(0,0) = 1.;
    P(1,1) = 1.;
  }

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const {
    b[0] += _visc* (U[0].x()*N.x()+U[0].y()*N.y());
    b[1] += _visc* (U[1].x()*N.x()+U[1].y()*N.y());
    
    double u = U[0].m();
    double v = U[1].m();
    double s = u/(u+_h);
    
    b[0] += N.m() * (-u*(1.-u) + s * v);
    b[1] += N.m() * (_k*_r* v - _k* s * v);
  }

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const {
    A(0,0) += _visc* (M.x()*N.x()+M.y()*N.y());
    A(1,1) += _visc* (M.x()*N.x()+M.y()*N.y());
    
    double u = U[0].m();
    double v = U[1].m();
    double MM = M.m()*N.m();
    
    double s = u/(u+_h);
    double t = _h/((u+_h)*(u+_h));
    
    A(0,0) += ( (-1.+2.*u) + t*v ) * MM;
    A(0,1) += s * MM;
    A(1,0) += -_k * t*v * MM;
    A(1,1) += ( _k*_r - _k*s ) * MM;
  }
};

class LocalInitialCondition : public InitialCondition
{
private:
  mutable double _us, _vs;
public:
  LocalInitialCondition(const LocalEquation* EQ) {
    _us = EQ->GetUs();
    _vs = EQ->GetVs();
  }
  
  std::string GetName() const { return "Local";}
  int GetNcomp() const {return 2;}  
  double operator()(int c, const Vertex2d& v) const {
  double x = v.x();
  double y = v.y();
  double eps1 = 2.e-7;
  double eps2 = 3.e-5;
  double eps3 = 1.2e-4;
  if(c==0)
    {
      double dist = - eps1*(x-0.1*y-225)*(x-0.1*y-675);
      return _us + dist;
    }
  else if(c==1)
    {
      return _vs - eps2*(x-450)-eps3*(y-150);
    }
  assert(0);
  }
};

class ProblemDescriptor : public ProblemDescriptorBase
{
public:
    
    std::string GetName() const {return "Local";}
    void BasicInit(const Gascoigne::ParamFile* pf) {
      GetEquationPointer() = new LocalEquation(GetParamFile());
      const LocalEquation* LEQ = dynamic_cast<const LocalEquation*>(GetEquation());
      GetInitialConditionPointer() = new LocalInitialCondition(LEQ);
    }
};

#endif
