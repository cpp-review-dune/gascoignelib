#include  "problem.h"
#include  "loop.h"
#include "domainfunctional.h"
#include  "filescanner.h"
#include  "stdloop.h"
#include "gascoignemesh2d.h"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"
#include  "simplematrix.h"
using namespace Gascoigne;
using namespace std;

extern double DELTAMIN;
/*---------------------------------------------------*/

class KineticEnergy  : public virtual DomainFunctional
{
protected:
  double rho;
  double Tref,Lref;
  mutable FemFunction  *H;

public:

  void SetFemData(FemData& q) const 
  {
    assert(q.find("H") != q.end() );
    H = &q["H"];
  }

  KineticEnergy(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("rho", &rho, 0.);
    DFH.insert("Tref",&Tref,0.0);
    DFH.insert("Lref",&Lref,0.0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
  }
  
  std::string GetName() const {return "Kinetic Energy";}
  double J(const FemFunction& U, const Vertex2d& v) const
  {
    if(v.x()>0.25 && v.y()>0.25)
      return rho * 0.5 * (*H)[0].m() * ( U[0].m()*Tref/Lref*U[0].m()*Tref/Lref + U[1].m()*Tref/Lref*U[1].m()*Tref/Lref);
    return 0.0;
    
  }
};
class Mass  : public virtual DomainFunctional
{
protected:
  double rho; 
 
  mutable FemFunction  *H;

public:
  void SetFemData(FemData& q) const 
  {
    assert(q.find("H") != q.end() );
    H = &q["H"];
  }

  Mass(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("rho", &rho, 0.);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    ExactValue() = 1404;
  }
  std::string GetName() const {return "Mass";}
  double J(const FemFunction& U, const Vertex2d& v) const
  {
    if ((*H)[1].m() > 0.15)
    return  (*H)[1].m();
    else return 0;
    
  }
};

class MeanSpeed  : public virtual DomainFunctional
{
protected:
  double rho; 
 
public:
  

  MeanSpeed(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    ExactValue() = 1404;
  }
  std::string GetName() const {return "MeanSpeed";}
  double J(const FemFunction& U, const Vertex2d& v) const
  {
   
    return  sqrt((U[0].m()*U[0].m())*(U[1].m()*U[1].m()));
  
    
  }
};


class MeanDeformation : public virtual DomainFunctional
{
protected:
  double rho; 
 
public:
  

  MeanDeformation(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    ExactValue() = 1404;
  }
  std::string GetName() const {return "MeanDeformation";}
  double J(const FemFunction& U, const Vertex2d& v) const
  {
   
      double div=U[0].x()+U[1].y();
      double shear=sqrt((U[0].x() - U[1].y())*(U[0].x() - U[1].y())+(U[0].y()-U[1].x())*(U[0].y()-U[1].x()));
    return  sqrt(div*div+shear*shear);
  
    
  }
};



class AvMass  : public virtual DomainFunctional
{
protected:
  double rho; 
  mutable FemFunction  *H;

public:
  void SetFemData(FemData& q) const 
  {
    assert(q.find("H") != q.end() );
    H = &q["H"];
  }

  AvMass(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("rho", &rho, 0.);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    ExactValue() = 1404;
  }
  std::string GetName() const {return "Mass";}
  double J(const FemFunction& U, const Vertex2d& v) const
  {
    if (v.x()>0.25 && v.x()< 0.5)
      return rho * (*H)[0].m();
  }
};


class Delta  : public virtual DomainFunctional
{
protected:
  double rho;
  double ellipse;
  double Tref;

public:
  
  Delta(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("rho", &rho, 0.);
    DFH.insert("ellipse",&ellipse,2.0);
    DFH.insert("Tref",&Tref,0.0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
  }
  std::string GetName() const {return "Mass";}
  double J(const FemFunction& U, const Vertex2d& v) const
  { 
    double dmin = DELTAMIN * Tref;
    double DELTAsquare =  dmin*dmin+
      (1.0 + pow(ellipse,-2.0)) * (U[0].x()*U[0].x() + U[1].y()*U[1].y())
      + pow(ellipse,-2.0) * pow(U[0].y() + U[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * U[0].x()*U[1].y();
    
    double DELTA = sqrt(DELTAsquare);
    
    return DELTA;
  }
};


ofstream ELLIPSE_OUT;

class Ellipse  : public virtual DomainFunctional
{
protected:
  double rho;
  double ellipse;
  double Tref;
  double Pstern;
  double C;
  mutable FemFunction  *H;
  
public:
  
  void SetFemData(FemData& q) const 
  {
    assert(q.find("H") != q.end() );
    H = &q["H"];
  }

  
  Ellipse(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("rho", &rho, 0.);
    DFH.insert("ellipse",&ellipse,2.0);
    DFH.insert("Tref",&Tref,0.0);
    DFH.insert("Pstern",&Pstern,0.0);
    DFH.insert("C",&C,20.0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
  }
  std::string GetName() const {return "Ellipse";}
  double J(const FemFunction& U, const Vertex2d& v) const
  { 

    double dmin = DELTAMIN * Tref;
    double DELTAsquare =  dmin*dmin+
      (1.0 + pow(ellipse,-2.0)) * (U[0].x()*U[0].x() + U[1].y()*U[1].y())
      + pow(ellipse,-2.0) * pow(U[0].y() + U[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * U[0].x()*U[1].y();
    
    double DELTA = sqrt(DELTAsquare);
 
    if (DELTAsquare<0)
      {
	cout<<DELTAsquare<<endl;
	abort();
      }
    double Pstar=Pstern*(*H)[0].m()*exp(-C*(1.0-(*H)[1].m()));
    double eta= Pstar/(2.0*DELTA);
    double zeta=eta*pow(ellipse,-2);


    double S1,S2,sqrS2,sqr1S2,sqr2S2,DIV,CONV;
    S1=(eta*(U[0].x()+U[1].y())-0.5*Pstar)/Pstar;
    sqrS2=sqrt(pow(U[0].x()-U[1].y(),2)+pow(U[0].y()+U[1].x(),2));
    S2=zeta*sqrS2/Pstar;



    DIV=U[0].x()+U[1].y();
   

    
     if (ELLIPSE_OUT.is_open())
      ELLIPSE_OUT << v.x() << "\t" << v.y() << "\t" << DELTA << "\t" << S1 << "\t" << S2 << "\t" << DELTA<< "\t" << endl;
    
    return sqrS2;
  }
};


#include  "domainfunctional.h"
#include  "residualfunctional.h" 
#include  "dirichletdatabycolor.h"

class Force1 : public virtual ResidualFunctional
{
public:
  Force1() : ResidualFunctional() 
  {
    __comps.push_back(0);
    __cols.insert(1);
    __scales.push_back(1.0);
    
    ExactValue() = 0.;
    __DD  = new DirichletDataByColor(GetComps(),GetColors(),GetScales());
  }

  std::string GetName() const { 
    return "Force 1";
  }
};

class Force2 : public virtual ResidualFunctional
{
public:
  Force2() : ResidualFunctional() 
  {
    __comps.push_back(1);
    __cols.insert(1);
    __scales.push_back(1.0);
    
    ExactValue() = 0.;
    __DD  = new DirichletDataByColor(GetComps(),GetColors(),GetScales());
  }

  std::string GetName() const { 
    return "Force 2";
  }
};



int main(int argc, char** argv)
{
  // parameter file as argument
  ParamFile paramfile("gascoigne.param");
  if(argc>=2) 
    paramfile.SetName(argv[1]);

  // define problem
  SeaIceProblem LPD;
  LPD.BasicInit(&paramfile);
  
  TGProblem TGPD;
  TGPD.BasicInit(&paramfile);
  
  OtherProblem OPD;
  OPD.BasicInit(&paramfile);
  

  ProblemContainer PC;
  PC.AddProblem("seaice", &LPD);
  PC.AddProblem("tg",    &TGPD);
  PC.AddProblem("other", &OPD);


  // output functionals (not in use now)
  FunctionalContainer FC;
  KineticEnergy K(&paramfile);
  Mass M(&paramfile);
  AvMass AvM(&paramfile);
  Delta D(&paramfile);
  Force1 F1;
  Force1 F2;
  //Ellipse EF(&paramfile);
  MeanSpeed MS(&paramfile);
  MeanDeformation MD(&paramfile);

<<<<<<< HEAD
  // FC.AddFunctional("0 Kinetic", &K);
  //  FC.AddFunctional("1 Mass", &M);
  // FC.AddFunctional("2 AvMass", &AvM);
  //  FC.AddFunctional("3 Delta", &D);
  // FC.AddFunctional("4 Ellipse", &EF);
=======
  FC.AddFunctional("0 Kinetic", &K);
  FC.AddFunctional("1 Mass", &M);
  FC.AddFunctional("2 AvMass", &AvM);
  FC.AddFunctional("3 Delta", &D);
  FC.AddFunctional("4 MeanSpeed", &MS);
  FC.AddFunctional("5 MeanDeformation", &MD);
>>>>>>> 57386a1770ecff28241d7293c1f24f9830f40954

  // loop for program control
  Loop loop;
  loop.BasicInit(&paramfile, &PC, &FC);

  // start solving
  loop.run("seaice");

  return 0;
}
