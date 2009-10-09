#ifndef  __localstokes_h
#define  __localstokes_h

#include  "stokeslps2d.h"
#include  "filescanner.h"
#include  "problemdescriptorbase.h"
#include  "dirichletdata.h"
#include  "dirichletdatabycolor.h"
#include  "periodicdata.h"

#include "residualfunctional.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

class LocalStokesDirichletData : public DirichletData 
{ 
public: 
  mutable double _d_speed;

  ~LocalStokesDirichletData(){ } 
  LocalStokesDirichletData(const ParamFile* paramfile){
    {
      DataFormatHandler DFH;
      DFH.insert("speed" , &_d_speed , 1.00) ;
      FileScanner FS(DFH); 
      FS.NoComplain(); 
      FS.readfile(paramfile,"Equation"); 
    }
  }

  std::string GetName() const {return "LocalStokesDirichletData";} 
 
  void operator()(DoubleVector& b, const Vertex2d& v, int i_color) const { 
    b.zero(); 
    
    if ( i_color==91){
      b[1] += _d_speed * 1.;
    }
  }
}; 

/* ----------------------------------------- */

class LocalStokesPeriodicData : public PeriodicData
{
  public:
    mutable double _d_speed;

    ~LocalStokesPeriodicData(){ }
    LocalStokesPeriodicData(const ParamFile* paramfile){}

    std::string GetName() const {return "LocalStokesPeriodicData";}

    void operator()(DoubleVector& b, const Vertex2d& v, int col) const {
      b.zero();
    }
    
};

/* ----------------------------------------- */

class LocalDragFunctional : public virtual Gascoigne::ResidualFunctional
{
  public:
    LocalDragFunctional() : ResidualFunctional()
  {
    __comps.push_back(1);
    __cols.insert(0);
    __scales.push_back(1);
    ExactValue() = 0.;

    __DD  = new Gascoigne::DirichletDataByColor(GetComps(),GetColors(),GetScales());
  }

  std::string GetName() const {
    return "LocalDrag";
  }
};


/* ----------------------------------------- */

class LocalLiftFunctional : public virtual Gascoigne::ResidualFunctional
{
  public:
    LocalLiftFunctional() : ResidualFunctional()
  {
    __comps.push_back(2);
    __cols.insert(0);
    __scales.push_back(1);
    ExactValue() = 0.;

    __DD  = new Gascoigne::DirichletDataByColor(GetComps(),GetColors(),GetScales());
  }

    std::string GetName() const {
      return "LocalLift";
    }
};

/* ----------------------------------------- */

class LocalStokesProblemDescriptor : public ProblemDescriptorBase
{
public:
    
    std::string GetName() const {return "LocalStokesProblemDescriptor";}
    void BasicInit(const ParamFile* pf) {
      GetParamFilePointer()        = pf;
      GetEquationPointer()         = new StokesLps2d(GetParamFile());

      GetDirichletDataPointer()    = new LocalStokesDirichletData (GetParamFile());
      GetPeriodicDataPointer()     = new LocalStokesPeriodicData  (GetParamFile());
      ProblemDescriptorBase::BasicInit(pf);
    }
};

#endif
