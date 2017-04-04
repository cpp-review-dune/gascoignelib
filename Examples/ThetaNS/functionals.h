/*----------------------------   functionals.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __functionals_H
#define __functionals_H
/*----------------------------   functionals.h     ---------------------------*/


#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainrighthandside.h"
#include  "boundaryrighthandside.h"
#include "weighteddiracrighthandside.h"
#include "filescanner.h"
#include  "domainfunctional.h"


extern double __TIME;
extern double __DT;

namespace Gascoigne
{
  static double F(double __nu, int c, Vertex2d v, double t)
  {
    return 0.0;
    
    if (c==0) return 0.0;
    if (c==2) return 0.0;
    return 1.0;
    
    double sc = 16.0 * __nu * M_PI*M_PI;
    double ex = t;
    double px = M_PI*v.x();
    double py = M_PI*v.y();
    double cc = static_cast<double> (c+1);
    return 0.125*sin(cc*px)*sin(cc*py) * (ex + sc * cc*cc * (1.0+ex));
  }


  class Pdrop : public WeightedPointFunctional
  {
  protected:
    
  public:
    Pdrop()
      {
	_weights.push_back(1.0);
	_weights.push_back(-1.0);
	_v2d.push_back(Vertex2d (0.25,0.2));
	_v2d.push_back(Vertex2d (0.15,0.2));
	_comps.push_back(0);
	_comps.push_back(0);
      }
    std::string GetName() const {return "pdrop";}
  };

  

  
  class VMEAN : public DomainFunctional
  {

  public:
    std::string GetName() const { return "VMEAN"; }
  VMEAN() : Functional() {}
    
    double J(const FemFunction& U, const Vertex2d& v) const
    {
      /* if (__TIME<=5.0) return 0; */
      /* if (__TIME>=6.0) return 0; */
      
      
      return (U[1].y()-U[2].x())*(U[1].y()-U[2].x());
    }
  };



  class GalerkinResidual : public DomainFunctional
  {

    double __nu,__lpsp0,__lpsv0;
    mutable FemFunction *dtU,*Z;
    mutable double __h;
  public:

    void SetCellSize(double h) const { __h = h;}
    void SetFemData(FemData& q) const
    {
      assert(q.find("dtU")!=q.end());
      assert(q.find("Z")!=q.end());
      dtU = &q["dtU"];
      Z   = &q["Z"];
    }

    std::string GetName() const { return "GalerkinResidual"; }
    GalerkinResidual(const ParamFile* pf)
      {
	DataFormatHandler DFH;
	DFH.insert("nu" ,    &__nu , 0.0);
	DFH.insert("lpsp" ,    &__lpsp0 , 0.0);
	DFH.insert("lpsv" ,    &__lpsv0 , 0.0);
	FileScanner FS(DFH, pf, "Equation");
	assert(__nu>0);
	assert(__lpsp0>0);
      }
    
    double J(const FemFunction& U, const Vertex2d& v) const
    {
      double res=0;
      double __lpsp = __lpsp0 / (__nu/ __h / __h + 1.0/__h);
      double __lpsv = __lpsv0 / (__nu/ __h / __h + 1.0/__h);



      
      // rhs
      for (int i=0;i<3;++i)
	res += F(__nu, i, v, __TIME) * (*Z)[i].m();


      
      
      //////////// form
      // time
      //      res -= 0.01*(*dtU)[0].m()*(*Z)[0].m();


      // stab
      res -= __lpsp*(U[0].x()*(*Z)[0].x()+U[0].y()*(*Z)[0].y());
      for (int i=0;i<2;++i)
	for (int j=0;j<2;++j)
	  for (int k=0;k<2;++k)
	    res -= __lpsv * U[j+1].m()*U[i+1][j+1] * U[k+1].m()*(*Z)[i+1][k+1];

      
      for (int i=1;i<3;++i) 
	res -= (*dtU)[i].m()*(*Z)[i].m();
      
      // tensor
      for (int i=1;i<3;++i) 
	res -= __nu * ( U[i].x()*(*Z)[i].x() + U[i].y()*(*Z)[i].y() );
      // convection
      res -= (U[1].m()*U[1].x()+U[2].m()*U[1].y())*(*Z)[1].m();
      res -= (U[1].m()*U[2].x()+U[2].m()*U[2].y())*(*Z)[2].m();

      /* // divergence */
      res -= (U[1].x()+U[2].y())*(*Z)[0].m();
      res -= -U[0].m()*(*Z)[1].x();
      res -= -U[0].m()*(*Z)[2].y();
      
      
      return res;
    }
  };

  
  class InitialResidual : public DomainFunctional
  {
    mutable FemFunction *Z;
  public:
    void SetFemData(FemData& q) const
    {
      assert(q.find("Z")!=q.end());
      Z   = &q["Z"];
    }

    std::string GetName() const { return "InitialResidual"; }
    
    double J(const FemFunction& U, const Vertex2d& v) const
    {
      double res=0;
      // rhs
      for (int i=0;i<3;++i)
	res -= U[i].m() * (*Z)[i].m();
      return res;
    }
  };



  
  class GalerkinAdjointResidual : public DomainFunctional
  {
    
    double __nu,__lpsp0,__lpsv0;
    mutable FemFunction *dtUH, *UH,*Z;
    mutable double __h;
    
  public:
    void SetCellSize(double h) const { __h = h;}
    void SetFemData(FemData& q) const
    {
      assert(q.find("UH")!=q.end());
      assert(q.find("dtUH")!=q.end());
      assert(q.find("Z")!=q.end());
      UH   = &q["UH"];
      dtUH = &q["dtUH"];
      Z    = &q["Z"];
    }

    std::string GetName() const { return "GalerkinResidual"; }
    GalerkinAdjointResidual(const ParamFile* pf)
      {
	DataFormatHandler DFH;
	DFH.insert("nu" ,    &__nu , 0.0);
	DFH.insert("lpsp" ,    &__lpsp0 , 0.0);
	DFH.insert("lpsv" ,    &__lpsv0 , 0.0);
	FileScanner FS(DFH, pf, "Equation");
	assert(__nu>0);
	assert(__lpsp0>0);
      }
    
    double J(const FemFunction& U, const Vertex2d& v) const
    {
      double res=0;
      double __lpsp = __lpsp0 / (__nu/ __h / __h + 1.0/__h);
      double __lpsv = __lpsv0 / (__nu/ __h / __h + 1.0/__h);
      
      // form
      // res -= 0.01*(*dtUH)[0].m()*(*Z)[0].m();

      // stab
      res -= __lpsp*((*UH)[0].x()*(*Z)[0].x()+(*UH)[0].y()*(*Z)[0].y());
      
      for (int i=0;i<2;++i)
	for (int j=0;j<2;++j)
	  for (int k=0;k<2;++k)
	    {
	      res -= __lpsv * __DT * (*UH)[j+1].m()*U[i+1][j+1] * U[k+1].m()*(*Z)[i+1][k+1];
	      res -= __lpsv * __DT * U[j+1].m()*(*UH)[i+1][j+1] * U[k+1].m()*(*Z)[i+1][k+1];
	      res -= __lpsv * __DT * U[j+1].m()*U[i+1][j+1] * (*UH)[k+1].m()*(*Z)[i+1][k+1];	
	    }

	
      for (int i=1;i<3;++i)
	res -= (*dtUH)[i].m()*(*Z)[i].m();
      for (int i=1;i<3;++i)
	res -= __nu * ( (*UH)[i].x()*(*Z)[i].x() + (*UH)[i].y()*(*Z)[i].y() );

      res -= (U[1].m()*(*UH)[1].x() + U[2].m()*(*UH)[1].y())*(*Z)[1].m();
      res -= (U[1].m()*(*UH)[2].x() + U[2].m()*(*UH)[2].y())*(*Z)[2].m();
      res -= ((*UH)[1].m()*U[1].x()+(*UH)[2].m()*U[1].y())*(*Z)[1].m();
      res -= ((*UH)[1].m()*U[2].x()+(*UH)[2].m()*U[2].y())*(*Z)[2].m();


      /* // divergence */
      res -= ((*UH)[1].x()+(*UH)[2].y())*(*Z)[0].m();
      res -= -((*UH)[0].m()*(*Z)[1].x());
      res -= -((*UH)[0].m()*(*Z)[2].y());
      
      return res;
    }
  };



  class VMEANRhs : public DomainRightHandSide
  {

  public:
    std::string GetName() const { return "VMEAN-Rhs"; }

    int GetNcomp() const {return 3;}
    
    mutable FemFunction* U;
          
    void SetFemData(FemData& q) const
    {
      assert(q.find("U")!=q.end());
      U = &q["U"];
    }

    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
    {
      /* if (__TIME<=5.0) return ; */
      /* if (__TIME>=6.0) return ; */

      b[1] +=  2.0*N.y()*((*U)[1].y()-(*U)[2].x());
      b[2] -= -2.0*N.x()*((*U)[1].y()-(*U)[2].x());
    }
  };
  

  

  
  class Pmean : public virtual BoundaryFunctional
  {
  protected:
    
  public:

    std::string GetName() const { return "PMean"; }
    
    Pmean(const ParamFile* pf)
      {
	/* DataFormatHandler DFH; */
	/* DFH.insert("nu" ,    &__nu , 0.0); */
	/* FileScanner FS(DFH, pf, "Equation"); */
	/* assert(__nu>0); */
      }
    
    double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const
    {
      if (color!=80) return 0;
      return U[0].m();
    }
    
  };

  class PmeanRhs : public virtual BoundaryRightHandSide
  {
  protected:
    
  public:
    std::string GetName() const { return "PMean-Rhs"; }
    int GetNcomp() const { return 3;} 

    PmeanRhs(const ParamFile* pf)
      {
	/* DataFormatHandler DFH; */
	/* DFH.insert("nu" ,    &__nu , 0.0); */
	/* FileScanner FS(DFH, pf, "Equation"); */
	/* assert(__nu>0); */
      }
    
    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v, const Vertex2d& n, int color) const
    {
      if (color!=80) return;
      b[0] += N.m();
    }
  };
  

  
  
  

class Bdrag : public virtual BoundaryFunctional
{
protected:
  double __nu;

public:
  std::string GetName() const { return "Bdrag"; }

  Bdrag(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("nu" ,    &__nu , 0.0);
    FileScanner FS(DFH, pf, "Equation");
    assert(__nu>0);
  }
    
  double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const
  {
    double __nu = 1.e-3;
    
    if (color!=80) return 0;
    
    Vertex2d r;
    r[0] = __nu * (n.x() * U[1].x() + n.y() * U[1].y()) - U[0].m() * n.x();
    r[1] = __nu * (n.x() * U[2].x() + n.y() * U[2].y()) - U[0].m() * n.y();
    return -20.0*r[0];
  }
};

  class LiftBRHS : public virtual BoundaryRightHandSide
  {
  protected:
    double __nu;
    
  public:
    std::string GetName() const { return "BLift-Rhs"; }
    int GetNcomp() const { return 3;} 
    
    LiftBRHS(const ParamFile* pf)
      {
	DataFormatHandler DFH;
	DFH.insert("nu" ,    &__nu , 0.0);
	FileScanner FS(DFH, pf, "Equation");
	assert(__nu>0);
      }
    
    void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v, const Vertex2d& n, int color) const;
  };
  
  
  class Blift : public virtual BoundaryFunctional
  {
  protected:
    double __nu;
    
  public:
    std::string GetName() const { return "Blift"; }
    
    Blift(const ParamFile* pf)
      {
	DataFormatHandler DFH;
	DFH.insert("nu" ,    &__nu , 0.0);
	FileScanner FS(DFH, pf, "Equation");
	assert(__nu>0);
      }
    
    double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const
    {
      double __nu = 1.e-3;
      
      if (color!=80) return 0;
      /* if (__TIME<=5.0) return 0; */
      /* if (__TIME>=6.0) return 0; */

      
      Vertex2d r;
      r[0] = __nu * (n.x() * U[1].x() + n.y() * U[1].y()) - U[0].m() * n.x();
      r[1] = __nu * (n.x() * U[2].x() + n.y() * U[2].y()) - U[0].m() * n.y();
      return -20.0*r[1];
    }
  };
  


class Drag : public virtual ResidualFunctional
{ 
  std::string GetName() const { return "drag"; } 
public:
  Drag()
  {
    __comps.push_back(1);
    __scales.push_back(20.0);
    __cols.insert(80);   
    __DD  = new DirichletDataByColor(GetComps(),GetColors(),GetScales());

  }
};
class Lift : public virtual ResidualFunctional
{ 
  std::string GetName() const { return "drag"; } 
public:
  Lift()
  {
    __comps.push_back(2);
    __scales.push_back(20.0);
    __cols.insert(80);   
    __DD  = new DirichletDataByColor(GetComps(),GetColors(),GetScales());
  }
};


}




/*----------------------------   functionals.h     ---------------------------*/
/* end of #ifndef __functionals_H */
#endif
/*----------------------------   functionals.h     ---------------------------*/
