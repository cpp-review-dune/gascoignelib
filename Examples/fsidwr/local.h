#ifndef  __LOCAL_h
#define  __LOCAL_h

#include  "problemdescriptorbase.h"
#include  "ale.h"
#include  "aleboundary.h"
#include  "ale_slow.h"
#include  "ale-dual.h"
#include  "dual-primal-equation.h"
//#include  "aleboundary.h"
#include  "zerodirichletdata.h"
#include  "dirichletdatabycolor.h"
#include  "constantrighthandside.h"
#include  "weighteddiracrighthandside.h"
#include  "residualfunctional.h"
#include  "domainfunctional.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
  

  class LocalDomainFunctionals_FlagForce : public virtual Gascoigne::ResidualFunctional
    {
    public:
      std::string _s_force_type;
    LocalDomainFunctionals_FlagForce(std::string s_force_type) : ResidualFunctional() 
	{
	  _s_force_type = s_force_type;
    

	  if (_s_force_type=="drag")
	    {
	      __comps.push_back(3);
	      __comps.push_back(3);
	      __cols.insert(80);
	      __cols.insert(81);
	      __scales.push_back(1.0);
	      __scales.push_back(1.0);
	    }
	  else if (_s_force_type=="drag-3d")
	    {
	      __comps.push_back(4);
	      __cols.insert(84);
	      __scales.push_back(1.0);
	    }
	  else if (_s_force_type=="lift-3d")
	    {
	      __comps.push_back(5);
	      __cols.insert(84);
	      __scales.push_back(1.0);
	    }
	  else if (_s_force_type=="lift")
	    {
	      __comps.push_back(4);
	      __comps.push_back(4);
	      __cols.insert(80);
	      __cols.insert(81);
	      __scales.push_back(1.0);
	      __scales.push_back(1.0);
	    }
	  else abort();
	  
	  ExactValue() = 0.;

	  __DD  = new Gascoigne::DirichletDataByColor(GetComps(),GetColors(),GetScales());
	}

      std::string GetName() const {
	return "LocalDomainFunctionals_FlagForce:"+_s_force_type;
      }
    };


  
  class DF : public DomainFunctional
  {
    Chi __chi;

  public:
    std::string GetName() const {return "DF";}
    
    DF(const ParamFile* pf)
      {
	std::string __solid_type;
	DataFormatHandler DFH;
	DFH.insert("solid_type" ,&__solid_type);
	FileScanner FS(DFH,pf,"Equation");
	__chi.BasicInit(__solid_type);
      }
    
    double J(const FemFunction& U, const Vertex2d& v) const 
    {
      int domain = __chi(v);
      if (domain==1)
	return  U[3].m() * U[3].m() + U[4].m() * U[4].m();
      return 0;
    }
    double J(const FemFunction& U, const Vertex3d& v) const 
    {
      int domain = __chi(v);
      if (domain==1)
	return  U[4].m() * U[4].m() + U[5].m() * U[5].m() + U[6].m() * U[6].m();
      return 0;
    }
  };
  

  template<int DIM>
    class DF_RHS : public DomainRightHandSide
    {
      Chi __chi;
      
      mutable FemFunction *__U;
      
    public:
      std::string GetName() const {return "DF-RHS";}
      
      int         GetNcomp() const { return 2*DIM+1; }
      void SetFemData(FemData& q) const 
      {
	if (q.find("U")==q.end())
	  {
	    std::cerr << " U nicht da!" << std::endl;
	    abort();
	  }
	__U = &q["U"];
      }
      
      DF_RHS(const ParamFile* pf)
	{
	  std::string __solid_type;
	  DataFormatHandler DFH;
	  DFH.insert("solid_type" ,&__solid_type);
	  FileScanner FS(DFH,pf,"Equation");
	  __chi.BasicInit(__solid_type);
	}
      
      void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
      {
	int domain = __chi(v);
	if (domain==1)
	{
	    b[3] += N.m() *2.0 * (*__U)[3].m();
	    b[4] += N.m() *2.0 * (*__U)[4].m();
	  }
      }
      void operator()(VectorIterator b, const TestFunction& N, const Vertex3d& v) const 
      {
	int domain = __chi(v);
	if (domain>0)
	  {
	    b[4] += N.m() *2.0 * (*__U)[4].m();
	    b[5] += N.m() *2.0 * (*__U)[5].m();
	    b[6] += N.m() *2.0 * (*__U)[6].m();
	  }
      }
    };


  
  
  class MyDD : public DirichletData
  {
    std::string __solid_type;
  public:
    
    std::string GetName() const {return "DD";}

    MyDD(const ParamFile* pf)
      {
	DataFormatHandler DFH;
	DFH.insert("solid_type" ,&__solid_type);
	FileScanner FS(DFH,pf,"Equation");
      }
    
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();
      return;
      
      if ((__solid_type=="benchmark")||(__solid_type=="squares"))
      	{
      	  if (col==0) b[1] = v.y() * (0.41-v.y()) / 0.205/0.205 * 0.3;
      	  return;
      	}
      else if (__solid_type=="simple")
      	{
      	  if (v.y()==3.0) b[1] = 1.0;
      	  return;
      	}
      else if (__solid_type=="art")
	{
	  if ((v.y()>0.5)&&(v.y()<1.5))
	    b[1] = (1.5-v.y())*(v.y()-0.5)*0.25;
	  
	}
      

    }
    void operator()(DoubleVector& b, const Vertex3d& v, int col) const
    {
      b.zero();
      if (__solid_type=="bench3d")
	{
	  if (v.x()==0)
	    b[1] = v.y() * (0.41-v.y()) * v.z() * (0.41-v.z())
	      * pow(0.205,-4.0) * 0.45;
	}
      else if (__solid_type=="new")
	{
	  if (v.x()==0.0)
	    b[1] = 
	      v.y()*(0.4-v.y()) / 0.2 / 0.2 *
	      (0.4-v.z()) * (v.z()+0.4) / 0.4 / 0.4
	      * 0.2;
	  return;
	}

    }
  };

  
  class ProblemDescriptor : public Gascoigne::ProblemDescriptorBase
    {
    public:
      
      std::string GetName() const {return "Local";}
      void BasicInit(const Gascoigne::ParamFile* pf)
      {
	GetEquationPointer()         = new Gascoigne::Ale<2>(pf);
	GetBoundaryEquationPointer() = new Gascoigne::AleBoundary<2>(pf);
	GetDirichletDataPointer()    = new Gascoigne::MyDD(pf);

	ProblemDescriptorBase::BasicInit(pf);
      }
    };

  class ProblemDescriptor3d : public Gascoigne::ProblemDescriptorBase
    {
    public:
      
      std::string GetName() const {return "Local";}
      void BasicInit(const Gascoigne::ParamFile* pf)
      {
	GetEquationPointer()         = new Gascoigne::Ale<3>(pf);
	GetBoundaryEquationPointer() = new Gascoigne::AleBoundary<3>(pf);
	GetDirichletDataPointer()    = new Gascoigne::MyDD(pf);
	
	ProblemDescriptorBase::BasicInit(pf);
      }
    };


  class DualProblemDescriptor : public Gascoigne::ProblemDescriptorBase
    {
      std::string __name;
      
    public:

      std::string GetName() const {return __name; }

      DualProblemDescriptor(const Gascoigne::ParamFile* pf,
			    const std::string name,
			    Application* RHS,
			    DirichletData* DD)	
	{
	  __name = name;
	  
	  GetEquationPointer() = new Gascoigne::DualPrimalEquation(new Ale<2>(pf));
	  //	  GetBoundaryEquationPointer() = new Gascoigne::DualPrimalBoundaryEquation(new AleBoundary<2>(pf));
	  GetRightHandSidePointer() = RHS;
	  GetDirichletDataPointer() = DD;
	  BasicInit(pf);
	}

      ~DualProblemDescriptor()
	{
	  RHS = NULL;
	  DD  = NULL;
	}

    };

    class DualProblemDescriptor3d : public Gascoigne::ProblemDescriptorBase
    {
      std::string __name;
      
    public:

      std::string GetName() const {return __name; }

      DualProblemDescriptor3d(const Gascoigne::ParamFile* pf,
			      const std::string name,
			      Application* RHS,
			      DirichletData* DD)	
	{
	  __name = name;
	  
	  GetEquationPointer() = new Gascoigne::DualPrimalEquation(new Ale<3>(pf));
	  //	  GetBoundaryEquationPointer() = new Gascoigne::DualPrimalBoundaryEquation(new AleBoundary<2>(pf));
	  GetRightHandSidePointer() = RHS;
	  GetDirichletDataPointer() = DD;
	  BasicInit(pf);
	}

      ~DualProblemDescriptor3d()
	{
	  RHS = NULL;
	  DD  = NULL;
	}

    };

}

#endif

