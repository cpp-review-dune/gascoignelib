/*----------------------------   dual-primal-equation.h     ---------------------------*/
/*      $Id: dual-primal-equation.h,v 1.3 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __dualprimalequation_H
#define __dualprimalequation_H
/*----------------------------   dual-primal-equation.h     ---------------------------*/


#include  "lpsequation.h"
#include  "chi.h"
#include  "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class DualPrimalEquation : public virtual LpsEquation
  {
  protected:
    mutable  LpsEquation*   __PE;
    mutable FemFunction*    __U;
    mutable EntryMatrix     __E;
    
    
  public:

    ~DualPrimalEquation()
      {
	assert(__PE);
	delete (__PE);
      }
  
  DualPrimalEquation(LpsEquation* PE) : __PE(PE)
    {
      assert(dynamic_cast<LpsEquation*> (__PE));

      __E.SetDimensionDof(1,1);
      __E.SetDimensionComp(GetNcomp(),GetNcomp());
      __E.resize();
      __E.SetDofIndex(0,0);
    }

    void SetFemData(FemData& q) const 
    {
      if (q.find("U")==q.end())
	{
	  std::cerr << " U nicht da!" << std::endl;
	  abort();
	}
      __U = &q["U"];
    }

      
    std::string GetName()  const { return "Primal Dual Equation" + __PE->GetName(); }
      
    int         GetNcomp() const { return __PE->GetNcomp(); }

    void point(double h, const FemFunction& Z, const Vertex<2>& v) const
    {
      __PE->point(h,(*__U),v);
    }
    void point(double h, const FemFunction& Z, const Vertex<3>& v) const
    {
      __PE->point(h,(*__U),v);
    }
    
    void Form(VectorIterator b, const FemFunction& Z, const TestFunction& N) const
    {
      for (int j=0;j<GetNcomp();++j)
	{
	  __E.zero();
	  __PE->Matrix(__E,(*__U),N , Z[j]);
	  for (int i=0;i<GetNcomp();++i)
	    b[i] += __E(j,i);
	}
    }
    
    void Matrix(EntryMatrix& A, const FemFunction& Z, const TestFunction& M, const TestFunction& N) const
    {
      abort();
    }
    
    void lpspoint(double h, const FemFunction& Z, const Vertex<2>& v) const
    { __PE->lpspoint(h,Z,v); }
    void lpspoint(double h, const FemFunction& Z, const Vertex<3>& v) const
    { __PE->lpspoint(h,Z,v); }
    
    void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& ZP, const TestFunction& N) const
    { }
    
    void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
    { }
  };




  class DualPrimalBoundaryEquation : public virtual BoundaryEquation
  {
  protected:
    mutable BoundaryEquation*   __PE;
    mutable FemFunction*        __U;
    mutable EntryMatrix         __E;
    
    
  public:

  DualPrimalBoundaryEquation(BoundaryEquation* PE) : __PE(PE)
    {
      assert(dynamic_cast<BoundaryEquation*> (__PE));
    
      __E.SetDimensionDof(1,1);
      __E.SetDimensionComp(GetNcomp(),GetNcomp());
      __E.resize();
      __E.SetDofIndex(0,0);
    }

    void SetFemData(FemData& q) const 
    {
      if (q.find("U")==q.end())
	{
	  std::cerr << " U nicht da!" << std::endl;
	  abort();
	}
      __U = &q["U"];
    }
  
  
    std::string GetName()  const { return "Primal Dual Boundary Equation" + __PE->GetName(); }
  
    int         GetNcomp() const { return __PE->GetNcomp(); }

    void pointboundary(double h, const FemFunction& Z, const Vertex<2>& v,
		       const Vertex<2>& n) const
    {
      __PE->pointboundary(h,(*__U),v,n);
    }
    void pointboundary(double h, const FemFunction& Z, const Vertex<3>& v,
		       const Vertex<3>& n) const
    {
      __PE->pointboundary(h,(*__U),v,n);
    }
    
    void Form(VectorIterator b, const FemFunction& Z, const TestFunction& N, int col) const 
    {
      for (int j=0;j<GetNcomp();++j)
	{
	  __E.zero();
	  __PE->Matrix(__E,(*__U),N , Z[j],col);
	  for (int i=0;i<GetNcomp();++i)
	    b[i] += __E(j,i);
	}
    }
    
    void Matrix(EntryMatrix& A, const FemFunction& Z, const TestFunction& M, const TestFunction& N, int col) const
    {
      abort();
    }
    
  };

}

/*----------------------------   dual-primal-equation.h     ---------------------------*/
/* end of #ifndef __dual-primal-equation_H */
#endif
/*----------------------------   dual-primal-equation.h     ---------------------------*/
