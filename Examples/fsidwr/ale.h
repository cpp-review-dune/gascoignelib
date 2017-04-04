/*----------------------------   ale.h     ---------------------------*/
/*      $Id: ale.h,v 1.5 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __ale_H
#define __ale_H
/*----------------------------   ale.h     ---------------------------*/


#include  "lpsequation.h"
#include  "boundaryequation.h"
#include  "chi.h"
#include  "paramfile.h"
#include  "alebase.h"


/*-----------------------------------------*/


extern bool __ADJOINT;
extern bool __ADJOINT1;
extern bool __ESTIMATE;


namespace Gascoigne
{
  

  template<int DIM>
    class Ale : public virtual LpsEquation, public AleBase<DIM>
    {
    protected:

      Chi   __chi;
      mutable int __domain;
      double __alpha_p0;
      mutable double __alpha_p, __h, __alpha_lps;
      double __mu_s,__nu_f, __rho_f,__rho_s, __lambda_s;
      double __lambda_m, __mu_m;
      std::string __extension;
      
      double __alpha_lps0;
      mutable Vertex<DIM> __v;
      

      double Ftilde(int i,int j) const
      { return AleBase<DIM>::__Ftilde(i,j); }
      double nablaV_Ftilde(int i,int j) const
      { return AleBase<DIM>::__nablaV_Ftilde(i,j); }
      double J() const
      { return AleBase<DIM>::__J;}
      
      
      
    public:


      
      Ale();
      Ale(const ParamFile* pf);
      
      std::string GetName()  const { return "FSI-Ale";}
      int         GetNcomp() const { return 2*DIM+1; }

      void point(double h, const FemFunction& U, const Vertex<DIM>& v) const;


      void lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const
      {
	__domain = __chi(v);
	__h = h;
	double V = 0;
	for (int i=0;i<DIM;++i)
	  V += U[1+i].m() * U[i+1].m();
	V = sqrt(V);
	__alpha_lps = __alpha_lps0 / (__nu_f/h/h + V / h);
	if (__ESTIMATE) __alpha_lps = 0.0;
      }
      
      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
      void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
      
      void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
      {
	if (__domain<0)
	  {
	    for (int i=0;i<DIM;++i)
	      b[0] += __alpha_lps/__rho_f * UP[0][i+1] * N[i+1];

	    /* double convN=0; */
	    /* for (int i=0;i<DIM;++i) */
	    /*   convN += U[i+1].m() * N[i+1]; */
	    
	    /* for (int i=0;i<DIM;++i) */
	    /*   for (int j=0;j<DIM;++j) */
	    /* 	b[i+1+DIM] += __alpha_lps * __rho_f * U[j+1].m() * UP[i+1][j+1] * convN; */
	  }
      }
      void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
      {	
	if (__domain<0)
	  {
	    for (int i=0;i<DIM;++i)
	      A(0,0) += __alpha_lps/__rho_f  * Mp[i+1] * Np[i+1];

	    /* double convN=0; */
	    /* for (int i=0;i<DIM;++i) */
	    /*   convN += U[i+1].m() * Np[i+1]; */
	    /* double convM=0; */
	    /* for (int i=0;i<DIM;++i) */
	    /*   convM += U[i+1].m() * Mp[i+1]; */

	    
	    /* for (int i=0;i<DIM;++i) */
	    /*   A(i+1+DIM,i+1) += __alpha_lps * __rho_f * convN * convM; */
	  }
	
      }

      
    };


  template<int DIM>
    class AleAdjoint : public virtual LpsEquation, public AleBase<DIM>
    {
    protected:

      Chi   __chi;
      mutable int __domain;
      double __alpha_p0;
      mutable double __alpha_p, __h, __alpha_lps;
      double __mu_s,__nu_f, __rho_f,__rho_s, __lambda_s;
      double __alpha_lps0;
      mutable Vertex<DIM> __v;
      
    public:


      
      AleAdjoint();
      AleAdjoint(const ParamFile* pf);
      
      std::string GetName()  const { return "FSI-Ale Adjoint";}
      int         GetNcomp() const { return 2*DIM+1; }

      void point(double h, const FemFunction& U, const Vertex<DIM>& v) const;


      void lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const
      {
	__domain = __chi(v);
	__h = h;
	__alpha_lps = __alpha_lps0 / (__nu_f/h/h + 0.2/h);
	if (__ESTIMATE) __alpha_lps = 0.0;
      }
      
      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
      void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const { abort(); }
      
      
      void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
      {
	if (__domain<0)
	  {
	    abort();
	    b[0] += __alpha_lps * (UP[0].x() * N.x() + UP[0].y() * N.y());
	  }
      }
      void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
      {
	abort();
	if (__domain<0)
	  A(0,0) += __alpha_lps  * (Mp.x() * Np.x() + Mp.y() * Np.y());
      }

      
    };
}

/*----------------------------   ale.h     ---------------------------*/
/* end of #ifndef __ale_H */
#endif
/*----------------------------   ale.h     ---------------------------*/
