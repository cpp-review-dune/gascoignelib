/*----------------------------   ale.h     ---------------------------*/
/*      $Id: ale_slow.h,v 1.3 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __ale_slow_H
#define __ale_slow_H
/*----------------------------   ale.h     ---------------------------*/


#include  "lpsequation.h"
#include  "chi.h"
#include  "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  template<int DIM>
    class AleSlow : public virtual LpsEquation
    {
    protected:
      /////////// FSI
      mutable double               __J;
      mutable nmatrix<double>      __F;
      mutable nmatrix<double>      __Ftilde;
      mutable fixarray<DIM,double> __conv_v;

      mutable nmatrix<double>      __nablaV_Ftilde;
      mutable double               __divergence_v;
      mutable double               __trace_E;
      ///////////

      
      double __rho_s,__rho_f, __nu_f,__mu_s,__lambda_s;
      double __alpha_p0,__alpha_u0,__alpha_v0;

      double         __delta_lap_p0, __delta_lps_p0, __delta_lps_v0;
      mutable double __delta_lap_p, __delta_lps_p, __delta_lps_v;
      

      mutable double __alpha_p,__alpha_u,__alpha_v;
      mutable int    __domain;
      

      Chi   __chi;
      
    public:
      
      AleSlow();
      AleSlow(const ParamFile* pf);
      
      std::string GetName()  const { return "FSI-Ale-Slow";}
      
      int         GetNcomp() const { return 1+DIM+DIM; }

      void compute_F(const FemFunction& U) const;
      void compute_Ftilde() const;
      void compute_convection(const FemFunction& U) const;




      double DV_conv_v(int i, int k, const FemFunction& U, const TestFunction& M) const;
      double DU_Ftilde(int i, int j, int k, const TestFunction& M) const;
      double DU_J(int k, const FemFunction& U, const TestFunction& M) const;
      double DU_trace_E(int k, const TestFunction& M) const;

      double DU_F(int i, int j, int k, const TestFunction& M) const;
    
      double DV_nablaV_Ftilde(int i, int j, int k, const TestFunction& M) const;
      double DU_nablaV_Ftilde(int i, int j, int k, const FemFunction& U, const TestFunction& M) const;
      
      double DU_divergence_v(int k, const FemFunction& U, const TestFunction& M) const;
      double DV_divergence_v(int k, const TestFunction& M) const;
      
      void point(double h, const FemFunction& U, const Vertex<DIM>& v) const;
      void point_solid(double h, const FemFunction& U, const Vertex<DIM>& v) const;
      void point_fluid(double h, const FemFunction& U, const Vertex<DIM>& v) const;
      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
      void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
      void MatrixFluid(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
      void MatrixSolid(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;


      /////////////////////////// LPS

      void lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const;
      void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
      void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;

      
    };
}

/*----------------------------   ale.h     ---------------------------*/
/* end of #ifndef __ale_H */
#endif
/*----------------------------   ale.h     ---------------------------*/
