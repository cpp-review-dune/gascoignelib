/*----------------------------   aleboundary.h     ---------------------------*/
/*      $Id: aleboundary.h,v 1.3 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __aleboundary_H
#define __aleboundary_H
/*----------------------------   aleboundary.h     ---------------------------*/



#include  "boundaryequation.h"
#include  "chi.h"
#include  "alebase.h"
#include  "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  template<int DIM>
    class AleBoundary : public virtual BoundaryEquation, public AleBase<DIM>
    {
    protected:
      
      Chi   __chi;
      
      mutable Vertex<DIM> __n,__v;
      mutable int __domain;
      double __p_left,__p_right;
      double __nu_f,__rho_f;
      mutable double __h;
      std::string __solid_type;
      
    public:
      
      AleBoundary() {abort();}
      AleBoundary(const ParamFile* pf);
      
      std::string GetName()  const { return "FSI-Ale Boundary";}
      
      int         GetNcomp() const { return 1+DIM+DIM; }

      void pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const;
      
      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const;
      void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N,int col) const;

    };



  template<int DIM>
    class AleBoundaryAdjoint : public virtual BoundaryEquation
    {
    protected:
      
      Chi   __chi;
      
      mutable Vertex<DIM> __n;
      mutable int __domain;
      double __nu_f;
      std::string __solid_type;      
    public:
      
      AleBoundaryAdjoint() {abort();}
      AleBoundaryAdjoint(const ParamFile* pf);
      
      std::string GetName()  const { return "FSI-Ale Boundary ADjoint";}
      
      int         GetNcomp() const { return 1+DIM+DIM; }

      void pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const;
      
      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const;
      void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N,int col) const;

    };
}



/*----------------------------   aleboundary.h     ---------------------------*/
/* end of #ifndef __aleboundary_H */
#endif
/*----------------------------   aleboundary.h     ---------------------------*/
