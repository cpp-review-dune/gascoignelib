#ifndef __gmres_h
#define __gmres_h
#include <stdio.h>
#include "info.h"

#include "nmatrix.h"
#include "nvector.h"
#include "mult.h"

namespace Gascoigne
{
static int             KMAX = 30-1;
static nmatrix<double> H(KMAX+1, KMAX);
static DoubleVector gmresgamma(KMAX+1), ci(KMAX), si(KMAX);
static int             firstvector;

/*
   Transformation of an upper Hessenberg matrix into
   triagonal structure by givens rotation of the last column
*/

inline void givens_rotation(DoubleVector& h, DoubleVector& b, 
			    DoubleVector& ci, DoubleVector& si, int col)
{
  for (int i=0 ; i<col ; i++)
    {
      double dummy = h[i];
      h[i]   =  ci[i]*dummy + si[i]*h[i+1];
      h[i+1] = -si[i]*dummy + ci[i]*h[i+1];
    }
  double r = 1./sqrt(h[col]*h[col] + h[col+1]*h[col+1]);
  si[col] = h[col+1] *r;
  ci[col] = h[col]   *r;
  h[col]  =  ci[col]*h[col] + si[col]*h[col+1];
  b[col+1]= -si[col]*b[col];
  b[col] *=  ci[col];
}

/*
   PDGMRES
   Restarted Preconditioned Direct Generalized Minimal Residual Method
*/    

template<class MATRIX, class VECTOR, class MEM, class INFO>
inline int gmres(MATRIX& A, VECTOR& x, VECTOR& b,
		      MEM& mem, INFO& info)
{
  int reached = 0;
  int kmax    = mem.n()-1;
  int j = 0;
  while (!reached)
    {
      info.iteration() = j*(kmax-1);
      reached = dgmres(A,x,b,mem,info);
      j++;
    }
  if (reached<0) return 1;
  return 0;
}

/*
   PDGMRES
   Preconditioned Direct Generalized Minimal Residual Method
*/    

template<class MATRIX, class VECTOR, class MEM, class INFO>
inline int dgmres(MATRIX& A, VECTOR& x, VECTOR& b,
		      MEM& mem, INFO& info)
{
  int k0   = info.iteration();
  int kmax = mem.n()-1;
  //nmatrix<double> H(kmax+1, kmax);
  //DoubleVector gamma(kmax+1), ci(kmax), si(kmax);

  //firstvector = 0;
  if (!firstvector)
    {
      H.zero();
      gmresgamma.zero();
      ci.zero();
      si.zero();
    }

  int i,k,reached=0, dim;
  int left_precondition = 1;

  VECTOR& v = mem[firstvector];
  VECTOR& p = mem[kmax];
  
  if (left_precondition)
  {
    A.residual(p,x,b);
    A.precondition(v,p);
  }
  else
  {
    A.residual(v,x,b);
  }
  
  double rho = sqrt(v*v);
  if(rho==0.) return 1;
  gmresgamma[0] = rho;
  
  v.equ(1./rho,v);
  
  for (k=firstvector ; k<kmax-1 && (!reached) ; k++)
  {
    VECTOR& vv = mem[k+1];
    if (left_precondition)
    {
      A.vmulteq(p, mem[k]);
      A.precondition(vv,p);
    }
    else
    {
      A.precondition(p,mem[k]);
      A.vmulteq(vv,p);
    }

    dim = k+1;

    /* Orthogonalization */
    
    DoubleVector h(kmax);
    for (i=0 ; i<dim ; i++)
    {
      h[i] = vv * mem[i];
      vv.add(-h[i],mem[i]);
    }
    double s = sqrt(vv*vv);
    h[k+1] = s;
    
    /* Re-orthogonalization */
    
    for (i=0 ; i<dim ; i++)
    {
      double htmp = vv * mem[i];
      h[i] += htmp;
      vv.add(-htmp,mem[i]);
    }
    s = sqrt(vv*vv);
    h[k+1] = s;
    
    vv.equ(1./s, vv);
    
    /*  Transformation into triagonal structure  */
    
    givens_rotation(h,gmresgamma,ci,si,k);

    /*  append vector on matrix  */

    for (i=0 ; i<dim ; i++)
      H(i,k) = h[i];

    /*  residual  */

    rho = fabs(gmresgamma[dim]);
    reached = info.check_residual(k0+k,rho);
  }
  firstvector = k;
  
  /*  Calculate solution  */  

  DoubleVector h(dim);
  nmatrix<double> H1(dim+1,dim);

  for (i=0 ; i<dim+1 ; i++)
    for (int j=0 ; j<dim ; j++)
      H1(i,j) = H(i,j);

  backward(h,H1,gmresgamma);

  if (left_precondition)
  {
    for (i=0 ; i<dim ; i++)
      x.add(h[i], mem[i]);
  }
  else
  {
    p = 0.;
    for (i=0 ; i<dim ; i++)
      p.add(h[i], mem[i]);
    A.precondition(v,p);
    x.add(1.,v);
  }
  
  return reached;
}
}

#endif
