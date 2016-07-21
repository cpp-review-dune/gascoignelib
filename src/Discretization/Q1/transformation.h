/**
*
* Copyright (C) 2004, 2005, 2010 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#ifndef __TransformationXX_h
#define __TransformationXX_h

#include  "nmatrix.h"
#include  "vertex.h"

/*-----------------------------------------------------*/

namespace Gascoigne
{
  template<int DIM,class BASE>
  class Transformation
  {
  protected:
    
    typedef nmatrix<double>    Matrix;
    
    BASE            B;
    mutable Matrix  X;
    mutable Matrix  dt, dti;
    
    inline void ComputeDT() const;
    
    // second derivatives tensor
    const nvector<Matrix>  ComputeDDT (const Vertex<DIM>& xi) const;
    
  public:
    
    Transformation();
    
    const Matrix& DT () const {return dt ;}
    const Matrix& DTI() const {return dti;}
    
    // inverse of second derivatives tensor
    const nvector<Matrix>  DDTI(const Vertex<DIM>& xi) const;
    
    inline double        J     () const;
    inline double        G     () const;
    inline Vertex<DIM>      x     () const;
    inline Vertex<DIM>      normal() const;
    inline void  init          (const Matrix& M) {X=M;}
    inline void  ReInit          (const Matrix& M) const {X=M;}
    inline void  point         (const Vertex<DIM>& xi) const;
    inline void  point_boundary(int ie, const Vertex<DIM-1>& s) const;
    inline void  GetCoordinates(Matrix& A) const { A.equ(1.,X);}
  };
  
  /*-----------------------------------------------------*/
  
  template<int DIM,class BASE>
  inline Transformation<DIM,BASE>::Transformation() : B()
  {
    int n = B.n();
    X.memory(DIM,n);
    dt .memory(DIM,DIM);
    dti.memory(DIM,DIM);
  }

/*-----------------------------------------------------*/

  template<int DIM,class BASE>
  inline Vertex<DIM>  Transformation<DIM,BASE>::x() const 
  {
    Vertex<DIM> xp;
    for(int i=0;i<B.n();i++)
      for (int d=0;d<DIM;++d)
	xp[d] += X(d,i) * B.phi(i);
    return xp;
  }

/*-----------------------------------------------------*/

  template<int DIM,class BASE>
  inline Vertex<DIM>  Transformation<DIM,BASE>::normal() const 
  {
    Vertex<DIM> xn;
    dti.mult(xn,*B.normal());
    double xx = sqrt(xn*xn);
    xn /= xx;
    return xn;
  }
  
/*-----------------------------------------------------*/

  template<int DIM,class BASE>
  inline void  Transformation<DIM,BASE>::ComputeDT() const
  {
    dt.zero();
    for(int i=0;i<B.n();i++)
      for (int d1=0;d1<DIM;++d1)
	for (int d2=0;d2<DIM;++d2)
	  dt(d1,d2) += X(d1,i) * B.Dphi(i,d2);
    

    for (int d1=0;d1<DIM;++d1)
      for (int d2=0;d2<DIM;++d2)
	dti(d1,d2) = dt(d2,d1);
    dti.gauss_jordan();
  }
   
  
  
  /*-----------------------------------------------------*/
  
  template<int DIM,class BASE>
  inline const nvector<nmatrix<double> > Transformation<DIM,BASE>::ComputeDDT(const Vertex<DIM>& xi) const
  {
    const_cast<BASE*> (&B)->point(xi);
    
    nvector<nmatrix<double> > ddt(DIM,nmatrix<double> (DIM,DIM));
    for (int i=0;i<DIM;++i) ddt[i].zero();
    
    // Warum so, wofuer wird das benötigt?
    std::cerr << "Ich verstehe die Aufgabe von ComputeDDT nicht. Das ist direkt aus der alten Version uebernommen." << std::endl;
    
    for (int i=0;i<B.n();++i)
      for (int d1=0;d1<DIM;++d1)
	for (int d2=0;d2<DIM;++d2)
	  ddt[i](d1,d2) += X(d1,i) * B.DDphi(d1, i , d2);
    return ddt;  

  const_cast<BASE*> (&B)->point(xi);


  // Alt 3D, 2D war entsprechend.

  // nvector<nmatrix<double> > ddt(3,nmatrix<double> (3,3));
  // for (int i=0;i<3;++i) ddt[i].zero();

  // for (int i=0;i<B.n();++i)
  //   {
  //     for (int j=0;j<3;++j)
  // 	{
  // 	  ddt[0](j,0) += X(j,i) * B.phi_xx(i);
  // 	  ddt[0](j,1) += X(j,i) * B.phi_xy(i);
  // 	  ddt[0](j,2) += X(j,i) * B.phi_xz(i);

  // 	  ddt[1](j,0) += X(j,i) * B.phi_xy(i);
  // 	  ddt[1](j,1) += X(j,i) * B.phi_yy(i);
  // 	  ddt[1](j,2) += X(j,i) * B.phi_yz(i);
	  
  // 	  ddt[2](j,0) += X(j,i) * B.phi_xz(i);
  // 	  ddt[2](j,1) += X(j,i) * B.phi_yz(i);
  // 	  ddt[2](j,2) += X(j,i) * B.phi_zz(i);
  // 	}
  //   }
  // return ddt;  

    
  }
  
  /*-----------------------------------------------------*/
  
  template<int DIM,class BASE>
  inline const nvector<nmatrix<double> > Transformation<DIM,BASE>::DDTI(const Vertex<DIM>& xi) const 
  {
    const nvector<nmatrix<double> >& ddt = ComputeDDT(xi);
    
    nvector<nmatrix<double> > ddti(DIM,nmatrix<double> (DIM,DIM));
    Matrix dti_ = dti;
    
    
    dti_.transpose();
    for (int i=0;i<DIM;++i)
      {
	nmatrix<double> tmp(DIM,DIM);
	
	ddti[i].zero();
	for (int d=0;d<DIM;++d)
	  ddti[i].add(-dti_(d,i),ddt[d]);
	
	ddti[i].mmult(tmp,dti_);
	dti_.mmult(ddti[i],tmp);
      }
    return ddti;  
  }
  
  /*-----------------------------------------------------*/
  
  template<int DIM,class BASE>
  inline void  Transformation<DIM,BASE>::point(const Vertex<DIM>& xi) const
  {
    B.point(xi);
    ComputeDT();
  }
  
/*-----------------------------------------------------*/
  
  template<int DIM,class BASE>
  inline void  Transformation<DIM,BASE>::point_boundary(int ie, const Vertex<DIM-1>& s) const 
  {
    B.point_boundary(ie,s);
    ComputeDT();
  }
  
  /*-----------------------------------------------------*/
  
  template<int DIM,class BASE>
  inline double  Transformation<DIM,BASE>::J() const  
  {
    double h = dt.det();
    if (h<0)
    {
      std::cout << "Transformation<DIM,BASE>::J()  h = " << h << std::endl;
    }
    return h;
  }
  
/*-----------------------------------------------------*/
  
  template<int DIM,class BASE>
  inline double  Transformation<DIM,BASE>::G() const  
  {
    if (DIM==2)
      {
	Vertex<DIM> xt;
	dt.mult(xt,*B.tangent());
	return xt.norm_l2();
      }
    else
      {
	double d1phi=0,d2phi=0,d12phi=0;
	const fixarray<2,int>& fc = *B.faces();
	for (int i=0;i<3;++i)
	  {
	    d1phi+=dt(i,fc[0])*dt(i,fc[0]);
	    d2phi+=dt(i,fc[1])*dt(i,fc[1]);
	    d12phi+=dt(i,fc[0])*dt(i,fc[1]);
	  }
	double h = d1phi*d2phi-d12phi*d12phi;
	assert(h>=0);
	return sqrt(h);
      }
    abort();
  }
}

#endif
