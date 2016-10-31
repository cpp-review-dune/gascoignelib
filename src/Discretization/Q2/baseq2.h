/**
*
* Copyright (C) 2004, 2005, 2011 by the Gascoigne 3D authors
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


#ifndef __baseq2_h
#define __baseq2_h

#include  <vector>
#include  <string>
#include  <utility>
#include  <cassert>
#include  "vertex.h"
#include  "numfixarray.h"
#include  "base.h"


namespace Gascoigne
{

  /////////////////////////////////////////////
  ///
  ///@brief
  ///  Basis on the reference element

  ///
  ///
  /////////////////////////////////////////////

  template<int DIM>
  class BaseQ2 : public Base<DIM>
  {
    
  protected:
    
    fixarray<3,double>        a,b,c;
    
    double psi   (int i, double x) const { return a[i] + b[i]*x + c[i]*x*x;}
    double psi_x (int i, double x) const { return b[i] + 2.*c[i]*x;       }
    double psi_xx(int i, double x) const { return 2.*c[i];       }
    
  public:
    
    BaseQ2()
    {
      Base<DIM>::N  .resize((DIM==2)?9:27);  
      Base<DIM>::DN .resize((DIM==2)?9:27);

      a[0] = 1.;  b[0] = -3.;  c[0] = 2.;
      a[1] = 0.;  b[1] =  4.;  c[1] = -4.;
      a[2] = 0.;  b[2] = -1.;  c[2] = 2.;
    }
    
    
    int  n() const {return (DIM==2)?9:27;}
    void point(const Vertex<DIM>& s) const;
    
    double phi   (int i) const {return Base<DIM>::N  [i];}
    double phi_x (int i) const {return Base<DIM>::DN [i].x();}
    double phi_y (int i) const {return Base<DIM>::DN [i].y();}
    double phi_z (int i) const {return Base<DIM>::DN [i].z();}
    double phi_xx(int i) const {
      std::cerr << "\"BaseQ23d::phi_xx\" not written!" << std::endl;
      abort();
    }
    double phi_yy(int i) const {
      std::cerr << "\"BaseQ23d::phi_yy\" not written!" << std::endl;
      abort();
    }
    double phi_zz(int i) const {
      std::cerr << "\"BaseQ23d::phi_zz\" not written!" << std::endl;
      abort();
    }
    double phi_xy(int i) const {
      std::cerr << "\"BaseQ23d::phi_xy\" not written!" << std::endl;
      abort();
    }
    double phi_xz(int i) const {
      std::cerr << "\"BaseQ23d::phi_xz\" not written!" << std::endl;
      abort();
    }
    double phi_yz(int i) const {
      std::cerr << "\"BaseQ23d::phi_yz\" not written!" << std::endl;
      abort();
    }
    
    double Dphi  (int i,int d) const  { return Base<DIM>::DN[i][d];  }
    double DDphi  (int i,int d1, int d2) const  
    { 
      std::cerr << "BaseQ1::DDphi not written!" << std::endl;
      abort();
    }

    const Vertex<DIM>&  phi_grad (int i) const {return Base<DIM>::DN [i];}
  };

  template<>
  inline void BaseQ2<2>::point(const Vertex<2>& s) const
  { 
    for(int i=0;i<9;i++) 
      { 
	int ix = i%3;
	int iy = i/3;
      
	Base<2>::N  [i]     = psi   (ix,s.x()) * psi   (iy,s.y());
	Base<2>::DN [i].x() = psi_x (ix,s.x()) * psi   (iy,s.y());
	Base<2>::DN [i].y() = psi   (ix,s.x()) * psi_x (iy,s.y());
      }
  }

  template<>
  inline void BaseQ2<3>::point(const Vertex<3>& s) const
  {
    for(int i=0;i<27;i++) 
      { 
	int ix = (i%9) % 3;
	int iy = (i%9) / 3;
	int iz = i/9;
	  
	double px = psi(ix,s.x());
	double py = psi(iy,s.y());
	double pz = psi(iz,s.z());
	  
	double dpx = psi_x(ix,s.x());
	double dpy = psi_x(iy,s.y());
	double dpz = psi_x(iz,s.z());
	  
	Base<3>::N  [i]     =  px *  py *  pz;
	Base<3>::DN [i].x() = dpx *  py *  pz;
	Base<3>::DN [i].y() =  px * dpy *  pz;
	Base<3>::DN [i].z() =  px *  py * dpz;
      }
    
  }
  
}


#endif
