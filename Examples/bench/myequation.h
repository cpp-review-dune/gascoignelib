/*----------------------------   equation.h     ---------------------------*/
/*      $Id: myequation.h,v 1.1 2009/10/26 08:42:26 richter Exp $                 */
#ifndef __equation_H
#define __equation_H
/*----------------------------   equation.h     ---------------------------*/


#define  __MyEquation_h

#include  "paramfile.h"
#include  "equation.h"
#include  "boundaryequation.h"
#include <fstream>

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

  
  class Field
  {
  public:
    int _NX,_NY;
    double _minx,_maxx,_miny,_maxy;
    vector<double> _DAT;

    void init(const string& fname, int NX, int NY,
	      double minx,double miny, double maxx,double maxy)
    {
      _NX = NX;
      _NY = NY;
      _minx=minx;_miny=miny;
      _maxx=maxx;_maxy=maxy;
      ifstream in(fname.c_str(), std::ifstream::binary);
      assert(in.is_open());

      _DAT.resize(_NX*_NY);
      for (int x=0;x<_NX;++x)      
	for (int y=0;y<_NY;++y)
	  in >> _DAT[y*_NX+x];
      string tut;
      in >> tut;
      cout << tut << endl;
      assert(tut=="x");
      in.close();
    }
    
    double operator()(double x,double y) const
    {
      double X = (x-_minx)/(_maxx-_minx);
      double Y = (y-_miny)/(_maxy-_miny);
      assert(X>=0);
      assert(X<1);
      assert(Y>=0);
      assert(Y<1);
      int iX = static_cast<int> (X*_NX);
      int iY = static_cast<int> (Y*_NY);
      assert(iX>=0);
      assert(iY>=0);
      assert(iX<_NX);
      assert(iY<_NY);
      return _DAT[iY*_NX+iX];
    }
  };
  
  
  

  
  class MyEquation : public virtual Equation, public virtual BoundaryEquation
  {
  protected:
    double alpha0,  tau0, Cstern, Pstern, Tref, Lref, MZ, ellipse, C, Cdw, rhow, theta_w,f;
    double deltamin, shock0;
    mutable double shock;
    mutable double alpha,tau;
    mutable double visc;
    mutable double uwx;
    mutable double uwy;
    
   // mutable double WWX,WWY;
    

   //mutable double eta,zeta,rho;
    mutable double rho;
    mutable FemFunction *oldu, *extu, *H;
    mutable double h_;

    mutable double gamma;
    double vin0;
    
    mutable double v_in, gamma0;
    mutable Vertex2d _n;

    Field Wx,Wy;
  public:

    void SetFemData(FemData& q) const 
    {
      assert(q.find("oldu") != q.end() );
      oldu = &q["oldu"];
      assert(q.find("extu") != q.end() );
      extu = &q["extu"];
      assert(q.find("H") != q.end() );
      H = &q["H"];
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    MyEquation() { abort(); }
    MyEquation(const ParamFile* pf);

    std::string GetName()  const { return "MyEquation";}

    int         GetNcomp() const {return 2;}

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;



  // Rand
  void pointboundary(double h, const FemFunction& U, const Vertex2d& v, const Vertex2d& n) const
  {
    gamma = gamma0 / h; 
    _n = n;
    v_in = vin0 * Tref/Lref;
  }

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
  {
    if (col == 3) // einstroemung
      {
	//b[0] += gamma * (U[0].m()-v_in) * N.m();
	//b[1] += gamma * U[1].m() * N.m();
      }

 else if (col==0)// Rand aussen
	{
	b[0] += gamma * U[0].m()* N.m();
	b[1] += gamma * U[1].m() * N.m();
      }
    else if (col == 4) 
      {
	b[0] += gamma * U[0].m()* N.m();
	b[1] += gamma * U[1].m() * N.m();
      }
    else if (col == 1) // oben/unten  ( v*n = 0) 

      {	b[0] += gamma * U[0].m()* N.m();
	b[1] += gamma * U[1].m() * N.m();
	
      }
    else if (col == 5) // insel links/rechts ( v*n = 0) 
      {
      	b[0] += gamma * U[0].m() * N.m();
	b[1] += gamma * U[1].m()* N.m();
      }
    else if (col == 6) // insel oben/unten ( v*n = 0) 
      {
	b[0] += gamma * U[0].m() * N.m();
	b[1] += gamma * U[1].m() * N.m();
      }
  }

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
  {
    if (col == 3) // einstroemung
      {
	//	A(0,0) += gamma * M.m() * N.m();
	//	A(1,1) += gamma * M.m() * N.m();
      }

    if (col == 0) // einstroemung
      {
		A(0,0) += gamma * M.m() * N.m();
		A(1,1) += gamma * M.m() * N.m();
         }

   else if (col == 4) // ausstroemung
      {

		A(0,0) += gamma * M.m() * N.m();
		A(1,1) += gamma * M.m() * N.m();
         }
    else if (col == 1) // oben/unten  ( v*n = 0) 
      {
		A(0,0) += gamma * M.m() * N.m();
		A(1,1) += gamma * M.m() * N.m();
      }
    else if (col == 5) // insel links/rechts ( v*n = 0) 
      {
      	A(0,0) += gamma * M.m() * N.m();
	A(1,1) += gamma * M.m() * N.m();
      }
    else if (col == 6) // insel oben/unten ( v*n = 0) 
      {
	A(0,0) += gamma * M.m() * N.m();
	A(1,1) += gamma * M.m()* N.m();
      }
  }



    
  };



  class TransportEquation : public virtual Equation
  {
  protected:
    double tau0;
    double shock0;
    mutable double shock;
    mutable double tau;
    mutable double Tref, Lref;

    mutable double h_;
    mutable double gamma, gamma0;
    mutable Vertex2d _n;

  public:

    void SetFemData(FemData& q) const 
    {
    }

    void point(double h, const FemFunction &U, const Vertex2d &v) const{};
    TransportEquation() { abort(); }
    TransportEquation(const ParamFile* pf){};

    std::string GetName()  const { return "Transport";}

    int         GetNcomp() const {return 2;}

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const{};
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const{};


};



 class OtherEquation : public virtual Equation
  {
  protected:
    mutable FemFunction  *UU;
    mutable FemFunction  *HH;
    mutable double rho;
        double alpha0,  tau0, Cstern, Pstern, Tref, Lref, MZ, ellipse, C, Cdw, rhow, theta_w,f, gamma;
    double deltamin, shock0;
    mutable double shock;
    mutable double alpha,tau;
    mutable double visc;
    mutable double uwx;
     mutable double uwy;

  public:

    void SetFemData(FemData& q) const 
    {
      assert(q.find("U") != q.end() );
      UU = &q["U"];
      assert(q.find("H") != q.end() );
      HH = &q["H"];
    }

    //    void point(double h, const FemFunction &U, const Vertex2d &v) const;
    OtherEquation() { abort(); }
    OtherEquation(const ParamFile* pf);

    std::string GetName()  const { return "Other";}

    int         GetNcomp() const {return 2;}

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
  };



  
}


/*----------------------------   equation.h     ---------------------------*/
/* end of #ifndef __equation_H */
#endif
/*----------------------------   equation.h     ---------------------------*/
