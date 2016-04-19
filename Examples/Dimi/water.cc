#include  "water.h"
#include  "filescanner.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  Water::Water(const ParamFile* paramfile) : BoundaryEquation()
  {
    if (1)
      {
	DataFormatHandler DFH;
	DFH.insert("alpha", &alpha0, 0.);
	DFH.insert("delta", &delta0, 0.);
	DFH.insert("visc", &visc, 0.);
	DFH.insert("v_in", &vin0, 0.);
	DFH.insert("gamma", &gamma0, 0.);
	DFH.insert("lps", &lps0, 0.);
	FileScanner FS(DFH);
	FS.NoComplain();
	FS.readfile(paramfile, "Water");
      }
    if (1)
      {
	DataFormatHandler DFH;
	DFH.insert("Tref", &Tref, 0.);
	DFH.insert("Lref", &Lref, 0.);
	FileScanner FS(DFH);
	FS.NoComplain();
	FS.readfile(paramfile, "Equation");
      }
  }
  
  
  void Water::point(double h, const FemFunction &U, const Vertex2d &v) const
  {
    // Stabilisierung muss auch auf Dimensionslose groessen transformiert werden
    alpha = alpha0 * h * h / visc;
    delta = delta0 * h;

    lps  = lps0*h;
  }
 
  /*-----------------------------------------*/

  void Water::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    // divergenz
    b[0] += (U[1].x()+U[2].y())*N.m();
    

    // nabla u, nabla phi
    for (int i=0;i<2;++i)
      {
	for (int j=0;j<2;++j)
	  b[i+1] += visc *  U[i+1][j+1] * N[j+1];
	b[i+1] -= U[0].m()*N[i+1];

	// v * nabla v_i
	for (int j=0;j<2;++j)
	  b[i+1] +=  U[j+1].m() * U[i+1][j+1] * N.m();
      }


    // stabilisierung
    for (int j=0;j<2;++j)
      b[0] += lps * U[0][j+1]*N[j+1];

    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
	for (int k=0;k<2;++k)
	  b[i+1] += lps * U[j+1].m()*U[i+1][j+1] * U[k+1].m()*N[k+1];

  }
  /*-----------------------------------------*/
  
  void Water::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  { 
    // divergenz
    A(0,1) += (M.x())*N.m();
    A(0,2) += (M.y())*N.m();
    

    // nabla u, nabla phi
    for (int i=0;i<2;++i)
      {
	for (int j=0;j<2;++j)
	  A(i+1,i+1) += visc *  M[j+1] * N[j+1];
	A(i+1,0) -= M.m()*N[i+1];
	
	// v * nabla v_i
	for (int j=0;j<2;++j)
	  {
	    A(i+1,j+1) +=  M.m() * U[i+1][j+1] * N.m();
	    A(i+1,i+1) +=  U[j+1].m() * M[j+1] * N.m();
	  }
      }


    /// Stabi
    for (int j=0;j<2;++j)
      A(0,0) += lps * M[j+1]*N[j+1];

    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
	for (int k=0;k<2;++k)
	  {
	    A(i+1,i+1) += lps * U[j+1].m()*M[j+1] * U[k+1].m()*N[k+1];
	    A(i+1,j+1) += lps * M.m()*U[i+1][j+1] * U[k+1].m()*N[k+1];
	    A(i+1,k+1) += lps * U[j+1].m()*U[i+1][j+1] * M.m()*N[k+1];
	  }

  }


  // Rand
  void Water::pointboundary(double h, const FemFunction& U, const Vertex2d& v, const Vertex2d& n) const
  {
    gamma = gamma0 / h; 
    _n = n;
    v_in = vin0 * Tref/Lref;
  }

  void Water::Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
  {
    if (col == 3) // einstroemung
      {
	b[1] += gamma * (U[1].m()-v_in) * N.m();
	b[2] += gamma * U[2].m() * N.m();
      }

 else if (col == 0) // oben/unten  ( v*n = 0) 
      {

	b[1] += gamma * (U[1].m()-v_in) * N.m();
	b[2] += gamma * U[2].m() * N.m();

      }
    else if (col == 4) // ausstroemung
      {}
    else if (col == 1) // oben/unten  ( v*n = 0) 
      {
	b[1] += gamma * U[1].m()*_n[0] * N.m();
	b[2] += gamma * U[2].m()*_n[1] * N.m();
      }
    else if (col == 5) // insel links/rechts ( v*n = 0) 
      {
      	b[1] += gamma * U[1].m() * N.m();
	b[2] += gamma * U[2].m() * N.m();
      }
    else if (col == 6) // insel oben/unten ( v*n = 0) 
      {
	b[1] += gamma * U[1].m()* N.m();
	b[2] += gamma * U[2].m()* N.m();
      }
  }

  void Water::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
  {
    if (col == 3) // einstroemung
      {
	A(1,1) += gamma * M.m() * N.m();
	A(2,2) += gamma * M.m() * N.m();
      }
 else if (col == 0) // oben/unten  ( v*n = 0) 
      {
	A(1,1) += gamma * M.m() * N.m();
	A(2,2) += gamma * M.m()* N.m();
      }
    else if (col == 4) // ausstroemung
      {}
    else if (col == 1) // oben/unten  ( v*n = 0) 
      {
	A(1,1) += gamma * M.m()*_n[0] * N.m();
	A(2,2) += gamma * M.m()*_n[1] * N.m();
      }
    else if (col == 5) // insel links/rechts ( v*n = 0) 
      {
      	A(1,1) += gamma * M.m() * N.m();
	A(2,2) += gamma * M.m()* N.m();
      }
    else if (col == 6) // insel oben/unten ( v*n = 0) 
      {
	A(1,1) += gamma * M.m() * N.m();
	A(2,2) += gamma * M.m() * N.m();
      }
  }

  // // Stabilisierung
  // void Water::lpspoint(double h, const FemFunction& U, const Vertex2d& v) const
  // {
  //   lps = lps0 *h;
  // }

  // void Water::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
  // {
  // }
  
  // void Water::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
  // {
  // }


}

/*-----------------------------------------*/
