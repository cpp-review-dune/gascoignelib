#include "navierstokes.h"
#include "filescanner.h"


// 0 = Theta
// 1 = BDF
int NSSCHEME;

namespace Gascoigne
{
  ////////////////////////////////////////////////// DATA

  void NavierStokesData::BasicInit(const ParamFile *pf)
  {
    DataFormatHandler DFH;
    DFH.insert("visc", &visc, 1.);
    DFH.insert("alpha", &alpha0, 0.);
    DFH.insert("dt", &dt, 0.);
    DFH.insert("theta", &theta, 0.);
    DFH.insert("fulltensor", &fulltensor, 0);
    FileScanner FS(DFH, pf, "Equation");
    assert(visc>0);
  }

  ////////////////////////////////////////////////// NavierStokes

  template <int DIM>
  void NavierStokes<DIM>::Form(VectorIterator b,
                               const FemFunction &U,
                               const TestFunction &N) const
  {
    for (int i = 0; i < DIM; ++i)
      {
	b[0] += U[i + 1][i + 1] * N.m();
	for (int j = 0; j < DIM; ++j)
	  {
	    b[i + 1] += data.visc * U[i + 1][j + 1] * N[j + 1];
	    if (data.fulltensor)
	      b[i + 1] += data.visc * U[j + 1][i + 1] * N[j + 1];

	
	    b[i + 1] += U[j + 1].m() * U[i + 1][j + 1] * N.m();
	  }

	b[i + 1] -= U[0].m() * N[i + 1];
      }
  }

  template <int DIM>
  void NavierStokes<DIM>::Matrix(EntryMatrix &A,
                                 const FemFunction &U,
                                 const TestFunction &M,
                                 const TestFunction &N) const
  {
	
    for (int i = 0; i < DIM; ++i)
      {
	A(0, i + 1) += M[i + 1] * N.m();
	for (int j = 0; j < DIM; ++j)
	  {
	    A(i + 1, i + 1) += data.visc * M[j + 1] * N[j + 1];
	    if (data.fulltensor)
	      A(i + 1, j + 1) += data.visc * M[i + 1] * N[j + 1];

	    A(i + 1, j + 1) += M.m() * U[i + 1][j + 1] * N.m();
	    A(i + 1, i + 1) += U[j + 1].m() * M[j + 1] * N.m();
	  }

	A(i + 1, 0) -= M.m() * N[i + 1];
      }


  }




  ////////////////////////////////////////////////// NavierStokesTangent

  template <int DIM>
  void NavierStokesTangent<DIM>::Form(VectorIterator b,
				      const FemFunction &W,
				      const TestFunction &N) const
  {
    for (int i = 0; i < DIM; ++i)
      {
	b[0] += W[i + 1][i + 1] * N.m();
	for (int j = 0; j < DIM; ++j)
	  {
	    b[i + 1] += data.visc * W[i + 1][j + 1] * N[j + 1];
	    if (data.fulltensor)
	      b[i + 1] += data.visc * W[j + 1][i + 1] * N[j + 1];

	
	    b[i + 1] += (*U)[j + 1].m() * W[i + 1][j + 1] * N.m();
	    b[i + 1] += W[j + 1].m() * (*U)[i + 1][j + 1] * N.m();
	  }

	b[i + 1] -= W[0].m() * N[i + 1];
      }
  }

  template <int DIM>
  void NavierStokesTangent<DIM>::Matrix(EntryMatrix &A,
					const FemFunction &W,
					const TestFunction &M,
					const TestFunction &N) const
  {
    for (int i = 0; i < DIM; ++i)
      {
	A(0,i+1) += M[i + 1] * N.m();
	for (int j = 0; j < DIM; ++j)
	  {
	    A(i+1,i+1) += data.visc * M[j + 1] * N[j + 1];
	    if (data.fulltensor)
	      A(i + 1,j+1) += data.visc * M[i + 1] * N[j + 1];

	
	    A(i + 1,i+1) += (*U)[j + 1].m() * M[j + 1] * N.m();
	    A(i + 1,j+1) += M.m() * (*U)[i + 1][j + 1] * N.m();
	  }

	A(i + 1,0) -= M.m() * N[i + 1];
      }
  }
  
  ////////////////////////////////////////////////// NavierStokesLps
  
  template<int DIM>
  void NavierStokesLps<DIM>::lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const
  {
    _h = h;
    _alpha = NavierStokes<DIM>::data.alpha0 /
      (NavierStokes<DIM>::data.visc/(_h*_h) + 0.3/_h);
  }
  
  template<int DIM>
  void NavierStokesLps<DIM>::StabForm(VectorIterator b,
				      const FemFunction &U,
				      const FemFunction &UP,
				      const TestFunction &NP) const
  {
    for (int i = 0; i < DIM; ++i)
      b[0] += _alpha * UP[0][i + 1] * NP[i + 1];
  }
  
  template<int DIM>
  void NavierStokesLps<DIM>::StabMatrix(EntryMatrix &A,
					const FemFunction &U,
					const TestFunction &Np,
					const TestFunction &Mp) const
  {
    for (int i = 0; i < DIM; ++i)
      A(0, 0) += _alpha * Mp[i + 1] * Np[i + 1];
  }

  ////////////////////////////////////////////////// NavierStokesTangeltLps
  
  template<int DIM>
  void NavierStokesTangentLps<DIM>::lpspoint(double h, const FemFunction &W, const Vertex<DIM> &v) const
  {
    _h = h;
    _alpha = NavierStokesTangent<DIM>::data.alpha0 /
      (NavierStokesTangent<DIM>::data.visc/(_h*_h) + 0.3/_h);
  }
  
  template<int DIM>
  void NavierStokesTangentLps<DIM>::StabForm(VectorIterator b,
					     const FemFunction &W,
					     const FemFunction &WP,
					     const TestFunction &NP) const
  {
    for (int i = 0; i < DIM; ++i)
      b[0] += _alpha * WP[0][i + 1] * NP[i + 1];
  }
  
  template<int DIM>
  void NavierStokesTangentLps<DIM>::StabMatrix(EntryMatrix &A,
					       const FemFunction &W,
					       const TestFunction &Np,
					       const TestFunction &Mp) const
  {
    for (int i = 0; i < DIM; ++i)
      A(0, 0) += _alpha * Mp[i + 1] * Np[i + 1];
  }
  
  ////////////////////////////////////////////////// Navier Stokes time

  template<int DIM>
  void
  NavierStokesLpsTime<DIM>::Form(VectorIterator b, const FemFunction &U, const TestFunction &N) const
  {      
    if (NSSCHEME == 0)
      {
	for (int i=0;i<DIM;++i)
	  {
	    b[0] += U[i + 1][i + 1] * N.m();
	  
	    b[i+1] += (U[i + 1].m() - (*OLD)[i+1].m()) / NavierStokes<DIM>::data.dt * N.m();
	  
	    for (int j = 0; j < DIM; ++j)
	      {
		b[i + 1] += NavierStokes<DIM>::data.theta * NavierStokes<DIM>::data.visc * U[i + 1][j + 1] * N[j + 1];
		if (NavierStokes<DIM>::data.fulltensor)
		  b[i + 1] += NavierStokes<DIM>::data.theta * NavierStokes<DIM>::data.visc * U[j + 1][i + 1] * N[j + 1];
	      
		b[i + 1] += NavierStokes<DIM>::data.theta * U[j + 1].m() * U[i + 1][j + 1] * N.m();
		b[i + 1] += (1.0-NavierStokes<DIM>::data.theta) * NavierStokes<DIM>::data.visc * (*OLD)[i + 1][j + 1] * N[j + 1];
		b[i + 1] += (1.0-NavierStokes<DIM>::data.theta) * (*OLD)[j + 1].m() * (*OLD)[i + 1][j + 1] * N.m();
	      }
	    b[i + 1] -= U[0].m() * N[i + 1];
	  }
      }
    else  if (NSSCHEME == 1)
      {
	// BDF
	for (int i=0;i<DIM;++i)
	  {
	    b[0] += U[i + 1][i + 1] * N.m();
	
	    // BDF 4:	    b[i+1] += (25.0/12.0*U[i + 1].m() +  (*OLD)[i+1].m()) / NavierStokes<DIM>::data.dt * N.m(); // in OLD steht die KOmbination der alten Werte 
	    // BDF 3: 
	    b[i+1] += (11.0/6.0*U[i + 1].m() +  (*OLD)[i+1].m()) / NavierStokes<DIM>::data.dt * N.m(); // in OLD steht die KOmbination der alten Werte 
	
	    for (int j = 0; j < DIM; ++j)
	      {
		b[i + 1] += NavierStokes<DIM>::data.visc * U[i + 1][j + 1] * N[j + 1];
		if (NavierStokes<DIM>::data.fulltensor)
		  b[i + 1] += NavierStokes<DIM>::data.visc * U[j + 1][i + 1] * N[j + 1];
	    
		b[i + 1] += U[j + 1].m() * U[i + 1][j + 1] * N.m();
	      }
	    b[i + 1] -= U[0].m() * N[i + 1];
	  }
      }
    else abort();
  }
  
  template<int DIM>
  void NavierStokesLpsTime<DIM>::Matrix(EntryMatrix &A,
					const FemFunction &U,
					const TestFunction &M,
					const TestFunction &N) const
  {
    if (NSSCHEME == 0)
      {
	// THETA
	for (int i = 0; i < DIM; ++i)
	  {
	    A(i+1,i+1) += M.m() * N.m() / NavierStokes<DIM>::data.dt;
	    A(0, i + 1) += M[i + 1] * N.m();
	    for (int j = 0; j < DIM; ++j)
	      {
		A(i + 1, i + 1) += NavierStokes<DIM>::data.theta * NavierStokes<DIM>::data.visc * M[j + 1] * N[j + 1];
		if (NavierStokes<DIM>::data.fulltensor)
		  A(i + 1, j + 1) += NavierStokes<DIM>::data.theta * NavierStokes<DIM>::data.visc * M[i + 1] * N[j + 1];
	    
		A(i + 1, j + 1) += NavierStokes<DIM>::data.theta * M.m() * U[i + 1][j + 1] * N.m();
		A(i + 1, i + 1) += NavierStokes<DIM>::data.theta * U[j + 1].m() * M[j + 1] * N.m();
	      }
	
	    A(i + 1, 0) -= M.m() * N[i + 1];
	  }
      }
    else if (NSSCHEME == 1)
      {
	// BDF-4
	for (int i = 0; i < DIM; ++i)
	  {
	    //	    A(i+1,i+1) += 25./12. * M.m() * N.m() / NavierStokes<DIM>::data.dt;
	    A(i+1,i+1) += 11./6. * M.m() * N.m() / NavierStokes<DIM>::data.dt;
	    A(0, i + 1) += M[i + 1] * N.m();
	    for (int j = 0; j < DIM; ++j)
	      {
		A(i + 1, i + 1) += NavierStokes<DIM>::data.visc * M[j + 1] * N[j + 1];
		if (NavierStokes<DIM>::data.fulltensor)
		  A(i + 1, j + 1) += NavierStokes<DIM>::data.visc * M[i + 1] * N[j + 1];
	    
		A(i + 1, j + 1) += M.m() * U[i + 1][j + 1] * N.m();
		A(i + 1, i + 1) += U[j + 1].m() * M[j + 1] * N.m();
	      }
	
	    A(i + 1, 0) -= M.m() * N[i + 1];
	  }
      }
  }

  ////////////////////////////////////////////////// BOUNDARY

  template<int DIM>
  void NavierStokesBoundary<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
  {
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	b[i+1] -= data.visc * U[j+1][i+1] * _n[j] * N.m();
  }
  template<int DIM>
  void NavierStokesBoundary<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
  {
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	A(i+1,j+1) -= data.visc * M[i+1] * _n[j] * N.m();
  }
  template<int DIM>
  void NavierStokesBoundary<DIM>::pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const
  {
    _n = n;
  }


  ////////////////////////////////////////////////// BOUNDARY Tangent

  template<int DIM>
  void NavierStokesTangentBoundary<DIM>::Form(VectorIterator b, const FemFunction& W, const TestFunction& N, int col) const
  {
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	b[i+1] -= data.visc * W[j+1][i+1] * _n[j] * N.m();
  }
  template<int DIM>
  void NavierStokesTangentBoundary<DIM>::Matrix(EntryMatrix& A, const FemFunction& W, const TestFunction& M, const TestFunction& N, int col) const
  {
    for (int i=0;i<DIM;++i)
      for (int j=0;j<DIM;++j)
	A(i+1,j+1) -= data.visc * M[i+1] * _n[j] * N.m();
  }
  template<int DIM>
  void NavierStokesTangentBoundary<DIM>::pointboundary(double h, const FemFunction& W, const Vertex<DIM>& v, const Vertex<DIM>& n) const
  {
    _n = n;
  }



  ////////////////////////////////// Boundary Time



  
  template<int DIM>
  void NavierStokesBoundaryTime<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
  {
    if (NSSCHEME == 0)
      {
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    {
	      b[i+1] -= data.theta * data.visc * U[j+1][i+1] * _n[j] * N.m();
	      b[i+1] -= (1.0-data.theta) * data.visc * (*OLD)[j+1][i+1] * _n[j] * N.m();
	    }
      }
    else if (NSSCHEME == 1)
      {
	// BDF-4
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    {
	      b[i+1] -= data.visc * U[j+1][i+1] * _n[j] * N.m();
	    }
      }
  }
  template<int DIM>
  void NavierStokesBoundaryTime<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
  {
    if (NSSCHEME == 0)
      {
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    A(i+1,j+1) -= data.theta * data.visc * M[i+1] * _n[j] * N.m();
      }
    else
      {
	// BDF-4
	for (int i=0;i<DIM;++i)
	  for (int j=0;j<DIM;++j)
	    A(i+1,j+1) -= data.visc * M[i+1] * _n[j] * N.m();
      }
  }
  template<int DIM>
  void NavierStokesBoundaryTime<DIM>::pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const
  {
    _n = n;
  }


  
  template class NavierStokes<2>;
  template class NavierStokes<3>;

  template class NavierStokesLps<2>;
  template class NavierStokesLps<3>;

  template class  NavierStokesLpsTime<2>;
  template class  NavierStokesLpsTime<3>;

  template class  NavierStokesBoundaryTime<2>;
  template class  NavierStokesBoundaryTime<3>;

  template class  NavierStokesBoundary<2>;
  template class  NavierStokesBoundary<3>;


  // Tangent

  template class NavierStokesTangent<2>;
  template class NavierStokesTangent<3>;

  template class NavierStokesTangentLps<2>;
  template class NavierStokesTangentLps<3>;

  template class  NavierStokesTangentBoundary<2>;
  template class  NavierStokesTangentBoundary<3>;


} // namespace Gascoigne
