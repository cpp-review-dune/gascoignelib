#include  "fsi.h"
#include  "filescanner.h"


using namespace std;

extern double THETA,TIME,DT;


/*-----------------------------------------*/

namespace Gascoigne
{
  FSI::FSI(const ParamFile* pf)
  {
    OLD = NULL;
    
    DataFormatHandler DFH;
    DFH.insert("rho_f" ,    &rho_f , 0.0);
    DFH.insert("nu_f" ,    &nu_f , 0.0);
    DFH.insert("lps" ,    &lps0 , 0.0);
    DFH.insert("P2min" ,    &P2min , 0.0);
    DFH.insert("P2max" ,    &P2max , 0.0);
    DFH.insert("P2period" ,    &P2period , 0.0);
    DFH.insert("P4" ,    &P4 , 0.0);
    DFH.insert("P7" ,    &P7 , 0.0);
    FileScanner FS(DFH, pf, "Equation");
    assert(rho_f>0);
    assert(nu_f>0);
    assert(lps0>0);
    
  }

  ////////////////////////////////////////////////// 



  void FSI::point(double h, const FemFunction& U, const Vertex<3>& v) const
  {
    // stabilization in fluid domain
    lps = lps0 * h * h / nu_f;
  }
  

  void FSI::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    for (int i=0;i<3;++i)
      {
	// time 
	b[i+1] += (U[i+1].m() - (*OLD)[i+1].m())/DT * N.m();
	
	// tensor
	for (int j=0;j<3;++j)
	  {
	    b[i+1] += THETA * rho_f * nu_f * U[i+1][j+1]*N[j+1];
	    b[i+1] += (1.0-THETA) * rho_f * nu_f * (*OLD)[i+1][j+1]*N[j+1];
	  }
	

	// conv
	for (int j=0;j<3;++j)
	  {
	    b[i+1] += THETA * rho_f * U[j+1].m() * U[i+1][j+1] * N.m();
	    b[i+1] += (1.0-THETA) *  rho_f * (*OLD)[j+1].m() * (*OLD)[i+1][j+1] * N.m();
	  }
	
	
	
	// pressure
	b[i+1] -= U[0].m() * N[i+1];

	// divergence
	b[0] += rho_f * U[i+1][i+1] * N.m();
      }
  }
  

  void FSI::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    for (int i=0;i<3;++i)
      {
	// time 
	A(i+1,i+1) += M.m() /DT * N.m();

	// tensor
	for (int j=0;j<3;++j)
	  A(i+1,i+1) += THETA * rho_f * nu_f * M[j+1]*N[j+1];
	
	// conv
	for (int j=0;j<3;++j)
	  {
	    A(i+1,j+1) += THETA * rho_f * M.m() * U[i+1][j+1] * N.m();
	    A(i+1,i+1) += THETA * rho_f * U[j+1].m() * M[j+1] * N.m();
	  }
	

	// pressure
	A(i+1,0) -= M.m() * N[i+1];
	
	// divergence
	A(0,i+1) += rho_f * M[i+1] * N.m();
      }
  }
  

  void FSI::Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
  {
    if (col==4) 
      for (int i=0;i<3;++i)
	b[i+1] += normal[i]*P4*N.m();
      
  }
  
  void FSI::Matrix(EntryMatrix& E, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const
  {
  }  
    
  void FSI::pointboundary(double h, const FemFunction& U, const Vertex3d& v, const Vertex3d& n) const
  {
    normal = n;

    // pulsierende bedingung an p2
    
    //P2 = P2min +  0.5 *  (P2max - P2min) * (1.0-cos(M_PI*2.0/P2period * TIME));
    
    
  }
  




  void FSI::lpspoint(double h, const FemFunction& U, const Vertex<3>& v) const
  {
    double vel = 1.0;
    
    lps = lps0 / (vel/h + nu_f/h/h);
  } 
  void FSI::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
  {
    for (int i=0;i<3;++i)
      b[0] += lps * UP[0][i+1] * N[i+1];

    double cN = 0;
    for (int i=0;i<3;++i)
      cN += U[i+1].m() * N[i+1];
    for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
	b[i+1] += lps * U[j+1].m() * UP[i+1][j+1] * cN;
  }
  
  void FSI::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
  {
    for (int i=0;i<3;++i)
      A(0,0) += lps *Mp[i+1] * Np[i+1];
    double cN = 0;
    for (int i=0;i<3;++i)
      cN += U[i+1].m() * Np[i+1];
    double cM = 0;
    for (int i=0;i<3;++i)
      cN += U[i+1].m() * Mp[i+1];
    for (int i=0;i<3;++i)
      A(i+1,i+1) += lps * cM * cN;
  }
  


}

/*-----------------------------------------*/
