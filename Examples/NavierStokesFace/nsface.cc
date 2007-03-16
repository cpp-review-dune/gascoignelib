#include "nsface.h"
#include "filescanner.h"


using namespace std;

namespace Gascoigne
{
  
  
  NavierStokesFace2d::NavierStokesFace2d(const ParamFile* pf) : FaceEquation()
  {
    DataFormatHandler DFH;
    DFH.insert("visc"   , &__visc , 1.);
    DFH.insert("alpha"  , &__alpha0,0.5);
    DFH.insert("delta"  , &__delta0, 0.);
    
    FileScanner FS(DFH, pf, "Equation");
  }

  // --------------------------------------------------

  void NavierStokesFace2d::point_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex2d& v,const Vertex2d& n) const 
  {
    __n=n;
    __alpha = __alpha0 / (1+__visc) * h*h;
  }
  
  // --------------------------------------------------

  void NavierStokesFace2d::FaceForm(VectorIterator b, const FemFunction& U1,const FemFunction& U2, const TestFunction& N1,const TestFunction& N2) const
  {
    double un=0,pn=0;
    for (int i=0;i<2;++i)
      {
	un += __n[i]*(U1[0][i+1]-U2[0][i+1]);
	pn += __n[i]*(N1[i+1]   -N2[i+1]);
      }
    
    b[0] += __alpha*un*pn;
  }
  
  
  void NavierStokesFace2d::FaceMatrix(EntryMatrix& A, const FemFunction& U1,const FemFunction& U2, const TestFunction& M1,const TestFunction& M2, const TestFunction& N1, const TestFunction& N2) const
  {
    double un=0,pn=0;
    for (int i=0;i<2;++i)
      {
	un += __n[i]*(M1[i+1]-M2[i+1]);
	pn += __n[i]*(N1[i+1]-N2[i+1]);
      }
    
    A(0,0) += __alpha*un*pn;
  }
  

}
