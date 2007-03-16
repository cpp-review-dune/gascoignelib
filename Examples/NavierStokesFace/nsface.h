/*----------------------------   nsface.h     ---------------------------*/
/*      $Id$                 */
#ifndef __nsface_H
#define __nsface_H
/*----------------------------   nsface.h     ---------------------------*/

#include "faceequation.h"
#include "paramfile.h"


using namespace std;

namespace Gascoigne
{
  
  class NavierStokesFace2d : public FaceEquation
    {
    protected:
      double __alpha0,__delta0,__visc;
      mutable double __alpha,__delta;
      
      mutable Vertex2d __n;
      
      
    public:
      
      NavierStokesFace2d(const ParamFile* pf);

      std::string GetName() const { return "NavierStokesFace2d";}
      int  GetNcomp() const { return 3; }
      
      void point_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex2d& v,const Vertex2d& n) const;
            
      void FaceForm(VectorIterator b, const FemFunction& U1,const FemFunction& U2, const TestFunction& N1,const TestFunction& N2) const;
      void FaceMatrix(EntryMatrix& A, const FemFunction& U1,const FemFunction& U2, const TestFunction& M1,const TestFunction& M2, const TestFunction& N1, const TestFunction& N2) const;
      
      
    };
  
}




/*----------------------------   nsface.h     ---------------------------*/
/* end of #ifndef __nsface_H */
#endif
/*----------------------------   nsface.h     ---------------------------*/
