/*----------------------------   faceequation.h     ---------------------------*/
/*      $Id$                 */
#ifndef __faceequation_H
#define __faceequation_H
/*----------------------------   faceequation.h     ---------------------------*/

#include  "entrymatrix.h"
#include  "vertex.h"
#include  "application.h"


/*-------------------------------------------------------------------------*/

namespace Gascoigne
{
  
  //////////////////////////////////////////////
  ///
  ///@brief
  /// Interface class for Equation

  ///
  ///
  //////////////////////////////////////////////

  class FaceEquation : public virtual Application
    {
    private:
      
    protected:
      
    public:
      //
      // Constructors
      //
      FaceEquation() : Application() {}
      virtual ~FaceEquation() {}


      virtual void point_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex2d& v,const Vertex2d& n) const {}
      virtual void point_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex3d& v,const Vertex3d& n) const {}
     
      virtual void pointmatrix_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex2d& v,const Vertex2d& n) const
	{
	  point_face(h,U1,U2,v,n);
	}
      virtual void pointmatrix_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex3d& v,const Vertex3d& n) const
	{
	  point_face(h,U1,U2,v,n);
	}
      
      virtual int GetNcomp() const=0;
      virtual void FaceForm(VectorIterator b, const FemFunction& U1, const FemFunction& U2, const TestFunction& N1, const TestFunction& N2) const=0;
      virtual void FaceMatrix(EntryMatrix& A, const FemFunction& U1, const FemFunction& U2, const TestFunction& M1, const TestFunction& M2, const TestFunction& N1, const TestFunction& N2) const=0;
  };
}

/*----------------------------   faceequation.h     ---------------------------*/
/* end of #ifndef __faceequation_H */
#endif
/*----------------------------   faceequation.h     ---------------------------*/
