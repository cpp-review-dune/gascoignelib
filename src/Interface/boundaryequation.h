#ifndef __BoundaryEquation_h
#define __BoundaryEquation_h

#include "application.h"
#include "vertex.h"
#include "entrymatrix.h"

namespace Gascoigne
{

/**********************************************************/

  class BoundaryEquation : public virtual Application
  {
    private:

    protected:

    public:
      BoundaryEquation() : Application() { }
      ~BoundaryEquation() { }

      virtual int GetNcomp() const=0;

      virtual void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const=0;
      virtual void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const=0;

      virtual void pointboundary(double h, const FemFunction& U, const Vertex2d& v, const Vertex2d& n) const {}
      virtual void pointboundary(double h, const FemFunction& U, const Vertex3d& v, const Vertex3d& n) const {}
      virtual void pointmatrixboundary(double h, const FemFunction& U, const Vertex2d& v, 
         const Vertex2d& n) const {
        pointboundary(h,U,v,n);
      }
      virtual void pointmatrixboundary(double h, const FemFunction& U, const Vertex3d& v, 
         const Vertex3d& n) const {
        pointboundary(h,U,v,n);
      }
  };

/**********************************************************/

}

#endif
