#ifndef __RobinData_h
#define __RobinData_h

#include "application.h"

namespace Gascoigne
{

/**********************************************************/

  class RobinData : public Application
  {
    private:

    protected:

    public:
      RobinData() : Application() { }
      ~RobinData() { }

      virtual void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const=0;
      virtual void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N, 
          int col) const=0;

      virtual void pointboundary(double h, const FemFunction& U, const Vertex2d& v, 
         const Vertex2d& n) const {}
      virtual void pointboundary(double h, const FemFunction& U, const Vertex3d& v, 
         const Vertex3d& n) const {}
      virtual void pointboundary(double h, const FemFunction& U, const FemData& Q, const Vertex2d& v, 
         const Vertex2d& n) const {
        pointboundary(h,U,v,n);
      }
      virtual void pointboundary(double h, const FemFunction& U, const FemData& Q, const Vertex3d& v, 
         const Vertex3d& n) const {
        pointboundary(h,U,v,n);
      }
  };

/**********************************************************/

}

#endif
