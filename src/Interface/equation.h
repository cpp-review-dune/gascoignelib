#ifndef __Equation_h
#define __Equation_h

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

  class MeshInterface;
  
  class Equation : public virtual Application
  {
    private:

    protected:

    public:
      //
      // Constructors
      //
      Equation() : Application() {}
      virtual ~Equation() {}

      virtual void OperatorStrong(DoubleVector& b, const FemFunction& U) const {
        std::cerr << "\"Equation::OperatorStrong\" not written!" << std::endl;
        abort();
      } 
      virtual void SetTimePattern(TimePattern& TP) const {
        std::cerr << "\"Equation::SetTimePattern\" not written!" << std::endl;
        abort();
      } 

      virtual void point(double h, const FemFunction& U, const Vertex2d& v) const {}
      virtual void point(double h, const FemFunction& U, const Vertex3d& v) const {}
     
      virtual void pointmatrix(double h, const FemFunction& U, const Vertex2d& v) const {
        point(h,U,v);
      }
      virtual void pointmatrix(double h, const FemFunction& U, const Vertex3d& v) const {
        point(h,U,v);
      }
     
      virtual int GetNcomp() const=0;
      virtual void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const=0;
      virtual void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const=0;
      // .cell is analogous to .point
      //virtual void cell(const MeshInterface* p_mesh, int cell_id, const LocalVector& U, const LocalNodeData& Q) const { }
  };
}

#endif
