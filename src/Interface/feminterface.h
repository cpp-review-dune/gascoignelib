#ifndef  __FemInterface_h
#define  __FemInterface_h

#include  "vertex.h"
#include  "nmatrix.h"
#include  "gascoigne.h"
#include  <string>

/*-----------------------------------------*/

namespace Gascoigne
{
  class FemInterface
  {
    private:

    protected:

    public:
      typedef nmatrix<double>    Matrix;

      FemInterface() {}
      virtual ~FemInterface() {}

      virtual std::string GetName() const=0;

      virtual int    n() const=0;
      virtual double J() const=0;
      virtual double G() const=0;

      virtual void x(Vertex2d& v) const {
        std::cerr << "\"FemInterface::x\" not written!" << std::endl;
        abort();
      }
      virtual void x(Vertex3d& v) const {
        std::cerr << "\"FemInterface::x\" not written!" << std::endl;
        abort();
      }

      virtual void normal(Vertex2d& v) const {
        std::cerr << "\"FemInterface::normal\" not written!" << std::endl;
        abort();
      }
      virtual void normal(Vertex3d& v) const {
        std::cerr << "\"FemInterface::normal\" not written!" << std::endl;
        abort();
      }

      virtual void point(const Vertex2d& v) const {
        std::cerr << "\"FemInterface::point\" not written!" << std::endl;
        abort();
      }
      virtual void point(const Vertex3d& v) const {
        std::cerr << "\"FemInterface::point\" not written!" << std::endl;
        abort();
      }

      virtual void  point_boundary(int ie, const Vertex1d& v) const {
        std::cerr << "\"FemInterface::point_boundary\" not written!" << std::endl;
        abort();
      }
      virtual void  point_boundary(int ie, const Vertex2d& v) const {
        std::cerr << "\"FemInterface::point_boundary\" not written!" << std::endl;
        abort();
      }

      virtual void ReInit(const Matrix& M) const=0;
      virtual void  init_test_functions(TestFunction& Phi, double w, int i) const=0;
  };
}

#endif
