#ifndef  __dirichletdata_h
#define  __dirichletdata_h

#include  "vertex.h"
#include  <set>
#include  "nvector.h"
#include  <string>
#include  "application.h"

//////////////////////////////////////////////
///
///@brief
/// Interface class for Dirichlet Boundary Conditions

/// void operator()(Vector& b, const Vertex2d& v, int col)
/// gets the coordinate v and color of boundarypart "col" and 
/// sets the values of b. b is a vector of length ncomp
///
//////////////////////////////////////////////


/*-----------------------------------------*/

namespace Gascoigne
{
  class DirichletData : public virtual Application
  {
    private:

    protected:

    public:
      DirichletData() : Application() {}
      virtual ~DirichletData() {}

      virtual void operator()(DoubleVector& b, const Vertex2d& v, int col) const {
        std::cerr << "\"DirichletData::operator()\" not written!" << std::endl;
        abort();
      }

      virtual void operator()(DoubleVector& b, const Vertex3d& v, int col) const {
        std::cerr << "\"DirichletData::operator()\" not written!" << std::endl;
        abort();
      }

      virtual std::set<int> preferred_colors()const {
        return std::set<int>();
      }
  };
}

#endif
