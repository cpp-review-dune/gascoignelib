#ifndef  __DomainFunctional_h
#define  __DomainFunctional_h


#include  "functional.h"
#include  "gascoigne.h"
#include  "vertex.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments DomainFunctional

  ///
  ///
  /////////////////////////////////////////////

  class DomainFunctional : public virtual Functional
  {
    private:

    protected:

    public:
      DomainFunctional() : Functional() {}
      virtual ~DomainFunctional() {}

      virtual double J(const FemFunction& U, const Vertex2d& v) const {
        std::cerr << "\"DomainFunctional::J\" not written" << std::endl; 
        abort();
      }

      virtual double J(const FemFunction& U, const Vertex3d& v) const {
        std::cerr << "\"DomainFunctional::J\" not written" << std::endl; 
        abort();
      }
  };
}

#endif
