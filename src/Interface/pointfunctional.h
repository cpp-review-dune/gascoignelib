#ifndef __PointFunctional_h
#define __PointFunctional_h

#include "functional.h"
#include  <vector>
#include  "vertex.h"

/**********************************************************/

namespace Gascoigne
{
  class PointFunctional : public virtual Functional
  {
    private:
      
    protected:
      std::vector<Vertex2d>  _v2d;
      std::vector<Vertex3d>  _v3d;

      std::vector<int>  _comps;

    public:
      PointFunctional() : Functional() {}
      virtual ~PointFunctional() {}
  
      virtual void BasicInit(const std::vector<Vertex2d>& v2d, const std::vector<int>& comps) {_v2d=v2d;_comps=comps;}
      virtual void BasicInit(const std::vector<Vertex3d>& v3d, const std::vector<int>& comps) {_v3d=v3d;_comps=comps;}
  
      virtual const std::vector<Vertex2d>& GetPoints2d() const { return _v2d;}
      virtual const std::vector<Vertex3d>& GetPoints3d() const { return _v3d;}
  
      virtual const std::vector<int>& GetComps() const { return _comps;}

      virtual double J(const std::vector<double>& u) const {
        std::cerr << "\"PointFunctional::J\" not written" << std::endl; 
        abort();
      }
  };
}

/**********************************************************/

#endif
