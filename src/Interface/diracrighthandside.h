#ifndef __DiracRightHandSide_h
#define __DiracRightHandSide_h

#include "righthandsidedata.h"
#include  <vector>
#include  "vertex.h"

/**********************************************************/
namespace Gascoigne
{
  class DiracRightHandSide : public virtual RightHandSideData
  {
  protected:
    std::vector<Vertex2d>  _v2d;
    std::vector<Vertex3d>  _v3d;
    
    std::vector<int>  _comps;
    
  public:
    DiracRightHandSide() : RightHandSideData() { }
    virtual ~DiracRightHandSide() { }

    virtual void BasicInit(const std::vector<Vertex2d>& v2d, const std::vector<int>& comps) {_v2d=v2d;_comps=comps;}
    virtual void BasicInit(const std::vector<Vertex3d>& v3d, const std::vector<int>& comps) {_v3d=v3d;_comps=comps;}
  
    virtual const std::vector<Vertex2d>& GetPoints2d() const { return _v2d;}
    virtual const std::vector<Vertex3d>& GetPoints3d() const { return _v3d;}
  
    virtual const std::vector<int>& GetComps()    const { return _comps;}

    virtual double operator()(int i, const Vertex2d& v) const {assert(0); return 0;}
    virtual double operator()(int i, const Vertex3d& v) const {assert(0); return 0;}

    virtual void operator()(int i, VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
      {
        b[_comps[i]] += N.m()* (*this)(i,v);
      }
    virtual void operator()(int i, VectorIterator b, const TestFunction& N, const Vertex3d& v) const 
      {
        b[_comps[i]] += N.m()* (*this)(i,v);
      }
};
}
/**********************************************************/

#endif
