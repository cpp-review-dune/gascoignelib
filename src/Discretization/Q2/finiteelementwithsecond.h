#ifndef  __FiniteElementWithSecond_h
#define  __FiniteElementWithSecond_h

#include  "finiteelement.h"

namespace Gascoigne
{
template<int DIM, int BDIM, class TRAFO, class BASE>
class FiniteElementWithSecond : public FiniteElement<DIM,BDIM,TRAFO,BASE>
{
  protected:
    
    typedef  FemInterface::Matrix   Matrix;
    
    mutable nvector<Matrix> hesse;
    
  public:

    void ComputeHesse(const Vertex2d& xi) const;
    void ComputeHesse(const Vertex3d& xi) const;

    std::string GetName() const {return "FiniteElementWithSecond";}
    
    FiniteElementWithSecond();
    
    void point(const Vertex<DIM>& v) const
    {
      FiniteElement<DIM,BDIM,TRAFO,BASE>::point(v);
      ComputeHesse(v);
    }

    void  init_test_functions(TestFunction& Phi, double w, int i) const
    {
      FiniteElement<DIM,BDIM,TRAFO,BASE>::init_test_functions(Phi,w,i);
      init_test_hesse(Phi,w,i);
    }
    void init_test_hesse(TestFunction& N, double w, int i) const;
};
}

#endif

