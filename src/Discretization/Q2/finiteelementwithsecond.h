#ifndef  __FiniteElementWithSecond_h
#define  __FiniteElementWithSecond_h

#include  "finiteelement.h"

template<int DIM, int BDIM, class TRAFO, class BASE>
class FiniteElementWithSecond : public Gascoigne::FiniteElement<DIM,BDIM,TRAFO,BASE>
{
  protected:
    
    typedef  Gascoigne::FemInterface::Matrix   Matrix;
    
    mutable Gascoigne::nvector<Matrix> hesse;
    
  public:

    void ComputeHesse(const Gascoigne::Vertex2d& xi) const;
    void ComputeHesse(const Gascoigne::Vertex3d& xi) const;

    std::string GetName() const {return "FiniteElementWithSecond";}
    
    FiniteElementWithSecond();
    
    void point(const Gascoigne::Vertex<DIM>& v) const
    {
      Gascoigne::FiniteElement<DIM,BDIM,TRAFO,BASE>::point(v);
      ComputeHesse(v);
    }

    void  init_test_functions(Gascoigne::TestFunction& Phi, double w, int i) const
    {
      Gascoigne::FiniteElement<DIM,BDIM,TRAFO,BASE>::init_test_functions(Phi,w,i);
      init_test_hesse(Phi,w,i);
    }
    void init_test_hesse(Gascoigne::TestFunction& N, double w, int i) const;
};

#endif

